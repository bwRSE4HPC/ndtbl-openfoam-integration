#ifndef TABLESOLVER_H_
#define TABLESOLVER_H_

#include "dictionary.H"
#include "error.H"
#include "fileName.H"
#include "List.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "wordList.H"

#include "ndtbl/ndtbl.hpp"

#include <array>
#include <cstdlib>
#include <cstddef>
#include <limits.h>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

namespace Foam
{
namespace combustionModels
{

namespace detail
{

/// Read a dictionary entry as an OpenFOAM fileName and return its path string.
inline std::string readPathEntry(const dictionary& dict, const word& key)
{
    fileName pathValue;
    dict.lookup(key) >> pathValue;
    return std::string(pathValue.c_str());
}

/// Read the ndtbl field name for an output, falling back to the output name.
inline std::string readFieldEntry(const dictionary& dict, const word& defaultName)
{
    word fieldName(defaultName);
    if (dict.found("field"))
    {
        dict.lookup("field") >> fieldName;
    }

    return std::string(fieldName.c_str());
}

/// Return the current host name for per-rank diagnostics.
inline std::string hostName()
{
    char name[256];
    name[0] = '\0';

    if (gethostname(name, sizeof(name)) != 0)
    {
        return "unknown";
    }

    name[sizeof(name) - 1] = '\0';
    return std::string(name);
}

/// Return the canonical path when available, otherwise preserve the input.
inline std::string realPath(const std::string& path)
{
    char resolved[PATH_MAX];
    if (realpath(path.c_str(), resolved) == nullptr)
    {
        return path;
    }

    return std::string(resolved);
}

} // End namespace detail

/// Loads ndtbl table files and evaluates configured output fields.
///
/// The solver is parameterized by the number of lookup coordinates and output
/// values. Its constructor reads the model coefficients, loads each unique
/// ndtbl path once per solver instance, and binds the requested output table
/// names to fields in the loaded table groups. A lookup evaluates every loaded
/// group at a query point and copies the bound field values into the
/// caller-provided output array.
template<std::size_t Dim>
class tableSolver
{
private:

    const std::size_t NOut;

    /// Interpolation schemes supported by the ndtbl runtime evaluator.
    enum class interpolationMethod
    {
        linear,
        cubic
    };

    /// Mapping from one requested output slot to one field in an ndtbl group.
    struct fieldBinding
    {
        /// Index in the caller-visible output value array.
        std::size_t outputIndex;

        /// Index of the corresponding field inside the ndtbl RuntimeFieldGroup.
        std::size_t fieldIndex;
    };

    /// One loaded ndtbl path and all output fields that are read from it.
    struct tableGroup
    {
        /// Path used to load the table group.
        std::string filePath;

        /// Runtime representation of the ndtbl axes and field data.
        ndtbl::RuntimeFieldGroup<Dim, scalar> group;

        /// Output-to-field bindings served by this loaded table group.
        std::vector<fieldBinding> bindings;
    };

    /// Unique ndtbl table groups needed by the configured output file paths.
    std::vector<tableGroup> groups_;

    /// Interpolation scheme used for all table evaluations.
    interpolationMethod interpolationMethod_;

    /// Return reusable per-thread storage for full-field interpolation results.
    static std::vector<scalar>& threadScratchValues()
    {
        static thread_local std::vector<scalar> scratchValues;
        return scratchValues;
    }

    /// Parse the interpolation scheme from coefficients and validate it.
    static interpolationMethod readInterpolationMethod(const dictionary& coeffs)
    {
        word method("linear");
        if (coeffs.found("interpolation"))
        {
            coeffs.lookup("interpolation") >> method;
        }

        if (method == "linear")
        {
            return interpolationMethod::linear;
        }

        if (method == "cubic")
        {
            return interpolationMethod::cubic;
        }

        FatalErrorInFunction
            << "Unknown ndtbl interpolation method " << method
            << ". Valid values are 'linear' and 'cubic'."
            << exit(FatalError);

        return interpolationMethod::linear;
    }

    /// Ensure that every axis has enough support points for cubic lookup.
    static void validateCubicSupport(const tableGroup& group)
    {
        const std::array<ndtbl::Axis, Dim> axes = group.group.axes();
        for (std::size_t axis = 0; axis < Dim; ++axis)
        {
            if (axes[axis].size() < 4)
            {
                FatalErrorInFunction
                    << "Cubic interpolation requires at least 4 support points "
                    << "on every axis; axis " << axis
                    << " in " << group.filePath
                    << " has " << axes[axis].size()
                    << exit(FatalError);
            }
        }
    }

public:

    /// Construct a solver from model coefficients and requested table outputs.
    ///
    /// The constructor verifies that the provided parameter and output lists
    /// match the template dimensions, reads the `ndtblTables` sub-dictionary,
    /// loads each distinct configured file path once per solver instance, and
    /// records which field in each loaded file supplies each output.
    /// Configuration, load, or binding failures are reported through OpenFOAM
    /// FatalError diagnostics.
    tableSolver(const dictionary& coeffs, const wordList& tableNames, const wordList& paramNames)
    :
        NOut(tableNames.size()),
        interpolationMethod_(readInterpolationMethod(coeffs))
    {
        if (paramNames.size() != static_cast<label>(Dim))
        {
            FatalErrorInFunction
                << "Received " << paramNames.size()
                << " lookup parameters but tableSolver expects " << Dim
                << exit(FatalError);
        }

        if (tableNames.size() != static_cast<label>(NOut))
        {
            FatalErrorInFunction
                << "Received " << tableNames.size()
                << " output table names but tableSolver expects " << NOut
                << exit(FatalError);
        }

        if (!coeffs.found("ndtblTables"))
        {
            FatalErrorInFunction
                << "Missing ndtblTables configuration in model coefficients"
                << exit(FatalError);
        }

        const dictionary& tableConfigs = coeffs.subDict("ndtblTables");
        std::map<std::string, std::size_t> groupIndexByPath;
        std::vector<bool> outputBound(NOut, false);

        forAll(tableNames, outputI)
        {
            const word& tableName = tableNames[outputI];
            const std::size_t outputIndex = static_cast<std::size_t>(outputI);

            if (!tableConfigs.found(tableName))
            {
                FatalErrorInFunction
                    << "Missing ndtblTables entry for output " << tableName
                    << exit(FatalError);
            }

            const dictionary& tableConfig = tableConfigs.subDict(tableName);
            const std::string filePath = detail::readPathEntry(tableConfig, "file");
            const std::string fieldName = detail::readFieldEntry(tableConfig, tableName);

            std::size_t groupIndex = 0;
            std::map<std::string, std::size_t>::const_iterator existing =
                groupIndexByPath.find(filePath);

            if (existing == groupIndexByPath.end())
            {
                try
                {
                    tableGroup newGroup;
                    newGroup.filePath = filePath;
                    newGroup.group = ndtbl::read_runtime_field_group<Dim, scalar>(filePath);
                    if (interpolationMethod_ == interpolationMethod::cubic)
                    {
                        validateCubicSupport(newGroup);
                    }
                    groups_.push_back(newGroup);
                    groupIndex = groups_.size() - 1;
                    groupIndexByPath[filePath] = groupIndex;
                }
                catch (const std::exception& error)
                {
                    FatalErrorInFunction
                        << "Failed to load ndtbl file " << filePath
                        << ": " << error.what()
                        << exit(FatalError);
                }
            }
            else
            {
                groupIndex = existing->second;
            }

            try
            {
                fieldBinding binding;
                binding.outputIndex = outputIndex;
                binding.fieldIndex = groups_[groupIndex].group.field_index(fieldName);
                groups_[groupIndex].bindings.push_back(binding);
                outputBound[outputIndex] = true;
            }
            catch (const std::exception& error)
            {
                FatalErrorInFunction
                    << "Failed to bind field " << fieldName
                    << " from " << filePath
                    << ": " << error.what()
                    << exit(FatalError);
            }
        }

        for (std::size_t outputI = 0; outputI < NOut; ++outputI)
        {
            if (!outputBound[outputI])
            {
                FatalErrorInFunction
                    << "Output index " << outputI
                    << " is not bound to an ndtbl field"
                    << exit(FatalError);
            }
        }
    }

    /// Evaluate the configured tables at a query point.
    ///
    /// The input array supplies the `Dim` lookup coordinates. The method
    /// evaluates each loaded ndtbl group with the configured interpolation
    /// scheme and writes the `NOut` bound field values into `values`.
    void lookup(const std::array<scalar, Dim>& x, std::vector<scalar>& values) const
    {
        std::array<double, Dim> query;
        for (std::size_t axis = 0; axis < Dim; ++axis)
        {
            query[axis] = x[axis];
        }

        for (std::size_t groupI = 0; groupI < groups_.size(); ++groupI)
        {
            const tableGroup& group = groups_[groupI];
            std::vector<scalar>& scratchValues = threadScratchValues();
            const std::size_t fieldCount = group.group.field_count();
            if (scratchValues.size() < fieldCount)
            {
                scratchValues.resize(fieldCount);
            }

            if (interpolationMethod_ == interpolationMethod::linear)
            {
                group.group.evaluate_all_linear_into(query, scratchValues.data());
            }
            else
            {
                group.group.evaluate_all_cubic_into(query, scratchValues.data());
            }

            for (std::size_t bindingI = 0; bindingI < group.bindings.size(); ++bindingI)
            {
                const fieldBinding& binding = group.bindings[bindingI];
                values[binding.outputIndex] = scratchValues[binding.fieldIndex];
            }
        }
    }

    /// Report mmap payload residency for every loaded table group.
    void reportResidency(const word& label, const bool rankDiagnostics = false) const
    {
        const scalar bytesToMiB = 1.0/(1024.0*1024.0);

        for (std::size_t groupI = 0; groupI < groups_.size(); ++groupI)
        {
            const tableGroup& group = groups_[groupI];
            ndtbl::residency_info info = { false, 0, 0, 0, 0, 0.0 };

            try
            {
                info = group.group.payload_residency();
            }
            catch (const std::exception& error)
            {
                WarningInFunction
                    << "Failed to query ndtbl residency for " << group.filePath
                    << ": " << error.what() << endl;
            }

            if (rankDiagnostics)
            {
                const std::string host = detail::hostName();
                std::ostringstream rankLine;
                rankLine
                    << "ndtbl residency rank [" << label << "] group " << groupI
                    << " rank " << Pstream::myProcNo()
                    << " host " << host
                    << " file \"" << group.filePath << "\""
                    << " realpath \"" << detail::realPath(group.filePath) << "\": ";

                if (info.available)
                {
                    rankLine
                        << "resident fraction = " << info.resident_fraction
                        << ", resident MiB = "
                        << static_cast<scalar>(info.resident_bytes)*bytesToMiB
                        << ", mapped MiB = "
                        << static_cast<scalar>(info.total_pages)
                         * static_cast<scalar>(info.page_size)
                         * bytesToMiB;
                    if (info.smaps_available)
                    {
                        rankLine
                            << ", smaps mapping MiB = "
                            << static_cast<scalar>(info.smaps_mapping_bytes)
                             * bytesToMiB
                            << ", smaps RSS MiB = "
                            << static_cast<scalar>(info.smaps_rss_bytes)
                             * bytesToMiB
                            << ", smaps PSS MiB = "
                            << static_cast<scalar>(info.smaps_pss_bytes)
                             * bytesToMiB
                            << ", smaps shared MiB = "
                            << static_cast<scalar>
                               (
                                   info.smaps_shared_clean_bytes
                                 + info.smaps_shared_dirty_bytes
                               )*bytesToMiB
                            << ", smaps private MiB = "
                            << static_cast<scalar>
                               (
                                   info.smaps_private_clean_bytes
                                 + info.smaps_private_dirty_bytes
                               )*bytesToMiB;
                    }
                    else
                    {
                        rankLine << ", smaps unavailable";
                    }
                }
                else
                {
                    rankLine
                        << "unavailable; enable ndtbl mmap diagnostics "
                        << "and mmap-backed reads";
                }

                List<string> rankLines(Pstream::nProcs());
                rankLines[Pstream::myProcNo()] = rankLine.str();
                Pstream::gatherList(rankLines);

                List<string> rankHosts(Pstream::nProcs());
                rankHosts[Pstream::myProcNo()] = host.c_str();
                Pstream::gatherList(rankHosts);

                List<scalar> rankSmapsAvailable(Pstream::nProcs(), scalar(0));
                rankSmapsAvailable[Pstream::myProcNo()] =
                    info.smaps_available ? 1 : 0;
                Pstream::gatherList(rankSmapsAvailable);

                List<scalar> rankSmapsRssMiB(Pstream::nProcs(), scalar(0));
                rankSmapsRssMiB[Pstream::myProcNo()] =
                    static_cast<scalar>(info.smaps_rss_bytes)*bytesToMiB;
                Pstream::gatherList(rankSmapsRssMiB);

                List<scalar> rankSmapsPssMiB(Pstream::nProcs(), scalar(0));
                rankSmapsPssMiB[Pstream::myProcNo()] =
                    static_cast<scalar>(info.smaps_pss_bytes)*bytesToMiB;
                Pstream::gatherList(rankSmapsPssMiB);

                if (Pstream::master())
                {
                    forAll(rankLines, rankI)
                    {
                        Info<< rankLines[rankI].c_str() << endl;
                    }

                    List<bool> hostPrinted(Pstream::nProcs(), false);
                    forAll(rankHosts, hostI)
                    {
                        if (hostPrinted[hostI])
                        {
                            continue;
                        }

                        Foam::label hostRankCount = 0;
                        scalar hostSmapsAvailable = 0;
                        scalar hostRssMiB = 0;
                        scalar hostPssMiB = 0;

                        forAll(rankHosts, rankI)
                        {
                            if (rankHosts[rankI] == rankHosts[hostI])
                            {
                                hostPrinted[rankI] = true;
                                ++hostRankCount;
                                hostSmapsAvailable += rankSmapsAvailable[rankI];
                                hostRssMiB += rankSmapsRssMiB[rankI];
                                hostPssMiB += rankSmapsPssMiB[rankI];
                            }
                        }

                        Info<< "ndtbl residency host [" << label << "] group "
                            << groupI
                            << " host " << rankHosts[hostI].c_str()
                            << ": ranks = " << hostRankCount;
                        if (hostSmapsAvailable == hostRankCount)
                        {
                            Info<< ", smaps RSS MiB sum = " << hostRssMiB
                                << ", smaps PSS MiB sum = " << hostPssMiB;
                        }
                        else
                        {
                            Info<< ", smaps available on "
                                << hostSmapsAvailable << " / "
                                << hostRankCount << " ranks";
                        }
                        Info<< endl;
                    }
                }
            }

            scalar rankCount = 1;
            scalar availableCount = info.available ? 1 : 0;
            scalar residentFractionSum = info.available
                ? static_cast<scalar>(info.resident_fraction)
                : 0;
            scalar residentFractionMin = info.available
                ? static_cast<scalar>(info.resident_fraction)
                : 1;
            scalar residentFractionMax = info.available
                ? static_cast<scalar>(info.resident_fraction)
                : 0;
            scalar residentBytes = info.available
                ? static_cast<scalar>(info.resident_bytes)
                : 0;
            scalar residentBytesMin = residentBytes;
            scalar residentBytesMax = residentBytes;
            scalar totalBytesMax = info.available
                ? static_cast<scalar>(info.total_pages) * static_cast<scalar>(info.page_size)
                : 0;

            Foam::reduce(rankCount, Foam::sumOp<scalar>());
            Foam::reduce(availableCount, Foam::sumOp<scalar>());
            Foam::reduce(residentFractionSum, Foam::sumOp<scalar>());
            Foam::reduce(residentFractionMin, Foam::minOp<scalar>());
            Foam::reduce(residentFractionMax, Foam::maxOp<scalar>());
            Foam::reduce(residentBytes, Foam::sumOp<scalar>());
            Foam::reduce(residentBytesMin, Foam::minOp<scalar>());
            Foam::reduce(residentBytesMax, Foam::maxOp<scalar>());
            Foam::reduce(totalBytesMax, Foam::maxOp<scalar>());

            if (Pstream::master())
            {
                Info<< "ndtbl residency [" << label << "] group " << groupI
                    << " file " << group.filePath << ": ";

                if (availableCount != rankCount)
                {
                    Info<< "available on " << availableCount << " / " << rankCount
                        << " ranks; enable ndtbl mmap diagnostics and mmap-backed reads"
                        << endl;
                    continue;
                }

                Info<< "resident fraction min/avg/max = "
                    << residentFractionMin << " / "
                    << residentFractionSum/rankCount << " / "
                    << residentFractionMax
                    << ", resident MiB min/avg/max = "
                    << residentBytesMin*bytesToMiB << " / "
                    << (residentBytes/rankCount)*bytesToMiB << " / "
                    << residentBytesMax*bytesToMiB
                    << ", mapped MiB max = " << totalBytesMax*bytesToMiB
                    << endl;
            }
        }
    }
};

} // End combustionModels namespace
} // End Foam namespace

#endif /* TABLESOLVER_H_ */
