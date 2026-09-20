// SPDX-FileCopyrightText: 2026 Thomas Isensee
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef TABLESOLVER_H_
#define TABLESOLVER_H_

#include "dictionary.H"
#include "error.H"
#include "fileName.H"
#include "List.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "Time.H"
#include "wordList.H"

#include "ndtbl/ndtbl.hpp"

#include <array>
#include <cstdlib>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits.h>
#include <limits>
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

/// Reject characters that would corrupt the one-record-per-line TSV format.
inline void validateTsvString(const char* field, const std::string& value)
{
    if (value.find_first_of("\t\r\n") != std::string::npos)
    {
        FatalErrorInFunction
            << "Cannot write ndtbl residency TSV field " << field
            << " because it contains a tab or newline: " << value
            << exit(FatalError);
    }
}

/// Append an available value or the TSV missing-value marker.
template<class Value>
inline void appendTsvOptional
(
    std::ostringstream& row,
    const bool available,
    const Value& value
)
{
    row << '\t';
    if (available)
    {
        row << value;
    }
    else
    {
        row << "nan";
    }
}

/// Build one version-1 rank/group residency TSV record.
inline std::string residencyTsvRow
(
    const std::size_t sample,
    const Time& runTime,
    const word& label,
    const std::size_t group,
    const std::string& host,
    const std::string& file,
    const std::string& resolvedFile,
    const ndtbl::residency_info& info
)
{
    validateTsvString("label", std::string(label.c_str()));
    validateTsvString("host", host);
    validateTsvString("file", file);
    validateTsvString("realpath", resolvedFile);

    std::ostringstream row;
    row << std::setprecision(std::numeric_limits<double>::max_digits10)
        << sample
        << '\t' << runTime.value()
        << '\t' << label
        << '\t' << runTime.elapsedClockTime()
        << '\t' << runTime.elapsedCpuTime()
        << '\t' << group
        << '\t' << Pstream::myProcNo()
        << '\t' << host
        << '\t' << file
        << '\t' << resolvedFile
        << '\t' << (info.available ? 1 : 0);

    appendTsvOptional(row, info.available, info.page_size);
    appendTsvOptional(row, info.available, info.total_pages);
    appendTsvOptional(row, info.available, info.resident_pages);
    appendTsvOptional(row, info.available, info.resident_bytes);
    appendTsvOptional(row, info.available, info.resident_fraction);

    row << '\t' << (info.smaps_available ? 1 : 0);
    appendTsvOptional(row, info.smaps_available, info.smaps_mapping_bytes);
    appendTsvOptional(row, info.smaps_available, info.smaps_rss_bytes);
    appendTsvOptional(row, info.smaps_available, info.smaps_pss_bytes);
    appendTsvOptional(row, info.smaps_available, info.smaps_shared_clean_bytes);
    appendTsvOptional(row, info.smaps_available, info.smaps_shared_dirty_bytes);
    appendTsvOptional(row, info.smaps_available, info.smaps_private_clean_bytes);
    appendTsvOptional(row, info.smaps_available, info.smaps_private_dirty_bytes);

    row << '\t' << (info.smaps_locked_available ? 1 : 0);
    appendTsvOptional
    (
        row,
        info.smaps_locked_available,
        info.smaps_locked_bytes
    );

    row << '\t' << (info.smaps_vm_flags_available ? 1 : 0);
    appendTsvOptional
    (
        row,
        info.smaps_vm_flags_available,
        info.smaps_lock_enabled ? 1 : 0
    );
    appendTsvOptional
    (
        row,
        info.smaps_vm_flags_available,
        info.smaps_lock_on_fault ? 1 : 0
    );

    row << '\t' << (info.process_vmlck_available ? 1 : 0);
    appendTsvOptional
    (
        row,
        info.process_vmlck_available,
        info.process_vmlck_bytes
    );

    row << '\n';
    return row.str();
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

    /// Index of the next residency report written by this solver instance.
    mutable std::size_t residencyReportSample_;

    /// Whether the current start-time TSV has received its version/header.
    mutable bool residencyTsvInitialized_;

    /// Write gathered version-1 TSV rows from the master process.
    void writeResidencyTsv(const Time& runTime, const std::string& rows) const
    {
        if (!Pstream::master())
        {
            return;
        }

        const word startTimeName =
            runTime.timeName(runTime.startTime().value());
        const fileName outputDir =
            runTime.globalPath()/"postProcessing"/"ndtblResidency"/startTimeName;
        const fileName outputPath = outputDir/"residency.tsv";

        if (!residencyTsvInitialized_ && !Foam::mkDir(outputDir))
        {
            FatalErrorInFunction
                << "Failed to create ndtbl residency output directory "
                << outputDir << exit(FatalError);
        }

        const std::ios_base::openmode mode =
            std::ios_base::out
          | (residencyTsvInitialized_
                ? std::ios_base::app
                : std::ios_base::trunc);
        std::ofstream output(outputPath.c_str(), mode);
        if (!output.is_open())
        {
            FatalErrorInFunction
                << "Failed to open ndtbl residency output file " << outputPath
                << exit(FatalError);
        }

        if (!residencyTsvInitialized_)
        {
            output
                << "# ndtbl_residency_tsv_version=1\n"
                << "sample\ttime\tlabel\telapsed_clock_s\telapsed_cpu_s"
                << "\tgroup\trank\thost\tfile\trealpath"
                << "\tmincore_available\tpage_size_bytes\ttotal_pages"
                << "\tresident_pages\tresident_bytes\tresident_fraction"
                << "\tsmaps_available\tsmaps_mapping_bytes"
                << "\tsmaps_rss_bytes\tsmaps_pss_bytes"
                << "\tsmaps_shared_clean_bytes\tsmaps_shared_dirty_bytes"
                << "\tsmaps_private_clean_bytes\tsmaps_private_dirty_bytes"
                << "\tsmaps_locked_available\tsmaps_locked_bytes"
                << "\tsmaps_vm_flags_available\tsmaps_lock_enabled"
                << "\tsmaps_lock_on_fault\tprocess_vmlck_available"
                << "\tprocess_vmlck_bytes\n";
        }

        output << rows;
        output.flush();
        if (!output.good())
        {
            FatalErrorInFunction
                << "Failed to write ndtbl residency output file " << outputPath
                << exit(FatalError);
        }

        output.close();
        if (output.fail())
        {
            FatalErrorInFunction
                << "Failed to close ndtbl residency output file " << outputPath
                << exit(FatalError);
        }

        residencyTsvInitialized_ = true;
    }

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
        interpolationMethod_(readInterpolationMethod(coeffs)),
        residencyReportSample_(0),
        residencyTsvInitialized_(false)
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

    /// Report mmap cache residency and Linux lock state for loaded table groups.
    void reportResidency
    (
        const word& label,
        const Time& runTime,
        const bool rankDiagnostics = false
    ) const
    {
        const scalar bytesToMiB = 1.0/(1024.0*1024.0);
        const std::string host = detail::hostName();
        std::ostringstream tsvRows;

        for (std::size_t groupI = 0; groupI < groups_.size(); ++groupI)
        {
            const tableGroup& group = groups_[groupI];
            ndtbl::residency_info info;

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

            const std::string resolvedFile = detail::realPath(group.filePath);
            const std::string localTsvRow = detail::residencyTsvRow
            (
                residencyReportSample_,
                runTime,
                label,
                groupI,
                host,
                group.filePath,
                resolvedFile,
                info
            );
            List<string> rankTsvRows(Pstream::nProcs());
            rankTsvRows[Pstream::myProcNo()] = localTsvRow;
            Pstream::gatherList(rankTsvRows);
            if (Pstream::master())
            {
                forAll(rankTsvRows, rankI)
                {
                    tsvRows << rankTsvRows[rankI].c_str();
                }
            }

            const bool onFault =
                info.smaps_vm_flags_available
             && info.smaps_lock_enabled
             && info.smaps_lock_on_fault;
            const bool eager =
                info.smaps_vm_flags_available
             && info.smaps_lock_enabled
             && !info.smaps_lock_on_fault;
            const bool unlocked =
                info.smaps_vm_flags_available
             && !info.smaps_lock_enabled
             && !info.smaps_lock_on_fault;
            const char* const lockPolicy = onFault
                ? "on-fault"
                : (eager ? "eager" : (unlocked ? "unlocked" : "unknown"));

            if (groupI == 0)
            {
                scalar rankCount = 1;
                scalar vmLckAvailableCount =
                    info.process_vmlck_available ? 1 : 0;
                scalar vmLckBytesSum = info.process_vmlck_available
                    ? static_cast<scalar>(info.process_vmlck_bytes)
                    : 0;
                scalar vmLckBytesMin = vmLckBytesSum;
                scalar vmLckBytesMax = vmLckBytesSum;

                Foam::reduce(rankCount, Foam::sumOp<scalar>());
                Foam::reduce(vmLckAvailableCount, Foam::sumOp<scalar>());
                Foam::reduce(vmLckBytesSum, Foam::sumOp<scalar>());
                Foam::reduce(vmLckBytesMin, Foam::minOp<scalar>());
                Foam::reduce(vmLckBytesMax, Foam::maxOp<scalar>());

                if (Pstream::master())
                {
                    Info<< "ndtbl memory [" << label << "]:" << endl;
                    if (vmLckAvailableCount == rankCount)
                    {
                        Info<< "  process/status: VmLck MiB min/avg/max = "
                            << vmLckBytesMin*bytesToMiB << " / "
                            << (vmLckBytesSum/rankCount)*bytesToMiB << " / "
                            << vmLckBytesMax*bytesToMiB << endl;
                    }
                    else
                    {
                        Info<< "  process/status: VmLck available on "
                            << vmLckAvailableCount << " / " << rankCount
                            << " ranks" << endl;
                    }
                }
            }

            if (rankDiagnostics)
            {
                if (groupI == 0)
                {
                    std::ostringstream processLine;
                    processLine
                        << "ndtbl memory rank [" << label << "] rank "
                        << Pstream::myProcNo() << " host " << host << ": ";
                    if (info.process_vmlck_available)
                    {
                        processLine
                            << "process/status VmLck MiB = "
                            << static_cast<scalar>(info.process_vmlck_bytes)
                             * bytesToMiB;
                    }
                    else
                    {
                        processLine << "process/status VmLck unavailable";
                    }

                    List<string> processLines(Pstream::nProcs());
                    processLines[Pstream::myProcNo()] = processLine.str();
                    Pstream::gatherList(processLines);
                    if (Pstream::master())
                    {
                        forAll(processLines, rankI)
                        {
                            Info<< processLines[rankI].c_str() << endl;
                        }
                    }
                }

                std::ostringstream rankLine;
                rankLine
                    << "ndtbl residency rank [" << label << "] group " << groupI
                    << " rank " << Pstream::myProcNo()
                    << " host " << host
                    << " file \"" << group.filePath << "\""
                    << " realpath \"" << resolvedFile << "\": ";

                if (info.available)
                {
                    rankLine
                        << "cache/mincore fraction = " << info.resident_fraction
                        << ", resident MiB = "
                        << static_cast<scalar>(info.resident_bytes)*bytesToMiB
                        << ", mapped MiB = "
                        << static_cast<scalar>(info.total_pages)
                         * static_cast<scalar>(info.page_size)
                         * bytesToMiB;
                    if (info.smaps_available)
                    {
                        rankLine
                            << "; mapping/smaps mapping MiB = "
                            << static_cast<scalar>(info.smaps_mapping_bytes)
                             * bytesToMiB
                            << ", RSS MiB = "
                            << static_cast<scalar>(info.smaps_rss_bytes)
                             * bytesToMiB
                            << ", PSS MiB = "
                            << static_cast<scalar>(info.smaps_pss_bytes)
                             * bytesToMiB
                            << ", shared MiB = "
                            << static_cast<scalar>
                               (
                                   info.smaps_shared_clean_bytes
                                 + info.smaps_shared_dirty_bytes
                               )*bytesToMiB
                            << ", private MiB = "
                            << static_cast<scalar>
                               (
                                   info.smaps_private_clean_bytes
                                 + info.smaps_private_dirty_bytes
                               )*bytesToMiB;
                        if (info.smaps_locked_available)
                        {
                            rankLine
                                << ", Locked MiB = "
                                << static_cast<scalar>(info.smaps_locked_bytes)
                                 * bytesToMiB;
                            if (info.smaps_mapping_bytes != 0)
                            {
                                rankLine
                                    << ", locked fraction = "
                                    << static_cast<scalar>
                                       (info.smaps_locked_bytes)
                                     / static_cast<scalar>
                                       (info.smaps_mapping_bytes);
                            }
                        }
                        else
                        {
                            rankLine << ", Locked unavailable";
                        }
                        rankLine << ", lock policy = " << lockPolicy;
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

                List<scalar> rankSmapsLockedAvailable
                (
                    Pstream::nProcs(),
                    scalar(0)
                );
                rankSmapsLockedAvailable[Pstream::myProcNo()] =
                    info.smaps_locked_available ? 1 : 0;
                Pstream::gatherList(rankSmapsLockedAvailable);

                List<scalar> rankSmapsLockedMiB(Pstream::nProcs(), scalar(0));
                rankSmapsLockedMiB[Pstream::myProcNo()] =
                    static_cast<scalar>(info.smaps_locked_bytes)*bytesToMiB;
                Pstream::gatherList(rankSmapsLockedMiB);

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
                        scalar hostSmapsLockedAvailable = 0;
                        scalar hostLockedMiB = 0;

                        forAll(rankHosts, rankI)
                        {
                            if (rankHosts[rankI] == rankHosts[hostI])
                            {
                                hostPrinted[rankI] = true;
                                ++hostRankCount;
                                hostSmapsAvailable += rankSmapsAvailable[rankI];
                                hostRssMiB += rankSmapsRssMiB[rankI];
                                hostPssMiB += rankSmapsPssMiB[rankI];
                                hostSmapsLockedAvailable +=
                                    rankSmapsLockedAvailable[rankI];
                                hostLockedMiB += rankSmapsLockedMiB[rankI];
                            }
                        }

                        Info<< "ndtbl residency host [" << label << "] group "
                            << groupI
                            << " host " << rankHosts[hostI].c_str()
                            << ": ranks = " << hostRankCount;
                        if (hostSmapsAvailable == hostRankCount)
                        {
                            Info<< ", smaps RSS MiB logical rank sum = "
                                << hostRssMiB
                                << ", smaps PSS MiB logical rank sum = "
                                << hostPssMiB;
                        }
                        else
                        {
                            Info<< ", smaps available on "
                                << hostSmapsAvailable << " / "
                                << hostRankCount << " ranks";
                        }
                        if (hostSmapsLockedAvailable == hostRankCount)
                        {
                            Info<< ", smaps Locked MiB logical rank sum = "
                                << hostLockedMiB;
                        }
                        else
                        {
                            Info<< ", smaps Locked available on "
                                << hostSmapsLockedAvailable << " / "
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
                ? static_cast<scalar>(info.total_pages)
                 * static_cast<scalar>(info.page_size)
                : 0;
            scalar smapsAvailableCount = info.smaps_available ? 1 : 0;
            scalar smapsRssBytesSum = info.smaps_available
                ? static_cast<scalar>(info.smaps_rss_bytes)
                : 0;
            scalar smapsRssBytesMin = smapsRssBytesSum;
            scalar smapsRssBytesMax = smapsRssBytesSum;
            scalar smapsMappingBytesMax = info.smaps_available
                ? static_cast<scalar>(info.smaps_mapping_bytes)
                : 0;
            scalar smapsLockedAvailableCount =
                info.smaps_locked_available ? 1 : 0;
            scalar smapsLockedBytesSum = info.smaps_locked_available
                ? static_cast<scalar>(info.smaps_locked_bytes)
                : 0;
            scalar smapsLockedBytesMin = smapsLockedBytesSum;
            scalar smapsLockedBytesMax = smapsLockedBytesSum;
            const scalar localLockedFraction =
                info.smaps_locked_available && info.smaps_mapping_bytes != 0
                ? static_cast<scalar>(info.smaps_locked_bytes)
                 / static_cast<scalar>(info.smaps_mapping_bytes)
                : 0;
            scalar smapsLockedFractionSum = localLockedFraction;
            scalar smapsLockedFractionMin = localLockedFraction;
            scalar smapsLockedFractionMax = localLockedFraction;
            scalar onFaultCount = onFault ? 1 : 0;
            scalar eagerCount = eager ? 1 : 0;
            scalar unlockedCount = unlocked ? 1 : 0;
            scalar unknownCount = onFault || eager || unlocked ? 0 : 1;

            Foam::reduce(rankCount, Foam::sumOp<scalar>());
            Foam::reduce(availableCount, Foam::sumOp<scalar>());
            Foam::reduce(residentFractionSum, Foam::sumOp<scalar>());
            Foam::reduce(residentFractionMin, Foam::minOp<scalar>());
            Foam::reduce(residentFractionMax, Foam::maxOp<scalar>());
            Foam::reduce(residentBytes, Foam::sumOp<scalar>());
            Foam::reduce(residentBytesMin, Foam::minOp<scalar>());
            Foam::reduce(residentBytesMax, Foam::maxOp<scalar>());
            Foam::reduce(totalBytesMax, Foam::maxOp<scalar>());
            Foam::reduce(smapsAvailableCount, Foam::sumOp<scalar>());
            Foam::reduce(smapsRssBytesSum, Foam::sumOp<scalar>());
            Foam::reduce(smapsRssBytesMin, Foam::minOp<scalar>());
            Foam::reduce(smapsRssBytesMax, Foam::maxOp<scalar>());
            Foam::reduce(smapsMappingBytesMax, Foam::maxOp<scalar>());
            Foam::reduce(smapsLockedAvailableCount, Foam::sumOp<scalar>());
            Foam::reduce(smapsLockedBytesSum, Foam::sumOp<scalar>());
            Foam::reduce(smapsLockedBytesMin, Foam::minOp<scalar>());
            Foam::reduce(smapsLockedBytesMax, Foam::maxOp<scalar>());
            Foam::reduce(smapsLockedFractionSum, Foam::sumOp<scalar>());
            Foam::reduce(smapsLockedFractionMin, Foam::minOp<scalar>());
            Foam::reduce(smapsLockedFractionMax, Foam::maxOp<scalar>());
            Foam::reduce(onFaultCount, Foam::sumOp<scalar>());
            Foam::reduce(eagerCount, Foam::sumOp<scalar>());
            Foam::reduce(unlockedCount, Foam::sumOp<scalar>());
            Foam::reduce(unknownCount, Foam::sumOp<scalar>());

            if (Pstream::master())
            {
                Info<< "ndtbl residency [" << label << "] group " << groupI
                    << " file \"" << group.filePath << "\":" << endl;

                if (availableCount != rankCount)
                {
                    Info<< "  cache/mincore: available on " << availableCount
                        << " / " << rankCount
                        << " ranks; enable ndtbl mmap diagnostics and "
                        << "mmap-backed reads" << endl;
                }
                else
                {
                    Info<< "  cache/mincore: resident MiB min/avg/max = "
                        << residentBytesMin*bytesToMiB << " / "
                        << (residentBytes/rankCount)*bytesToMiB << " / "
                        << residentBytesMax*bytesToMiB
                        << "; fraction min/avg/max = "
                        << residentFractionMin << " / "
                        << residentFractionSum/rankCount << " / "
                        << residentFractionMax
                        << "; mapped MiB max = "
                        << totalBytesMax*bytesToMiB << endl;
                }

                if (smapsAvailableCount != rankCount)
                {
                    Info<< "  mapping/smaps: available on "
                        << smapsAvailableCount << " / " << rankCount
                        << " ranks" << endl;
                }
                else
                {
                    Info<< "  mapping/smaps: RSS MiB min/avg/max = "
                        << smapsRssBytesMin*bytesToMiB << " / "
                        << (smapsRssBytesSum/rankCount)*bytesToMiB << " / "
                        << smapsRssBytesMax*bytesToMiB;
                    if (smapsLockedAvailableCount == rankCount)
                    {
                        Info<< "; Locked MiB min/avg/max = "
                            << smapsLockedBytesMin*bytesToMiB << " / "
                            << (smapsLockedBytesSum/rankCount)*bytesToMiB
                            << " / " << smapsLockedBytesMax*bytesToMiB
                            << "; locked fraction min/avg/max = "
                            << smapsLockedFractionMin << " / "
                            << smapsLockedFractionSum/rankCount << " / "
                            << smapsLockedFractionMax;
                    }
                    else
                    {
                        Info<< "; Locked available on "
                            << smapsLockedAvailableCount << " / "
                            << rankCount << " ranks";
                    }
                    Info<< "; mapping MiB max = "
                        << smapsMappingBytesMax*bytesToMiB << endl;
                }

                Info<< "  lock policy: on-fault = " << onFaultCount << " / "
                    << rankCount << " ranks; eager = " << eagerCount << " / "
                    << rankCount << "; unlocked = " << unlockedCount << " / "
                    << rankCount << "; unknown = " << unknownCount << " / "
                    << rankCount << endl;
            }
        }

        writeResidencyTsv(runTime, tsvRows.str());
        ++residencyReportSample_;
    }
};

} // End combustionModels namespace
} // End Foam namespace

#endif /* TABLESOLVER_H_ */
