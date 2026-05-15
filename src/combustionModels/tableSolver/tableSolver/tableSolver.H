#ifndef TABLESOLVER_H_
#define TABLESOLVER_H_

#include "dictionary.H"
#include "error.H"
#include "fileName.H"
#include "wordList.H"

#include "ndtbl/ndtbl.hpp"

#include <array>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

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

} // End namespace detail

/// Loads ndtbl table files and evaluates configured output fields.
///
/// The solver is parameterized by the number of lookup coordinates and output
/// values. Its constructor reads the model coefficients, loads each unique
/// ndtbl path once per solver instance, and binds the requested output table
/// names to fields in the loaded table groups. A lookup evaluates every loaded
/// group at a query point and copies the bound field values into the
/// caller-provided output array.
template<std::size_t Dim, std::size_t NOut>
class tableSolver
{
private:

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
        std::array<bool, NOut> outputBound = {};

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
    void lookup(const std::array<scalar, Dim>& x, std::array<scalar, NOut>& values) const
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
};

} // End combustionModels namespace
} // End Foam namespace

#endif /* TABLESOLVER_H_ */
