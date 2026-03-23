/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::combustionModels::tableSolver

Description
    Interface between ndtbl table look-up and the combustion model.

\*---------------------------------------------------------------------------*/

#ifndef TABLESOLVER_H_
#define TABLESOLVER_H_

#include "dictionary.H"
#include "error.H"
#include "fileName.H"
#include "wordList.H"

#include "ndtbl/io.hpp"

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

inline std::string readPathEntry(const dictionary& dict, const word& key)
{
    fileName pathValue;
    dict.lookup(key) >> pathValue;
    return std::string(pathValue.c_str());
}

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

template<std::size_t Dim, std::size_t NOut>
class tableSolver
{
private:

    struct fieldBinding
    {
        std::size_t outputIndex;
        std::size_t fieldIndex;
    };

    struct tableGroup
    {
        std::string filePath;
        ndtbl::RuntimeFieldGroup<Dim> group;
        std::vector<fieldBinding> bindings;
        mutable std::vector<double> scratchValues;
    };

    std::vector<tableGroup> groups_;

public:

    tableSolver(const dictionary& coeffs, const wordList& tableNames, const wordList& paramNames)
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

        forAll(tableNames, outputI)
        {
            const word& tableName = tableNames[outputI];

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
                    newGroup.group = ndtbl::read_group<Dim>(filePath);
                    newGroup.scratchValues.resize(newGroup.group.field_count(), 0.0);
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
                binding.outputIndex = static_cast<std::size_t>(outputI);
                binding.fieldIndex = groups_[groupIndex].group.field_index(fieldName);
                groups_[groupIndex].bindings.push_back(binding);
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
    }

    void interpolate(const std::array<scalar, Dim>& x, std::array<scalar, NOut>& values) const
    {
        values.fill(scalar(0));

        std::array<double, Dim> query;
        for (std::size_t axis = 0; axis < Dim; ++axis)
        {
            query[axis] = x[axis];
        }

        for (std::size_t groupI = 0; groupI < groups_.size(); ++groupI)
        {
            const tableGroup& group = groups_[groupI];
            try
            {
                group.group.evaluate_all_into(query, group.scratchValues.data());

                for (std::size_t bindingI = 0; bindingI < group.bindings.size(); ++bindingI)
                {
                    const fieldBinding& binding = group.bindings[bindingI];
                    values[binding.outputIndex] =
                        static_cast<scalar>(group.scratchValues[binding.fieldIndex]);
                }
            }
            catch (const std::exception& error)
            {
                FatalErrorInFunction
                    << "Failed to evaluate ndtbl file " << group.filePath
                    << ": " << error.what()
                    << exit(FatalError);
            }
        }
    }
};

} // End combustionModels namespace
} // End Foam namespace

#endif /* TABLESOLVER_H_ */
