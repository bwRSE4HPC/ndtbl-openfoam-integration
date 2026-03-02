/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "tableSolver.H"
#include "Pstream.H"

namespace Foam
{
namespace combustionModels
{

defineTypeNameAndDebug(tableSolver, 0);
defineRunTimeSelectionTable(tableSolver, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tableSolver::tableSolver(const fvMesh& mesh, const wordList& tableNames, const wordList& paramNames)
:
    tableNames_(tableNames),
    paramNames_(paramNames),
    paramSize_(tableNames_.size()),
    tables_(tableNames_.size())
{
    forAll(tableNames_, i)
    {
        tableNames_[i] = tableNames_[i] + "_table";
        tables_.set(i, multidimensionalTable::New(mesh, tableNames_[i], paramNames_));

        paramSize_[i] = tables_[i].tableSize();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tableSolver::~tableSolver()
{}

// * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * //

List<int> tableSolver::upperBounds(const scalarList& x, const label& i) const
{
    List<int> ub(x.size(), 0);
    int tlb;

    forAll(ub, j)
    {
        tlb = min(x[j] * (paramSize_[i][j] -1), (paramSize_[i][j] -2));
        ub[j] = tlb + 1;
    }

    return ub;
}

scalarList tableSolver::position(const List<int>& ub, const scalarList& x, const label& i) const
{
    scalarList pos(x.size(), 0.);

    forAll(ub, j)
    {
        pos[j] =  x[j] * (paramSize_[i][j] - 1) - (ub[j] - 1);
    }

    return pos;
}

scalar tableSolver::interpolate(const List<int>& ub , const scalarList& pos, const label& i) const
{
    return tables_[i].interpolate(ub, pos);
}

int tableSolver::sizeTableNames() const
{
    return tableNames_.size();
}

} // End combustionModels namespace
} // End Foam namespace
