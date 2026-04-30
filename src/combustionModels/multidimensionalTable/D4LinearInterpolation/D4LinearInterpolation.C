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


#include "../D4LinearInterpolation/D4LinearInterpolation.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(D4LinearInterpolation, 0);
addToRunTimeSelectionTable(multidimensionalTable, D4LinearInterpolation, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

D4LinearInterpolation::D4LinearInterpolation(const fvMesh& mesh, const word& tableName)
:
    multidimensionalTable(mesh, tableName),
    tableValues_()
{
    IOdictionary currentTable
    (
       IOobject
       (
          tableName,
          mesh.time().constant(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
       )
    );
    tableValues_ = currentTable.lookup<List<List<List<scalarList> > > >(tableName);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

D4LinearInterpolation::~D4LinearInterpolation()
{}

// * * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //

inline scalar D4LinearInterpolation::interpolate(const List<int>& ub, const scalarList& pos) const
{
    scalar c000 = tableValues_[ub[0]-1][ub[1]-1][ub[2]-1][ub[3]-1]*(1-pos[0]) + tableValues_[ub[0]][ub[1]-1][ub[2]-1][ub[3]-1]*pos[0];
    scalar c100 = tableValues_[ub[0] -1][ub[1]][ub[2] -1][ub[3] -1]*(1-pos[0]) + tableValues_[ub[0]][ub[1]][ub[2] -1][ub[3]-1]*pos[0];
    scalar c010 = tableValues_[ub[0] -1][ub[1] -1][ub[2]][ub[3] -1]*(1-pos[0]) + tableValues_[ub[0]][ub[1] -1][ub[2]][ub[3]-1]*pos[0];
    scalar c001 = tableValues_[ub[0] -1][ub[1] -1][ub[2] -1][ub[3]]*(1-pos[0]) + tableValues_[ub[0]][ub[1] -1][ub[2] -1][ub[3]]*pos[0];
    scalar c110 = tableValues_[ub[0] -1][ub[1]][ub[2]][ub[3] -1]*(1-pos[0]) + tableValues_[ub[0]][ub[1]][ub[2]][ub[3]-1]*pos[0];
    scalar c101 = tableValues_[ub[0] -1][ub[1]][ub[2] -1][ub[3]]*(1-pos[0]) + tableValues_[ub[0]][ub[1]][ub[2] -1][ub[3]]*pos[0];
    scalar c011 = tableValues_[ub[0] -1][ub[1] -1][ub[2]][ub[3]]*(1-pos[0]) + tableValues_[ub[0]][ub[1] -1][ub[2]][ub[3]]*pos[0];
    scalar c111 = tableValues_[ub[0] -1][ub[1]][ub[2]][ub[3]]*(1-pos[0]) + tableValues_[ub[0]][ub[1]][ub[2]][ub[3]]*pos[0];

    scalar c00 = c000*(1-pos[1]) + c100*pos[1];
    scalar c10 = c010*(1-pos[1]) + c110*pos[1];
    scalar c01 = c001*(1-pos[1]) + c101*pos[1];
    scalar c11 = c011*(1-pos[1]) + c111*pos[1];

    scalar c0 = c00*(1-pos[2]) + c10*pos[2];
    scalar c1 = c01*(1-pos[2]) + c11*pos[2];

    return c0*(1-pos[3]) + c1*pos[3];
}

inline List<int> D4LinearInterpolation::tableSize() const
{
    List<int> tableSize(4);

    tableSize[0] = tableValues_.size();
    tableSize[1] = tableValues_[0].size();
    tableSize[2] = tableValues_[0][0].size();
    tableSize[3] = tableValues_[0][0][0].size();

    return tableSize;
}

} // End Foam namespace
