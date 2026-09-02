// SPDX-License-Identifier: GPL-3.0-or-later

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

\*---------------------------------------------------------------------------*/

#include "fgmModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::fgmModel::fgmModel
(
        const word& modelType,
        const fluidReactionThermo& thermo,
        const compressibleMomentumTransportModel& turb,
        const fluidReactionThermophysicalTransportModel& trans,
        const word& combustionProperties
)
:
    combustionModel(modelType, thermo, turb, trans, combustionProperties)
{}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::combustionModels::fgmModel::~fgmModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::combustionModels::fgmModel::read()
{
    if (combustionModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
