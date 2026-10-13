// SPDX-License-Identifier: GPL-3.0-or-later

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
