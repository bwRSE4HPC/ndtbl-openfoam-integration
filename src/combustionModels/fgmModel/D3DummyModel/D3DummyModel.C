/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "D3DummyModel.H"
#include "addToRunTimeSelectionTable.H"

#include <array>

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(D3DummyModel, 0);
    addToRunTimeSelectionTable(combustionModel, D3DummyModel, dictionary);
}
}

// * * * * * * * * * * * * *  Private member functions  * * * * * * * * * * * //

void Foam::combustionModels::D3DummyModel::tableLookup()
{
    // Look-up parameter fields
    const scalarField& Param1Cells = Param1_.internalField();
    const scalarField& Param2Cells = Param2_.internalField();
    const scalarField& Param3Cells = Param3_.internalField();

    // Tabulated parameter fields
    scalarField& Table1Cells = Table1_.primitiveFieldRef();
    scalarField& Table2Cells = Table2_.primitiveFieldRef();
    scalarField& Table3Cells = Table3_.primitiveFieldRef();

    // Local (in cell) look-up parameters
    std::array<scalar, 3> x = { scalar(0), scalar(0), scalar(0) };
    std::vector<scalar> values = { scalar(0), scalar(0), scalar(0) };

    // For internal cells
    forAll(Param1Cells, cellI)
    {
        // Store local values of parameters
        x[0] = Param1Cells[cellI];
        x[1] = Param2Cells[cellI];
        x[2] = Param3Cells[cellI];

        solver_.lookup(x, values);
        Table1Cells[cellI] = values[0];
        Table2Cells[cellI] = values[1];
        Table3Cells[cellI] = values[2];
    }

    // For bundary faces
    forAll(Param1_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pParam1 = Param1_.boundaryField()[patchi];
        const fvPatchScalarField& pParam2 = Param2_.boundaryField()[patchi];
        const fvPatchScalarField& pParam3 = Param3_.boundaryField()[patchi];

        fvPatchScalarField& pTable1 = Table1_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pTable2 = Table2_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pTable3 = Table3_.boundaryFieldRef()[patchi];

        forAll(pParam1 , facei)
        {
            // Store local values of parameters
            x[0] = pParam1[facei];
            x[1] = pParam2[facei];
            x[2] = pParam3[facei];

            solver_.lookup(x, values);
            pTable1[facei] = values[0];
            pTable2[facei] = values[1];
            pTable3[facei] = values[2];
        }
    }
}

Foam::hashedWordList Foam::combustionModels::D3DummyModel::parameters()
{
    hashedWordList paramNames;
    paramNames.append("Param1");
    paramNames.append("Param2");
    paramNames.append("Param3");

    return paramNames;
}

Foam::hashedWordList Foam::combustionModels::D3DummyModel::tables()
{
    hashedWordList tableNames;

    tableNames.append("Table1");
    tableNames.append("Table2");
    tableNames.append("Table3");

    return tableNames;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::D3DummyModel::D3DummyModel
(
    const word& modelType,
    const fluidReactionThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const fluidReactionThermophysicalTransportModel& trans,
    const word& combustionProperties
)
:
    fgmModel(modelType, thermo, turb, trans, combustionProperties),
    gammaParam_("gammaParam", dimensionSet(0, 2, -1, 0, 0), this->coeffs().lookup<scalar>("gammaParam")),
    rhoParam_("rhoParam", dimensionSet(1, -3, 0, 0, 0), this->coeffs().lookup<scalar>("rhoParam")),
    Param1_
    (
        IOobject
        (
            "Param1",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    Param2_
    (
        IOobject
        (
            "Param2",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    Param3_
    (
        IOobject
        (
            "Param3",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    Table1_
    (
        IOobject
        (
            "Table1",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    Table2_
    (
        IOobject
        (
            "Table2",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    Table3_
    (
        IOobject
        (
            "Table3",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    tableInitializationStartCpuTime_(this->mesh().time().elapsedCpuTime()),
    ndtblDiagnostics_(false),
    solver_(tableSolver<3>(this->coeffs(), tables(), parameters()))
{
    if (this->coeffs().found("ndtblDiagnostics"))
    {
        this->coeffs().lookup("ndtblDiagnostics") >> ndtblDiagnostics_;
    }

    Info<< "Total table loading/initialization time = "
        << this->mesh().time().elapsedCpuTime() - tableInitializationStartCpuTime_
        << " s" << endl;
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::combustionModels::D3DummyModel::~D3DummyModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::combustionModels::D3DummyModel::correct()
{
    // Access to phi for transport equations
    const surfaceScalarField& phi = mesh_.lookupObject<surfaceScalarField>("phi");

    // Transport equations for tabulation parameters
    {
        volScalarField& Param_ = Param1_;
        #include "../ParamEqn.H"
    }
    {
        volScalarField& Param_ = Param2_;
        #include "../ParamEqn.H"
    }
    {
        volScalarField& Param_ = Param3_;
        #include "../ParamEqn.H"
    }

    // Perform table lookup
    tableLookup();
}

void Foam::combustionModels::D3DummyModel::reportTableResidency
(
    const word& label
) const
{
    if (!ndtblDiagnostics_)
    {
        return;
    }

    solver_.reportResidency(label);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
