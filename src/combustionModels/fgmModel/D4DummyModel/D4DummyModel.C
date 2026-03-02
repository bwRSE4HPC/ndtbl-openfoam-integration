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

#include "D4DummyModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(D4DummyModel, 0);
    addToRunTimeSelectionTable(combustionModel, D4DummyModel, dictionary);
}
}

// * * * * * * * * * * * * *  Private member functions  * * * * * * * * * * * //

void Foam::combustionModels::D4DummyModel::tableLookup()
{
    // Look-up parameter fields
    const scalarField& Param1Cells = Param1_.internalField();
    const scalarField& Param2Cells = Param2_.internalField();
    const scalarField& Param3Cells = Param3_.internalField();
    const scalarField& Param4Cells = Param4_.internalField();

    // Tabulated parameter fields
    scalarField& Table1Cells = Table1_.primitiveFieldRef();
    scalarField& Table2Cells = Table2_.primitiveFieldRef();
    scalarField& Table3Cells = Table3_.primitiveFieldRef();
    scalarField& Table4Cells = Table4_.primitiveFieldRef();

    // Local (in cell) look-up parameters
    scalarList x(4); // Local parameter value
    List<int> ub(4); // Upper bound index for linear interpolation
    scalarList pos(4); // Interpolation weight

    // For internal cells
    forAll(Param1Cells, cellI)
    {
        // Store local values of parameters
        x[0] = Param1Cells[cellI];
        x[1] = Param2Cells[cellI];
        x[2] = Param3Cells[cellI];
        x[3] = Param4Cells[cellI];

        // Determine upper bound indices and weights for interpolation
        ub = solver_.upperBounds(x, solver_.sizeTableNames() - 1);
        pos = solver_.position(ub, x, solver_.sizeTableNames() - 1);

        // Perform interpolation
        Table1Cells[cellI] = solver_.interpolate(ub, pos, (solver_.sizeTableNames() - 4));
        Table2Cells[cellI] = solver_.interpolate(ub, pos, (solver_.sizeTableNames() - 3));
        Table3Cells[cellI] = solver_.interpolate(ub, pos, (solver_.sizeTableNames() - 2));
        Table4Cells[cellI] = solver_.interpolate(ub, pos, (solver_.sizeTableNames() - 1));
    }

    // For bundary faces
    forAll(Param1_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pParam1 = Param1_.boundaryField()[patchi];
        const fvPatchScalarField& pParam2 = Param2_.boundaryField()[patchi];
        const fvPatchScalarField& pParam3 = Param3_.boundaryField()[patchi];
        const fvPatchScalarField& pParam4 = Param4_.boundaryField()[patchi];

        fvPatchScalarField& pTable1 = Table1_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pTable2 = Table2_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pTable3 = Table3_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pTable4 = Table4_.boundaryFieldRef()[patchi];

        forAll(pParam1 , facei)
        {
            // Store local values of parameters
            x[0] = pParam1[facei];
            x[1] = pParam2[facei];
            x[2] = pParam3[facei];
            x[3] = pParam4[facei];

            // Determine upper bound indices and weights for interpolation
            ub = solver_.upperBounds(x, solver_.sizeTableNames() - 1);
            pos = solver_.position(ub, x, solver_.sizeTableNames() - 1);

            // Perform interpolation
            pTable1[facei] = solver_.interpolate(ub, pos, (solver_.sizeTableNames() - 4));
            pTable2[facei] = solver_.interpolate(ub, pos, (solver_.sizeTableNames() - 3));
            pTable3[facei] = solver_.interpolate(ub, pos, (solver_.sizeTableNames() - 2));
            pTable4[facei] = solver_.interpolate(ub, pos, (solver_.sizeTableNames() - 1));
        }
    }
}

Foam::hashedWordList Foam::combustionModels::D4DummyModel::parameters()
{
    hashedWordList paramNames;
    paramNames.append("Param1");
    paramNames.append("Param2");
    paramNames.append("Param3");
    paramNames.append("Param4");

    return paramNames;
}

Foam::hashedWordList Foam::combustionModels::D4DummyModel::tables()
{
    hashedWordList tableNames;

    tableNames.append("Table1");
    tableNames.append("Table2");
    tableNames.append("Table3");
    tableNames.append("Table4");

    return tableNames;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::D4DummyModel::D4DummyModel
(
    const word& modelType,
    const fluidReactionThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const fluidReactionThermophysicalTransportModel& trans,
    const word& combustionProperties
)
:
    fgmModel(modelType, thermo, turb, trans, combustionProperties),
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
    Param4_
    (
        IOobject
        (
            "Param4",
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
    Table4_
    (
        IOobject
        (
            "Table4",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    solver_(tableSolver(this->mesh(), tables(), parameters()))
{}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::combustionModels::D4DummyModel::~D4DummyModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::combustionModels::D4DummyModel::correct()
{
    // Perform table lookup
    tableLookup();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
