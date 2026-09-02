// SPDX-License-Identifier: GPL-3.0-or-later

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

Class
    Foam::combustionModels::D3DummyModel

Description
    Base class for D3DummyModel.

SourceFiles
    D3DummyModel.C

\*---------------------------------------------------------------------------*/

#ifndef D3DummyModel_H
#define D3DummyModel_H

#include "fgmModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class D3DummyModel Declaration
\*---------------------------------------------------------------------------*/

class D3DummyModel
:
    public fgmModel
{

private:

    // Private Member

        //- Disallow copy construct
        D3DummyModel(const D3DummyModel&);

        //- Disallow default bitwise assignment
        void operator=(const D3DummyModel&);

        //- Diffusion coefficient
        dimensionedScalar gammaParam_;

        //- Density for parameter transport
        dimensionedScalar rhoParam_;

        //- Look-up parameter
        volScalarField Param1_;
        volScalarField Param2_;
        volScalarField Param3_;

        //- Tabulated fields
        volScalarField Table1_;
        volScalarField Table2_;
        volScalarField Table3_;

        //- CPU time before table loading and initialization
        scalar tableInitializationStartCpuTime_;

        //- Enable ndtbl residency reporting
        bool ndtblDiagnostics_;

        //- Enable per-rank ndtbl residency reporting
        bool ndtblRankDiagnostics_;

        //- Tables
        tableSolver<3> solver_;

    // Member Functions

        //- Perform table lookup
        virtual void tableLookup();

        //- Return parameter list
        virtual hashedWordList parameters();

        //- Return list of tables
        virtual hashedWordList tables();

public:

    // Constructors

        //- Runtime type information
        TypeName("D3DummyModel");

        //- Construct from components
        D3DummyModel
        (
            const word& modelType,
            const fluidReactionThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const fluidReactionThermophysicalTransportModel& trans,
            const word& combustionProperties
        );

        //- Destructor
        virtual ~D3DummyModel();

    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Report ndtbl residency diagnostics
        virtual void reportTableResidency(const word& label) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
