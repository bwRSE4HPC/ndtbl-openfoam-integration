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
    Foam::combustionModels::D4DummyModel

Description
    Base class for D4DummyModel.

SourceFiles
    D4DummyModel.C

\*---------------------------------------------------------------------------*/

#ifndef D4DummyModel_H
#define D4DummyModel_H

#include "fgmModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class D4DummyModel Declaration
\*---------------------------------------------------------------------------*/

class D4DummyModel
:
    public fgmModel
{

private:

    // Private Member

        //- Disallow copy construct
        D4DummyModel(const D4DummyModel&);

        //- Disallow default bitwise assignment
        void operator=(const D4DummyModel&);

        //- Diffusion coefficient
        dimensionedScalar gammaParam_;

        //- Density for parameter transport
        dimensionedScalar rhoParam_;

        //- Look-up parameter
        volScalarField Param1_;
        volScalarField Param2_;
        volScalarField Param3_;
        volScalarField Param4_;

        //- Tabulated fields
        volScalarField Table1_;
        volScalarField Table2_;
        volScalarField Table3_;
        volScalarField Table4_;

        //- CPU time before table loading and initialization
        scalar tableInitializationStartCpuTime_;

        //- Tables
        tableSolver solver_;

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
        TypeName("D4DummyModel");

        //- Construct from components
        D4DummyModel
        (
            const word& modelType,
            const fluidReactionThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const fluidReactionThermophysicalTransportModel& trans,
            const word& combustionProperties
        );

        //- Destructor
        virtual ~D4DummyModel();

    // Member Functions

        //- Correct combustion rate
        virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
