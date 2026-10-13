// SPDX-License-Identifier: GPL-3.0-or-later

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
