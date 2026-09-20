// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef fgmModel_H
#define fgmModel_H

#include "combustionModel.H"
#include "basicSpecieMixture.H"

#include "hashedWordList.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "tableSolver.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class fgmModel Declaration
\*---------------------------------------------------------------------------*/

class fgmModel
:
    public combustionModel
{

private:

    // Private Member

        //- Disallow copy construct
        fgmModel(const fgmModel&);

        //- Disallow default bitwise assignment
        void operator=(const fgmModel&);

protected:

public:

    // Constructors

        //- Construct from components
        fgmModel
        (
            const word& modelType,
            const fluidReactionThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const fluidReactionThermophysicalTransportModel& trans,
            const word& combustionProperties
        );

        //- Destructor
        virtual ~fgmModel();

    // Member Functions

        //- Correct combustion rate
        virtual void correct() = 0;

    // I-O
        //- Update properties
        virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
