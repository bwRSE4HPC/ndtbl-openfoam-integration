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

Class
    Foam::combustionModels::fgmModel

Description
    Base class for fgmModel.

SourceFiles
    fgmModel.C

\*---------------------------------------------------------------------------*/

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
