/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Application
    bwRSE4HPCFoam

Description
    Dummy solver to evaluate memory requirement and run-time for table look-up

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidReactionThermo.H"
#include "combustionModel.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidReactionThermophysicalTransportModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createDyMControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createRhoUfIfPresent.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    const scalar timeLoopStartCpuTime = runTime.elapsedCpuTime();

    while (pimple.run(runTime))
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
            reaction->correct();
        }

        runTime.write();
    }

    Info<< "Total time loop runtime = "
        << runTime.elapsedCpuTime() - timeLoopStartCpuTime
        << " s" << endl;

    // Evaluate total memory usage
    std::ifstream status("/proc/self/status");
    std::string line;
    double memMB = 0.0;

    while (std::getline(status, line))
    {
        if (line.substr(0,6) == "VmRSS:")
        {
            std::istringstream iss(line);
            std::string key, unit;
            double value;
            iss >> key >> value >> unit;
            memMB = value / 1024.0;
            break;
        }
    }

    Foam::scalar localMB = static_cast<Foam::scalar>(memMB);
    Foam::scalar totalMB = localMB;
    label commID = 0;
    label returnValue = 0;

    // Reduce across all processes
    Foam::reduce(localMB, Foam::sumOp<Foam::scalar>(), 0, commID, returnValue);
    totalMB = localMB;

    if (Foam::Pstream::master())
    {
        Info << "Total memory across all processes: " << totalMB << " MB" << endl;
    }

    Info<< "Execution time end of simulation = " << runTime.elapsedCpuTime() << " s" << endl;

    return 0;
}


// ************************************************************************* //
