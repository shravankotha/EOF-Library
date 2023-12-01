/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    mhdInterFoam = Elmer (EM) + interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

-------------------------------------------------------------------------------

Original interFoam solver is part of OpenFOAM

Modified by: Juris Vencels

\*---------------------------------------------------------------------------*/
// shravan - The original interFoam solver details can be found at: https://openfoamwiki.net/index.php/InterFoam

#include "fvCFD.H"	// shravan - This file brings in the most fundamental tools for performing any finite volume calculation. This file in fact just includes a bunch of other files, each of which represents a building block of the edifice of the finite volume technique. see https://openfoamwiki.net/index.php/OpenFOAM_guide/Introduction_to_OpenFOAM_Programming,_A_Walk_Through_reactingFOAM
#include "dynamicFvMesh.H"	// shravan - Abstract base class for geometry and/or topology changing fvMesh. see https://cfd.direct/openfoam/free-software/dynamic-meshes/. The settings to control the mesh are taken from dynamicMeshDict file in levitation/constant/ folder for ex.
#include "CMULES.H"	// shravan - Solve a convective-only transport equation using an explicit universal multi-dimensional limiter to correct an implicit conservative bounded obtained using rigorously bounded schemes such as Euler-implicit in time upwind in space
#include "EulerDdtScheme.H"	// shravan - Basic first-order Euler implicit/explicit ddt using only the current and previous time-step values. https://www.openfoam.com/documentation/guides/latest/doc/guide-schemes-time.html
#include "localEulerDdtScheme.H"	// shravan -     Local time-step first-order Euler implicit/explicit ddt. The reciprocal of the local time-step field is looked-up from the database. This scheme should only be used for steady-state computations using transient codes where local time-stepping is preferably to under-relaxation for transport consistency reasons. 
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"	// shravan - An immiscible incompressible two-phase mixture transport model
#include "turbulentTransportModel.H"		// shravan -     Typedefs for turbulence, RAS and LES models for incompressible flow based on the standard laminar transport package.
#include "pimpleControl.H"	// shravan -     PIMPLE control class to supply convergence information/checks for the PIMPLE loop. May also be used to for PISO-based algorithms as PISO controls are a sub-set of PIMPLE controls.
#include "fvOptions.H"
#include "CorrectPhi.H"	// shravan - Flux correction functions to ensure continuity.
#include "fvcSmooth.H"	// shravan - Provides functions smooth spread and sweep which use the FaceCellWave algorithm to smooth and redistribute the first field argument.
#include "Elmer.H"	// shravan - Coupling with FEM open source multiphysical simulation software Elmer. This is the only header file that is not part of interFoam
					// shravan - Elmer.H and Elmer.C files have all the methods used to exchange information between Elmer and openfoam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])	// shravan - argc - number of arguments enetered from the command line and argv[] is the streing containing the arguments
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"	// shravan - DyM indicates that the mesh is dynamic.
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"	// shravan - Creates and initialises the velocity field Uf if required.

    turbulence->validate();		// shravan - same as (*turbulence).validate(). i.e. dereference the object pointed to by turbulence and get its method called validate(). See https://github.com/OpenFOAM/OpenFOAM-6/blob/master/src/TurbulenceModels/turbulenceModels/turbulenceModel.C

    if (!LTS)	// shravan - LTS==Local time stepping
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    // Send fields to Elmer
    Elmer<dynamicFvMesh> sending(mesh,1); // 1=send, -1=receive	// shravan - sending is an object that is an instance of template class Foam::Elmer. since, Foam::Elmer has a constructor method called Elmer, it is executed when sending(mesh,1) is instantiated here. This class has a method called sendStatus as defined in Elmer.C file
    sending.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    elcond = alpha1 * elcond_melt;	// shravan - elcond = electrical conducitivty of the melt
    sending.sendScalar(elcond);	// shravan - electrical conductivity sent from openFoam and received by Elmer

    // Receive fields from Elmer
    Elmer<dynamicFvMesh> receiving(mesh,-1); // 1=send, -1=receive
    receiving.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    receiving.recvVector(JxB_recv);	// shravan - magnetic force is received by openFoam

    while (runTime.run())
    {
        JxB = JxB_recv*alpha1;	// shravan - alpha1 is the volume fraction of melt in any given cell. So multiplying JxB_recv with alpha1 gives the force acting on the melt in that cell.

        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Check whether we need to update electromagnetic stuff with Elmer
        scalar maxRelDiff_local = (max(mag(alpha_old - alpha1))).value();	// shravan - alpha is the fluid fraction parameter in VoF method. See equations/ur notes. Look at the max of abs. values of diff between old and new alpha values over the entire domain and decide

        bool doElmer = false;
        if(maxRelDiff_local>maxRelDiff && (maxRelDiff<SMALL || maxRelDiff+SMALL<=1.0)) {	// shravan - maxRelDiff is given in const/transportProperties file
            doElmer = true;
        }

        if(doElmer || !runTime.run()) {
            alpha_old = alpha1;

            // Send fields to Elmer
            sending.sendStatus(runTime.run());	// shravan - here Elmer does EM simulation
            elcond = alpha1 * elcond_melt;	// shravan - sending the timestep status, electrical conductivity
            sending.sendScalar(elcond);

            // Receive fields form Elmer
            receiving.sendStatus(runTime.run());	// shravan - receive the data from Elmer i.e. the magnetic body force JxB_recv
            receiving.recvVector(JxB_recv);
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
