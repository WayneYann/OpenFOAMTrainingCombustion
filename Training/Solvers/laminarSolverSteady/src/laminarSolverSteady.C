/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
\*-----------------------------------------------------------------------*/

// OpenSMOKE++
#include "OpenSMOKEpp"
#include "maps/Maps_CHEMKIN"
#include "linearModelChemistry.H"

// OpenFOAM
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "simpleControl.H"
#include "interpolation.H"
#include "fvOptions.H"
#include "fvcSmooth.H"
#include "pressureControl.H"
#include "sparkModel.H"

int main(int argc, char *argv[])
{
	// OpenFOAM stuff
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	#include "createControl.H"
	#include "readGravitationalAcceleration.H"

	// OpenSMOKE++
	#include "createBasicFields.H"
	#include "createChemicalFields.H"
	#include "createSourceFields.H"
	#include "transportProperties.H"
	#include "createAdditionalFields.H"

	// Linear model for reacting source term
	linearModelChemistry chemistry(thermoMap, kineticsMap);

	// OpenFOAM stuff
        #include "createMRF.H"
	#include "createFvOptions.H"
	#include "initContinuityErrs.H"

	dimensionedScalar initialMass = fvc::domainIntegrate(rho);

	Info<< "\nStarting time loop\n" << endl;
	while (simple.loop())
        {
		Info<< "Time = " << runTime.timeName() << nl << endl;

		double t0 = runTime.value() - runTime.deltaT().value();
		double tf = runTime.value();

		// Pressure-velocity SIMPLE corrector
		{
			#include "UEqn.H"

			#include "transportProperties.H"
			#include "chemistry.H"

			#include "YEqn.H"
			#include "TEqn.H" 

			if (simple.consistent())
				#include "pcEqn.H"
			else
				#include "pEqn.H"
		}

		runTime.write();
		
         	Info 	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
              		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
              		<< nl << endl;
    	}

    	Info<< "End\n" << endl;

    	return 0;
}


// ************************************************************************* //
