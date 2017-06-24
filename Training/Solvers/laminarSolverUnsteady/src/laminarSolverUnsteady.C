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

// Batch reactor ODE system
#include "batchAdiabaticOdeSystem.H"

// OpenFOAM
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
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
	#include "createTimeControls.H"
	#include "initContinuityErrs.H"
	#include "createMRF.H"
	#include "createFvOptions.H"
	#include "readGravitationalAcceleration.H"

	// OpenSMOKE++
	#include "createBasicFields.H"
	#include "createChemicalFields.H"
	#include "transportProperties.H"
	#include "createAdditionalFields.H"

	// OpenFOAM stuff
	#include "compressibleCourantNo.H"
	#include "setInitialDeltaT.H"

	Info<< "\nStarting time loop\n" << endl;
	while (runTime.run())
	{
		#include "readTimeControls.H"
		#include "compressibleCourantNo.H"
		#include "setDeltaT.H"

		runTime++;
		Info<< "Time = " << runTime.timeName() << nl << endl;

		if (pimple.nCorrPIMPLE() <= 1)
		{
		    #include "rhoEqn.H"
		}

		// Pressure-velocity PIMPLE corrector loop
		while (pimple.loop())
		{
			// Transport
			#include "UEqn.H"
			#include "YEqn.H"
			#include "TEqn.H"

			// Chemistry
			#include "chemistry.H"
			
			// Update transport properties
			#include "transportProperties.H"

			// Pressure corrector loop
			while (pimple.correct())
			{
				if (pimple.consistent())
				{
					#include "pcEqn.H"
				}
				else
				{
					#include "pEqn.H"
				}
			}
		}

		runTime.write();

		Info 	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;
	}

	Info<< "End\n" << endl;

	return 0;
}

