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
	argList::validArgs.append("ODESolver");

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

//	#include "createAdditionalFields.H"

//	argList::validArgs.append("ODESolver");
//	argList args(argc, argv);

	// 
/*

	// Operating conditions
	const double T = 1000.;
	const double P = 101325.;
	Eigen::VectorXd xStart(thermoMap.NumberOfSpecies());
	xStart.setZero();
	xStart(thermoMap.IndexOfSpecies("H2")-1)  = 0.11;
	xStart(thermoMap.IndexOfSpecies("CO")-1)  = 0.11;
	xStart(thermoMap.IndexOfSpecies("O2")-1)  = 0.17;
	xStart(thermoMap.IndexOfSpecies("N2")-1)  = 0.61;

	// Create the ODE system as object of type batchOdeSystem
	batchOdeSystem batch(thermoMap, kineticsMap, T,P);

	// Create dictionary and add the odeSolver name
	dictionary dict;
	dict.add("solver", args[1]);

	// Create the selected ODE system solver
	autoPtr<ODESolver> odeSolver = ODESolver::New(batch, dict);


	// Initialize the ODE system fields (concentrations in kmol/m3)
	const double cTot = P/(PhysicalConstants::R_J_kmol*T);
	scalarField cStart(thermoMap.NumberOfSpecies());
	for(unsigned int i=0;i<thermoMap.NumberOfSpecies();i++)
		cStart[i] = cTot*xStart[i];

	// ODE integration parameters
	const label n = 1000;		// number of steps (used only for writing output)
	scalar tStart = 0.;		// start time (in s)
    	scalar tEnd = 1e-3; 		// end time (in s)
	scalar dtStart = 1e-8;		// initial time step (in s)
	scalar dt = (tEnd-tStart)/n;	// time step (in s)

	// Required to store dcdt
	scalarField dcStart(thermoMap.NumberOfSpecies());

	// Output file
	std::ofstream fOutput("Solution.out", std::ios::out);
	fOutput.setf(std::ios::scientific);
	fOutput.setf(std::ios::left);

	fOutput << std::setw(16) << "time(1)";
 	for(unsigned int i=0;i<thermoMap.NumberOfSpecies();i++)
	{
		std::string label = thermoMap.NamesOfSpecies()[i] + "(" + std::to_string(i+2) + ")";
		fOutput << std::setw(16) << label;
	}
	fOutput << std::endl;

	// Integration loop
	for (label i=0; i<n; i++)
	{
		std::cout << std::scientific << tStart << std::endl;
		
		fOutput << std::setw(16) << tStart;
 		for(unsigned int i=0;i<thermoMap.NumberOfSpecies();i++)
			fOutput << std::setw(16) << cStart[i];
		fOutput << std::endl;

		batch.derivatives(tStart, cStart, dcStart);
		odeSolver->solve(tStart, tStart + dt, cStart, dtStart);
		tStart += dt;
	}

	fOutput.close();
*/
	return 0;
}

