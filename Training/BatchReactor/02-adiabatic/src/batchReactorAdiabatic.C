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

#include "batchAdiabaticOdeSystem.H"

int main(int argc, char *argv[])
{

	argList::validArgs.append("ODESolver");
	argList args(argc, argv);

	boost::filesystem::path file_path = "../../../../PreProcessing/POLIMI_H2CO_1412/kinetics-POLIMI_H2CO_1412/kinetics.xml";

	// Open XML file containing the thermodynamic data
	rapidxml::xml_document<> doc;
	std::vector<char> xml_string;
	OpenSMOKE::OpenInputFileXML(doc, xml_string, file_path);

	OpenSMOKE::ThermodynamicsMap_CHEMKIN thermoMap(doc);
	OpenSMOKE::KineticsMap_CHEMKIN kineticsMap(thermoMap, doc);  

	// Operating conditions
	double T = 1000.;
	double P = 101325.;
	Eigen::VectorXd moleFractions(thermoMap.NumberOfSpecies());
	moleFractions.setZero();
	moleFractions(thermoMap.IndexOfSpecies("H2")-1)  = 0.11;
	moleFractions(thermoMap.IndexOfSpecies("CO")-1)  = 0.11;
	moleFractions(thermoMap.IndexOfSpecies("O2")-1)  = 0.17;
	moleFractions(thermoMap.IndexOfSpecies("N2")-1)  = 0.61;

	// Create the ODE system as object of type batchOdeSystem
	batchAdiabaticOdeSystem batch(thermoMap, kineticsMap);

	// Create dictionary and add the odeSolver name
	dictionary dict;
	dict.add("solver", args[1]);

	// Create the selected ODE system solver
	autoPtr<ODESolver> odeSolver = ODESolver::New(batch, dict);


	// Initialize the ODE system fields (concentrations in kmol/m3)
	double cTot = P/(PhysicalConstants::R_J_kmol*T);
	scalarField c(thermoMap.NumberOfSpecies());
	for(unsigned int i=0;i<thermoMap.NumberOfSpecies();i++)
		c[i] = cTot*moleFractions(i);

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
	fOutput << std::setw(16) << "T(2)";
	fOutput << std::setw(16) << "P(3)";
 	for(unsigned int i=0;i<thermoMap.NumberOfSpecies();i++)
	{
		std::string label = thermoMap.NamesOfSpecies()[i] + "(" + std::to_string(i+4) + ")";
		fOutput << std::setw(16) << label;
	}
	fOutput << std::endl;

	// Integration loop
	for (label i=0; i<n; i++)
	{
		// Info
		std::cout << std::setw(16) << std::scientific << tStart;
		std::cout << std::setw(16) << std::scientific << T;
		std::cout << std::setw(16) << std::scientific << P;
		std::cout << std::endl;
		
		// Print on file
		fOutput << std::setw(16) << tStart;
		fOutput << std::setw(16) << T;
		fOutput << std::setw(16) << P;
 		for(unsigned int i=0;i<thermoMap.NumberOfSpecies();i++)
			fOutput << std::setw(16) << c[i];
		fOutput << std::endl;

		// From concentrations to mass fractions
		cTot = std::accumulate(c.begin(), c.end(), 0.);
		for(unsigned int i=0;i<thermoMap.NumberOfSpecies();i++)
			moleFractions(i) = c[i]/cTot;

		// Enthalpy
		thermoMap.SetTemperature(T);
		thermoMap.SetPressure(P);
		const double H = thermoMap.hMolar_Mixture_From_MoleFractions(moleFractions.data());
		batch.setEnthalpy(H);

		// Solve ODE system
		batch.setInitialTemperature(T);
		batch.setInitialPressure(P);
		batch.derivatives(tStart, c, dcStart);
		odeSolver->solve(tStart, tStart + dt, c, dtStart);
		tStart += dt;

		// Temperature
		cTot = std::accumulate(c.begin(), c.end(), 0.);
		for(unsigned int i=0;i<thermoMap.NumberOfSpecies();i++)
			moleFractions(i) = c[i]/cTot;
		T = thermoMap.GetTemperatureFromEnthalpyAndMoleFractions(H, P, moleFractions.data(), T);
		P = cTot*(PhysicalConstants::R_J_kmol*T);	
	}

	fOutput.close();

	return 0;
}

