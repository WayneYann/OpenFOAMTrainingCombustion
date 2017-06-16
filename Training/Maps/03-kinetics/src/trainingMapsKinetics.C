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

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

int main(int argc, char *argv[])
{
	// OpenSMOKE++ Maps	
	OpenSMOKE::ThermodynamicsMap_CHEMKIN* 	thermoMap;
	OpenSMOKE::KineticsMap_CHEMKIN*		kineticsMap;
	

	// Import map from preprocessed XML file
	{
		boost::filesystem::path file_path = "../../../../PreProcessing/POLIMI_H2_1412/kinetics-POLIMI_H2_1412/kinetics.xml";

		// Open XML file containing the thermodynamic data
		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc, xml_string, file_path);

		// Import
		thermoMap = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
		kineticsMap = new OpenSMOKE::KineticsMap_CHEMKIN(*thermoMap, doc);  
	}

	
	// Mixture (molar basis): H2(18%), H(1%), OH(1%), O2(10%), N2(70%)  
	Eigen::VectorXd x(thermoMap->NumberOfSpecies());
	x.setZero();
	x(thermoMap->IndexOfSpecies("H2")-1) = 0.18;
	x(thermoMap->IndexOfSpecies("H")-1)  = 0.01;
        x(thermoMap->IndexOfSpecies("OH")-1) = 0.01;
	x(thermoMap->IndexOfSpecies("O2")-1) = 0.10;
	x(thermoMap->IndexOfSpecies("N2")-1) = 0.70;

	// Operating conditions
	const double T = 1000.;
	const double P = 101325.;

	// Set maps
	thermoMap->SetTemperature(T);
	thermoMap->SetPressure(P);
	kineticsMap->SetTemperature(T);
	kineticsMap->SetPressure(P);

	// Concentrations (in kmol/m3)
	const double cTot = P/(PhysicalConstants::R_J_kmol*T);
	Eigen::VectorXd c(thermoMap->NumberOfSpecies());
	c = cTot*x;
	
	// Basic functions
	{
		// Reaction rates (in kmol/m3/s)
		Eigen::VectorXd r(kineticsMap->NumberOfReactions());
		kineticsMap->ReactionRates(c.data());
		kineticsMap->GiveMeReactionRates(r.data());

		// Forward and backward reaction rates (the sum is equal to the reaction rate) (in kmol/m3/s)
		Eigen::VectorXd rf(kineticsMap->NumberOfReactions());
		kineticsMap->GetForwardReactionRates(rf.data());
		Eigen::VectorXd rb(kineticsMap->NumberOfReactions());
		kineticsMap->GetBackwardReactionRates(rb.data());

		// Formation rates (in kmol/m3/s)
		Eigen::VectorXd R(thermoMap->NumberOfSpecies());
		kineticsMap->FormationRates(R.data());

		// Production and destruction rates (the sum is equal to the formation rate) (in kmol/m3/s)
		Eigen::VectorXd RP(thermoMap->NumberOfSpecies());
		Eigen::VectorXd RD(thermoMap->NumberOfSpecies());
		kineticsMap->ProductionAndDestructionRates(RP.data(), RD.data());

		// Derivative of formation rates with respect to concentrations (in 1/s)
		Eigen::MatrixXd dR_over_dC(thermoMap->NumberOfSpecies(), thermoMap->NumberOfSpecies());
		kineticsMap->DerivativesOfFormationRates(c.data(), &dR_over_dC);

		// Heat release (in J/m3/s)
		const double Q = kineticsMap->HeatRelease(R.data());


		// Print
		std::cout << std::endl;
		std::cout << "Formation rates (in kmol/m3/s)" << std::endl;
		for (unsigned int i=0;i<thermoMap->NumberOfSpecies();i++)
			std::cout 	<< std::setw(16) << std::left << thermoMap->NamesOfSpecies()[i]
					<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << R(i) 
					<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << RP(i)
					<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << RD(i)
					<< std::endl;

		// Print
		std::cout << std::endl;
		std::cout << "Reaction rates (in kmol/m3/s)" << std::endl;
		for (unsigned int i=0;i<kineticsMap->NumberOfReactions();i++)
			std::cout 	<< std::setw(5)  << std::left << i
					<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << r(i) 
					<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << rf(i)
					<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << rb(i)
					<< std::endl;

		// Print
		std::cout << std::endl;
		std::cout << "Derivative of formation rates with respect to concentrations (in 1/s)" << std::endl;
		for (unsigned int i=0;i<thermoMap->NumberOfSpecies();i++)
		{
			for (unsigned int j=0;j<thermoMap->NumberOfSpecies();j++)
				std::cout << std::setw(11) << std::right << std::scientific << std::setprecision(3) << dR_over_dC(i,j);
			std::cout << std::endl;
		} 
		
	}
}

