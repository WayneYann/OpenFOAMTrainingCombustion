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
	OpenSMOKE::ThermodynamicsMap_CHEMKIN* 		thermoMap;
	OpenSMOKE::TransportPropertiesMap_CHEMKIN*	tranMap;
	

	// Import map from preprocessed XML file
	{
		boost::filesystem::path file_path = "../../../PreProcessing/POLIMI_H2_1412/kinetics-POLIMI_H2_1412/kinetics.xml";

		// Open XML file containing the thermodynamic data
		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc, xml_string, file_path);

		// Import
		thermoMap = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
		tranMap = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc); 
	}

	// Transport properties
	{
		// Mixture (molar basis): H2(18%), O2(10%), H(1%), OH(1%), N2(70%)  
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
		tranMap->SetTemperature(T);
		tranMap->SetPressure(P);

		// Mixture averaged properties
		{
			// Molecular weight (in kg/kmol)
			const double mw = thermoMap->MolecularWeight_From_MoleFractions(x.data());

			// Dynamic viscosity (in kg/m/s)
			const double eta = tranMap->DynamicViscosity(x.data());

			// Thermal conductivity (in W/m/K)
			const double lambda = tranMap->ThermalConductivity(x.data());

			// Planck mean absorption coefficient (in 1/m)
			const double kPlanck = tranMap->kPlanckMix(x.data());

			// Molecular diffusion coefficients (in m2/s)
			Eigen::VectorXd GammaMix(thermoMap->NumberOfSpecies());
			tranMap->MassDiffusionCoefficients(GammaMix.data(), x.data());

			// Thermal diffusion ratios (i.e. Soret effect)
			Eigen::VectorXd TetaMix(thermoMap->NumberOfSpecies());
			tranMap->ThermalDiffusionRatios(TetaMix.data(), x.data());

			// Thermal diffusion coefficients (i.e. Soret effect) (in m2/s)
			Eigen::VectorXd GammaSoretMix(thermoMap->NumberOfSpecies());
			for (unsigned int i=0;i<thermoMap->NumberOfSpecies();i++)
				GammaSoretMix(i) = GammaMix(i)*TetaMix(i)*thermoMap->MW(i)/mw;
					

			// Print
			std::cout << std::endl;
			std::cout << "Transport properties (mixture)" << std::endl;
			std::cout << " * eta(kg/m/s):   " << std::setw(16) << std::scientific << eta << std::endl;
			std::cout << " * lambda(W/m/K): " << std::setw(16) << std::scientific << lambda << std::endl;
			std::cout << " * kPlanck(1/m):  " << std::setw(16) << std::scientific << kPlanck << std::endl;

			// Print
			std::cout << std::endl;
			std::cout << "Diffusion coefficients" << std::endl;
			for (unsigned int i=0;i<thermoMap->NumberOfSpecies();i++)
				std::cout 	<< std::setw(16) << std::left << thermoMap->NamesOfSpecies()[i]
						<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << GammaMix(i) 
						<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << TetaMix(i)
						<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << GammaSoretMix(i)
						<< std::endl;			
		}
	}
}

