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
	// OpenSMOKE++ Thermodynamic Map	
	OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermoMap;

	// Import map from preprocessed XML file
	{
		boost::filesystem::path file_path = "../../../../PreProcessing/POLIMI_H2_1412/kinetics-POLIMI_H2_1412/kinetics.xml";

		// Open XML file containing the thermodynamic data
		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc, xml_string, file_path);

		// Import
		thermoMap = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
	}

	// Number of species
	std::cout << "Species: " << thermoMap->NumberOfSpecies() << std::endl;

	// List of species and molecular weights (in kg/kmol)
	std::cout << "List of species and molecular weights (in kg/kmol)" << std::endl;
	for (unsigned int i=0;i<thermoMap->NumberOfSpecies();i++)
		std::cout 	<< std::setw(16) << std::left << thermoMap->NamesOfSpecies()[i]
				<< std::setw(12) << std::right << std::fixed << std::setprecision(4) << thermoMap->MW(i)
				<< std::endl;

	// Thermodynamic properties
	{
		// Mixture (molar basis): H2(10%), O2(20%), N2(70%)  
		Eigen::VectorXd x(thermoMap->NumberOfSpecies());
		x.setZero();
		x(thermoMap->IndexOfSpecies("H2")-1) = 0.10;
		x(thermoMap->IndexOfSpecies("O2")-1) = 0.20;
		x(thermoMap->IndexOfSpecies("N2")-1) = 0.70;

		// Operating conditions
		const double T = 1000.;
		const double P = 101325.;

		// Set thermodynamic map
		thermoMap->SetTemperature(T);
		thermoMap->SetPressure(P);

		// Mixture averaged properties
		{
			// Molecular weight (in kg/kmol)
			const double mw = thermoMap->MolecularWeight_From_MoleFractions(x.data());

			// Constant pressure specific heat (in J/kmol/K)
			const double cp = thermoMap->cpMolar_Mixture_From_MoleFractions(x.data());

			// Enthalpy (in J/kmol)
			const double h = thermoMap->hMolar_Mixture_From_MoleFractions(x.data());

			// Entropy (in J/kmol/K)
			const double s = thermoMap->sMolar_Mixture_From_MoleFractions(x.data());

			// Internal energy (in J/kmol)
			const double u = thermoMap->uMolar_Mixture_From_MoleFractions(x.data());

			// Gibb's free energy (in J/kmol)
			const double g = thermoMap->gMolar_Mixture_From_MoleFractions(x.data());

			// Helmotz's free energy (in J/kmol)
			const double a = thermoMap->aMolar_Mixture_From_MoleFractions(x.data());

			// Print
			std::cout << std::endl;
			std::cout << "Thermodynamic properties of mixture (mass units)" << std::endl;
			std::cout << " * MW: " << std::setw(16) << std::scientific << mw << std::endl;
			std::cout << " * cp: " << std::setw(16) << std::scientific << cp/mw << std::endl;
			std::cout << " * h:  " << std::setw(16) << std::scientific << h/mw  << std::endl;
			std::cout << " * s:  " << std::setw(16) << std::scientific << s/mw  << std::endl;
			std::cout << " * u:  " << std::setw(16) << std::scientific << u/mw  << std::endl;
			std::cout << " * g:  " << std::setw(16) << std::scientific << g/mw  << std::endl;
			std::cout << " * a:  " << std::setw(16) << std::scientific << a/mw  << std::endl;
		}

		// Species properties
		{

			// Constant pressure specific heat (in J/kmol/K)
			Eigen::VectorXd cp_species(thermoMap->NumberOfSpecies());
			thermoMap->cpMolar_Species(cp_species.data());
	
			// Enthalpy (in J/kmol)
			Eigen::VectorXd h_species(thermoMap->NumberOfSpecies());
			thermoMap->hMolar_Species(h_species.data());

			// Entropy (in J/kmol/K)
			Eigen::VectorXd s_species(thermoMap->NumberOfSpecies());
			thermoMap->sMolar_Species(s_species.data());

			// Internal energy (in J/kmol)
			Eigen::VectorXd u_species(thermoMap->NumberOfSpecies());
			thermoMap->uMolar_Species(u_species.data());

			// Gibb's free energy (in J/kmol)
			Eigen::VectorXd g_species(thermoMap->NumberOfSpecies());
			thermoMap->gMolar_Species(g_species.data());

			// Helmotz's free energy (in J/kmol)
			Eigen::VectorXd a_species(thermoMap->NumberOfSpecies());
			thermoMap->aMolar_Species(a_species.data());

			// Print
			std::cout << std::endl;
			std::cout << "Thermodynamic properties of single species (molar units)" << std::endl;
			for (unsigned int i=0;i<thermoMap->NumberOfSpecies();i++)
				std::cout 	<< std::setw(16) << std::left << thermoMap->NamesOfSpecies()[i]
						<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << cp_species(i) 
						<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << h_species(i) 
						<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << s_species(i) 
						<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << u_species(i) 
						<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << g_species(i) 
						<< std::setw(16) << std::right << std::scientific << std::setprecision(6) << a_species(i) 
						<< std::endl;
		}

		// Temperature from enthalpy
		{
			// Enthalpy (in J/kmol)
			double h = thermoMap->hMolar_Mixture_From_MoleFractions(x.data());
			
			// Small perturbation (+5%)
			h *= 1.05;

			// Temperature from enthalpy
			const double TFirstGuess = 300.;
			const double TStar = thermoMap->GetTemperatureFromEnthalpyAndMoleFractions(h, P, x.data(), TFirstGuess);

			// Print
			std::cout << std::endl;
			std::cout << "Temperature from enthalpy" << std::endl;
			std::cout << " * H(J/kmol): " << std::setw(16) << std::scientific << h << std::endl;
			std::cout << " * T(K):      " << std::setw(16) << std::scientific << TStar << std::endl;
		}
	}

	// Using different types of vectors
	{
		// STL vectors
		{
			// Mixture (molar basis): H2(10%), O2(20%), N2(70%)  
			std::vector<double> x(thermoMap->NumberOfSpecies());
			std::fill(x.begin(), x.end(), 0.);
			x[thermoMap->IndexOfSpecies("H2")-1] = 0.10;
			x[thermoMap->IndexOfSpecies("O2")-1] = 0.20;
			x[thermoMap->IndexOfSpecies("N2")-1] = 0.70;

			// Operating conditions
			const double T = 1000.;
			const double P = 101325.;

			// Set thermodynamic map
			thermoMap->SetTemperature(T);
			thermoMap->SetPressure(P);

			// Molecular weight (in kg/kmol)
			const double mw = thermoMap->MolecularWeight_From_MoleFractions(x.data());

			// Print
			std::cout << "MW from STL vector: " << mw << std::endl;
		}

		// OpenSMOKE++ Vectors
		// Warning: they start from 1 (FORTRAN style)
		{
			// Mixture (molar basis): H2(10%), O2(20%), N2(70%)  
			OpenSMOKE::OpenSMOKEVectorDouble x(thermoMap->NumberOfSpecies());
			x = 0.;
			x[thermoMap->IndexOfSpecies("H2")] = 0.10;
			x[thermoMap->IndexOfSpecies("O2")] = 0.20;
			x[thermoMap->IndexOfSpecies("N2")] = 0.70;

			// Operating conditions
			const double T = 1000.;
			const double P = 101325.;

			// Set thermodynamic map
			thermoMap->SetTemperature(T);
			thermoMap->SetPressure(P);

			// Molecular weight (in kg/kmol)
			const double mw = thermoMap->MolecularWeight_From_MoleFractions(x.GetHandle());

			// Print
			std::cout << "MW from OpenSMOKE++ vector: " << mw << std::endl;
		}
	}
}

