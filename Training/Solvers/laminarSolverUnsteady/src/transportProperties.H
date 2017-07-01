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

{
	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info<< "Properties evaluation... " ;

	Eigen::VectorXd GammaMixVector(thermoMap.NumberOfSpecies());
	Eigen::VectorXd massFractions(thermoMap.NumberOfSpecies());
	Eigen::VectorXd moleFractions(thermoMap.NumberOfSpecies());

	// Internal fields
	{
		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField(); 
		scalarField& psiCells =  psi.ref();
		scalarField& cpCells = cp.ref();
		scalarField& etaCells  =  eta.ref();
		scalarField& lambdaCells = lambda.ref();
	
		forAll(TCells, celli)
		{
			thermoMap.SetPressure(pCells[celli]);
			thermoMap.SetTemperature(TCells[celli]);
			
			transportMap.SetPressure(pCells[celli]);
			transportMap.SetTemperature(TCells[celli]);
				
			// From mass fractions to mole fractions
			double mw;
			for(int i=0;i<thermoMap.NumberOfSpecies();i++)
				massFractions(i) = Y[i].internalField()[celli];
			thermoMap.MoleFractions_From_MassFractions(moleFractions.data(), mw, massFractions.data());

			// Psi variable [s2/m2]
			psiCells[celli]  = mw/(PhysicalConstants::R_J_kmol*TCells[celli]);

			// Constant pressure specific heat [J/kg/K]
			cpCells[celli] = thermoMap.cpMolar_Mixture_From_MoleFractions(moleFractions.data())/mw;

			// Dynamic viscosity [kg/m/s]
			etaCells[celli] = transportMap.DynamicViscosity(moleFractions.data());

			// Thermal conductivity [W/m/K]
			lambdaCells[celli] = transportMap.ThermalConductivity(moleFractions.data());
			
			// Diffusion coefficients [m2/s]
			transportMap.MassDiffusionCoefficients(GammaMixVector.data(), moleFractions.data());
			for(int i=0;i<thermoMap.NumberOfSpecies();i++)
				GammaMix[i].ref()[celli] = GammaMixVector(i);
		}
    	}

	// Boundaries
	forAll(T.boundaryField(), patchi)
	{
		fvPatchScalarField& pp = p.boundaryFieldRef()[patchi];
		fvPatchScalarField& pT = T.boundaryFieldRef()[patchi];
		fvPatchScalarField& ppsi  = psi.boundaryFieldRef()[patchi];
		fvPatchScalarField& pcp = cp.boundaryFieldRef()[patchi];
		fvPatchScalarField& peta   = eta.boundaryFieldRef()[patchi];
		fvPatchScalarField& plambda = lambda.boundaryFieldRef()[patchi];	
        
		forAll(pT, facei)
		{
			thermoMap.SetPressure(pp[facei]);
			thermoMap.SetTemperature(pT[facei]);

			transportMap.SetPressure(pp[facei]);
			transportMap.SetTemperature(pT[facei]);
			
			// From mass fractions to mole fractions
			double mw;
			for(int i=0;i<thermoMap.NumberOfSpecies();i++)
				massFractions(i) = Y[i].boundaryField()[patchi][facei];
			thermoMap.MoleFractions_From_MassFractions(moleFractions.data(), mw, massFractions.data());

			// Psi variable [s2/m2]
			ppsi[facei]  = mw/(PhysicalConstants::R_J_kmol*pT[facei]);

			// Constant pressure specific heat [J/kg/K]
			pcp[facei] = thermoMap.cpMolar_Mixture_From_MoleFractions(moleFractions.data())/mw;

			// Dynamic viscosity [kg/m/s]
			peta[facei] = transportMap.DynamicViscosity(moleFractions.data());

			// Thermal conductivity [W/m/K]
			plambda[facei] = transportMap.ThermalConductivity(moleFractions.data());
				
			// Diffusion coefficients [m2/s]
			transportMap.MassDiffusionCoefficients(GammaMixVector.data(), moleFractions.data());
			for(int i=0;i<thermoMap.NumberOfSpecies();i++)
				GammaMix[i].boundaryFieldRef()[patchi][facei] = GammaMixVector(i);
		}
	}

	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}
