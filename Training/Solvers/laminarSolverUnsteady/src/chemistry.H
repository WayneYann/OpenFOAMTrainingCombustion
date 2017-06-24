{
	Info <<" * Solving homogeneous chemistry (OpenSMOKE solver)... " <<endl;

	scalar tStart = runTime.value() - runTime.deltaT().value();
	scalar tEnd = runTime.value();
	scalar dtStart = 1e-8;

	scalarField& TCells = T.ref();
	scalarField& pCells = p.ref();

	const unsigned int NC  = thermoMap.NumberOfSpecies();

	Eigen::VectorXd massFractions(NC);
	Eigen::VectorXd moleFractions(NC);
	scalarField c(NC);
	scalarField dcStart(NC);

	// Create dictionary and add the odeSolver name
	dictionary dict;
	dict.add("solver", args[1]);

	// Create the ODE system as object of type batchOdeSystem
	batchAdiabaticOdeSystem batch(thermoMap, kineticsMap);

	// Create the selected ODE system solver
	autoPtr<ODESolver> odeSolver = ODESolver::New(batch, dict);
			
	unsigned int counter = 0;
			
	double tStartCpu = OpenSMOKE::OpenSMOKEGetCpuTime();
	forAll(TCells, celli)
	{
		// From mass fractions to mole fractions
		double mw;
		for(unsigned int i=0;i<NC;i++)
			massFractions(i) = Y[i].internalField()[celli];
		thermoMap.MoleFractions_From_MassFractions(moleFractions.data(), mw, massFractions.data());

		// Concentrations [kmol/m3]
		double cTot = pCells[celli]/(PhysicalConstants::R_J_kmol*TCells[celli]);
		for(unsigned int i=0;i<thermoMap.NumberOfSpecies();i++)
			c[i] = cTot*moleFractions(i);

		// Enthalpy
		thermoMap.SetTemperature(TCells[celli]);
		thermoMap.SetPressure(pCells[celli]);
		const double H = thermoMap.hMolar_Mixture_From_MoleFractions(moleFractions.data());
		batch.setEnthalpy(H);

		// Solve ODE system
		batch.setTemperature(TCells[celli]);
		batch.setPressure(pCells[celli]);
		batch.derivatives(tStart, c, dcStart);
		odeSolver->solve(tStart, tEnd, c, dtStart);

		// From concentrations to mass fractions
		cTot = std::accumulate(c.begin(), c.end(), 0.);
		for(unsigned int i=0;i<NC;i++)
			moleFractions(i) = c[i]/cTot;
		thermoMap.MassFractions_From_MoleFractions(massFractions.data(), mw, moleFractions.data());
		for(unsigned int i=0;i<NC;i++)
			Y[i].ref()[celli] = massFractions(i);

		// Temperature
		TCells[celli] = thermoMap.GetTemperatureFromEnthalpyAndMoleFractions(H, pCells[celli], moleFractions.data(), TCells[celli]);

		// Counter
		if (counter%(int(0.20*mesh.nCells())+1) == 0)
			Info <<"   Done: " << counter << "/" << mesh.nCells() << endl;

		counter++;
	}
	
	double tEndCpu = OpenSMOKE::OpenSMOKEGetCpuTime();
			
	Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEndCpu-tStartCpu)/double(mesh.nCells())*1000. << " ms per cell)" << endl;

}
