{
	fvScalarMatrix TEqn
	(
		  cp*fvm::ddt(rho,T)
		+ cp*fvm::div(phi, T)		    
		 ==
		  fvm::laplacian(lambda,T)
		+ cp*fvOptions(rho, T)
	);

	TEqn.relax();
	fvOptions.constrain(TEqn);
	TEqn.solve(mesh.solver("T"));
	fvOptions.correct(T);
		

	Info<< " * T gas min/max (after transport)   = " << min(T).value() << ", " << max(T).value() << endl;

	// Spark
	spark.SetSpark(runTime.value(), T);
}
