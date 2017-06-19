{
	fvScalarMatrix TEqn
	(
		cp*fvm::div(phi, T)		    
		 ==
		fvm::laplacian(lambda,T)
		+ sourceExplicit[Y.size()] 
                + fvm::Sp(sourceImplicit[Y.size()],T)
		+ cp*fvOptions(rho, T)
	);

	TEqn.relax();
	fvOptions.constrain(TEqn);
	TEqn.solve(mesh.solver("T"));
	fvOptions.correct(T);

	T.max(260.);
	T.min(5000.);
	Info<< " * T gas min/max (after transport)   = " << min(T).value() << ", " << max(T).value() << endl;

	// Spark
	spark.SetSpark(runTime.value(), T);
}
