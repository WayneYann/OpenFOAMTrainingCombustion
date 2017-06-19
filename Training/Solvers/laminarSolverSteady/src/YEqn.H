// Discretization schemes
tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi)")
    )
);

// Transport equations for species
{
    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
    
    volScalarField Yt = 0.0*Y[0];

    for (label i=0; i<Y.size(); i++)
    {
        if (i != inertIndex)
        {
           	volScalarField& Yi = Y[i];
	    	volScalarField& GammaMixi = GammaMix[i];	 
		
		fvScalarMatrix YiEqn
		(
		    	mvConvection->fvmDiv(phi, Yi)
		      - fvm::laplacian(rho*GammaMixi, Yi) 
	                ==  
		        sourceExplicit[i] 
		      + fvm::Sp(sourceImplicit[i],Yi)						
		      +	fvOptions(rho, Yi)
		);

		// Solve
		YiEqn.relax();
		fvOptions.constrain(YiEqn);
		YiEqn.solve(mesh.solver("Yi"));
		fvOptions.correct(Yi);	
		
		// Sum of mass fractions
	    	Yi.max(0.0);
	   	Yt += Yi;
        }
    }

    Info << "Inert species is " << Y[inertIndex].name() << " with local index equal to " << inertIndex << endl;
    Y[inertIndex] = scalar(1.0) - Yt;
    Y[inertIndex].max(0.0);

    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
	
    Info << "Transport equations of species solved in " << tEnd - tStart << " s " << endl;
}
