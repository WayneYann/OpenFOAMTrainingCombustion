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
