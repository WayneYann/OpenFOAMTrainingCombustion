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

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
	fvm::ddt(rho, U) 
	+ fvm::div(phi, U)
	+ MRF.DDt(rho, U)
	- fvm::laplacian(eta, U)
	- fvc::div(eta*dev2(fvc::grad(U)().T()))
	==
          rho*g
	+ fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
	solve(UEqn == -fvc::grad(p));

	fvOptions.correct(U);
	K = 0.5*magSqr(U);
}

Info << " Umin: " << min(U).value() << " Umax: " << max(U).value() << endl;
