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
	Info <<" * Updating chemistry (implicit)... " <<endl;

	double tStartCpu = OpenSMOKE::OpenSMOKEGetCpuTime();	
	{
		const unsigned int NC = thermoMap.NumberOfSpecies();

		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField(); 

		Eigen::VectorXd J(NC+1);
		Eigen::VectorXd Source(NC+1);	
		Eigen::VectorXd y(NC+1);

		forAll(TCells, celli)
		{
			for(int i=0;i<NC;i++)
				y(i) = Y[i].internalField()[celli];
			y(NC) = TCells[celli];

			chemistry.reactionSourceTerms(thermoMap, kineticsMap, y, pCells[celli], Source);
			chemistry.reactionJacobian( thermoMap, kineticsMap, y, pCells[celli], J );

			for(int i=0;i<NC+1;i++)
				sourceImplicit[i].ref()[celli] = J(i);

			for(int i=0;i<NC+1;i++)
				sourceExplicit[i].ref()[celli] = Source(i) - J(i)*y(i);
		}
	}
	double tEndCpu = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEndCpu - tStartCpu << " s  (" << (tEndCpu-tStartCpu)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}
