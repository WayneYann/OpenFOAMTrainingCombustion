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

PtrList<volScalarField> sourceImplicit(thermoMap.NumberOfSpecies()+1);
PtrList<volScalarField> sourceExplicit(thermoMap.NumberOfSpecies()+1);

for (int i=0;i<thermoMap.NumberOfSpecies();i++)
{
	sourceImplicit.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"SI_" + thermoMap.NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("SI", dimensionSet(1, -3, -1, 0, 0), 0.0)
		)
	);

	sourceExplicit.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"SE_" + thermoMap.NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("SE", dimensionSet(1, -3, -1, 0, 0), 0.0)
		)
	);
}

sourceImplicit.set
(
	thermoMap.NumberOfSpecies(),
	new volScalarField
	(
		IOobject
		(
			"SI_T",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("SIT", dimensionSet(1, -1, -3, -1, 0), 0.0)
	)
);

sourceExplicit.set
(
	thermoMap.NumberOfSpecies(),
	new volScalarField
	(
		IOobject
		(
			"SE_T",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("SET", dimensionSet(1, -1, -3, 0, 0), 0.0)
	)
);
