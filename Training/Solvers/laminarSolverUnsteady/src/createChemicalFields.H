Info<< "Reading Solver options\n" << endl;
IOdictionary solverOptions
(
	IOobject
	(
		"solverOptions",
		U.time().constant(),
		U.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

const dictionary& kineticsDictionary = solverOptions.subDict("Kinetics");
Foam::string kinetics_folder = kineticsDictionary.lookup("folder");

// Open XML file containing the thermodynamic data
rapidxml::xml_document<> doc;
std::vector<char> xml_string;
OpenSMOKE::OpenInputFileXML(doc, xml_string, kinetics_folder + "/kinetics.xml");

OpenSMOKE::ThermodynamicsMap_CHEMKIN 		thermoMap(doc);
OpenSMOKE::KineticsMap_CHEMKIN 			kineticsMap(thermoMap, doc);
OpenSMOKE::TransportPropertiesMap_CHEMKIN	transportMap(doc);

// Inert species
word inertSpecies(kineticsDictionary.lookup("inertSpecies"));
label inertIndex = thermoMap.IndexOfSpecies(inertSpecies)-1;

// Detect spark
const dictionary& sparkDictionary = solverOptions.subDict("Spark");
sparkModel spark;
spark.Read(sparkDictionary);


Info << endl;
Info << " * Reading chemical species mole fractions or mass fractions..." << endl;

PtrList<volScalarField> Y(thermoMap.NumberOfSpecies());
PtrList<volScalarField> GammaMix(thermoMap.NumberOfSpecies());


// Loop over all the species in the kinetic mechanism
for(int i=0;i<thermoMap.NumberOfSpecies();i++)
{
	Info << "   Reading species: " << thermoMap.NamesOfSpecies()[i] << " > ";

	IOobject header
	(
		thermoMap.NamesOfSpecies()[i],
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);

	// Check if field exists and can be read
	if (header.headerOk())
	{
		Info << thermoMap.NamesOfSpecies()[i] << endl;
		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					thermoMap.NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);
	}
	else
	{
		Info << "Ydefault" << endl;

		volScalarField Ydefault
		(
			IOobject
			(
				"Ydefault",
				mesh.time().timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			mesh
		);

		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					thermoMap.NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				Ydefault
			)
		);
	}

	// Molecular diffusivities of species
	GammaMix.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"GammaMix_"+ thermoMap.NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("GammaMix", dimensionSet(0, 2, -1, 0, 0), 1e-16)
		)
	);
}

