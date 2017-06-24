//- Creating field rho [kg/m3]
Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p*psi
);

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, simple.dict());
mesh.setFluxRequired(p.name());
dimensionedScalar initialMass = fvc::domainIntegrate(rho);

// Kinetic energy
Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

// Interpolation scheme
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(Y, i)
{
	fields.add(Y[i]);
}

// Spark
spark.Setup(mesh);
