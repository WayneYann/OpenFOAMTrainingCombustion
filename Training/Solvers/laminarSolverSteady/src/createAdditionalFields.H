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

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);

dimensionedScalar rhoMax(simple.dict().lookup("rhoMax"));
dimensionedScalar rhoMin(simple.dict().lookup("rhoMin"));

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
