Info<< "Reading solid thermophysical properties\n" << endl;
autoPtr<fluidReactionThermo> psThermo(fluidReactionThermo::New(mesh,"solid"));
fluidReactionThermo& thermo_solid = psThermo();
thermo_solid.validate(args.executable(), "h.solid", "e.solid");

basicSpecieMixture& composition_solid = thermo_solid.composition();
PtrList<volScalarField>& Y_solid = composition_solid.Y();

volScalarField rho_solid
(
    IOobject
    (
        "rhos",
        runTime.timeName(),
        mesh
    ),
    thermo_solid.rho()
);

Info<< "Reading field alpha\n" << endl;
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


//- Create dummy fields for solid
volVectorField U_solid
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("U",dimVelocity,vector::zero)
);

surfaceScalarField phi_solid(fvc::flux(U_solid));

autoPtr<compressible::momentumTransportModel> turbulence_solid
(
    compressible::momentumTransportModel::New
    (
        rho_solid,
        U_solid,
        phi_solid,
        thermo_solid
    )
);

autoPtr<fluidReactionThermophysicalTransportModel>
thermophysicalTransport_solid
(
    fluidReactionThermophysicalTransportModel::New
    (
        turbulence_solid(),
        thermo_solid
    )
);

Info<< "Creating solid reaction model\n" << endl;
autoPtr<combustionModel> reaction_solid(combustionModel::New(thermo_solid, turbulence_solid()));