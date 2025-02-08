Info<< "Reading thermophysical properties\n" << endl;
autoPtr<fluidReactionThermo> pThermo_solid(fluidReactionThermo::New(mesh,"solid"));
fluidReactionThermo& thermo_solid = pThermo_solid();
thermo_solid.validate(args.executable(), "e", "solid");

basicSpecieMixture& composition_solid = thermo_solid.composition();
PtrList<volScalarField>& Y_solid = composition_solid.Y();

volScalarField rho_solid
(
    IOobject
    (
        "rho_solid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE        
    ),
    thermo_solid.rho()
);


//- Dummy solid fields
volVectorField U_solid
(
    IOobject
    (
        "U_solid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("U",dimVelocity,vector::zero)
);

surfaceScalarField phi_solid("phi_solid",fvc::flux(U_solid));

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

Info<< "Creating solid thermophysical transport model\n" << endl;
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


// IOdictionary combustionDict
// (
//     IOobject
//     (
//         "combustionProperties.solid",
//         runTime.constant(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::NO_WRITE
//     )
// );

//wordList speciesToTransferNames(combustionDict.lookup("solidToGasTransferSpecies"));

List<std::pair<bool,label>> speciesTransferList(Y_solid.size());

forAll(Y_solid, i)
{
    speciesTransferList[i].first = false;
    speciesTransferList[i].second = -1;

    forAll(Y,igas)
    {
        if(Y[igas].name() == Y_solid[i].name())
        {
            speciesTransferList[i].first = true;
            speciesTransferList[i].second = igas;
        }
    }
     
}


// forAll(Y_solid, i)
// {
//     speciesTransferList[i].first = false;
//     speciesTransferList[i].second = -1;

//     forAll(speciesToTransferNames,namei)
//     {
//         if(speciesToTransferNames[namei] == Y_solid[i].name())
//         {
//             forAll(Y,igas)
//             {
//                 if(Y[igas].name() == Y_solid[i].name())
//                 {
//                     speciesTransferList[i].first = true;
//                     speciesTransferList[i].second = igas;
//                 }
//             }

//             if(speciesTransferList[i].first == false)
//             {
//                 FatalErrorInFunction<< "Species to transfer should be in the gas phase too!\n"
//                                     << abort(FatalError);
//             }
//         }
//     }
     
// }
