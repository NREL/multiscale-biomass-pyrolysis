{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        )
      + thermophysicalTransport->divq(he)
     ==
        reaction->Qdot()
      + fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}

//! There are issues with solving the equation for the solid
//! Most of all the interphase transfer
//- Solve equation for solid
// {
//     volScalarField& he = thermo_solid.he();

//     fvScalarMatrix EEqn
//     (
//         fvm::ddt(rho_solid, he)
//       + thermophysicalTransport_solid->divq(he)
//      ==
//         reaction_solid->Qdot()
//       + fvModels.source(rho_solid, he)
//     );

//     EEqn.relax();

//     fvConstraints.constrain(EEqn);

//     EEqn.solve();

//     fvConstraints.constrain(he);

//     thermo_solid.correct();

//     Info<< "min/max(T.solid) = "
//         << min(T_solid).value() << ", " << max(T_solid).value() << endl;
// }

//! Instead of solving for he, just assume thermal equilibrium
T_solid = T;
