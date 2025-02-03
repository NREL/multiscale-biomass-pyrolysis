tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

reaction->correct();

forAll(Y, i)
{
    if (composition.solve(i))
    {
        volScalarField& Yi = Y[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(phi, Yi)
          + thermophysicalTransport->divj(Yi)
         ==
            reaction->R(Yi)
          + fvModels.source(rho, Yi)
        );

        YiEqn.relax();

        fvConstraints.constrain(YiEqn);

        YiEqn.solve("Yi");

        fvConstraints.constrain(Yi);
    }
}

composition.normalise();

//- Solve solid reactions
reaction_solid->correct();

forAll(Y_solid,i)
{
    if(composition.solve(i))
    {

        volScalarField& Yi = Y_solid[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(alpha*rho_solid, Yi)
         ==
            reaction_solid->R(Yi)
          + fvModels.source(rho_solid, Yi)
        );

        YiEqn.relax();

        fvConstraints.constrain(YiEqn);

        YiEqn.solve("Yi");

        fvConstraints.constrain(Yi);        
    }
}

composition_solid.normalise();

//- Recompute alpha
{
    volScalarField frac_vol_lost = (thermo_solid.rho() - rho_solid)/rho_solid;

    alpha = alpha*(1.0 - frac_vol_lost);

    rho_solid = thermo_solid.rho();
}
