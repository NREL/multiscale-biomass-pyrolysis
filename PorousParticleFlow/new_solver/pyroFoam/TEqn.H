	fvScalarMatrix TEqn
		(
		 fvm::ddt(T)
		 + fvm::div(phi, T)
		 - fvm::laplacian(vfrac*sthdiff+(1-vfrac)*gthdiff, T)
		 ==
		 fvModels.source(T)
		);

	TEqn.relax();
	fvConstraints.constrain(TEqn);
	TEqn.solve();
	fvConstraints.constrain(T);
