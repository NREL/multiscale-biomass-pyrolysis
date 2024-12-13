	fvScalarMatrix TEqn
		(
		 fvm::ddt(T)
		 + fvm::div(phi, T)
		 - fvm::laplacian(thdiff, T)
		 ==
		 fvModels.source(T)
		);

	TEqn.relax();
	fvConstraints.constrain(TEqn);
	TEqn.solve();
	fvConstraints.constrain(T);
