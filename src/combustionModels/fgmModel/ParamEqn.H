//- Passive scalar transport for Param_
fvScalarMatrix ParamEqn
(
      fvm::ddt(rhoParam_, Param_)
    + fvm::div(phi, Param_)
    - fvm::laplacian(rhoParam_*gammaParam_, Param_)
);

ParamEqn.solve();
