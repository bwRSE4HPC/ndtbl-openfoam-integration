//- Passive scalar transport for Param1
fvScalarMatrix Param1Eqn
(
      fvm::ddt(rhoParam_, Param1_)
    + fvm::div(phi, Param1_)
     - fvm::laplacian(rhoParam_*gammaParam_, Param1_)
);

Param1Eqn.solve();

//- Passive scalar transport for Param2
fvScalarMatrix Param2Eqn
(
      fvm::ddt(rhoParam_, Param2_)
    + fvm::div(phi, Param2_)
     - fvm::laplacian(rhoParam_*gammaParam_, Param2_)
);

Param2Eqn.solve();

//- Passive scalar transport for Param3
fvScalarMatrix Param3Eqn
(
      fvm::ddt(rhoParam_, Param3_)
    + fvm::div(phi, Param3_)
     - fvm::laplacian(rhoParam_*gammaParam_, Param3_)
);

Param3Eqn.solve();

//- Passive scalar transport for Param4
fvScalarMatrix Param4Eqn
(
      fvm::ddt(rhoParam_, Param4_)
    + fvm::div(phi, Param4_)
     - fvm::laplacian(rhoParam_*gammaParam_, Param4_)
);

Param4Eqn.solve();
