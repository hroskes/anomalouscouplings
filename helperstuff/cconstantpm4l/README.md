This is a standalone piece to calculate the c constants for D_bkg_m4l (the D_bkg_kin part is from Ulascan).

Procedure:
1. print the supermela probabilities
2. change these lines to from 105-140:
   https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement/blob/71f2f5458464e45ad85a9386be2c1247332476da/MELA/src/SuperMELA.cc#L57
   https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement/blob/71f2f5458464e45ad85a9386be2c1247332476da/MELA/src/SuperMELA.cc#L147-L149
3. scram b and print supermela probabilities again
4. get the ratio
