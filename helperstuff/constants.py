from math import sqrt

#https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement/blob/master/MELA/src/Mela.cc#L307
CJLSTg4decay_mix = 2.521
CJLSTg2decay_mix = 1.638
CJLSTg1prime2decay_mix = 12046.01

#https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement/blob/master/MELA/src/Mela.cc#L630
CJLSTg4decay_pure = [sqrt(7.0), sqrt(7.0), sqrt(6.0)]#0 is 4mu, 1 is 4e, 2 is 2e2mu
CJLSTg2decay_pure = [sqrt(2.3), sqrt(2.3), sqrt(2.1)]
CJLSTg1prime2decay_pure = 12046.01#?

g4decay = 2.55052
g2decay = 1.65684
g1prime2decay_gen = -12100.42
g1prime2decay_reco = 12100.42

#https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWG/Higgs_XSBR_YR4_update.xlsx
SMXSggH = (44.14       #'YR4 SM 13TeV'!B24   (ggH cross section, m=125)
             *1000)    #                     (pb to fb)
SMBR2L2l = (5.897E-05  #'YR4 SM BR'!CO25     (2e2mu BR, m=125)
             *3)       #                     (include 2e2tau, 2mu2tau)
SMXS2L2l = SMXSggH * SMBR2L2l

JHUXS2L2la1   = 7.1517173;      JHUXS2L2la1err   = 0.23044650E-03
JHUXS2L2la2   = 2.5849908;      JHUXS2L2la2err   = 0.77294379E-04
JHUXS2L2la3   = 1.0954288;      JHUXS2L2la3err   = 0.43443941E-04
JHUXS2L2lL1   = 0.48754258E-07; JHUXS2L2lL1err   = 0.15112345E-11
JHUXS2L2la1a2 = 26.073377;      JHUXS2L2la1a2err = 0.74026916E-03
JHUXS2L2la1a3 = 14.278034;      JHUXS2L2la1a3err = 0.45318122E-03
JHUXS2L2la1L1 = 0.23727827;     JHUXS2L2la1L1err = 0.15982277E-04    #using g1prime2 = -12100.42

if __name__ == "__main__":
    print "All of the following should be 0:"
    print "a1XS - g2**2*a2XS          = {}%".format((JHUXS2L2la1 - g2decay**2           * JHUXS2L2la2       ) / JHUXS2L2la1 * 100)
    print "a1XS - g4**2*a3XS          = {}%".format((JHUXS2L2la1 - g4decay**2           * JHUXS2L2la3       ) / JHUXS2L2la1 * 100)
    print "a1XS - g1prime2**2*L1XS    = {}%".format((JHUXS2L2la1 - g1prime2decay_gen**2 * JHUXS2L2lL1       ) / JHUXS2L2la1 * 100)
    print "a1XS + g4**2*a3XS - a1a3XS = {}%".format((JHUXS2L2la1 + g4decay**2 * JHUXS2L2la3 - JHUXS2L2la1a3 ) / JHUXS2L2la1 * 100)
    print "2*a1XS - a1a3XS            = {}%".format((2*JHUXS2L2la1 - JHUXS2L2la1a3                          ) / (2*JHUXS2L2la1) * 100)

#Set them to exactly 0
values = [JHUXS2L2la1, JHUXS2L2la2*g2decay**2, JHUXS2L2la3*g4decay**2, JHUXS2L2lL1*g1prime2decay_gen**2]
errors = [JHUXS2L2la1err, JHUXS2L2la2err*g2decay**2, JHUXS2L2la3err*g4decay**2, JHUXS2L2lL1err*g1prime2decay_gen**2]

JHUXS2L2la1, JHUXS2L2la1err = sum(value/error**2 for value, error in zip(values, errors)) / sum(1/error**2 for error in errors), sum(1/error**2 for error in errors)**-.5
JHUXS2L2la2, JHUXS2L2la2err = JHUXS2L2la1 * g2decay**2, JHUXS2L2la1err * g2decay**2
JHUXS2L2la3, JHUXS2L2la3err = JHUXS2L2la1 * g4decay**2, JHUXS2L2la1err * g4decay**2
JHUXS2L2lL1, JHUXS2L2lL1err = JHUXS2L2la1 * g1prime2decay_gen**2, JHUXS2L2la1err * g1prime2decay_gen**2

JHUXS2L2la1a3, JHUXS2L2la1a3err = JHUXS2L2la1*2, JHUXS2L2la1err*2

del values, errors
