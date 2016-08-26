from math import sqrt

#https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement/blob/c6d45de/MELA/src/Mela.cc#L307
CJLSTg4decay_mix = 2.521
CJLSTg2decay_mix = 1.638
CJLSTg1prime2decay_mix = 12046.01

#https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement/blob/c6d45de/MELA/src/Mela.cc#L630
CJLSTg4decay_pure = {13*13*13*13: sqrt(7.0), 11*11*11*11: sqrt(7.0), 11*11*13*13: sqrt(6.0)}
CJLSTg2decay_pure = {13*13*13*13: sqrt(2.3), 11*11*11*11: sqrt(2.3), 11*11*13*13: sqrt(2.1)}
CJLSTg1prime2decay_pure = 12046.01#?

g2decay = 1.663195
g4decay = 2.55502
g1prime2decay_gen = -12110.20   #for the sample
g1prime2decay_reco = 12110.20   #for discriminants

g2VBF = 0.271965
g4VBF = 0.297979
g1prime2VBF_gen = -2158.21
g1prime2VBF_reco = 2158.21

#https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWG/Higgs_XSBR_YR4_update.xlsx
SMXSggH = (44.14       #'YR4 SM 13TeV'!B24   (ggH cross section, m=125)
             *1000)    #                     (pb to fb)
SMBR2L2l = (5.897E-05  #'YR4 SM BR'!CO25     (2e2mu BR, m=125)
             *3)       #                     (include 2e2tau, 2mu2tau)
SMXSVBF = (3.782E+00   #'YR4 SM 13TeV'!B24   (VBF cross section, m=125)
             *1000)    #                     (pb to fb)
SMXSggH2L2l = SMXSggH * SMBR2L2l
SMXSVBF2L2l = SMXSVBF * SMBR2L2l


JHUXSggH2L2la1   = 7.1517173;      JHUXSggH2L2la1err   = 0.23044650E-03
JHUXSggH2L2la2   = 2.5849908;      JHUXSggH2L2la2err   = 0.77294379E-04
JHUXSggH2L2la3   = 1.0954288;      JHUXSggH2L2la3err   = 0.43443941E-04
JHUXSggH2L2lL1   = 0.48754258E-07; JHUXSggH2L2lL1err   = 0.15112345E-11
JHUXSggH2L2la1a2 = 26.073377;      JHUXSggH2L2la1a2err = 0.74026916E-03
JHUXSggH2L2la1a3 = 14.278034;      JHUXSggH2L2la1a3err = 0.45318122E-03
JHUXSggH2L2la1L1 = 0.23727827;     JHUXSggH2L2la1L1err = 0.15982277E-04    #using g1prime2 = -12100.42

JHUXSVBF2L2la1   = 968.674284006;     JHUXSVBF2L2la1err   = 0.075115702763
JHUXSVBF2L2la2   = 13102.7106117;     JHUXSVBF2L2la2err   = 0.522399748272
JHUXSVBF2L2la3   = 10909.5390002;     JHUXSVBF2L2la3err   = 0.50975030067
JHUXSVBF2L2lL1   = 0.000208309799883; JHUXSVBF2L2lL1err   = 1.24640942579e-08
JHUXSVBF2L2la1a2 = 2207.72848655;     JHUXSVBF2L2la1a2err = 0.126379327428
JHUXSVBF2L2la1a3 = 1937.20646111;     JHUXSVBF2L2la1a3err = 0.122617320785
JHUXSVBF2L2la1L1 = 2861.21349769;     JHUXSVBF2L2la1L1err = 0.0771278408768

if __name__ == "__main__":
    print "All of the following should be 0:"
    print
    print "  decay:"
    print "    a1XS - g2**2*a2XS          = {}%".format((JHUXSggH2L2la1 - g2decay**2           * JHUXSggH2L2la2          ) / JHUXSggH2L2la1 * 100)
    print "    a1XS - g4**2*a3XS          = {}%".format((JHUXSggH2L2la1 - g4decay**2           * JHUXSggH2L2la3          ) / JHUXSggH2L2la1 * 100)
    print "    a1XS - g1prime2**2*L1XS    = {}%".format((JHUXSggH2L2la1 - g1prime2decay_gen**2 * JHUXSggH2L2lL1          ) / JHUXSggH2L2la1 * 100)
    print "    a1XS + g4**2*a3XS - a1a3XS = {}%".format((JHUXSggH2L2la1 + g4decay**2 * JHUXSggH2L2la3 - JHUXSggH2L2la1a3 ) / JHUXSggH2L2la1 * 100)
    print "    2*a1XS - a1a3XS            = {}%".format((2*JHUXSggH2L2la1 - JHUXSggH2L2la1a3                             ) / (2*JHUXSggH2L2la1) * 100)
    print
    print "  VBF:"
    print "    a1XS - g2**2*a2XS          = {}%".format((JHUXSVBF2L2la1 - g2VBF**2           * JHUXSVBF2L2la2          ) / JHUXSVBF2L2la1 * 100)
    print "    a1XS - g4**2*a3XS          = {}%".format((JHUXSVBF2L2la1 - g4VBF**2           * JHUXSVBF2L2la3          ) / JHUXSVBF2L2la1 * 100)
    print "    a1XS - g1prime2**2*L1XS    = {}%".format((JHUXSVBF2L2la1 - g1prime2VBF_gen**2 * JHUXSVBF2L2lL1          ) / JHUXSVBF2L2la1 * 100)
    print "    a1XS + g4**2*a3XS - a1a3XS = {}%".format((JHUXSVBF2L2la1 + g4VBF**2 * JHUXSVBF2L2la3 - JHUXSVBF2L2la1a3 ) / JHUXSVBF2L2la1 * 100)
    print "    2*a1XS - a1a3XS            = {}%".format((2*JHUXSVBF2L2la1 - JHUXSVBF2L2la1a3                           ) / (2*JHUXSVBF2L2la1) * 100)

#Set them to exactly 0

#decay

values = [JHUXSggH2L2la1, JHUXSggH2L2la2*g2decay**2, JHUXSggH2L2la3*g4decay**2, JHUXSggH2L2lL1*g1prime2decay_gen**2]
errors = [JHUXSggH2L2la1err, JHUXSggH2L2la2err*g2decay**2, JHUXSggH2L2la3err*g4decay**2, JHUXSggH2L2lL1err*g1prime2decay_gen**2]

JHUXSggH2L2la1, JHUXSggH2L2la1err = sum(value/error**2 for value, error in zip(values, errors)) / sum(1/error**2 for error in errors), sum(1/error**2 for error in errors)**-.5
JHUXSggH2L2la2, JHUXSggH2L2la2err = JHUXSggH2L2la1 * g2decay**2, JHUXSggH2L2la1err * g2decay**2
JHUXSggH2L2la3, JHUXSggH2L2la3err = JHUXSggH2L2la1 * g4decay**2, JHUXSggH2L2la1err * g4decay**2
JHUXSggH2L2lL1, JHUXSggH2L2lL1err = JHUXSggH2L2la1 * g1prime2decay_gen**2, JHUXSggH2L2la1err * g1prime2decay_gen**2

JHUXSggH2L2la1a3, JHUXSggH2L2la1a3err = JHUXSggH2L2la1*2, JHUXSggH2L2la1err*2

#VBF

values = [JHUXSVBF2L2la1, JHUXSVBF2L2la2*g2VBF**2, JHUXSVBF2L2la3*g4VBF**2, JHUXSVBF2L2lL1*g1prime2VBF_gen**2]
errors = [JHUXSVBF2L2la1err, JHUXSVBF2L2la2err*g2VBF**2, JHUXSVBF2L2la3err*g4VBF**2, JHUXSVBF2L2lL1err*g1prime2VBF_gen**2]

JHUXSVBF2L2la1, JHUXSVBF2L2la1err = sum(value/error**2 for value, error in zip(values, errors)) / sum(1/error**2 for error in errors), sum(1/error**2 for error in errors)**-.5
JHUXSVBF2L2la2, JHUXSVBF2L2la2err = JHUXSVBF2L2la1 * g2VBF**2, JHUXSVBF2L2la1err * g2VBF**2
JHUXSVBF2L2la3, JHUXSVBF2L2la3err = JHUXSVBF2L2la1 * g4VBF**2, JHUXSVBF2L2la1err * g4VBF**2
JHUXSVBF2L2lL1, JHUXSVBF2L2lL1err = JHUXSVBF2L2la1 * g1prime2VBF_gen**2, JHUXSVBF2L2la1err * g1prime2VBF_gen**2

JHUXSVBF2L2la1a3, JHUXSVBF2L2la1a3err = JHUXSVBF2L2la1*2, JHUXSVBF2L2la1err*2

del values, errors
