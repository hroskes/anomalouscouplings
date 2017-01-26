#!/usr/bin/env python
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
ghzgs1prime2decay_gen = -7613.351302119843
ghzgs1prime2decay_reco = 7613.351302119843

g2VBF = 0.271965
g4VBF = 0.297979
g1prime2VBF_gen = -2158.21
g1prime2VBF_reco = 2158.21
ghzgs1prime2VBF_gen = -4091.051456694223
ghzgs1prime2VBF_reco = 4091.051456694223

g2ZH = 0.112481
g4ZH = 0.144057
g1prime2ZH_gen = -517.788
g1prime2ZH_reco = 517.788
ghzgs1prime2ZH_gen = -642.9534550379002
ghzgs1prime2ZH_reco = 642.9534550379002

g2WH = 0.0998956
g4WH = 0.1236136
g1prime2WH_gen = -525.274
g1prime2WH_reco = 525.274
ghzgs1prime2WH_gen = -1
ghzgs1prime2WH_reco = 1

ghg4HJJ = 1.0062
kappa_tilde_ttH = 1.6

#https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWG/Higgs_XSBR_YR4_update.xlsx
SMXSggH  = (44.14      #'YR4 SM 13TeV'!B24   (ggH cross section, m=125)
             *1000)    #                     (pb to fb)
SMBR2L2l = (5.897E-05  #'YR4 SM BR'!CO25     (2e2mu BR, m=125)
             *3)       #                     (include 2e2tau, 2mu2tau)
SMXSVBF  = (3.782E+00  #'YR4 SM 13TeV'!B24   (VBF cross section, m=125)
             *1000)    #                     (pb to fb)
SMXSWH   = (1.373E+00  #'YR4 SM 13TeV'!R24   (WH cross section, m=125)
             *1000)    #                     (pb to fb)
SMXSWpH  = (8.400E-01  #'YR4 SM 13TeV'!X24   (W+H cross section, m=125)
             *1000)    #                     (pb to fb)
SMXSWmH  = (5.328E-01  #'YR4 SM 13TeV'!X24   (W-H cross section, m=125)
             *1000)    #                     (pb to fb)
SMXSZH   = (8.839E-01  #'YR4 SM 13TeV'!AB24  (ZH cross section, m=125)
             *1000)    #                     (pb to fb)
SMXSttH  = (5.071E-01  #'YR4 SM 13TeV'!AK24  (ttH cross section, m=125)
             *1000)    #                     (pb to fb)

SMXSggH2L2l = SMXSggH * SMBR2L2l
SMXSVBF2L2l = SMXSVBF * SMBR2L2l
SMXSZH2L2l = SMXSZH * SMBR2L2l
SMXSWH2L2l = SMXSWH * SMBR2L2l
SMXSWpH2L2l = SMXSWpH * SMBR2L2l
SMXSWmH2L2l = SMXSWmH * SMBR2L2l
SMXSttH2L2l = SMXSttH * SMBR2L2l

#############
from datetime import date
if date.today() > date(2017, 1, 26):
    raise ValueError("Update mix L1Zg xsecs!")
#############

JHUXSggH2L2la1           = 7.1517173;      JHUXSggH2L2la1err            = 0.23044650E-03
JHUXSggH2L2la2           = 2.5849908;      JHUXSggH2L2la2err            = 0.77294379E-04
JHUXSggH2L2la3           = 1.0954288;      JHUXSggH2L2la3err            = 0.43443941E-04
JHUXSggH2L2lL1           = 0.48754258E-07; JHUXSggH2L2lL1err            = 0.15112345E-11
JHUXSggH2L2la1a2         = 26.073377;      JHUXSggH2L2la1a2err          = 0.74026916E-03
JHUXSggH2L2la1a3         = 14.278034;      JHUXSggH2L2la1a3err          = 0.45318122E-03
JHUXSggH2L2la1L1         = 0.23727827;     JHUXSggH2L2la1L1err          = 0.15982277E-04    #using g1prime2 = -12100.42
JHUXSggH2L2la1_photoncut = 7.1542865;      JHUXSggH2L2la1err_photoncut  = 0.27583899E-03    #This is bigger than the uncut xsec.  Fix below.
JHUXSggH2L2lL1Zg         = 0.12338393E-06; JHUXSggH2L2lL1Zgerr          = 0.58002047E-11
JHUXSggH2L2la1L1Zg       = 13.169996;      JHUXSggH2L2la1L1Zgerr        = 0.18226318E-01

JHUXSVBFa1           = 968.674284006;      JHUXSVBFa1err            = 0.075115702763
JHUXSVBFa2           = 13102.7106117;      JHUXSVBFa2err            = 0.522399748272
JHUXSVBFa3           = 10909.5390002;      JHUXSVBFa3err            = 0.50975030067
JHUXSVBFL1           = 2.083097999e-4;     JHUXSVBFL1err            = 1.24640942579e-08
JHUXSVBFa1a2         = 2207.72848655;      JHUXSVBFa1a2err          = 0.126379327428
JHUXSVBFa1a3         = 1937.20646111;      JHUXSVBFa1a3err          = 0.122617320785
JHUXSVBFa1L1         = 2861.21349769;      JHUXSVBFa1L1err          = 0.0771278408768
JHUXSVBFa1_photoncut = 834.24595;          JHUXSVBFa1err_photoncut  = 0.41631690
JHUXSVBFL1Zg         = 0.49845301E-04;     JHUXSVBFL1Zgerr          = 0.21806623E-07
JHUXSVBFa1L1Zg       = 1394.3489;          JHUXSVBFa1L1Zgerr        = 10.828923

JHUXSZHa1           = 9022.36;        JHUXSZHa1err            = 1.17
JHUXSZHa2           = 713123;         JHUXSZHa2err            = 103
JHUXSZHa3           = 434763.7;       JHUXSZHa3err            = 62.2
JHUXSZHL1           = 33652.46e-6;    JHUXSZHL1err            = 4.19e-6
JHUXSZHa1a2         = 4258.966;       JHUXSZHa1a2err          = 0.783
JHUXSZHa1a3         = 18040.66;       JHUXSZHa1a3err          = 2.73
JHUXSZHa1L1         = 6852.307;       JHUXSZHa1L1err          = 0.929
JHUXSZHa1_photoncut = 1437150.9;      JHUXSZHa1err_photoncut  = 289.72492
JHUXSZHL1Zg         = 3.4592597;      JHUXSZHL1Zgerr          = 0.56602360E-03
JHUXSZHa1L1Zg       = 3479155.1;      JHUXSZHa1L1Zgerr        = 6901.4960

JHUXSWHa1           = 30998.54;       JHUXSWHa1err           = 2.50
JHUXSWHa2           = 3106339;        JHUXSWHa2err           = 308
JHUXSWHa3           = 2028656;        JHUXSWHa3err           = 191
JHUXSWHL1           = 11234.91e-5;    JHUXSWHL1err           = 1.10e-5
JHUXSWHa1a2         = 16486.68;       JHUXSWHa1a2err         = 2.01
JHUXSWHa1a3         = 62001.57;       JHUXSWHa1a3err         = 5.54
JHUXSWHa1L1         = 25302.37;       JHUXSWHa1L1err         = 2.67
JHUXSWHa1_photoncut = None;           JHUXSWHa1err_photoncut = None
JHUXSWHL1Zg         = 0;              JHUXSWHL1Zgerr         = 0
JHUXSWHa1L1Zg       = None;           JHUXSWHa1errL1Zg       = None

#VH cross sections changed somewhere between c85a387eaf3a8a92f5893e5293ed3c3d36107e16 and fbf449150f4df49f66b21c1638adb02b68a308d0
#correct for that
def fixVH():
    VHxsecratio = 158.49737883504775; VHxsecratioerr = 0.39826108447619935
    for V in "Z", "W":
        for coupling in "a1", "a2", "a3", "L1", "a1a2", "a1a3", "a1L1":
            name = "JHUXS{}H{}".format(V, coupling); errorname = name+"err"
            newvalue = globals()[name] * VHxsecratio
            newerror = sqrt(globals()[name] * VHxsecratioerr + VHxsecratio * globals()[errorname])
            globals()[name] = newvalue
            globals()[errorname] = newerror
fixVH(); del fixVH
JHUXSWHa1_photoncut = JHUXSWHa1;      JHUXSWHa1err_photoncut = JHUXSWHa1err
JHUXSWHa1L1Zg       = JHUXSWHa1;      JHUXSWHa1errL1Zg       = JHUXSWHa1err


JHUXSHJJa2       = 14583.61;       JHUXSHJJa2err       = 0.94
JHUXSHJJa3       = 14397.13;       JHUXSHJJa3err       = 0.97
JHUXSHJJa2a3     = 29169.2;        JHUXSHJJa2a3err     = 2.1

JHUXSttHkappa    = 0.912135589;    JHUXSttHkappaerr    = 0.00143032
JHUXSttHkappatilde = 0.35609194;   JHUXSttHkappatildeerr = 0.000492662
JHUXSttHkappakappatilde = 1.8231162489;   JHUXSttHkappakappatildeerr = 0.00254131

if JHUXSggH2L2la1_photoncut > JHUXSggH2L2la1: JHUXSggH2L2la1_photoncut = JHUXSggH2L2la1; JHUXSggH2L2la1err_photoncut = JHUXSggH2L2la1err
if JHUXSVBFa1_photoncut > JHUXSVBFa1: JHUXSVBFa1_photoncut = JHUXSVBFa1; JHUXSVBFa1err_photoncut = JHUXSVBFa1err
if JHUXSZHa1_photoncut > JHUXSZHa1: JHUXSZHa1_photoncut = JHUXSZHa1; JHUXSZHa1err_photoncut = JHUXSZHa1err

if __name__ == "__main__":
    print "All of the following should be 0:"
    print
    print "  decay:"
    print "    a1XS - g2**2*a2XS             = {}%".format((JHUXSggH2L2la1           - g2decay**2               * JHUXSggH2L2la2  ) / JHUXSggH2L2la1           * 100)
    print "    a1XS - g4**2*a3XS             = {}%".format((JHUXSggH2L2la1           - g4decay**2               * JHUXSggH2L2la3  ) / JHUXSggH2L2la1           * 100)
    print "    a1XS - g1prime2**2*L1XS       = {}%".format((JHUXSggH2L2la1           - g1prime2decay_gen**2     * JHUXSggH2L2lL1  ) / JHUXSggH2L2la1           * 100)
    print "    a1XS - ghzgs1prime2**2*L1ZgXS = {}%".format((JHUXSggH2L2la1_photoncut - ghzgs1prime2decay_gen**2 * JHUXSggH2L2lL1Zg) / JHUXSggH2L2la1_photoncut * 100)
    print "    a1XS + g4**2*a3XS - a1a3XS    = {}%".format((JHUXSggH2L2la1 + g4decay**2 * JHUXSggH2L2la3 - JHUXSggH2L2la1a3 ) / JHUXSggH2L2la1 * 100)
    print "    2*a1XS - a1a3XS               = {}%".format((2*JHUXSggH2L2la1 - JHUXSggH2L2la1a3                             ) / (2*JHUXSggH2L2la1) * 100)
    print
    print "  VBF:"
    print "    a1XS - g2**2*a2XS             = {}%".format((JHUXSVBFa1           - g2VBF**2               * JHUXSVBFa2  ) / JHUXSVBFa1           * 100)
    print "    a1XS - g4**2*a3XS             = {}%".format((JHUXSVBFa1           - g4VBF**2               * JHUXSVBFa3  ) / JHUXSVBFa1           * 100)
    print "    a1XS - g1prime2**2*L1XS       = {}%".format((JHUXSVBFa1           - g1prime2VBF_gen**2     * JHUXSVBFL1  ) / JHUXSVBFa1           * 100)
    print "    a1XS - ghzgs1prime2**2*L1ZgXS = {}%".format((JHUXSVBFa1_photoncut - ghzgs1prime2VBF_gen**2 * JHUXSVBFL1Zg) / JHUXSVBFa1_photoncut * 100)
    print "    a1XS + g4**2*a3XS - a1a3XS    = {}%".format((JHUXSVBFa1 + g4VBF**2 * JHUXSVBFa3 - JHUXSVBFa1a3 ) / JHUXSVBFa1 * 100)
    print "    2*a1XS - a1a3XS               = {}%".format((2*JHUXSVBFa1 - JHUXSVBFa1a3                       ) / (2*JHUXSVBFa1) * 100)
    print "    W+H + W-H - WH                = {}%".format((SMXSWpH + SMXSWmH - SMXSWH) / SMXSWH)
    print
    print "  ZH:"
    print "    a1XS - g2**2*a2XS             = {}%".format((JHUXSZHa1           - g2ZH**2               * JHUXSZHa2  ) / JHUXSZHa1 * 100)
    print "    a1XS - g4**2*a3XS             = {}%".format((JHUXSZHa1           - g4ZH**2               * JHUXSZHa3  ) / JHUXSZHa1 * 100)
    print "    a1XS - g1prime2**2*L1XS       = {}%".format((JHUXSZHa1           - g1prime2ZH_gen**2     * JHUXSZHL1  ) / JHUXSZHa1 * 100)
    print "    a1XS - ghzgs1prime2**2*L1ZgXS = {}%".format((JHUXSZHa1_photoncut - ghzgs1prime2ZH_gen**2 * JHUXSZHL1Zg) / JHUXSZHa1 * 100)
    print "    a1XS + g4**2*a3XS - a1a3XS    = {}%".format((JHUXSZHa1 + g4ZH**2 * JHUXSZHa3 - JHUXSZHa1a3 ) / JHUXSZHa1 * 100)
    print "    2*a1XS - a1a3XS               = {}%".format((2*JHUXSZHa1 - JHUXSZHa1a3                     ) / (2*JHUXSZHa1) * 100)
    print
    print "  WH:"
    print "    a1XS - g2**2*a2XS          = {}%".format((JHUXSWHa1 - g2WH**2           * JHUXSWHa2     ) / JHUXSWHa1 * 100)
    print "    a1XS - g4**2*a3XS          = {}%".format((JHUXSWHa1 - g4WH**2           * JHUXSWHa3     ) / JHUXSWHa1 * 100)
    print "    a1XS - g1prime2**2*L1XS    = {}%".format((JHUXSWHa1 - g1prime2WH_gen**2 * JHUXSWHL1     ) / JHUXSWHa1 * 100)
    print "    a1XS + g4**2*a3XS - a1a3XS = {}%".format((JHUXSWHa1 + g4WH**2 * JHUXSWHa3 - JHUXSWHa1a3 ) / JHUXSWHa1 * 100)
    print "    2*a1XS - a1a3XS            = {}%".format((2*JHUXSWHa1 - JHUXSWHa1a3                     ) / (2*JHUXSWHa1) * 100)
    print "    WpHXS + WmHXS - WHXS       = {}%".format((SMXSWpH2L2l + SMXSWmH2L2l - SMXSWH2L2l        ) / SMXSWH2L2l * 100)
    print
    print "  HJJ:"
    print "    a2XS - g4**2*a3XS          = {}%".format((JHUXSHJJa2 - ghg4HJJ**2              * JHUXSHJJa3   ) / JHUXSHJJa2 * 100)
    print "    a2XS + g4**2*a3XS - a2a3XS = {}%".format((JHUXSHJJa2 + ghg4HJJ**2 * JHUXSHJJa3 - JHUXSHJJa2a3 ) / JHUXSHJJa2 * 100)
    print "    2*a2XS - a2a3XS            = {}%".format((2*JHUXSHJJa2 - JHUXSHJJa2a3                         ) / (2*JHUXSHJJa2) * 100)
    print
    print "  ttH:"
    print "    kappaXS - kappa_tilde**2*kappatildeXS                     = {}%".format((JHUXSttHkappa - kappa_tilde_ttH**2                      * JHUXSttHkappatilde      ) / JHUXSttHkappa * 100)
    print "    kappaXS + kappa_tilde**2*kappatildeXS - kappakappatildeXS = {}%".format((JHUXSttHkappa + kappa_tilde_ttH**2 * JHUXSttHkappatilde - JHUXSttHkappakappatilde ) / JHUXSttHkappa * 100)
    print "    2*kappaXS - kappakappatildeXS                             = {}%".format((2*JHUXSttHkappa - JHUXSttHkappakappatilde                                         ) / (2*JHUXSttHkappa) * 100)

#Set them to exactly 0

#decay

photoncutsame = (JHUXSggH2L2la1_photoncut == JHUXSggH2L2la1)
assert photoncutsame
values = [JHUXSggH2L2la1, JHUXSggH2L2la2*g2decay**2, JHUXSggH2L2la3*g4decay**2, JHUXSggH2L2lL1*g1prime2decay_gen**2, JHUXSggH2L2lL1Zg*ghzgs1prime2decay_gen**2]
errors = [JHUXSggH2L2la1err, JHUXSggH2L2la2err*g2decay**2, JHUXSggH2L2la3err*g4decay**2, JHUXSggH2L2lL1err*g1prime2decay_gen**2, JHUXSggH2L2lL1Zgerr*ghzgs1prime2decay_gen**2]

JHUXSggH2L2la1, JHUXSggH2L2la1err = sum(value/error**2 for value, error in zip(values, errors)) / sum(1/error**2 for error in errors), sum(1/error**2 for error in errors)**-.5
JHUXSggH2L2la2, JHUXSggH2L2la2err = JHUXSggH2L2la1 / g2decay**2, JHUXSggH2L2la1err / g2decay**2
JHUXSggH2L2la3, JHUXSggH2L2la3err = JHUXSggH2L2la1 / g4decay**2, JHUXSggH2L2la1err / g4decay**2
JHUXSggH2L2lL1, JHUXSggH2L2lL1err = JHUXSggH2L2la1 / g1prime2decay_gen**2, JHUXSggH2L2la1err / g1prime2decay_gen**2
JHUXSggH2L2lL1Zg, JHUXSggH2L2lL1Zgerr = JHUXSggH2L2la1 / ghzgs1prime2decay_gen**2, JHUXSggH2L2la1err / ghzgs1prime2decay_gen**2

JHUXSggH2L2la1a3, JHUXSggH2L2la1a3err = JHUXSggH2L2la1*2, JHUXSggH2L2la1err*2

if photoncutsame:
    JHUXSggH2L2la1_photoncut = JHUXSggH2L2la1
    JHUXSggH2L2la1err_photoncut = JHUXSggH2L2la1err

#VBF

photoncutsame = (JHUXSVBFa1_photoncut == JHUXSVBFa1)
assert not photoncutsame

values = [JHUXSVBFa1, JHUXSVBFa2*g2VBF**2, JHUXSVBFa3*g4VBF**2, JHUXSVBFL1*g1prime2VBF_gen**2]
errors = [JHUXSVBFa1err, JHUXSVBFa2err*g2VBF**2, JHUXSVBFa3err*g4VBF**2, JHUXSVBFL1err*g1prime2VBF_gen**2]

JHUXSVBFa1, JHUXSVBFa1err = sum(value/error**2 for value, error in zip(values, errors)) / sum(1/error**2 for error in errors), sum(1/error**2 for error in errors)**-.5
JHUXSVBFa2, JHUXSVBFa2err = JHUXSVBFa1 / g2VBF**2, JHUXSVBFa1err / g2VBF**2
JHUXSVBFa3, JHUXSVBFa3err = JHUXSVBFa1 / g4VBF**2, JHUXSVBFa1err / g4VBF**2
JHUXSVBFL1, JHUXSVBFL1err = JHUXSVBFa1 / g1prime2VBF_gen**2, JHUXSVBFa1err / g1prime2VBF_gen**2

JHUXSVBFa1a3, JHUXSVBFa1a3err = JHUXSVBFa1*2, JHUXSVBFa1err*2

values = [JHUXSVBFa1_photoncut, JHUXSVBFL1Zg*ghzgs1prime2VBF_gen**2]
errors = [JHUXSVBFa1err_photoncut, JHUXSVBFL1Zgerr*ghzgs1prime2VBF_gen**2]

JHUXSVBFa1_photoncut, JHUXSVBFa1err_photonuct = sum(value/error**2 for value, error in zip(values, errors)) / sum(1/error**2 for error in errors), sum(1/error**2 for error in errors)**-.5
JHUXSVBFL1Zg, JHUXSVBFL1Zgerr = JHUXSVBFa1_photoncut / ghzgs1prime2VBF_gen**2, JHUXSVBFa1err_photoncut / ghzgs1prime2VBF_gen**2

if photoncutsame:
    JHUXSVBFa1_photoncut = JHUXSVBFa1
    JHUXSVBFa1err_photoncut = JHUXSVBFa1err

#ZH

photoncutsame = (JHUXSZHa1_photoncut == JHUXSZHa1)
assert photoncutsame

values = [JHUXSZHa1, JHUXSZHa2*g2ZH**2, JHUXSZHa3*g4ZH**2, JHUXSZHL1*g1prime2ZH_gen**2, JHUXSZHL1Zg*ghzgs1prime2ZH_gen**2]
errors = [JHUXSZHa1err, JHUXSZHa2err*g2ZH**2, JHUXSZHa3err*g4ZH**2, JHUXSZHL1Zgerr*g1prime2ZH_gen**2, JHUXSZHL1Zgerr*ghzgs1prime2ZH_gen**2]

JHUXSZHa1, JHUXSZHa1err = sum(value/error**2 for value, error in zip(values, errors)) / sum(1/error**2 for error in errors), sum(1/error**2 for error in errors)**-.5
JHUXSZHa2, JHUXSZHa2err = JHUXSZHa1 / g2ZH**2, JHUXSZHa1err / g2ZH**2
JHUXSZHa3, JHUXSZHa3err = JHUXSZHa1 / g4ZH**2, JHUXSZHa1err / g4ZH**2
JHUXSZHL1, JHUXSZHL1err = JHUXSZHa1 / g1prime2ZH_gen**2, JHUXSZHa1err / g1prime2ZH_gen**2
JHUXSZHL1Zg, JHUXSZHL1Zgerr = JHUXSZHa1 / ghzgs1prime2ZH_gen**2, JHUXSZHa1err / ghzgs1prime2ZH_gen**2

JHUXSZHa1a3, JHUXSZHa1a3err = JHUXSZHa1*2, JHUXSZHa1err*2

if photoncutsame:
    JHUXSZHa1_photoncut = JHUXSZHa1
    JHUXSZHa1err_photoncut = JHUXSZHa1err
del photoncutsame

#WH

values = [JHUXSWHa1, JHUXSWHa2*g2WH**2, JHUXSWHa3*g4WH**2, JHUXSWHL1*g1prime2WH_gen**2]
errors = [JHUXSWHa1err, JHUXSWHa2err*g2WH**2, JHUXSWHa3err*g4WH**2, JHUXSWHL1err*g1prime2WH_gen**2]

JHUXSWHa1, JHUXSWHa1err = sum(value/error**2 for value, error in zip(values, errors)) / sum(1/error**2 for error in errors), sum(1/error**2 for error in errors)**-.5
JHUXSWHa2, JHUXSWHa2err = JHUXSWHa1 / g2WH**2, JHUXSWHa1err / g2WH**2
JHUXSWHa3, JHUXSWHa3err = JHUXSWHa1 / g4WH**2, JHUXSWHa1err / g4WH**2
JHUXSWHL1, JHUXSWHL1err = JHUXSWHa1 / g1prime2WH_gen**2, JHUXSWHa1err / g1prime2WH_gen**2
JHUXSWHL1Zg, JHUXSWHL1Zgerr = 0, 0

JHUXSWHa1a3, JHUXSWHa1a3err = JHUXSWHa1*2, JHUXSWHa1err*2

#HJJ

values = [JHUXSHJJa2, JHUXSHJJa3*ghg4HJJ**2]
errors = [JHUXSHJJa2err, JHUXSHJJa3err*ghg4HJJ**2]

JHUXSHJJa2, JHUXSHJJa2err = sum(value/error**2 for value, error in zip(values, errors)) / sum(1/error**2 for error in errors), sum(1/error**2 for error in errors)**-.5
JHUXSHJJa3, JHUXSHJJa3err = JHUXSHJJa2 / ghg4HJJ**2, JHUXSHJJa2err / ghg4HJJ**2

JHUXSHJJa2a3, JHUXSHJJa2a3err = JHUXSHJJa2*2, JHUXSHJJa2err*2

del values, errors

#ttH

values = [JHUXSttHkappa, JHUXSttHkappatilde*kappa_tilde_ttH**2]
errors = [JHUXSttHkappaerr, JHUXSttHkappatildeerr*kappa_tilde_ttH**2]

JHUXSttHkappa, JHUXSttHkappaerr = sum(value/error**2 for value, error in zip(values, errors)) / sum(1/error**2 for error in errors), sum(1/error**2 for error in errors)**-.5
JHUXSttHkappatilde, JHUXSttHkappatildeerr = JHUXSttHkappa / kappa_tilde_ttH**2, JHUXSttHkappaerr / kappa_tilde_ttH**2

JHUXSttHkappakappatilde, JHUXSttHkappakappatildeerr = JHUXSttHkappa*2, JHUXSttHkappaerr*2

del values, errors

#define interference xsecs instead of mixture xsecs

JHUXSggH2L2la1a2,   JHUXSggH2L2la1a2err   = JHUXSggH2L2la1a2   - 2*JHUXSggH2L2la1, sqrt(JHUXSggH2L2la1a2err  **2 + 4*JHUXSggH2L2la1err**2)
JHUXSggH2L2la1a3,   JHUXSggH2L2la1a3err   = JHUXSggH2L2la1a3   - 2*JHUXSggH2L2la1, 0
JHUXSggH2L2la1L1,   JHUXSggH2L2la1L1err   = JHUXSggH2L2la1L1   - 2*JHUXSggH2L2la1, sqrt(JHUXSggH2L2la1L1err  **2 + 4*JHUXSggH2L2la1err**2)
JHUXSggH2L2la1L1Zg, JHUXSggH2L2la1L1Zgerr = JHUXSggH2L2la1L1Zg - 2*JHUXSggH2L2la1, sqrt(JHUXSggH2L2la1L1Zgerr**2 + 4*JHUXSggH2L2la1err**2)

JHUXSVBFa1a2,   JHUXSVBFa1a2err   = JHUXSVBFa1a2   - 2*JHUXSVBFa1, sqrt(JHUXSVBFa1a2err  **2 + 4*JHUXSVBFa1err**2)
JHUXSVBFa1a3,   JHUXSVBFa1a3err   = JHUXSVBFa1a3   - 2*JHUXSVBFa1, 0
JHUXSVBFa1L1,   JHUXSVBFa1L1err   = JHUXSVBFa1L1   - 2*JHUXSVBFa1, sqrt(JHUXSVBFa1L1err  **2 + 4*JHUXSVBFa1err**2)
JHUXSVBFa1L1Zg, JHUXSVBFa1L1Zgerr = JHUXSVBFa1L1Zg - 2*JHUXSVBFa1, sqrt(JHUXSVBFa1L1Zgerr**2 + 4*JHUXSVBFa1err**2)

JHUXSZHa1a2,   JHUXSZHa1a2err   = JHUXSZHa1a2   - 2*JHUXSZHa1, sqrt(JHUXSZHa1a2err  **2 + 4*JHUXSZHa1err**2)
JHUXSZHa1a3,   JHUXSZHa1a3err   = JHUXSZHa1a3   - 2*JHUXSZHa1, 0
JHUXSZHa1L1,   JHUXSZHa1L1err   = JHUXSZHa1L1   - 2*JHUXSZHa1, sqrt(JHUXSZHa1L1err  **2 + 4*JHUXSZHa1err**2)
JHUXSZHa1L1Zg, JHUXSZHa1L1Zgerr = JHUXSZHa1L1Zg - 2*JHUXSZHa1, sqrt(JHUXSZHa1L1Zgerr**2 + 4*JHUXSZHa1err**2)

JHUXSWHa1a2,   JHUXSWHa1a2err   = JHUXSWHa1a2   - 2*JHUXSWHa1, sqrt(JHUXSWHa1a2err  **2 + 4*JHUXSWHa1err**2)
JHUXSWHa1a3,   JHUXSWHa1a3err   = JHUXSWHa1a3   - 2*JHUXSWHa1, 0
JHUXSWHa1L1,   JHUXSWHa1L1err   = JHUXSWHa1L1   - 2*JHUXSWHa1, sqrt(JHUXSWHa1L1err  **2 + 4*JHUXSWHa1err**2)
JHUXSWHa1L1Zg, JHUXSWHa1L1Zgerr = 0, 0

JHUXSHJJa2a3, JHUXSHJJa2a3err = JHUXSHJJa2a3 - 2*JHUXSHJJa2, 0
JHUXSttHkappakappatilde, JHUXSttHkappakappatildeerr = JHUXSttHkappakappatilde - 2*JHUXSttHkappa, 0

g2VH = sqrt((JHUXSZHa1 + JHUXSWHa1) / (JHUXSZHa2 + JHUXSWHa2))
g4VH = sqrt((JHUXSZHa1 + JHUXSWHa1) / (JHUXSZHa3 + JHUXSWHa3))
g1prime2VH_gen = -sqrt((JHUXSZHa1 + JHUXSWHa1) / (JHUXSZHL1 + JHUXSWHL1))
g1prime2VH_reco = -g1prime2VH_gen
ghzgs1prime2VH_gen = -sqrt((JHUXSZHa1 + JHUXSWHa1) / (JHUXSZHL1Zg + JHUXSWHL1Zg))
ghzgs1prime2VH_reco = -ghzgs1prime2VH_gen

#defined this way, just make sure
assert (JHUXSggH2L2la1a3 == JHUXSVBFa1a3 == JHUXSZHa1a3 == JHUXSWHa1a3 == JHUXSHJJa2a3 == JHUXSttHkappakappatilde
              == JHUXSWHL1Zg == JHUXSWHa1L1Zg == 0)
for k, v in globals().items():
    if "__" in k: continue
    assert v is not None, k
del k, v
