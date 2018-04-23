#!/usr/bin/env python
"""
Conventions for the constants here:

The top few constants are directly from JHUGen mod_Parameters.

g2HZZ and similar: these are from
https://twiki.cern.ch/twiki/bin/view/CMS/Run2MCProductionforHiggsProperties
They are sqrt(sigma1/sigmai) for the given process,
calculated by JHUGen in 2015. They were used for all MC samples,
and are also used to calculate fai and discriminants.
fai is particularly important, since these numbers are published
in HIG-17-011 (at least for decay; for production they appear in
published plots for fa3, and supplemental plots for the others).
They are basically a convention at this point and should not be changed.

SMXS are from YR4.  It's documented exactly how to
get them from the excel spreadsheet.

JHUXS - these are measured from JHUGen.  Currently, if you take
sqrt(sigma1/sigmai) you should get something reasonably
consistent with gidecay (or giVBF or whatever), but that's
not guaranteed to be the case.  In particular, when we
change PDFs that will no longer be true for production cross
sections.

The JHUXS for SM are calculated for g1=1, and the pure anomalous
JHUXS are calculated for gi=1.  For the mixture xsecs,
READ CAREFULLY:

The numbers written explicitly in this file are calculated with
g1=1, gi=gidecay or giVBF or ....  This is to maximize statistics.
(Similarly, for the mixtures between two anomalous couplings, they
are calculated with gi=gi, gj=gj)
Later in this file, they are processed further: first, we subtract
(JHUXSa1 + gi**2 * JHUXSai) [or similar for the aiaj xsecs].
We are then left with pure inteference.

The numbers are then divided by gi (or gi*gj).

At this point, if __name__ == "__main__" we print some things
that are expected to be ~0.  Check these if you make any changes
in constants.  They include sigma1 - sigmai*gi**2, which should
be removed when we change pdfs, as mentioned above.

Then, some of these are defined to be exactly 0.  For example
any interference between scalar and pseudoscalar, without cuts.

This is what you get when you import these constants and use them
in python.
"""

try:
  import uncertainties
except ImportError:
  raise ImportError("Need to install uncertainties!")

from uncertainties import ufloat
from uncertainties.umath import sqrt

M_Z = 91.1876
Ga_Z = 2.4952
aL = -0.53762
aR = 0.46238
e = 0.8431872482432357  # = cL_lep = cR_lep from mod_Parameters
L1 = 10000.

g2HZZ = 1.65684
g4HZZ = 2.55052
g1prime2HZZ = -12100.42   #for the sample
ghzgs1prime2HZZ_gen = -7613.351302119843
eLHZZ = sqrt(7.2310297E+00 / 1.4347981E+01)
eRHZZ = sqrt(7.2310297E+00 / 1.3952140E+00)

g2VBF = 0.27196538
g4VBF = 0.297979018705
g1prime2VBF = -2158.21307286
ghzgs1prime2VBF_gen = -4091.051456694223

g2ZH = 0.112481
g4ZH = 0.144057
g1prime2ZH = -517.788
ghzgs1prime2ZH_gen = -642.9534550379002

g2WH = 0.0998956
g4WH = 0.1236136
g1prime2WH = -525.274
ghzgs1prime2WH_gen = -1000

g2VH = 0.10430356645812816 #sqrt((JHUXSZHa1 + JHUXSWHa1*normalize_WH_to_ZH) / (JHUXSZHa2 + JHUXSWHa2*normalize_WH_to_ZH))
g4VH = 0.13053750671388425 #sqrt((JHUXSZHa1 + JHUXSWHa1*normalize_WH_to_ZH) / (JHUXSZHa3 + JHUXSWHa3*normalize_WH_to_ZH))
g1prime2VH = -522.3034453633128 #-sqrt((JHUXSZHa1 + JHUXSWHa1*normalize_WH_to_ZH) / (JHUXSZHL1 + JHUXSWHL1*normalize_WH_to_ZH))
ghzgs1prime2VH_gen = -1027.387141119873 #-sqrt((JHUXSZHa1 + JHUXSWHa1*normalize_WH_to_ZH) / (JHUXSZHL1Zg + JHUXSWHL1Zg*normalize_WH_to_ZH))
nominal_normalize_WH_to_ZH = 0.15070409765374365

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

JHUXSggH2L2la1             = ufloat(7.1517173,      0.23044650E-03)
JHUXSggH2L2la2             = ufloat(2.5849908,      0.77294379E-04)
JHUXSggH2L2la3             = ufloat(1.0954288,      0.43443941E-04)
JHUXSggH2L2lL1             = ufloat(0.48754258E-07, 0.15112345E-11)
JHUXSggH2L2lL1Zg           = ufloat(0.12338393E-06, 0.58002047E-11)
JHUXSggH2L2la1a2           = ufloat(26.073377,      0.74026916E-03)
JHUXSggH2L2la1a3           = ufloat(14.278034,      0.45318122E-03)
JHUXSggH2L2la1L1           = ufloat(0.23727827,     0.15982277E-04)    #using g1prime2 = -12100.42
JHUXSggH2L2la1L1Zg         = ufloat(13.145473,      0.53550812E-03)
JHUXSggH2L2la2a3           = ufloat(14.3016,        0.000642865)
JHUXSggH2L2la2L1           = ufloat(2.33083,        0.000141056)
JHUXSggH2L2la2L1Zg         = ufloat(2.33096,        0.000140816)
JHUXSggH2L2la3L1           = ufloat(14.3016,        0.000521711)
JHUXSggH2L2la3L1Zg         = ufloat(14.3019,        0.000521397)
JHUXSggH2L2lL1L1Zg         = ufloat(1.6472524E+01,  2.2689281E-02)

JHUXSggH2L2leL = 1.4347981E+01
JHUXSggH2L2leR = 1.3952140E+01

JHUXSVBFa1                 = ufloat(968.674284006,      0.075115702763)
JHUXSVBFa2                 = ufloat(13102.7106117,      0.522399748272)
JHUXSVBFa3                 = ufloat(10909.5390002,      0.50975030067)
JHUXSVBFL1                 = ufloat(2.083097999e-4,     1.24640942579e-08)
JHUXSVBFL1Zg               = ufloat(0.49845301E-04,     0.21806623E-07)
JHUXSVBFa1a2               = ufloat(2207.72848655,      0.126379327428)
JHUXSVBFa1a3               = ufloat(1937.20646111,      0.122617320785)
JHUXSVBFa1L1               = ufloat(2861.21349769,      0.0771278408768)
JHUXSVBFa1L1Zg             = ufloat(1410.5494,          0.68216715)
JHUXSVBFa2a3               = ufloat(1936.8416,          0.55683416)
JHUXSVBFa2L1               = ufloat(2507.0486,          0.79805153)
JHUXSVBFa2L1Zg             = ufloat(2433.0553,          0.78466670)
JHUXSVBFa3L1               = ufloat(1939.6777,          0.72505911)
JHUXSVBFa3L1Zg             = ufloat(1865.6689,          0.71037371)
JHUXSVBFL1L1Zg             = ufloat(916.21310,          0.34424016)

JHUXSZHa1                  = ufloat(9022.36,        1.17)
JHUXSZHa2                  = ufloat(713123,         103)
JHUXSZHa3                  = ufloat(434763.7,       62.2)
JHUXSZHL1                  = ufloat(33652.46e-6,    4.19e-6)
JHUXSZHL1Zg                = ufloat(3.4592597,      0.56602360E-03)
JHUXSZHa1a2                = ufloat(4258.966,       0.783)
JHUXSZHa1a3                = ufloat(18040.66,       2.73)
JHUXSZHa1L1                = ufloat(6852.307,       0.929)
JHUXSZHa1L1Zg              = ufloat(3479412.2,      486.4)
JHUXSZHa2a3                = ufloat(2.87412e+06,    166.887)
JHUXSZHa2L1                = ufloat(4.48796e+06,    215.167)
JHUXSZHa2L1Zg              = ufloat(4.48622e+06,    1502.75)
JHUXSZHa3L1                = ufloat(2.87379e+06,    207.548)
JHUXSZHa3L1Zg              = ufloat(2.8737e+06,     614.514)
JHUXSZHL1L1Zg              = ufloat(1.43739e+06,    239.856)

JHUXSWHa1                  = ufloat(30998.54,       2.50)
JHUXSWHa2                  = ufloat(3106339,        308)
JHUXSWHa3                  = ufloat(2028656,        191)
JHUXSWHL1                  = ufloat(11234.91e-5,    1.10e-5)
JHUXSWHL1Zg                = ufloat(0,              0)
JHUXSWHa1a2                = ufloat(16486.68,       2.01)
JHUXSWHa1a3                = ufloat(62001.57,       5.54)
JHUXSWHa1L1                = ufloat(25302.37,       2.67)
JHUXSWHa2a3                = ufloat(2.96308e+07,    2132.86)
JHUXSWHa2L1                = ufloat(4.47663e+07,    4985.14)
JHUXSWHa3L1                = ufloat(2.96285e+07,    3423.37)

#VH cross sections changed somewhere between c85a387eaf3a8a92f5893e5293ed3c3d36107e16 and fbf449150f4df49f66b21c1638adb02b68a308d0
#correct for that
def fixVH():
    ZHxsecratio = ufloat(158.49737883504775, 0.39826108447619935)
    WHxsecratio = ZHxsecratio*3
    for V in "Z", "W":
        VHxsecratio = locals()[V+"Hxsecratio"]
        for coupling in "a1", "a2", "a3", "L1", "a1a2", "a1a3", "a1L1":
            name = "JHUXS{}H{}".format(V, coupling)
            newvalue = globals()[name] * VHxsecratio
            globals()[name] = newvalue
fixVH(); del fixVH
JHUXSWHa1L1Zg       = JHUXSWHa1
JHUXSWHa2L1Zg       = g2WH**2 * JHUXSWHa2
JHUXSWHa3L1Zg       = g4WH**2 * JHUXSWHa3
JHUXSWHL1L1Zg       = g1prime2WH**2 * JHUXSWHL1

JHUXSHJJa2       = ufloat(14583.61,       0.94)
JHUXSHJJa3       = ufloat(14397.13,       0.97)
JHUXSHJJa2a3     = ufloat(29169.2,        2.1)

JHUXSttHkappa    = ufloat(0.912135589,    0.00143032)
JHUXSttHkappatilde = ufloat(0.35609194,   0.000492662)
JHUXSttHkappakappatilde = ufloat(1.8231162489,   0.00254131)


#Subtract the pure component from the interference, then divide by (gi*gj)
JHUXSggH2L2la1a2   = (JHUXSggH2L2la1a2   -                        JHUXSggH2L2la1 - g2HZZ              **2 * JHUXSggH2L2la2  ) / (g2HZZ                                      )
JHUXSggH2L2la1a3   = (JHUXSggH2L2la1a3   -                        JHUXSggH2L2la1 - g4HZZ              **2 * JHUXSggH2L2la3  ) / (g4HZZ                                      )
JHUXSggH2L2la1L1   = (JHUXSggH2L2la1L1   -                        JHUXSggH2L2la1 - g1prime2HZZ    **2 * JHUXSggH2L2lL1  ) / (g1prime2HZZ                            )
JHUXSggH2L2la1L1Zg = (JHUXSggH2L2la1L1Zg -                        JHUXSggH2L2la1 - ghzgs1prime2HZZ_gen**2 * JHUXSggH2L2lL1Zg) / (ghzgs1prime2HZZ_gen                        )
JHUXSggH2L2la2a3   = (JHUXSggH2L2la2a3   - g2HZZ          **2 * JHUXSggH2L2la2 - g4HZZ              **2 * JHUXSggH2L2la3  ) / (g2HZZ               * g4HZZ              )
JHUXSggH2L2la2L1   = (JHUXSggH2L2la2L1   - g2HZZ          **2 * JHUXSggH2L2la2 - g1prime2HZZ    **2 * JHUXSggH2L2lL1  ) / (g2HZZ               * g1prime2HZZ    )
JHUXSggH2L2la2L1Zg = (JHUXSggH2L2la2L1Zg - g2HZZ          **2 * JHUXSggH2L2la2 - ghzgs1prime2HZZ_gen**2 * JHUXSggH2L2lL1Zg) / (g2HZZ               * ghzgs1prime2HZZ_gen)
JHUXSggH2L2la3L1   = (JHUXSggH2L2la3L1   - g4HZZ          **2 * JHUXSggH2L2la3 - g1prime2HZZ    **2 * JHUXSggH2L2lL1  ) / (g4HZZ               * g1prime2HZZ    )
JHUXSggH2L2la3L1Zg = (JHUXSggH2L2la3L1Zg - g4HZZ          **2 * JHUXSggH2L2la3 - ghzgs1prime2HZZ_gen**2 * JHUXSggH2L2lL1Zg) / (g4HZZ               * ghzgs1prime2HZZ_gen)
JHUXSggH2L2lL1L1Zg = (JHUXSggH2L2lL1L1Zg - g1prime2HZZ**2 * JHUXSggH2L2lL1 - ghzgs1prime2HZZ_gen**2 * JHUXSggH2L2lL1Zg) / (g1prime2HZZ     * ghzgs1prime2HZZ_gen)

JHUXSVBFa1a2       = (JHUXSVBFa1a2       -                        JHUXSVBFa1     - g2VBF                **2 * JHUXSVBFa2      ) / (g2VBF                                        )
JHUXSVBFa1a3       = (JHUXSVBFa1a3       -                        JHUXSVBFa1     - g4VBF                **2 * JHUXSVBFa3      ) / (g4VBF                                        )
JHUXSVBFa1L1       = (JHUXSVBFa1L1       -                        JHUXSVBFa1     - g1prime2VBF      **2 * JHUXSVBFL1      ) / (g1prime2VBF                              )
JHUXSVBFa1L1Zg     = (JHUXSVBFa1L1Zg     -                        JHUXSVBFa1     - ghzgs1prime2VBF_gen  **2 * JHUXSVBFL1Zg    ) / (ghzgs1prime2VBF_gen                          )
JHUXSVBFa2a3       = (JHUXSVBFa2a3       - g2VBF            **2 * JHUXSVBFa2     - g4VBF                **2 * JHUXSVBFa3      ) / (g2VBF                 * g4VBF                )
JHUXSVBFa2L1       = (JHUXSVBFa2L1       - g2VBF            **2 * JHUXSVBFa2     - g1prime2VBF      **2 * JHUXSVBFL1      ) / (g2VBF                 * g1prime2VBF      )
JHUXSVBFa2L1Zg     = (JHUXSVBFa2L1Zg     - g2VBF            **2 * JHUXSVBFa2     - ghzgs1prime2VBF_gen  **2 * JHUXSVBFL1Zg    ) / (g2VBF                 * ghzgs1prime2VBF_gen  )
JHUXSVBFa3L1       = (JHUXSVBFa3L1       - g4VBF            **2 * JHUXSVBFa3     - g1prime2VBF      **2 * JHUXSVBFL1      ) / (g4VBF                 * g1prime2VBF      )
JHUXSVBFa3L1Zg     = (JHUXSVBFa3L1Zg     - g4VBF            **2 * JHUXSVBFa3     - ghzgs1prime2VBF_gen  **2 * JHUXSVBFL1Zg    ) / (g4VBF                 * ghzgs1prime2VBF_gen  )
JHUXSVBFL1L1Zg     = (JHUXSVBFL1L1Zg     - g1prime2VBF  **2 * JHUXSVBFL1     - ghzgs1prime2VBF_gen  **2 * JHUXSVBFL1Zg    ) / (g1prime2VBF       * ghzgs1prime2VBF_gen  )

JHUXSZHa1a2        = (JHUXSZHa1a2        -                        JHUXSZHa1      - g2ZH                 **2 * JHUXSZHa2       ) / (g2ZH                                         )
JHUXSZHa1a3        = (JHUXSZHa1a3        -                        JHUXSZHa1      - g4ZH                 **2 * JHUXSZHa3       ) / (g4ZH                                         )
JHUXSZHa1L1        = (JHUXSZHa1L1        -                        JHUXSZHa1      - g1prime2ZH       **2 * JHUXSZHL1       ) / (g1prime2ZH                               )
JHUXSZHa1L1Zg      = (JHUXSZHa1L1Zg      -                        JHUXSZHa1      - ghzgs1prime2ZH_gen   **2 * JHUXSZHL1Zg     ) / (ghzgs1prime2ZH_gen                           )
JHUXSZHa2a3        = (JHUXSZHa2a3        - g2ZH             **2 * JHUXSZHa2      - g4ZH                 **2 * JHUXSZHa3       ) / (g2ZH                  * g4ZH                 )
JHUXSZHa2L1        = (JHUXSZHa2L1        - g2ZH             **2 * JHUXSZHa2      - g1prime2ZH       **2 * JHUXSZHL1       ) / (g2ZH                  * g1prime2ZH       )
JHUXSZHa2L1Zg      = (JHUXSZHa2L1Zg      - g2ZH             **2 * JHUXSZHa2      - ghzgs1prime2ZH_gen   **2 * JHUXSZHL1Zg     ) / (g2ZH                  * ghzgs1prime2ZH_gen   )
JHUXSZHa3L1        = (JHUXSZHa3L1        - g4ZH             **2 * JHUXSZHa3      - g1prime2ZH       **2 * JHUXSZHL1       ) / (g4ZH                  * g1prime2ZH       )
JHUXSZHa3L1Zg      = (JHUXSZHa3L1Zg      - g4ZH             **2 * JHUXSZHa3      - ghzgs1prime2ZH_gen   **2 * JHUXSZHL1Zg     ) / (g4ZH                  * ghzgs1prime2ZH_gen   )
JHUXSZHL1L1Zg      = (JHUXSZHL1L1Zg      - g1prime2ZH   **2 * JHUXSZHL1      - ghzgs1prime2ZH_gen   **2 * JHUXSZHL1Zg     ) / (g1prime2ZH        * ghzgs1prime2ZH_gen   )

JHUXSWHa1a2        = (JHUXSWHa1a2        -                        JHUXSWHa1      - g2WH                 **2 * JHUXSWHa2       ) / (g2WH                                         )
JHUXSWHa1a3        = (JHUXSWHa1a3        -                        JHUXSWHa1      - g4WH                 **2 * JHUXSWHa3       ) / (g4WH                                         )
JHUXSWHa1L1        = (JHUXSWHa1L1        -                        JHUXSWHa1      - g1prime2WH       **2 * JHUXSWHL1       ) / (g1prime2WH                               )
JHUXSWHa1L1Zg      = 0
JHUXSWHa2a3        = (JHUXSWHa2a3        - g2WH             **2 * JHUXSWHa2      - g4WH                 **2 * JHUXSWHa3       ) / (g2WH                  * g4WH                 )
JHUXSWHa2L1        = (JHUXSWHa2L1        - g2WH             **2 * JHUXSWHa2      - g1prime2WH       **2 * JHUXSWHL1       ) / (g2WH                  * g1prime2WH       )
JHUXSWHa2L1Zg      = 0
JHUXSWHa3L1        = (JHUXSWHa3L1        - g4WH             **2 * JHUXSWHa3      - g1prime2WH       **2 * JHUXSWHL1       ) / (g4WH                  * g1prime2WH       )
JHUXSWHa3L1Zg      = 0
JHUXSWHL1L1Zg      = 0

JHUXSHJJa2a3       = (JHUXSHJJa2a3       -                        JHUXSHJJa2     - ghg4HJJ              **2 * JHUXSHJJa3      ) / (ghg4HJJ                                      )

JHUXSttHkappakappatilde = (JHUXSttHkappakappatilde - JHUXSttHkappa - kappa_tilde_ttH**2 * JHUXSttHkappatilde) / kappa_tilde_ttH

normalize_WH_to_ZH = SMXSWH / JHUXSWHa1 / (SMXSZH / JHUXSZHa1)

if __name__ == "__main__":
    print "All of the following should be 0:"
    print
    print "  decay:"
    print "    a1XS -           g2**2*    a2XS = {:%}".format((JHUXSggH2L2la1 - g2HZZ**2               * JHUXSggH2L2la2    ) / JHUXSggH2L2la1)
    print "    a1XS -           g4**2*    a3XS = {:%}".format((JHUXSggH2L2la1 - g4HZZ**2               * JHUXSggH2L2la3    ) / JHUXSggH2L2la1)
    print "    a1XS -     g1prime2**2*    L1XS = {:%}".format((JHUXSggH2L2la1 - g1prime2HZZ**2     * JHUXSggH2L2lL1    ) / JHUXSggH2L2la1)
    print "    a1XS - ghzgs1prime2**2*  L1ZgXS = {:%}".format((JHUXSggH2L2la1 - ghzgs1prime2HZZ_gen**2 * JHUXSggH2L2lL1Zg  ) / JHUXSggH2L2la1)
    print "                        g4*  a1a3XS = {:%}".format((                 g4HZZ                  * JHUXSggH2L2la1a3  ) / JHUXSggH2L2la1)
    print "                     g2*g4*  a2a3XS = {:%}".format((                 g2HZZ*g4HZZ          * JHUXSggH2L2la2a3  ) / JHUXSggH2L2la1)
    print "               g1prime2*g4*  a3L1XS = {:%}".format((                g1prime2HZZ*g4HZZ * JHUXSggH2L2la3L1  ) / JHUXSggH2L2la1)
    print "           ghzgs1prime2*g4*a3L1ZgXS = {:%}".format((            ghzgs1prime2HZZ_gen*g4HZZ * JHUXSggH2L2la3L1Zg) / JHUXSggH2L2la1)
    print
    print "  VBF:"
    print "    a1XS -           g2**2*    a2XS = {:%}".format((JHUXSVBFa1     - g2VBF**2                 * JHUXSVBFa2        ) / JHUXSVBFa1    )
    print "    a1XS -           g4**2*    a3XS = {:%}".format((JHUXSVBFa1     - g4VBF**2                 * JHUXSVBFa3        ) / JHUXSVBFa1    )
    print "    a1XS -     g1prime2**2*    L1XS = {:%}".format((JHUXSVBFa1     - g1prime2VBF**2       * JHUXSVBFL1        ) / JHUXSVBFa1    )
#   print "    a1XS - ghzgs1prime2**2*  L1ZgXS = {:%}".format((JHUXSVBFa1     - ghzgs1prime2VBF_gen**2   * JHUXSVBFL1Zg      ) / JHUXSVBFa1    )
    print "                        g4*  a1a3XS = {:%}".format((                 g4VBF                    * JHUXSVBFa1a3      ) / JHUXSVBFa1    )
    print "                     g2*g4*  a2a3XS = {:%}".format((                 g2VBF  *g4VBF            * JHUXSVBFa2a3      ) / JHUXSVBFa1    )
    print "               g1prime2*g4*  a3L1XS = {:%}".format((                g1prime2VBF  *g4VBF   * JHUXSVBFa3L1      ) / JHUXSVBFa1    )
#   print "           ghzgs1prime2*g4*a3L1ZgXS = {:%}".format((            ghzgs1prime2VBF_gen  *g4VBF   * JHUXSVBFa3L1Zg    ) / JHUXSVBFa1    )
    print
    print "  ZH:"
    print "    a1XS -           g2**2*    a2XS = {:%}".format((JHUXSZHa1      - g2ZH**2                  * JHUXSZHa2         ) / JHUXSZHa1     )
    print "    a1XS -           g4**2*    a3XS = {:%}".format((JHUXSZHa1      - g4ZH**2                  * JHUXSZHa3         ) / JHUXSZHa1     )
    print "    a1XS -     g1prime2**2*    L1XS = {:%}".format((JHUXSZHa1      - g1prime2ZH**2        * JHUXSZHL1         ) / JHUXSZHa1     )
    print "    a1XS - ghzgs1prime2**2*  L1ZgXS = {:%}".format((JHUXSZHa1      - ghzgs1prime2ZH_gen**2    * JHUXSZHL1Zg       ) / JHUXSZHa1     )
    print "                        g4*  a1a3XS = {:%}".format((                 g4ZH                     * JHUXSZHa1a3       ) / JHUXSZHa1     )
    print "                     g2*g4*  a2a3XS = {:%}".format((                 g2ZH   *g4ZH             * JHUXSZHa2a3       ) / JHUXSZHa1     )
    print "               g1prime2*g4*  a3L1XS = {:%}".format((                g1prime2ZH   *g4ZH    * JHUXSZHa3L1       ) / JHUXSZHa1     )
    print "           ghzgs1prime2*g4*a3L1ZgXS = {:%}".format((            ghzgs1prime2ZH_gen   *g4ZH    * JHUXSZHa3L1Zg     ) / JHUXSZHa1     )
    print
    print "  WH:"
    print "    a1XS -           g2**2*    a2XS = {:%}".format((JHUXSWHa1      - g2WH**2                  * JHUXSWHa2         ) / JHUXSWHa1     )
    print "    a1XS -           g4**2*    a3XS = {:%}".format((JHUXSWHa1      - g4WH**2                  * JHUXSWHa3         ) / JHUXSWHa1     )
    print "    a1XS -     g1prime2**2*    L1XS = {:%}".format((JHUXSWHa1      - g1prime2WH**2        * JHUXSWHL1         ) / JHUXSWHa1     )
#   print "    a1XS - ghzgs1prime2**2*  L1ZgXS = {:%}".format((JHUXSWHa1      - ghzgs1prime2WH_gen**2    * JHUXSWHL1Zg       ) / JHUXSWHa1     )
    print "                        g4*  a1a3XS = {:%}".format((                 g4WH                     * JHUXSWHa1a3       ) / JHUXSWHa1     )
    print "                     g2*g4*  a2a3XS = {:%}".format((                 g2WH   *g4WH             * JHUXSWHa2a3       ) / JHUXSWHa1     )
    print "               g1prime2*g4*  a3L1XS = {:%}".format((                g1prime2WH   *g4WH    * JHUXSWHa3L1       ) / JHUXSWHa1     )
#   print "           ghzgs1prime2*g4*a3L1ZgXS = {:%}".format((            ghzgs1prime2WH_gen   *g4WH    * JHUXSWHa3L1Zg     ) / JHUXSWHa1     )
    print "    WpHXS + WmHXS - WHXS            = {:%}".format((SMXSWpH2L2l + SMXSWmH2L2l - SMXSWH2L2l                        ) / SMXSWH2L2l    )
    print
    print "  HJJ:"
    print "    a2XS -           g4**2*    a3XS = {:%}".format((JHUXSHJJa2     - ghg4HJJ**2               * JHUXSHJJa3        ) / JHUXSHJJa2    )
    print "                        g4*  a2a3XS = {:%}".format((                 ghg4HJJ                  * JHUXSHJJa2a3      ) / JHUXSHJJa2    )
    print
    print "  ttH:"
    kt = kappa_tilde_ttH
    JHUXSttHk = JHUXSttHkappa
    JHUXSttHkt = JHUXSttHkappatilde
    JHUXSttHkkt = JHUXSttHkappakappatilde
    print "     kXS -           k~**2*    k~XS = {:%}".format((JHUXSttHk      - kt**2                    * JHUXSttHkt        ) / JHUXSttHk     )
    print "                        k~*   kk~XS = {:%}".format((                 kt                       * JHUXSttHkkt       ) / JHUXSttHk     )
    print
    print "  VH:"
    print "    a1XS -           g2**2*    a2XS  = {:%}".format((JHUXSZHa1 + JHUXSWHa1*normalize_WH_to_ZH - g2VH              **2 * (JHUXSZHa2   + JHUXSWHa2  *normalize_WH_to_ZH)) / (JHUXSZHa1 + JHUXSWHa1*normalize_WH_to_ZH))
    print "    a1XS -           g4**2*    a3XS  = {:%}".format((JHUXSZHa1 + JHUXSWHa1*normalize_WH_to_ZH - g4VH              **2 * (JHUXSZHa3   + JHUXSWHa3  *normalize_WH_to_ZH)) / (JHUXSZHa1 + JHUXSWHa1*normalize_WH_to_ZH))
    print "    a1XS -     g1prime2**2*    L1XS  = {:%}".format((JHUXSZHa1 + JHUXSWHa1*normalize_WH_to_ZH - g1prime2VH    **2 * (JHUXSZHL1   + JHUXSWHL1  *normalize_WH_to_ZH)) / (JHUXSZHa1 + JHUXSWHa1*normalize_WH_to_ZH))
    print "    a1XS - ghzgs1prime2**2*  L1ZgXS  = {:%}".format((JHUXSZHa1 + JHUXSWHa1*normalize_WH_to_ZH - ghzgs1prime2VH_gen**2 * (JHUXSZHL1Zg + JHUXSWHL1Zg*normalize_WH_to_ZH)) / (JHUXSZHa1 + JHUXSWHa1*normalize_WH_to_ZH))
    del kt, JHUXSttHk, JHUXSttHkt, JHUXSttHkkt

#Set them to exactly 0

JHUXSggH2L2la1a3   = \
JHUXSggH2L2la2a3   = \
JHUXSggH2L2la3L1   = \
JHUXSggH2L2la3L1Zg = 0

JHUXSVBFa1a3   = \
JHUXSVBFa2a3   = \
JHUXSVBFa3L1   = \
JHUXSVBFa3L1Zg = 0

JHUXSZHa1a3   = \
JHUXSZHa2a3   = \
JHUXSZHa3L1   = \
JHUXSZHa3L1Zg = 0

JHUXSWHa1a3   = \
JHUXSWHa2a3   = \
JHUXSWHa3L1   = \
JHUXSWHa3L1Zg = 0

JHUXSHJJa2a3 = 0

JHUXSttHkappakappatilde = 0

#defined this way, just make sure
for _ in """
  JHUXSggH2L2la1a3 JHUXSggH2L2la2a3 JHUXSggH2L2la3L1 JHUXSggH2L2la3L1Zg
  JHUXSVBFa1a3 JHUXSVBFa2a3 JHUXSVBFa3L1 JHUXSVBFa3L1Zg
  JHUXSZHa1a3 JHUXSZHa2a3 JHUXSZHa3L1 JHUXSZHa3L1Zg
  JHUXSWHa1a3 JHUXSWHa2a3 JHUXSWHa3L1
  JHUXSHJJa2a3 JHUXSttHkappakappatilde
  JHUXSWHL1Zg JHUXSWHa1L1Zg JHUXSWHa2L1Zg JHUXSWHa3L1Zg JHUXSWHL1L1Zg
""".split():
  assert globals()[_] == 0, (_, globals()[_])
for k, v in globals().items():
    if "__" in k: continue
    assert v is not None, k
del k, v, _
