#!/usr/bin/env python
"""
Conventions for the constants here:

The top few constants are directly from JHUGen mod_Parameters.

g2decay and similar: these are from
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

from math import sqrt

M_Z = 91.1876
Ga_Z = 2.4952
aL = -0.53762
aR = 0.46238
e = 0.8431872482432357  # = cL_lep = cR_lep from mod_Parameters
L1 = 10000.

g2decay = 1.663195
g4decay = 2.55502
g1prime2decay_gen = -12110.20   #for the sample
g1prime2decay_reco = 12110.20   #for discriminants
ghzgs1prime2decay_gen = -7613.351302119843
ghzgs1prime2decay_reco = 7613.351302119843
eLdecay = sqrt(7.2310297E+00 / 1.4347981E+01)
eRdecay = sqrt(7.2310297E+00 / 1.3952140E+00)

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

JHUXSggH2L2la1             = 7.1517173;      JHUXSggH2L2la1err             = 0.23044650E-03
JHUXSggH2L2la2             = 2.5849908;      JHUXSggH2L2la2err             = 0.77294379E-04
JHUXSggH2L2la3             = 1.0954288;      JHUXSggH2L2la3err             = 0.43443941E-04
JHUXSggH2L2lL1             = 0.48754258E-07; JHUXSggH2L2lL1err             = 0.15112345E-11
JHUXSggH2L2lL1Zg           = 0.12338393E-06; JHUXSggH2L2lL1Zgerr           = 0.58002047E-11
JHUXSggH2L2la1a2           = 26.073377;      JHUXSggH2L2la1a2err           = 0.74026916E-03
JHUXSggH2L2la1a3           = 14.278034;      JHUXSggH2L2la1a3err           = 0.45318122E-03
JHUXSggH2L2la1L1           = 0.23727827;     JHUXSggH2L2la1L1err           = 0.15982277E-04    #using g1prime2 = -12100.42
JHUXSggH2L2la1L1Zg         = 13.145473;      JHUXSggH2L2la1L1Zgerr         = 0.53550812E-03
JHUXSggH2L2la2a3           = 14.3016;        JHUXSggH2L2la2a3err           = 0.000642865
JHUXSggH2L2la2L1           = 2.33083;        JHUXSggH2L2la2L1err           = 0.000141056
JHUXSggH2L2la2L1Zg         = 2.33096;        JHUXSggH2L2la2L1Zgerr         = 0.000140816
JHUXSggH2L2la3L1           = 14.3016;        JHUXSggH2L2la3L1err           = 0.000521711
JHUXSggH2L2la3L1Zg         = 14.3019;        JHUXSggH2L2la3L1Zgerr         = 0.000521397
JHUXSggH2L2lL1L1Zg         = 1.6472524E+01;  JHUXSggH2L2lL1L1Zgerr         = 2.2689281E-02

JHUXSggH2L2leL = 1.4347981E+01
JHUXSggH2L2leR = 1.3952140E+01

JHUXSVBFa1             = 968.674284006;      JHUXSVBFa1err             = 0.075115702763
JHUXSVBFa2             = 13102.7106117;      JHUXSVBFa2err             = 0.522399748272
JHUXSVBFa3             = 10909.5390002;      JHUXSVBFa3err             = 0.50975030067
JHUXSVBFL1             = 2.083097999e-4;     JHUXSVBFL1err             = 1.24640942579e-08
JHUXSVBFL1Zg           = 0.49845301E-04;     JHUXSVBFL1Zgerr           = 0.21806623E-07
JHUXSVBFa1a2           = 2207.72848655;      JHUXSVBFa1a2err           = 0.126379327428
JHUXSVBFa1a3           = 1937.20646111;      JHUXSVBFa1a3err           = 0.122617320785
JHUXSVBFa1L1           = 2861.21349769;      JHUXSVBFa1L1err           = 0.0771278408768
JHUXSVBFa1L1Zg         = 1410.5494;          JHUXSVBFa1L1Zgerr         = 0.68216715
JHUXSVBFa2a3           = 1936.8416;          JHUXSVBFa2a3err           = 0.55683416
JHUXSVBFa2L1           = 2507.0486;          JHUXSVBFa2L1err           = 0.79805153
JHUXSVBFa2L1Zg         = 2433.0553;          JHUXSVBFa2L1Zgerr         = 0.78466670
JHUXSVBFa3L1           = 1939.6777;          JHUXSVBFa3L1err           = 0.72505911
JHUXSVBFa3L1Zg         = 1865.6689;          JHUXSVBFa3L1Zgerr         = 0.71037371
JHUXSVBFL1L1Zg         = 916.21310;          JHUXSVBFL1L1Zgerr         = 0.34424016

JHUXSZHa1             = 9022.36;        JHUXSZHa1err             = 1.17
JHUXSZHa2             = 713123;         JHUXSZHa2err             = 103
JHUXSZHa3             = 434763.7;       JHUXSZHa3err             = 62.2
JHUXSZHL1             = 33652.46e-6;    JHUXSZHL1err             = 4.19e-6
JHUXSZHL1Zg           = 3.4592597;      JHUXSZHL1Zgerr           = 0.56602360E-03
JHUXSZHa1a2           = 4258.966;       JHUXSZHa1a2err           = 0.783
JHUXSZHa1a3           = 18040.66;       JHUXSZHa1a3err           = 2.73
JHUXSZHa1L1           = 6852.307;       JHUXSZHa1L1err           = 0.929
JHUXSZHa1L1Zg         = 3479412.2;      JHUXSZHa1L1Zgerr         = 486.4
JHUXSZHa2a3           = 2.87412e+06;    JHUXSZHa2a3err           = 166.887
JHUXSZHa2L1           = 4.48796e+06;    JHUXSZHa2L1err           = 215.167
JHUXSZHa2L1Zg         = 4.48622e+06;    JHUXSZHa2L1Zgerr         = 1502.75
JHUXSZHa3L1           = 2.87379e+06;    JHUXSZHa3L1err           = 207.548
JHUXSZHa3L1Zg         = 2.8737e+06;     JHUXSZHa3L1Zgerr         = 614.514
JHUXSZHL1L1Zg         = 1.43739e+06;    JHUXSZHL1L1Zgerr         = 239.856

JHUXSWHa1           = 30998.54;       JHUXSWHa1err           = 2.50
JHUXSWHa2           = 3106339;        JHUXSWHa2err           = 308
JHUXSWHa3           = 2028656;        JHUXSWHa3err           = 191
JHUXSWHL1           = 11234.91e-5;    JHUXSWHL1err           = 1.10e-5
JHUXSWHL1Zg         = 0;              JHUXSWHL1Zgerr         = 0
JHUXSWHa1a2         = 16486.68;       JHUXSWHa1a2err         = 2.01
JHUXSWHa1a3         = 62001.57;       JHUXSWHa1a3err         = 5.54
JHUXSWHa1L1         = 25302.37;       JHUXSWHa1L1err         = 2.67
JHUXSWHa2a3         = 2.96308e+07;    JHUXSWHa2a3err         = 2132.86
JHUXSWHa2L1         = 4.47663e+07;    JHUXSWHa2L1err         = 4985.14
JHUXSWHa3L1         = 2.96285e+07;    JHUXSWHa3L1err         = 3423.37

#VH cross sections changed somewhere between c85a387eaf3a8a92f5893e5293ed3c3d36107e16 and fbf449150f4df49f66b21c1638adb02b68a308d0
#correct for that
def fixVH():
    ZHxsecratio = 158.49737883504775; ZHxsecratioerr = 0.39826108447619935
    WHxsecratio = ZHxsecratio*3;      WHxsecratioerr = ZHxsecratioerr*3
    for V in "Z", "W":
        VHxsecratio, VHxsecratioerr = locals()[V+"Hxsecratio"], locals()[V+"Hxsecratioerr"]
        for coupling in "a1", "a2", "a3", "L1", "a1a2", "a1a3", "a1L1":
            name = "JHUXS{}H{}".format(V, coupling); errorname = name+"err"
            newvalue = globals()[name] * VHxsecratio
            newerror = sqrt(globals()[name] * VHxsecratioerr + VHxsecratio * globals()[errorname])
            globals()[name] = newvalue
            globals()[errorname] = newerror
fixVH(); del fixVH
JHUXSWHa1L1Zg       = JHUXSWHa1;      JHUXSWHa1L1Zgerr       = JHUXSWHa1err
JHUXSWHa2L1Zg       = JHUXSWHa1;      JHUXSWHa2L1Zgerr       = JHUXSWHa1err
JHUXSWHa3L1Zg       = JHUXSWHa1;      JHUXSWHa3L1Zgerr       = JHUXSWHa1err
JHUXSWHL1L1Zg       = JHUXSWHa1;      JHUXSWHL1L1Zgerr       = JHUXSWHa1err

JHUXSHJJa2       = 14583.61;       JHUXSHJJa2err       = 0.94
JHUXSHJJa3       = 14397.13;       JHUXSHJJa3err       = 0.97
JHUXSHJJa2a3     = 29169.2;        JHUXSHJJa2a3err     = 2.1

JHUXSttHkappa    = 0.912135589;    JHUXSttHkappaerr    = 0.00143032
JHUXSttHkappatilde = 0.35609194;   JHUXSttHkappatildeerr = 0.000492662
JHUXSttHkappakappatilde = 1.8231162489;   JHUXSttHkappakappatildeerr = 0.00254131


#Subtract the pure component from the interference, then divide by (gi*gj)
JHUXSggH2L2la1a2   = (JHUXSggH2L2la1a2   -                        JHUXSggH2L2la1 - g2decay              **2 * JHUXSggH2L2la2  ) / (g2decay                                      )
JHUXSggH2L2la1a3   = (JHUXSggH2L2la1a3   -                        JHUXSggH2L2la1 - g4decay              **2 * JHUXSggH2L2la3  ) / (g4decay                                      )
JHUXSggH2L2la1L1   = (JHUXSggH2L2la1L1   -                        JHUXSggH2L2la1 - g1prime2decay_gen    **2 * JHUXSggH2L2lL1  ) / (g1prime2decay_gen                            )
JHUXSggH2L2la1L1Zg = (JHUXSggH2L2la1L1Zg -                        JHUXSggH2L2la1 - ghzgs1prime2decay_gen**2 * JHUXSggH2L2lL1Zg) / (ghzgs1prime2decay_gen                        )
JHUXSggH2L2la2a3   = (JHUXSggH2L2la2a3   - g2decay          **2 * JHUXSggH2L2la2 - g4decay              **2 * JHUXSggH2L2la3  ) / (g2decay               * g4decay              )
JHUXSggH2L2la2L1   = (JHUXSggH2L2la2L1   - g2decay          **2 * JHUXSggH2L2la2 - g1prime2decay_gen    **2 * JHUXSggH2L2lL1  ) / (g2decay               * g1prime2decay_gen    )
JHUXSggH2L2la2L1Zg = (JHUXSggH2L2la2L1Zg - g2decay          **2 * JHUXSggH2L2la2 - ghzgs1prime2decay_gen**2 * JHUXSggH2L2lL1Zg) / (g2decay               * ghzgs1prime2decay_gen)
JHUXSggH2L2la3L1   = (JHUXSggH2L2la3L1   - g4decay          **2 * JHUXSggH2L2la3 - g1prime2decay_gen    **2 * JHUXSggH2L2lL1  ) / (g4decay               * g1prime2decay_gen    )
JHUXSggH2L2la3L1Zg = (JHUXSggH2L2la3L1Zg - g4decay          **2 * JHUXSggH2L2la3 - ghzgs1prime2decay_gen**2 * JHUXSggH2L2lL1Zg) / (g4decay               * ghzgs1prime2decay_gen)
JHUXSggH2L2lL1L1Zg = (JHUXSggH2L2lL1L1Zg - g1prime2decay_gen**2 * JHUXSggH2L2lL1 - ghzgs1prime2decay_gen**2 * JHUXSggH2L2lL1Zg) / (g1prime2decay_gen     * ghzgs1prime2decay_gen)

JHUXSVBFa1a2       = (JHUXSVBFa1a2       -                        JHUXSVBFa1     - g2VBF                **2 * JHUXSVBFa2      ) / (g2VBF                                        )
JHUXSVBFa1a3       = (JHUXSVBFa1a3       -                        JHUXSVBFa1     - g4VBF                **2 * JHUXSVBFa3      ) / (g4VBF                                        )
JHUXSVBFa1L1       = (JHUXSVBFa1L1       -                        JHUXSVBFa1     - g1prime2VBF_gen      **2 * JHUXSVBFL1      ) / (g1prime2VBF_gen                              )
JHUXSVBFa1L1Zg     = (JHUXSVBFa1L1Zg     -                        JHUXSVBFa1     - ghzgs1prime2VBF_gen  **2 * JHUXSVBFL1Zg    ) / (ghzgs1prime2VBF_gen                          )
JHUXSVBFa2a3       = (JHUXSVBFa2a3       - g2VBF            **2 * JHUXSVBFa2     - g4VBF                **2 * JHUXSVBFa3      ) / (g2VBF                 * g4VBF                )
JHUXSVBFa2L1       = (JHUXSVBFa2L1       - g2VBF            **2 * JHUXSVBFa2     - g1prime2VBF_gen      **2 * JHUXSVBFL1      ) / (g2VBF                 * g1prime2VBF_gen      )
JHUXSVBFa2L1Zg     = (JHUXSVBFa2L1Zg     - g2VBF            **2 * JHUXSVBFa2     - ghzgs1prime2VBF_gen  **2 * JHUXSVBFL1Zg    ) / (g2VBF                 * ghzgs1prime2VBF_gen  )
JHUXSVBFa3L1       = (JHUXSVBFa3L1       - g4VBF            **2 * JHUXSVBFa3     - g1prime2VBF_gen      **2 * JHUXSVBFL1      ) / (g4VBF                 * g1prime2VBF_gen      )
JHUXSVBFa3L1Zg     = (JHUXSVBFa3L1Zg     - g4VBF            **2 * JHUXSVBFa3     - ghzgs1prime2VBF_gen  **2 * JHUXSVBFL1Zg    ) / (g4VBF                 * ghzgs1prime2VBF_gen  )
JHUXSVBFL1L1Zg     = (JHUXSVBFL1L1Zg     - g1prime2VBF_gen  **2 * JHUXSVBFL1     - ghzgs1prime2VBF_gen  **2 * JHUXSVBFL1Zg    ) / (g1prime2VBF_gen       * ghzgs1prime2VBF_gen  )

JHUXSZHa1a2        = (JHUXSZHa1a2        -                        JHUXSZHa1      - g2ZH                 **2 * JHUXSZHa2       ) / (g2ZH                                         )
JHUXSZHa1a3        = (JHUXSZHa1a3        -                        JHUXSZHa1      - g4ZH                 **2 * JHUXSZHa3       ) / (g4ZH                                         )
JHUXSZHa1L1        = (JHUXSZHa1L1        -                        JHUXSZHa1      - g1prime2ZH_gen       **2 * JHUXSZHL1       ) / (g1prime2ZH_gen                               )
JHUXSZHa1L1Zg      = (JHUXSZHa1L1Zg      -                        JHUXSZHa1      - ghzgs1prime2ZH_gen   **2 * JHUXSZHL1Zg     ) / (ghzgs1prime2ZH_gen                           )
JHUXSZHa2a3        = (JHUXSZHa2a3        - g2ZH             **2 * JHUXSZHa2      - g4ZH                 **2 * JHUXSZHa3       ) / (g2ZH                  * g4ZH                 )
JHUXSZHa2L1        = (JHUXSZHa2L1        - g2ZH             **2 * JHUXSZHa2      - g1prime2ZH_gen       **2 * JHUXSZHL1       ) / (g2ZH                  * g1prime2ZH_gen       )
JHUXSZHa2L1Zg      = (JHUXSZHa2L1Zg      - g2ZH             **2 * JHUXSZHa2      - ghzgs1prime2ZH_gen   **2 * JHUXSZHL1Zg     ) / (g2ZH                  * ghzgs1prime2ZH_gen   )
JHUXSZHa3L1        = (JHUXSZHa3L1        - g4ZH             **2 * JHUXSZHa3      - g1prime2ZH_gen       **2 * JHUXSZHL1       ) / (g4ZH                  * g1prime2ZH_gen       )
JHUXSZHa3L1Zg      = (JHUXSZHa3L1Zg      - g4ZH             **2 * JHUXSZHa3      - ghzgs1prime2ZH_gen   **2 * JHUXSZHL1Zg     ) / (g4ZH                  * ghzgs1prime2ZH_gen   )
JHUXSZHL1L1Zg      = (JHUXSZHL1L1Zg      - g1prime2ZH_gen   **2 * JHUXSZHL1      - ghzgs1prime2ZH_gen   **2 * JHUXSZHL1Zg     ) / (g1prime2ZH_gen        * ghzgs1prime2ZH_gen   )

JHUXSWHa1a2        = (JHUXSWHa1a2        -                        JHUXSWHa1      - g2WH                 **2 * JHUXSWHa2       ) / (g2WH                                         )
JHUXSWHa1a3        = (JHUXSWHa1a3        -                        JHUXSWHa1      - g4WH                 **2 * JHUXSWHa3       ) / (g4WH                                         )
JHUXSWHa1L1        = (JHUXSWHa1L1        -                        JHUXSWHa1      - g1prime2WH_gen       **2 * JHUXSWHL1       ) / (g1prime2WH_gen                               )
JHUXSWHa1L1Zg      = (JHUXSWHa1L1Zg      -                        JHUXSWHa1      - ghzgs1prime2WH_gen   **2 * JHUXSWHL1Zg     ) / (ghzgs1prime2WH_gen                           )
JHUXSWHa2a3        = (JHUXSWHa2a3        - g2WH             **2 * JHUXSWHa2      - g4WH                 **2 * JHUXSWHa3       ) / (g2WH                  * g4WH                 )
JHUXSWHa2L1        = (JHUXSWHa2L1        - g2WH             **2 * JHUXSWHa2      - g1prime2WH_gen       **2 * JHUXSWHL1       ) / (g2WH                  * g1prime2WH_gen       )
JHUXSWHa2L1Zg      = (JHUXSWHa2L1Zg      - g2WH             **2 * JHUXSWHa2      - ghzgs1prime2WH_gen   **2 * JHUXSWHL1Zg     ) / (g2WH                  * ghzgs1prime2WH_gen   )
JHUXSWHa3L1        = (JHUXSWHa3L1        - g4WH             **2 * JHUXSWHa3      - g1prime2WH_gen       **2 * JHUXSWHL1       ) / (g4WH                  * g1prime2WH_gen       )
JHUXSWHa3L1Zg      = (JHUXSWHa3L1Zg      - g4WH             **2 * JHUXSWHa3      - ghzgs1prime2WH_gen   **2 * JHUXSWHL1Zg     ) / (g4WH                  * ghzgs1prime2WH_gen   )
JHUXSWHL1L1Zg      = (JHUXSWHL1L1Zg      - g1prime2WH_gen   **2 * JHUXSWHL1      - ghzgs1prime2WH_gen   **2 * JHUXSWHL1Zg     ) / (g1prime2WH_gen        * ghzgs1prime2WH_gen   )

JHUXSHJJa2a3       = (JHUXSHJJa2a3       -                        JHUXSHJJa2     - ghg4HJJ              **2 * JHUXSHJJa3      ) / (ghg4HJJ                                      )

JHUXSttHkappakappatilde = (JHUXSttHkappakappatilde - JHUXSttHkappa - kappa_tilde_ttH**2 * JHUXSttHkappatilde) / kappa_tilde_ttH

if __name__ == "__main__":
    print "All of the following should be 0:"
    print
    print "  decay:"
    print "    a1XS -           g2**2*    a2XS = {}%".format((JHUXSggH2L2la1 - g2decay**2               * JHUXSggH2L2la2    ) / JHUXSggH2L2la1 * 100)
    print "    a1XS -           g4**2*    a3XS = {}%".format((JHUXSggH2L2la1 - g4decay**2               * JHUXSggH2L2la3    ) / JHUXSggH2L2la1 * 100)
    print "    a1XS -     g1prime2**2*    L1XS = {}%".format((JHUXSggH2L2la1 - g1prime2decay_gen**2     * JHUXSggH2L2lL1    ) / JHUXSggH2L2la1 * 100)
    print "    a1XS - ghzgs1prime2**2*  L1ZgXS = {}%".format((JHUXSggH2L2la1 - ghzgs1prime2decay_gen**2 * JHUXSggH2L2lL1Zg  ) / JHUXSggH2L2la1 * 100)
    print "                        g4*  a1a3XS = {}%".format((                 g4decay                  * JHUXSggH2L2la1a3  ) / JHUXSggH2L2la1 * 100)
    print "                     g2*g4*  a2a3XS = {}%".format((                 g2decay*g4decay          * JHUXSggH2L2la2a3  ) / JHUXSggH2L2la1 * 100)
    print "               g1prime2*g4*  a3L1XS = {}%".format((                g1prime2decay_gen*g4decay * JHUXSggH2L2la3L1  ) / JHUXSggH2L2la1 * 100)
    print "           ghzgs1prime2*g4*a3L1ZgXS = {}%".format((            ghzgs1prime2decay_gen*g4decay * JHUXSggH2L2la3L1Zg) / JHUXSggH2L2la1 * 100)
    print
    print "  VBF:"
    print "    a1XS -           g2**2*    a2XS = {}%".format((JHUXSVBFa1     - g2VBF**2                 * JHUXSVBFa2        ) / JHUXSVBFa1     * 100)
    print "    a1XS -           g4**2*    a3XS = {}%".format((JHUXSVBFa1     - g4VBF**2                 * JHUXSVBFa3        ) / JHUXSVBFa1     * 100)
    print "    a1XS -     g1prime2**2*    L1XS = {}%".format((JHUXSVBFa1     - g1prime2VBF_gen**2       * JHUXSVBFL1        ) / JHUXSVBFa1     * 100)
#   print "    a1XS - ghzgs1prime2**2*  L1ZgXS = {}%".format((JHUXSVBFa1     - ghzgs1prime2VBF_gen**2   * JHUXSVBFL1Zg      ) / JHUXSVBFa1     * 100)
    print "                        g4*  a1a3XS = {}%".format((                 g4VBF                    * JHUXSVBFa1a3      ) / JHUXSVBFa1     * 100)
    print "                     g2*g4*  a2a3XS = {}%".format((                 g2VBF  *g4VBF            * JHUXSVBFa2a3      ) / JHUXSVBFa1     * 100)
    print "               g1prime2*g4*  a3L1XS = {}%".format((                g1prime2VBF_gen  *g4VBF   * JHUXSVBFa3L1      ) / JHUXSVBFa1     * 100)
#   print "           ghzgs1prime2*g4*a3L1ZgXS = {}%".format((            ghzgs1prime2VBF_gen  *g4VBF   * JHUXSVBFa3L1Zg    ) / JHUXSVBFa1     * 100)
    print
    print "  ZH:"
    print "    a1XS -           g2**2*    a2XS = {}%".format((JHUXSZHa1      - g2ZH**2                  * JHUXSZHa2         ) / JHUXSZHa1      * 100)
    print "    a1XS -           g4**2*    a3XS = {}%".format((JHUXSZHa1      - g4ZH**2                  * JHUXSZHa3         ) / JHUXSZHa1      * 100)
    print "    a1XS -     g1prime2**2*    L1XS = {}%".format((JHUXSZHa1      - g1prime2ZH_gen**2        * JHUXSZHL1         ) / JHUXSZHa1      * 100)
    print "    a1XS - ghzgs1prime2**2*  L1ZgXS = {}%".format((JHUXSZHa1      - ghzgs1prime2ZH_gen**2    * JHUXSZHL1Zg       ) / JHUXSZHa1      * 100)
    print "                        g4*  a1a3XS = {}%".format((                 g4ZH                     * JHUXSZHa1a3       ) / JHUXSZHa1      * 100)
    print "                     g2*g4*  a2a3XS = {}%".format((                 g2ZH   *g4ZH             * JHUXSZHa2a3       ) / JHUXSZHa1      * 100)
    print "               g1prime2*g4*  a3L1XS = {}%".format((                g1prime2ZH_gen   *g4ZH    * JHUXSZHa3L1       ) / JHUXSZHa1      * 100)
    print "           ghzgs1prime2*g4*a3L1ZgXS = {}%".format((            ghzgs1prime2ZH_gen   *g4ZH    * JHUXSZHa3L1Zg     ) / JHUXSZHa1      * 100)
    print
    print "  WH:"
    print "    a1XS -           g2**2*    a2XS = {}%".format((JHUXSWHa1      - g2WH**2                  * JHUXSWHa2         ) / JHUXSWHa1      * 100)
    print "    a1XS -           g4**2*    a3XS = {}%".format((JHUXSWHa1      - g4WH**2                  * JHUXSWHa3         ) / JHUXSWHa1      * 100)
    print "    a1XS -     g1prime2**2*    L1XS = {}%".format((JHUXSWHa1      - g1prime2WH_gen**2        * JHUXSWHL1         ) / JHUXSWHa1      * 100)
#   print "    a1XS - ghzgs1prime2**2*  L1ZgXS = {}%".format((JHUXSWHa1      - ghzgs1prime2WH_gen**2    * JHUXSWHL1Zg       ) / JHUXSWHa1      * 100)
    print "                        g4*  a1a3XS = {}%".format((                 g4WH                     * JHUXSWHa1a3       ) / JHUXSWHa1      * 100)
    print "                     g2*g4*  a2a3XS = {}%".format((                 g2WH   *g4WH             * JHUXSWHa2a3       ) / JHUXSWHa1      * 100)
    print "               g1prime2*g4*  a3L1XS = {}%".format((                g1prime2WH_gen   *g4WH    * JHUXSWHa3L1       ) / JHUXSWHa1      * 100)
#   print "           ghzgs1prime2*g4*a3L1ZgXS = {}%".format((            ghzgs1prime2WH_gen   *g4WH    * JHUXSWHa3L1Zg     ) / JHUXSWHa1      * 100)
    print "    WpHXS + WmHXS - WHXS            = {}%".format((SMXSWpH2L2l + SMXSWmH2L2l - SMXSWH2L2l                        ) / SMXSWH2L2l     * 100)
    print
    print "  HJJ:"
    print "    a2XS -           g4**2*    a3XS = {}%".format((JHUXSHJJa2     - ghg4HJJ**2               * JHUXSHJJa3        ) / JHUXSHJJa2     * 100)
    print "                        g4*  a2a3XS = {}%".format((                 ghg4HJJ                  * JHUXSHJJa2a3      ) / JHUXSHJJa2     * 100)
    print
    print "  ttH:"
    kt = kappa_tilde_ttH
    JHUXSttHk = JHUXSttHkappa
    JHUXSttHkt = JHUXSttHkappatilde
    JHUXSttHkkt = JHUXSttHkappakappatilde
    print "     kXS -           k~**2*    k~XS = {}%".format((JHUXSttHk      - kt**2                    * JHUXSttHkt        ) / JHUXSttHk      * 100)
    print "                        k~*   kk~XS = {}%".format((                 kt                       * JHUXSttHkkt       ) / JHUXSttHk      * 100)
    del kt, JHUXSttHk, JHUXSttHkt, JHUXSttHkkt

#Set them to exactly 0

JHUXSggH2L2la1a3   = JHUXSggH2L2la1a3err   = \
JHUXSggH2L2la2a3,  = JHUXSggH2L2la2a3err   = \
JHUXSggH2L2la3L1,  = JHUXSggH2L2la3L1err   = \
JHUXSggH2L2la3L1Zg = JHUXSggH2L2la3L1Zgerr = 0

JHUXSVBFa1a3   = JHUXSVBFa1a3err   = \
JHUXSVBFa2a3   = JHUXSVBFa2a3err   = \
JHUXSVBFa3L1   = JHUXSVBFa3L1err   = \
JHUXSVBFa3L1Zg = JHUXSVBFa3L1Zgerr = 0

JHUXSZHa1a3   = JHUXSZHa1a3err   = \
JHUXSZHa2a3   = JHUXSZHa2a3err   = \
JHUXSZHa3L1   = JHUXSZHa3L1err   = \
JHUXSZHa3L1Zg = JHUXSZHa3L1Zgerr = 0

JHUXSWHa1a3   = JHUXSWHa1a3err   = \
JHUXSWHa2a3   = JHUXSWHa2a3err   = \
JHUXSWHa3L1   = JHUXSWHa3L1err   = \
JHUXSWHa3L1Zg = JHUXSWHa3L1Zgerr = 0

JHUXSHJJa2a3 = JHUXSHJJa2a3err = 0

JHUXSttHkappakappatilde = JHUXSttHkappakappatildeerr = 0

g2VH = sqrt((JHUXSZHa1 + JHUXSWHa1) / (JHUXSZHa2 + JHUXSWHa2))
g4VH = sqrt((JHUXSZHa1 + JHUXSWHa1) / (JHUXSZHa3 + JHUXSWHa3))
g1prime2VH_gen = -sqrt((JHUXSZHa1 + JHUXSWHa1) / (JHUXSZHL1 + JHUXSWHL1))
g1prime2VH_reco = -g1prime2VH_gen
ghzgs1prime2VH_gen = -sqrt((JHUXSZHa1 + JHUXSWHa1) / (JHUXSZHL1Zg + JHUXSWHL1Zg))
ghzgs1prime2VH_reco = -ghzgs1prime2VH_gen

#defined this way, just make sure
for _ in """
  JHUXSggH2L2la1a3 JHUXSggH2L2la2a3 JHUXSggH2L2la3L1 JHUXSggH2L2la1a3_photoncut JHUXSggH2L2la2a3_photoncut JHUXSggH2L2la3L1_photoncut JHUXSggH2L2la3L1Zg
  JHUXSVBFa1a3 JHUXSVBFa2a3 JHUXSVBFa3L1 JHUXSVBFa1a3_photoncut JHUXSVBFa2a3_photoncut JHUXSVBFa3L1_photoncut JHUXSVBFa3L1Zg
  JHUXSZHa1a3 JHUXSZHa2a3 JHUXSZHa3L1 JHUXSZHa1a3_photoncut JHUXSZHa2a3_photoncut JHUXSZHa3L1_photoncut JHUXSZHa3L1Zg
  JHUXSWHa1a3 JHUXSWHa2a3 JHUXSWHa3L1
  JHUXSHJJa2a3 JHUXSttHkappakappatilde
  JHUXSWHL1Zg JHUXSWHa1L1Zg JHUXSWHa2L1Zg JHUXSWHa3L1Zg JHUXSWHL1L1Zg
""".split():
  assert globals()[_] == 0, (_, globals()[_])
for k, v in globals().items():
    if "__" in k: continue
    assert v is not None, k
del k, v, _
