from collections import namedtuple
from math import pi
import re

from config import defaultnbins
from utilities import rreplace

class Discriminant(namedtuple("Discriminant", "name title bins min max identifier formula")):
    def __new__(cls, name, title, bins, min, max, identifier=None, formula=None):
        if identifier is None: identifier = name
        if formula is None: formula = name
        min = float(min)
        max = float(max)
        return super(Discriminant, cls).__new__(cls, name=name, title=title, bins=bins, min=min, max=max, identifier=identifier, formula=formula)

decaydiscriminants = [
    Discriminant("D_bkg", "D_{bkg}", defaultnbins, 0, 1),
    Discriminant("D_bkg_ResUp", "D_{bkg}^{ResUp}", defaultnbins, 0, 1),
    Discriminant("D_bkg_ResDown", "D_{bkg}^{ResDown}", defaultnbins, 0, 1),
    Discriminant("D_bkg_ScaleUp", "D_{bkg}^{ScaleUp}", defaultnbins, 0, 1),
    Discriminant("D_bkg_ScaleDown", "D_{bkg}^{ScaleDown}", defaultnbins, 0, 1),

    Discriminant("D_0minus_decay", "D_{0-}^{dec}", defaultnbins, 0, 1),
    Discriminant("D_CP_decay", "D_{CP}^{dec}", defaultnbins, -0.5, 0.5),
    Discriminant("D_0hplus_decay", "D_{0h+}^{dec}", defaultnbins, 0, 1),
    Discriminant("D_int_decay", "D_{int}^{dec}", defaultnbins, 0, 1),
    Discriminant("D_L1_decay", "D_{#Lambda1}^{dec}", defaultnbins, .3, 1),
    Discriminant("D_L1int_decay", "D_{#Lambda1int}^{dec}", defaultnbins, -1, -.75),
    Discriminant("D_L1Zg_decay", "D_{#Lambda1}^{Z#gamma,dec}", defaultnbins, .25, 1),
    Discriminant("D_L1Zgint_decay", "D_{#Lambda1Z#gammaint}^{dec}", defaultnbins, -.25, .25),
    Discriminant("D_L1L1Zg_decay", "D_{#Lambda1}^{ZZ/Z#gamma,dec}", defaultnbins, 0, 1),
    Discriminant("D_L1L1Zgint_decay", "D_{#Lambda1int}^{ZZ/Z#gamma,dec}", defaultnbins, -.3, .25),

    Discriminant("D_CP_decay_new", "D_{CP}^{dec}", defaultnbins, -1, 1),
    Discriminant("D_int_decay_new", "D_{int}^{dec}", defaultnbins, -1, 1),

    Discriminant("D_eL_decay", "D_{#epsilonL}^{dec}", defaultnbins, .12, .5),
    Discriminant("D_eR_decay", "D_{#epsilonL}^{dec}", defaultnbins, .12, .5),
    Discriminant("D_eLeR_decay", "D_{#epsilonL/#epsilonR}^{dec}", defaultnbins, .42, .59),
    Discriminant("D_eLint_decay", "D_{#epsilonLint}^{dec}", defaultnbins, .55, 1),
    Discriminant("D_eRint_decay", "D_{#epsilonLint}^{dec}", defaultnbins, -1, -.5),
    Discriminant("D_eLeRint_decay", "D_{#epsilonL/#epsilonRint}^{dec}", defaultnbins, -.8, .2),

    Discriminant("Z1Mass", "m_{1}", defaultnbins, 20, 100),
    Discriminant("Z2Mass", "m_{2}", defaultnbins, 0, 62),
    Discriminant("Phi", "#Phi", defaultnbins, -pi, pi),
    Discriminant("phistarZ1", "dummy", 1, -pi, pi),
    Discriminant("phistarZ2", "dummy", 1, -pi, pi),

    Discriminant("D_4couplings_decay_raw", "D_{4}^{dec,raw}", 162, 0, 162),
    Discriminant("D_4couplings_decay", "D_{4}^{dec}", 162-50, 0, 162-50),
    Discriminant("D_CP_decay", "D_{CP}^{dec}", 2, -1, 1, identifier="D_CP_decay_2bins"),

    Discriminant("D_4couplings_photons_decay_raw", "D_{4#gamma}^{dec,raw}", 1296, 0, 1296),

    Discriminant("D_bkg", "D_{bkg}", 3, -0.3, 1.2, identifier="D_bkg_3bins"),   #set min and max to have boundaries at 0.2 and 0.7
    Discriminant("D_bkg_ResUp", "D_{bkg}^{ResUp}", 3, -0.3, 1.2, identifier="D_bkg_ResUp_3bins"),
    Discriminant("D_bkg_ResDown", "D_{bkg}^{ResDown}", 3, -0.3, 1.2, identifier="D_bkg_ResDown_3bins"),
    Discriminant("D_bkg_ScaleUp", "D_{bkg}^{ScaleUp}", 3, -0.3, 1.2, identifier="D_bkg_ScaleUp_3bins"),
    Discriminant("D_bkg_ScaleDown", "D_{bkg}^{ScaleDown}", 3, -0.3, 1.2, identifier="D_bkg_ScaleDown_3bins"),

    Discriminant("D_0minus_decay", "D_{0-}^{dec}", 3, 0, 1, identifier="D_0minus_decay_3bins"),
    Discriminant("D_0hplus_decay", "D_{0h+}^{dec}", 3, .3, .9, identifier="D_0hplus_decay_3bins", formula="max(min(D_0hplus_decay, .8), .4)"),
    Discriminant("D_L1_decay", "D_{#Lambda1}^{dec}", 3, .3, 1.05, identifier="D_L1_decay_3bins", formula="max(D_L1_decay, .4)"),
    Discriminant("D_L1Zg_decay", "D_{#Lambda1}^{Z#gamma,dec}", 3, .25, .7, identifier="D_L1Zg_decay_3bins", formula="max(min(D_L1Zg_decay, .6), .3)"),
    Discriminant("D_int_decay", "D_{int}^{dec}", 2, 0, 1.6, identifier="D_int_decay_2bins"),

    Discriminant("D_bkg", "D_{bkg}", 10, 0, 1, identifier="D_bkg_10bins"),
    Discriminant("D_0minus_decay", "D_{0-}^{dec}", 10, 0, 1, identifier="D_0minus_decay_10bins"),
    Discriminant("D_0hplus_decay", "D_{0h+}^{dec}", 10, 0, 1, identifier="D_0hplus_decay_10bins"),
    Discriminant("D_L1_decay", "D_{#Lambda1}^{dec}", 10, 0, 1, identifier="D_L1_decay_10bins"),
    Discriminant("D_L1Zg_decay", "D_{#Lambda1}^{Z#gamma,dec}", 10, 0, 1, identifier="D_L1Zg_decay_10bins"),
    Discriminant("D_CP_decay", "D_{CP}^{dec}", 10, -1, 1, identifier="D_CP_decay_10bins"),
    Discriminant("D_int_decay", "D_{int}^{dec}", 10, -1, 1, identifier="D_int_decay_10bins"),
    Discriminant("D_CP_decay_new", "D_{CP}^{dec}", 10, -1, 1, identifier="D_CP_decay_new_10bins"),
    Discriminant("D_int_decay_new", "D_{int}^{dec}", 10, -1, 1, identifier="D_int_decay_new_10bins"),

    Discriminant("D_bkg", "D_{bkg}", 20, 0, 1, identifier="D_bkg_20bins"),
    Discriminant("D_bkg_ResUp", "D_{bkg}^{ResUp}", 20, 0, 1, identifier="D_bkg_ResUp_20bins"),
    Discriminant("D_bkg_ResDown", "D_{bkg}^{ResDown}", 20, 0, 1, identifier="D_bkg_ResDown_20bins"),
    Discriminant("D_bkg_ScaleUp", "D_{bkg}^{ScaleUp}", 20, 0, 1, identifier="D_bkg_ScaleUp_20bins"),
    Discriminant("D_bkg_ScaleDown", "D_{bkg}^{ScaleDown}", 20, 0, 1, identifier="D_bkg_ScaleDown_20bins"),

    Discriminant("D_0minus_decay", "D_{0-}^{dec}", 20, 0, 1, identifier="D_0minus_decay_20bins"),
    Discriminant("D_0hplus_decay", "D_{0h+}^{dec}", 20, 0, 1, identifier="D_0hplus_decay_20bins"),
    Discriminant("D_L1_decay", "D_{#Lambda1}^{dec}", 20, 0, 1, identifier="D_L1_decay_20bins"),
    Discriminant("D_L1Zg_decay", "D_{#Lambda1}^{Z#gamma,dec}", 20, 0, 1, identifier="D_L1Zg_decay_20bins"),
    Discriminant("D_CP_decay_new", "D_{CP}^{dec}", 20, -1, 1, identifier="D_CP_decay_new_20bins"),
    Discriminant("D_int_decay_new", "D_{int}^{dec}", 20, -1, 1, identifier="D_int_decay_new_20bins"),

    Discriminant("D_bkg", "D_{bkg}", 30, 0, 1, identifier="D_bkg_30bins"),
    Discriminant("D_0minus_decay", "D_{0-}^{dec}", 30, 0, 1, identifier="D_0minus_decay_30bins"),
    Discriminant("D_0hplus_decay", "D_{0h+}^{dec}", 30, 0, 1, identifier="D_0hplus_decay_30bins"),
    Discriminant("D_L1_decay", "D_{#Lambda1}^{dec}", 30, 0, 1, identifier="D_L1_decay_30bins"),
    Discriminant("D_L1Zg_decay", "D_{#Lambda1}^{Z#gamma,dec}", 30, 0, 1, identifier="D_L1Zg_decay_30bins"),
    Discriminant("D_CP_decay_new", "D_{CP}^{dec}", 30, -1, 1, identifier="D_CP_decay_new_30bins"),
    Discriminant("D_int_decay_new", "D_{int}^{dec}", 30, -1, 1, identifier="D_int_decay_new_30bins"),

    Discriminant("ZZPt", "p_{T}^{4l}", 6, 100, 700, identifier="ZZPt_boosted"),
    Discriminant("ZZPt", "p_{T}^{4l}", 3, 0, 180, identifier="ZZPt_VBF1jtagged"),
    Discriminant("ZZPt", "p_{T}^{4l}", 4, 0, 400, identifier="ZZPt_VHLepttagged"),
]
jetdiscriminants = [
    Discriminant("D_0minus_VBF", "D_{0-}^{VBF}", defaultnbins, 0, 1),
    Discriminant("D_CP_VBF", "D_{CP}^{VBF}", defaultnbins, -1, 1),
    Discriminant("D_0hplus_VBF", "D_{0h+}^{VBF}", defaultnbins, 0, 1),
    Discriminant("D_int_VBF", "D_{int}^{VBF}", defaultnbins, -1, 1),
    Discriminant("D_L1_VBF", "D_{#Lambda1}^{VBF}", defaultnbins, 0, 1),
    Discriminant("D_L1int_VBF", "D_{#Lambda1int}^{VBF}", defaultnbins, -1, 1),
    Discriminant("D_L1Zg_VBF", "D_{#Lambda1Z#gamma}^{VBF}", defaultnbins, 0, 1),
    Discriminant("D_L1Zgint_VBF", "D_{#Lambda1Z#gammaint}^{VBF}", defaultnbins, -.35, 0),
    Discriminant("D_0minus_VBFdecay", "D_{0-}^{VBF+dec}", defaultnbins, 0, 1),
    Discriminant("D_0hplus_VBFdecay", "D_{0h+}^{VBF+dec}", defaultnbins, 0, 1),
    Discriminant("D_L1_VBFdecay", "D_{#Lambda1}^{VBF+dec}", defaultnbins, 0, 1),
    Discriminant("D_L1Zg_VBFdecay", "D_{#Lambda1}^{Z#gamma,VBF+dec}", defaultnbins, 0, 1),

    Discriminant("D_4couplings_VBFdecay_raw", "D_{4}^{VBF+dec,raw}", 162, 0, 162),
    Discriminant("D_4couplings_VBFdecay", "D_{4}^{VBF+dec}", 162-86, 0, 162-86),

    Discriminant("D_CP_VBF", "D_{CP}^{dec}", 2, -1, 1, identifier="D_CP_VBF_2bins"),

    Discriminant("D_4couplings_photons_VBF_raw", "D_{4#gamma}^{VBF,raw}", 1296, 0, 1296),
    Discriminant("D_4couplings_photons_VBFdecay_raw", "D_{4#gamma}^{VBF+dec,raw}", 1296, 0, 1296),

    Discriminant("D_bkg_VBFdecay", "D_{bkg}", 10, 0, 1, identifier="D_bkg_VBFdecay_10bins"),
    Discriminant("D_bkg_VBFdecay_ResUp", "D_{bkg}^{ResUp}", 10, 0, 1, identifier="D_bkg_VBFdecay_ResUp_10bins"),
    Discriminant("D_bkg_VBFdecay_ResDown", "D_{bkg}^{ResDown}", 10, 0, 1, identifier="D_bkg_VBFdecay_ResDown_10bins"),
    Discriminant("D_bkg_VBFdecay_ScaleUp", "D_{bkg}^{ScaleUp}", 10, 0, 1, identifier="D_bkg_VBFdecay_ScaleUp_10bins"),
    Discriminant("D_bkg_VBFdecay_ScaleDown", "D_{bkg}^{ScaleDown}", 10, 0, 1, identifier="D_bkg_VBFdecay_ScaleDown_10bins"),

    Discriminant("D_bkg_VBFdecay", "D_{bkg}", 3, -0.3, 1.2, identifier="D_bkg_VBFdecay_3bins"),   #set min and max to have boundaries at 0.2 and 0.7
    Discriminant("D_bkg_VBFdecay_ResUp", "D_{bkg}^{ResUp}", 3, -0.3, 1.2, identifier="D_bkg_VBFdecay_ResUp_3bins"),
    Discriminant("D_bkg_VBFdecay_ResDown", "D_{bkg}^{ResDown}", 3, -0.3, 1.2, identifier="D_bkg_VBFdecay_ResDown_3bins"),
    Discriminant("D_bkg_VBFdecay_ScaleUp", "D_{bkg}^{ScaleUp}", 3, -0.3, 1.2, identifier="D_bkg_VBFdecay_ScaleUp_3bins"),
    Discriminant("D_bkg_VBFdecay_ScaleDown", "D_{bkg}^{ScaleDown}", 3, -0.3, 1.2, identifier="D_bkg_VBFdecay_ScaleDown_3bins"),

    Discriminant("D_0minus_VBFdecay", "D_{0-}^{VBF+dec}", 3, -0.7, 1.7, identifier="D_0minus_VBFdecay_3bins"),  #set min and max to have boundaries at .1 and .9
    Discriminant("D_0hplus_VBFdecay", "D_{0-}^{VBF+dec}", 3, -0.7, 1.7, identifier="D_0hplus_VBFdecay_3bins"),  #set min and max to have boundaries at .1 and .9
    Discriminant("D_L1_VBFdecay", "D_{0-}^{VBF+dec}", 3, -0.7, 1.7, identifier="D_L1_VBFdecay_3bins"),  #set min and max to have boundaries at .1 and .9
    Discriminant("D_L1Zg_VBFdecay", "D_{0-}^{VBF+dec}", 3, -0.6, 1.5, identifier="D_L1Zg_VBFdecay_3bins"),  #set min and max to have boundaries at .1 and .8
    Discriminant("D_CP_VBF_new", "D_{CP}^{VBF}", 2, -1, 1, identifier="D_CP_VBF_new_2bins"),
    Discriminant("D_int_VBF_new", "D_{int}^{VBF}", 2, -1, 1, identifier="D_int_VBF_new_2bins"),

    Discriminant("D_0minus_VBFdecay", "D_{0-}^{VBF+dec}", 10, 0, 1, identifier="D_0minus_VBFdecay_10bins"),
    Discriminant("D_0hplus_VBFdecay", "D_{0h+}^{VBF+dec}", 10, 0, 1, identifier="D_0hplus_VBFdecay_10bins"),
    Discriminant("D_L1_VBFdecay", "D_{#Lambda1}^{VBF+dec}", 10, 0, 1, identifier="D_L1_VBFdecay_10bins"),
    Discriminant("D_L1Zg_VBFdecay", "D_{#Lambda1}^{Z#gamma,VBF+dec}", 10, 0, 1, identifier="D_L1Zg_VBFdecay_10bins"),
    Discriminant("D_CP_VBF", "D_{CP}^{VBF}", 10, -1, 1, identifier="D_CP_VBF_10bins"),
    Discriminant("D_int_VBF", "D_{int}^{VBF}", 10, -1, 1, identifier="D_int_VBF_10bins"),
    Discriminant("D_CP_VBF_new", "D_{CP}^{VBF}", 10, -1, 1, identifier="D_CP_VBF_new_10bins"),
    Discriminant("D_int_VBF_new", "D_{int}^{VBF}", 10, -1, 1, identifier="D_int_VBF_new_10bins"),

    Discriminant("D_bkg_VBFdecay", "D_{bkg}", 20, 0, 1, identifier="D_bkg_VBFdecay_20bins"),
    Discriminant("D_0minus_VBFdecay", "D_{0-}^{VBF+dec}", 20, 0, 1, identifier="D_0minus_VBFdecay_20bins"),
    Discriminant("D_0hplus_VBFdecay", "D_{0h+}^{VBF+dec}", 20, 0, 1, identifier="D_0hplus_VBFdecay_20bins"),
    Discriminant("D_L1_VBFdecay", "D_{#Lambda1}^{VBF+dec}", 20, 0, 1, identifier="D_L1_VBFdecay_20bins"),
    Discriminant("D_L1Zg_VBFdecay", "D_{#Lambda1}^{Z#gamma,VBF+dec}", 20, 0, 1, identifier="D_L1Zg_VBFdecay_20bins"),
    Discriminant("D_CP_VBF_new", "D_{CP}^{VBF}", 20, -1, 1, identifier="D_CP_VBF_new_20bins"),
    Discriminant("D_int_VBF_new", "D_{int}^{VBF}", 20, -1, 1, identifier="D_int_VBF_new_20bins"),

    Discriminant("D_0minus_HadVH", "D_{0-}^{VH}", defaultnbins, 0, 1),
    Discriminant("D_CP_HadVH", "D_{CP}^{VH}", defaultnbins, -.35, .35),
    Discriminant("D_0hplus_HadVH", "D_{0h+}^{VH}", defaultnbins, 0, 1),
    Discriminant("D_int_HadVH", "D_{int}^{VH}", defaultnbins, -1, .05),
    Discriminant("D_L1_HadVH", "D_{#Lambda1}^{VH}", defaultnbins, 0, 1),
    Discriminant("D_L1int_HadVH", "D_{#Lambda1int}^{VH}", defaultnbins, -1, 0),
    Discriminant("D_L1Zg_HadVH", "D_{#Lambda1}^{Z#gamma,VH}", defaultnbins, 0, 1),
    Discriminant("D_L1Zgint_HadVH", "D_{#Lambda1Z#gammaint}^{VH}", defaultnbins, 0, .4),
    Discriminant("D_0minus_HadVHdecay", "D_{0-}^{VH+dec}", defaultnbins, 0, 1),
    Discriminant("D_0hplus_HadVHdecay", "D_{0h+}^{VH+dec}", defaultnbins, 0, 1),
    Discriminant("D_L1_HadVHdecay", "D_{#Lambda1}^{VH+dec}", defaultnbins, 0, 1),
    Discriminant("D_L1Zg_HadVHdecay", "D_{#Lambda1}^{Z#gamma,VH+dec}", defaultnbins, 0, 1),

    Discriminant("D_4couplings_HadVHdecay_raw", "D_{4}^{VH+dec,raw}", 162, 0, 162),
    Discriminant("D_4couplings_HadVHdecay", "D_{4}^{VH+dec}", 162-73, 0, 162-73),

    Discriminant("D_CP_HadVH", "D_{CP}^{dec}", 2, -1, 1, identifier="D_CP_HadVH_2bins"),

    Discriminant("D_4couplings_photons_HadVH_raw", "D_{4#gamma}^{VH,raw}", 36, 0, 36),
    Discriminant("D_4couplings_photons_HadVHdecay_raw", "D_{4#gamma}^{VH+dec,raw}", 36, 0, 36),

    Discriminant("D_bkg_HadVHdecay", "D_{bkg}", 10, 0, 1, identifier="D_bkg_HadVHdecay_10bins"),
    Discriminant("D_bkg_HadVHdecay_ResUp", "D_{bkg}^{ResUp}", 10, 0, 1, identifier="D_bkg_HadVHdecay_ResUp_10bins"),
    Discriminant("D_bkg_HadVHdecay_ResDown", "D_{bkg}^{ResDown}", 10, 0, 1, identifier="D_bkg_HadVHdecay_ResDown_10bins"),
    Discriminant("D_bkg_HadVHdecay_ScaleUp", "D_{bkg}^{ScaleUp}", 10, 0, 1, identifier="D_bkg_HadVHdecay_ScaleUp_10bins"),
    Discriminant("D_bkg_HadVHdecay_ScaleDown", "D_{bkg}^{ScaleDown}", 10, 0, 1, identifier="D_bkg_HadVHdecay_ScaleDown_10bins"),

    Discriminant("D_bkg_HadVHdecay", "D_{bkg}", 3, -0.4, 1.4, identifier="D_bkg_HadVHdecay_3bins"),   #set min and max to have boundaries at 0.2 and 0.8
    Discriminant("D_bkg_HadVHdecay_ResUp", "D_{bkg}^{ResUp}", 3, -0.4, 1.4, identifier="D_bkg_HadVHdecay_ResUp_3bins"),
    Discriminant("D_bkg_HadVHdecay_ResDown", "D_{bkg}^{ResDown}", 3, -0.4, 1.4, identifier="D_bkg_HadVHdecay_ResDown_3bins"),
    Discriminant("D_bkg_HadVHdecay_ScaleUp", "D_{bkg}^{ScaleUp}", 3, -0.4, 1.4, identifier="D_bkg_HadVHdecay_ScaleUp_3bins"),
    Discriminant("D_bkg_HadVHdecay_ScaleDown", "D_{bkg}^{ScaleDown}", 3, -0.4, 1.4, identifier="D_bkg_HadVHdecay_ScaleDown_3bins"),

    Discriminant("D_0minus_HadVHdecay", "D_{0-}^{VH+dec}", 3, -0.4, 1.4, identifier="D_0minus_HadVHdecay_3bins"),   #min and max to have boundaries at .2 and .8
    Discriminant("D_0hplus_HadVHdecay", "D_{0h+}^{VH+dec}", 3, 0, 1, identifier="D_0hplus_HadVHdecay_3bins"),
    Discriminant("D_L1_HadVHdecay", "D_{#Lambda1}^{VH+dec}", 3, 0, 1, identifier="D_L1_HadVHdecay_3bins"),
    Discriminant("D_L1Zg_HadVHdecay", "D_{#Lambda1}^{Z#gamma,VH+dec}", 3, -0.7, 1.7, identifier="D_L1Zg_HadVHdecay_3bins"),  #set min and max to have boundaries at .1 and .9
    Discriminant("D_int_HadVH_new", "D_{int}^{VH+dec}", 2, -2.2, 1, identifier="D_int_HadVH_new_2bins"),  #set min and max to have the boundary at -0.6

    Discriminant("D_0minus_HadVHdecay", "D_{0-}^{VH+dec}", 10, 0, 1, identifier="D_0minus_HadVHdecay_10bins"),
    Discriminant("D_0hplus_HadVHdecay", "D_{0h+}^{VH+dec}", 10, 0, 1, identifier="D_0hplus_HadVHdecay_10bins"),
    Discriminant("D_L1_HadVHdecay", "D_{#Lambda1}^{VH+dec}", 10, 0, 1, identifier="D_L1_HadVHdecay_10bins"),
    Discriminant("D_L1Zg_HadVHdecay", "D_{#Lambda1}^{Z#gamma,VH+dec}", 10, 0, 1, identifier="D_L1Zg_HadVHdecay_10bins"),
    Discriminant("D_CP_HadVH", "D_{CP}^{VH}", 10, -1, 1, identifier="D_CP_HadVH_10bins"),
    Discriminant("D_int_HadVH", "D_{int}^{VH}", 10, -1, 1, identifier="D_int_HadVH_10bins"),
    Discriminant("D_CP_HadVH_new", "D_{CP}^{VH}", 10, -1, 1, identifier="D_CP_HadVH_new_10bins"),
    Discriminant("D_int_HadVH_new", "D_{int}^{VH}", 10, -1, 1, identifier="D_int_HadVH_new_10bins"),

    Discriminant("D_bkg_HadVHdecay", "D_{bkg}", 20, 0, 1, identifier="D_bkg_HadVHdecay_20bins"),
    Discriminant("D_0minus_HadVHdecay", "D_{0-}^{VH+dec}", 20, 0, 1, identifier="D_0minus_HadVHdecay_20bins"),
    Discriminant("D_0hplus_HadVHdecay", "D_{0h+}^{VH+dec}", 20, 0, 1, identifier="D_0hplus_HadVHdecay_20bins"),
    Discriminant("D_L1_HadVHdecay", "D_{#Lambda1}^{VH+dec}", 20, 0, 1, identifier="D_L1_HadVHdecay_20bins"),
    Discriminant("D_L1Zg_HadVHdecay", "D_{#Lambda1}^{Z#gamma,VH+dec}", 20, 0, 1, identifier="D_L1Zg_HadVHdecay_20bins"),
    Discriminant("D_CP_HadVH_new", "D_{CP}^{VH}", 20, -1, 1, identifier="D_CP_HadVH_new_20bins"),
    Discriminant("D_int_HadVH_new", "D_{int}^{VH}", 20, -1, 1, identifier="D_int_HadVH_new_20bins"),

    Discriminant("D_STXS_stage1p1", "D_{STXS}^{1.1}", 22, 0, 22),
    Discriminant("D_STXS_stage1p1", "D_{STXS}^{1.1}", 15, 0, 15, identifier="D_STXS_stage1p1_untagged", formula="D_STXS_stage1p1 - (D_STXS_stage1p1 >= 18) * 7"),
    Discriminant("D_STXS_stage1p1", "D_{STXS}^{1.1}", 5, 11, 16, identifier="D_STXS_stage1p1_VBF"),
    Discriminant("D_STXS_stage1p1", "D_{STXS}^{1.1}", 2, 16, 18, identifier="D_STXS_stage1p1_HadVH"),
]
categorydiscriminants = [
    Discriminant("D_2jet_0plus", "D_{2jet}^{VBF, 0+}", defaultnbins, 0, 1),
    Discriminant("D_2jet_0minus", "D_{2jet}^{VBF, 0-}", defaultnbins, 0, 1),
    Discriminant("D_2jet_a2", "D_{2jet}^{VBF, 0h+}", defaultnbins, 0, 1),
    Discriminant("D_2jet_L1", "D_{2jet}^{VBF, #Lambda1}", defaultnbins, 0, 1),
    Discriminant("D_2jet_L1Zg", "D_{2jet}^{VBF, #Lambda1Z#gamma}", defaultnbins, 0, 1),
    Discriminant("D_HadZH_0plus", "D_{2jet}^{ZH, 0+}", defaultnbins, 0, 1),
    Discriminant("D_HadZH_0minus", "D_{2jet}^{ZH, 0-}", defaultnbins, 0, 1),
    Discriminant("D_HadZH_a2", "D_{2jet}^{ZH, 0h+}", defaultnbins, 0, 1),
    Discriminant("D_HadZH_L1", "D_{2jet}^{ZH, #Lambda1}", defaultnbins, 0, 1),
    Discriminant("D_HadZH_L1Zg", "D_{2jet}^{ZH, #Lambda1Z#gamma}", defaultnbins, 0, 1),
    Discriminant("D_HadWH_0plus", "D_{2jet}^{WH, 0+}", defaultnbins, 0, 1),
    Discriminant("D_HadWH_0minus", "D_{2jet}^{WH, 0-}", defaultnbins, 0, 1),
    Discriminant("D_HadWH_a2", "D_{2jet}^{WH, 0h+}", defaultnbins, 0, 1),
    Discriminant("D_HadWH_L1", "D_{2jet}^{WH, #Lambda1}", defaultnbins, 0, 1),
    Discriminant("D_HadWH_L1Zg", "D_{2jet}^{WH, #Lambda1Z#gamma}", defaultnbins, 0, 1),
]

discriminants = decaydiscriminants + jetdiscriminants + categorydiscriminants
discriminants += list(
    Discriminant(name+"_"+JEC, rreplace(title, "}", ", "+JEC, 1), nbins, min, max, identifier+"_"+JEC, re.sub(r"(\b\w*\b)", r"\1"+"_"+JEC, formula))
        for name, title, nbins, min, max, identifier, formula in jetdiscriminants
        for JEC in ("JECUp", "JECDn", "JESUp", "JESDn")
)
if len(discriminants) != len({d.identifier for d in discriminants}):
    raise ValueError("Multiple discriminants have the same identifier")

discriminants = {d.identifier: d for d in discriminants}

del decaydiscriminants, jetdiscriminants, categorydiscriminants

otherplottablethings = {
    d.identifier: d for d in [
        Discriminant("DiJetMass", "m_{JJ}", 50, 0, 150),
    ]
}

def discriminant(identifier):
    if isinstance(identifier, Discriminant): return identifier
    if identifier in discriminants:
        return discriminants[identifier]
    if identifier in otherplottablethings:
        return otherplottablethings[identifier]
    raise KeyError("Unknown discriminant {}".format(identifier))

assert not set(otherplottablethings) & set(discriminants), set(otherplottablethings) & set(discriminants)
