from collections import namedtuple
from math import pi

from config import defaultnbins
from utilities import rreplace

class Discriminant(namedtuple("Discriminant", "name title bins min max identifier")):
    def __new__(cls, name, title, bins, min, max, identifier=None):
        if identifier is None: identifier = name
        return super(Discriminant, cls).__new__(cls, name=name, title=title, bins=bins, min=min, max=max, identifier=identifier)

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
    Discriminant("D_L1int_decay", "D_{#Lambda1int}^{dec}", defaultnbins, .75, 1),
    Discriminant("D_L1Zg_decay", "D_{#Lambda1}^{Z#gamma,dec}", defaultnbins, .25, 1),
    Discriminant("D_L1Zgint_decay", "D_{#Lambda1Z#gammaint}^{dec}", defaultnbins, -.25, .25),
    Discriminant("D_L1L1Zg_decay", "D_{#Lambda1}^{ZZ/Z#gamma,dec}", defaultnbins, 0, 1),
    Discriminant("D_L1L1Zgint_decay", "D_{#Lambda1int}^{ZZ/Z#gamma,dec}", defaultnbins, -.3, .25),

    Discriminant("D_eL_decay", "D_{#epsilonL}^{dec}", defaultnbins, .12, .5),
    Discriminant("D_eR_decay", "D_{#epsilonL}^{dec}", defaultnbins, .12, .5),
    Discriminant("D_eLeR_decay", "D_{#epsilonL/#epsilonR}^{dec}", defaultnbins, .42, .59),
    Discriminant("D_eLint_decay", "D_{#epsilonLint}^{dec}", defaultnbins, .55, 1),
    Discriminant("D_eRint_decay", "D_{#epsilonLint}^{dec}", defaultnbins, -1, -.5),
    Discriminant("D_eLeRint_decay", "D_{#epsilonL/#epsilonRint}^{dec}", defaultnbins, -.8, .2),

    Discriminant("Z1Mass", "m_{1}", defaultnbins, 20, 100),
    Discriminant("Z2Mass", "m_{2}", defaultnbins, 0, 62),
    Discriminant("Phi", "#Phi", defaultnbins, -pi, pi),
    Discriminant("phistarZ2", "dummy", 1, -pi, pi),

    Discriminant("D_STXS_stage0", "D_{STXS0}", 2, 0, 2),

    Discriminant("D_4couplings_decay_raw", "D_{4}^{dec,raw}", 162, 0, 162),
    Discriminant("D_4couplings_decay", "D_{4}^{dec}", 162-40, 0, 162-40),
    Discriminant("D_CP_decay", "D_{CP}^{dec}", 2, -0.5, 0.5, identifier="D_CP_decay_2bins"),

    Discriminant("D_bkg", "D_{bkg}", 10, 0, 1, identifier="D_bkg_10bins"),
    Discriminant("D_bkg", "D_{bkg}", 20, 0, 1, identifier="D_bkg_20bins"),
    Discriminant("D_0minus_decay", "D_{0-}^{dec}", 20, 0, 1, identifier="D_0minus_decay_20bins"),
    Discriminant("D_CP_decay", "D_{CP}^{dec}", 20, -0.5, 0.5, identifier="D_CP_decay_20bins"),
    Discriminant("D_0hplus_decay", "D_{0h+}^{dec}", 20, 0, 1, identifier="D_0hplus_decay_20bins"),
    Discriminant("D_int_decay", "D_{int}^{dec}", 20, 0, 1, identifier="D_int_decay_20bins"),
    Discriminant("D_L1_decay", "D_{#Lambda1}^{dec}", 20, .3, 1, identifier="D_L1_decay_20bins"),
    Discriminant("D_L1int_decay", "D_{#Lambda1int}^{dec}", 20, .75, 1, identifier="D_L1int_decay_20bins"),
    Discriminant("D_L1Zg_decay", "D_{#Lambda1}^{Z#gamma,dec}", 20, .25, 1, identifier="D_L1Zg_decay_20bins"),
    Discriminant("D_L1Zgint_decay", "D_{#Lambda1Z#gammaint}^{dec}", 20, -.25, .25, identifier="D_L1Zgint_decay_20bins"),
]
VBFdiscriminants = [
    Discriminant("D_0minus_VBF", "D_{0-}^{VBF}", defaultnbins, 0, 1),
    Discriminant("D_CP_VBF", "D_{CP}^{VBF}", defaultnbins, -1, 1),
    Discriminant("D_0hplus_VBF", "D_{0h+}^{VBF}", defaultnbins, 0, 1),
    Discriminant("D_int_VBF", "D_{int}^{VBF}", defaultnbins, -1, 1),
    Discriminant("D_L1_VBF", "D_{#Lambda1}^{VBF}", defaultnbins, 0, 1),
    Discriminant("D_L1int_VBF", "D_{#Lambda1int}^{VBF}", defaultnbins, -1, 1),
    Discriminant("D_L1Zg_VBF", "D_{#Lambda1Z#gamma}^{VBF}", defaultnbins, 0, 1),
    Discriminant("D_L1Zgint_VBF", "D_{#Lambda1Z#gammaint}^{VBF}", defaultnbins, 0, .35),
    Discriminant("D_0minus_VBFdecay", "D_{0-}^{VBF+dec}", defaultnbins, 0, 1),
    Discriminant("D_0hplus_VBFdecay", "D_{0h+}^{VBF+dec}", defaultnbins, 0, 1),
    Discriminant("D_L1_VBFdecay", "D_{#Lambda1}^{VBF+dec}", defaultnbins, 0, 1),
    Discriminant("D_L1Zg_VBFdecay", "D_{#Lambda1}^{Z#gamma,VBF+dec}", defaultnbins, 0, 1),

    Discriminant("D_STXS_ggH_stage1", "D_{STXS1}^{ggH}", 12, 0, 12),
    Discriminant("D_STXS_VBF_stage1", "D_{STXS1}^{VBF}", 6, 0, 6),

    Discriminant("D_4couplings_VBFdecay_raw", "D_{4}^{VBF+dec,raw}", 162, 0, 162),
    Discriminant("D_4couplings_VBFdecay", "D_{4}^{VBF+dec}", 162-80, 0, 162-80),

    Discriminant("D_CP_VBF", "D_{CP}^{dec}", 2, -0.5, 0.5, identifier="D_CP_VBF_2bins"),

    Discriminant("D_0minus_VBFdecay", "D_{0-}^{VBF+dec}", 20, 0, 1, identifier="D_0minus_VBFdecay_20bins"),
    Discriminant("D_0hplus_VBFdecay", "D_{0h+}^{VBF+dec}", 20, 0, 1, identifier="D_0hplus_VBFdecay_20bins"),
    Discriminant("D_L1_VBFdecay", "D_{#Lambda1}^{VBF+dec}", 20, 0, 1, identifier="D_L1_VBFdecay_20bins"),
    Discriminant("D_L1Zg_VBFdecay", "D_{#Lambda1}^{Z#gamma,VBF+dec}", 20, 0, 1, identifier="D_L1Zg_VBFdecay_20bins"),
    Discriminant("D_CP_VBF", "D_{CP}^{VBF}", 20, -1, 1, identifier="D_CP_VBF_20bins"),
    Discriminant("D_int_VBF", "D_{int}^{VBF}", 20, -1, 1, identifier="D_int_VBF_20bins"),
]
VHdiscriminants = [
    Discriminant("D_0minus_HadVH", "D_{0-}^{VH}", defaultnbins, 0, 1),
    Discriminant("D_CP_HadVH", "D_{CP}^{VH}", defaultnbins, -.35, .35),
    Discriminant("D_0hplus_HadVH", "D_{0h+}^{VH}", defaultnbins, 0, 1),
    Discriminant("D_int_HadVH", "D_{int}^{VH}", defaultnbins, -1, .05),
    Discriminant("D_L1_HadVH", "D_{#Lambda1}^{VH}", defaultnbins, 0, 1),
    Discriminant("D_L1int_HadVH", "D_{#Lambda1int}^{VH}", defaultnbins, 0, 1),
    Discriminant("D_L1Zg_HadVH", "D_{#Lambda1}^{Z#gamma,VH}", defaultnbins, 0, 1),
    Discriminant("D_L1Zgint_HadVH", "D_{#Lambda1Z#gammaint}^{VH}", defaultnbins, -.4, 0),
    Discriminant("D_0minus_HadVHdecay", "D_{0-}^{VH+dec}", defaultnbins, 0, 1),
    Discriminant("D_0hplus_HadVHdecay", "D_{0h+}^{VH+dec}", defaultnbins, 0, 1),
    Discriminant("D_L1_HadVHdecay", "D_{#Lambda1}^{VH+dec}", defaultnbins, 0, 1),
    Discriminant("D_L1Zg_HadVHdecay", "D_{#Lambda1}^{Z#gamma,VH+dec}", defaultnbins, 0, 1),

    Discriminant("D_4couplings_HadVHdecay_raw", "D_{4}^{VH+dec,raw}", 162, 0, 162),
    Discriminant("D_4couplings_HadVHdecay", "D_{4}^{VH+dec}", 162-67, 0, 162-67),

    Discriminant("D_CP_HadVH", "D_{CP}^{dec}", 2, -0.5, 0.5, identifier="D_CP_HadVH_2bins"),

    Discriminant("D_0minus_HadVHdecay", "D_{0-}^{VH+dec}", 20, 0, 1, identifier="D_0minus_HadVHdecay_20bins"),
    Discriminant("D_0hplus_HadVHdecay", "D_{0h+}^{VH+dec}", 20, 0, 1, identifier="D_0hplus_HadVHdecay_20bins"),
    Discriminant("D_L1_HadVHdecay", "D_{#Lambda1}^{VH+dec}", 20, 0, 1, identifier="D_L1_HadVHdecay_20bins"),
    Discriminant("D_L1Zg_HadVHdecay", "D_{#Lambda1}^{Z#gamma,VH+dec}", 20, 0, 1, identifier="D_L1Zg_HadVHdecay_20bins"),
    Discriminant("D_CP_HadVH", "D_{CP}^{VH}", 20, -.35, .35, identifier="D_CP_HadVH_20bins"),
    Discriminant("D_int_HadVH", "D_{int}^{VH}", 20, -1, .05, identifier="D_int_HadVH_20bins"),
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

discriminants = decaydiscriminants + VBFdiscriminants + VHdiscriminants + categorydiscriminants
discriminants += list(
    Discriminant(name+"_"+JEC, rreplace(title, "}", ", "+JEC, 1), nbins, min, max, identifier+"_"+JEC)
        for name, title, nbins, min, max, identifier in VBFdiscriminants+VHdiscriminants
        for JEC in ("JECUp", "JECDn")
)
if len(discriminants) != len({d.identifier for d in discriminants}):
    raise ValueError("Multiple discriminants have the same identifier")

discriminants = {d.identifier: d for d in discriminants}

del decaydiscriminants, VBFdiscriminants, VHdiscriminants, categorydiscriminants

otherplottablethings = {
    d.identifier: d for d in [
        Discriminant("ZZPt", "p_{T}^{ZZ}", defaultnbins, 0, 500),
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
