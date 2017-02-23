from collections import namedtuple

from config import defaultnbins
from utilities import rreplace

Discriminant = namedtuple("Discriminant", "name title bins min max")

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
    Discriminant("D_L1_decay", "D_{#Lambda1}^{dec}", defaultnbins, 0, 1),
    Discriminant("D_L1int_decay", "D_{#Lambda1int}^{dec}", defaultnbins, 0, 1),
    Discriminant("D_L1Zg_decay", "D_{#Lambda1Z#gamma}^{dec}", defaultnbins, 0, 1),
    Discriminant("D_L1Zgint_decay", "D_{#Lambda1Z#gammaint}^{dec}", defaultnbins, -.3, .25),
]
VBFdiscriminants = [
    Discriminant("D_0minus_VBF", "D_{0-}^{VBF}", defaultnbins, 0, 1),
    Discriminant("D_CP_VBF", "D_{CP}^{VBF}", defaultnbins, -1, 1),
    Discriminant("D_0hplus_VBF", "D_{0h+}^{VBF}", defaultnbins, 0, 1),
    Discriminant("D_int_VBF", "D_{int}^{VBF}", defaultnbins, -1, 1),
    Discriminant("D_L1_VBF", "D_{#Lambda1}^{VBF}", defaultnbins, 0, 1),
    Discriminant("D_L1int_VBF", "D_{#Lambda1int}^{VBF}", defaultnbins, -1, 1),
    Discriminant("D_L1Zg_VBF", "D_{#LambdaZ#gamma1}^{VBF}", defaultnbins, 0, 1),
    Discriminant("D_L1Zgint_VBF", "D_{#Lambda1Z#gammaint}^{VBF}", defaultnbins, 0, .35),
    Discriminant("D_0minus_VBFdecay", "D_{0-}^{VBFdec}", defaultnbins, 0, 1),
    Discriminant("D_0hplus_VBFdecay", "D_{0h+}^{VBFdec}", defaultnbins, 0, 1),
    Discriminant("D_L1_VBFdecay", "D_{#Lambda1}^{VBFdec}", defaultnbins, 0, 1),
    Discriminant("D_L1Zg_VBFdecay", "D_{#Lambda1Z#gamma}^{VBFdec}", defaultnbins, 0, 1),
]
VHdiscriminants = [
    Discriminant("D_0minus_HadVH", "D_{0-}^{VH}", defaultnbins, 0, 1),
    Discriminant("D_CP_HadVH", "D_{CP}^{VH}", defaultnbins, -.35, .35),
    Discriminant("D_0hplus_HadVH", "D_{0h+}^{VH}", defaultnbins, 0, 1),
    Discriminant("D_int_HadVH", "D_{int}^{VH}", defaultnbins, -1, .05),
    Discriminant("D_L1_HadVH", "D_{#Lambda1}^{VH}", defaultnbins, 0, 1),
    Discriminant("D_L1int_HadVH", "D_{#Lambda1int}^{VH}", defaultnbins, 0, 1),
    Discriminant("D_L1Zg_HadVH", "D_{#Lambda1Z#gamma}^{VH}", defaultnbins, 0, 1),
    Discriminant("D_L1Zgint_HadVH", "D_{#Lambda1Z#gammaint}^{VH}", defaultnbins, -.4, 0),
    Discriminant("D_0minus_HadVHdecay", "D_{0-}^{VHdec}", defaultnbins, 0, 1),
    Discriminant("D_0hplus_HadVHdecay", "D_{0h+}^{VHdec}", defaultnbins, 0, 1),
    Discriminant("D_L1_HadVHdecay", "D_{#Lambda1}^{VHdec}", defaultnbins, 0, 1),
    Discriminant("D_L1Zg_HadVHdecay", "D_{#Lambda1Z#gamma}^{VHdec}", defaultnbins, 0, 1),
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
    Discriminant(name+"_"+JEC, rreplace(title, "}", ", "+JEC, 1), nbins, min, max)
        for name, title, nbins, min, max in VBFdiscriminants+VHdiscriminants
        for JEC in ("JECUp", "JECDn")
)
if len(discriminants) != len({d.name for d in discriminants}):
    raise ValueError("Multiple discriminants have the same name")

discriminants = {d.name: d for d in discriminants}

del decaydiscriminants, VBFdiscriminants, VHdiscriminants, categorydiscriminants

otherplottablethings = {
    d.name: d for d in [
        Discriminant("ZZPt", "p_{T}^{ZZ}", defaultnbins, 0, 500),
        Discriminant("DiJetMass", "m_{JJ}", 50, 0, 150),
    ]
}

def discriminant(name):
    if isinstance(name, Discriminant): return name
    if name in discriminants:
        return discriminants[name]
    if name in otherplottablethings:
        return otherplottablethings[name]
    raise KeyError("Unknown discriminants {}".format(name))

assert not any(_ in otherplottablethings for _ in discriminants)
