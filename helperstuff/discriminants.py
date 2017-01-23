from collections import namedtuple
from config import defaultnbins

Discriminant = namedtuple("Discriminant", "name title bins min max")

minmax_g1jgik = {
    ("VBF", "g4", 1, ""): (-3, 3),
    ("VBF", "g4", 2, ""): (0, 5),
    ("VBF", "g4", 3, ""): (-4, 4),
    ("VBF", "g2", 1, ""): (-2, 3),
    ("VBF", "g2", 2, ""): (-1, 5),
    ("VBF", "g2", 3, ""): (-4, 4),
    ("VBF", "g1prime2", 1, ""): (-1, 4),
    ("VBF", "g1prime2", 2, ""): (-1, 2),
    ("VBF", "g1prime2", 3, ""): (-3, 3),

    ("VBF", "g4", 1, "_prime"): (-1, 1),
    ("VBF", "g4", 2, "_prime"): (0, 1),
    ("VBF", "g4", 3, "_prime"): (-.5, .5),
    ("VBF", "g2", 1, "_prime"): (-.5, 1.25),
    ("VBF", "g2", 2, "_prime"): (-.5, 1.5),
    ("VBF", "g2", 3, "_prime"): (-1, 1.25),
    ("VBF", "g1prime2", 1, "_prime"): (-.3, 1),
    ("VBF", "g1prime2", 2, "_prime"): (-.5, 1),
    ("VBF", "g1prime2", 3, "_prime"): (-.7, 1),

    ("ZH", "g4", 1, ""): (-1, 1),  #long tails go into over/underflow
    ("ZH", "g4", 2, ""): (0, 5),
    ("ZH", "g4", 3, ""): (-1, 1),
    ("ZH", "g2", 1, ""): (-3, 4),
    ("ZH", "g2", 2, ""): (-1, 2),   #long tails
    ("ZH", "g2", 3, ""): (-2, 2),
    ("ZH", "g1prime2", 1, ""): (0, 4),
    ("ZH", "g1prime2", 2, ""): (0, 5),
    ("ZH", "g1prime2", 3, ""): (0, 3),

    ("ZH", "g4", 1, "_prime"): (-.5, .5),
    ("ZH", "g4", 2, "_prime"): (-.05, 1),
    ("ZH", "g4", 3, "_prime"): (-.4, .4),  #long tails
    ("ZH", "g2", 1, "_prime"): (-1, 1),
    ("ZH", "g2", 2, "_prime"): (-.5, 1),
    ("ZH", "g2", 3, "_prime"): (-1, 1),
    ("ZH", "g1prime2", 1, "_prime"): (0, 1.3),
    ("ZH", "g1prime2", 2, "_prime"): (0, 1.5),
    ("ZH", "g1prime2", 3, "_prime"): (0, 1.3),

    ("WH", "g4", 1, ""): (-1, 1),
    ("WH", "g4", 2, ""): (0, 5),   #long tails
    ("WH", "g4", 3, ""): (-.6, .6),
    ("WH", "g2", 1, ""): (-1, 4),
    ("WH", "g2", 2, ""): (-1, 2),   #long tails
    ("WH", "g2", 3, ""): (-2, 2),
    ("WH", "g1prime2", 1, ""): (0, 4),
    ("WH", "g1prime2", 2, ""): (0, 5),
    ("WH", "g1prime2", 3, ""): (0, 3),

    ("WH", "g4", 1, "_prime"): (-.4, .4),
    ("WH", "g4", 2, "_prime"): (0, 1),
    ("WH", "g4", 3, "_prime"): (-.4, .4),
    ("WH", "g2", 1, "_prime"): (-.5, 1),
    ("WH", "g2", 2, "_prime"): (-.5, 1),
    ("WH", "g2", 3, "_prime"): (-1, 1),
    ("WH", "g1prime2", 1, "_prime"): (0, 1.3),
    ("WH", "g1prime2", 2, "_prime"): (0, 1.5),
    ("WH", "g1prime2", 3, "_prime"): (0, 1.3),
}

discriminants = {
    d.name: d for d in [
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
        Discriminant("D_0minus_VBF", "D_{0-}^{VBF}", defaultnbins, 0, 1),
        Discriminant("D_CP_VBF", "D_{CP}^{VBF}", defaultnbins, -1, 1),
        Discriminant("D_0hplus_VBF", "D_{0h+}^{VBF}", defaultnbins, 0, 1),
        Discriminant("D_int_VBF", "D_{int}^{VBF}", defaultnbins, -1, 1),
        Discriminant("D_L1_VBF", "D_{#Lambda1}^{VBF}", defaultnbins, 0, 1),
        Discriminant("D_L1int_VBF", "D_{#Lambda1int}^{VBF}", defaultnbins, -1, 1),
        Discriminant("D_L1Zg_VBF", "D_{#Lambda1}^{VBF}", defaultnbins, 0, 1),
        Discriminant("D_L1Zgint_VBF", "D_{#Lambda1int}^{VBF}", defaultnbins, 0, .35),
        Discriminant("D_0minus_VBFdecay", "D_{0-}^{VBFdec}", defaultnbins, 0, 1),
        Discriminant("D_0hplus_VBFdecay", "D_{0h+}^{VBFdec}", defaultnbins, 0, 1),
        Discriminant("D_L1_VBFdecay", "D_{#Lambda1}^{VBFdec}", defaultnbins, 0, 1),
        Discriminant("D_L1Zg_VBFdecay", "D_{#Lambda1Z#gamma}^{VBFdec}", defaultnbins, 0, 1),
#    ] + [
#        Discriminant(
#                     "D_g1{}_{}{}_VBFdecay{}".format(i, gj, 4-i, prime),
#                     "D{}^[VBFdec]_[g_[1]^[{}]g_[{}]^[{}]]".format("'" if prime else "", i, gj, 4-i).replace("[", "{").replace("]", "}"),
#                     defaultnbins,
#                     *minmax_g1jgik["VBF", gj, i, prime]
#                    )
#            for prime in ("", "_prime")
#            for gj in ("g4", "g2", "g1prime2")
#            for i in range(1, 4)
    ] + [
        Discriminant("D_0minus_HadZH", "D_{0-}^{ZH}", defaultnbins, 0, 1),
        Discriminant("D_CP_HadZH", "D_{CP}^{ZH}", defaultnbins, -.4, .4),
        Discriminant("D_0hplus_HadZH", "D_{0h+}^{ZH}", defaultnbins, 0, 1),
        Discriminant("D_int_HadZH", "D_{int}^{ZH}", defaultnbins, -1, .05),
        Discriminant("D_L1_HadZH", "D_{#Lambda1}^{ZH}", defaultnbins, 0, 1),
        Discriminant("D_L1int_HadZH", "D_{#Lambda1int}^{ZH}", defaultnbins, 0, 1),
        Discriminant("D_L1Zg_HadZH", "D_{#Lambda1Z#gamma}^{ZH}", defaultnbins, 0, 1),
        Discriminant("D_L1Zgint_HadZH", "D_{#Lambda1Z#gammaint}^{ZH}", defaultnbins, -1, 1),
        Discriminant("D_0minus_HadZHdecay", "D_{0-}^{ZHdec}", defaultnbins, 0, 1),
        Discriminant("D_0hplus_HadZHdecay", "D_{0h+}^{ZHdec}", defaultnbins, 0, 1),
        Discriminant("D_L1_HadZHdecay", "D_{#Lambda1}^{ZHdec}", defaultnbins, 0, 1),
        Discriminant("D_L1Zg_HadZHdecay", "D_{#Lambda1Z#gamma}^{ZHdec}", defaultnbins, 0, 1),
#    ] + [
#        Discriminant(
#                     "D_g1{}_{}{}_HadZHdecay{}".format(i, gj, 4-i, prime),
#                     "D{}^[ZHdec]_[g_[1]^[{}]g_[{}]^[{}]]".format("'" if prime else "", i, gj, 4-i).replace("[", "{").replace("]", "}"),
#                     defaultnbins,
#                     *minmax_g1jgik["ZH", gj, i, prime]
#                    )
#            for prime in ("", "_prime")
#            for gj in ("g4", "g2", "g1prime2")
#            for i in range(1, 4)
    ] + [
        Discriminant("D_0minus_HadWH", "D_{0-}^{WH}", defaultnbins, 0, 1),
        Discriminant("D_CP_HadWH", "D_{CP}^{WH}", defaultnbins, -.1, .1),
        Discriminant("D_0hplus_HadWH", "D_{0h+}^{WH}", defaultnbins, 0, 1),
        Discriminant("D_int_HadWH", "D_{int}^{WH}", defaultnbins, -1, .05),
        Discriminant("D_L1_HadWH", "D_{#Lambda1}^{WH}", defaultnbins, 0, 1),
        Discriminant("D_L1int_HadWH", "D_{#Lambda1int}^{WH}", defaultnbins, 0, 1),
        Discriminant("D_0minus_HadWHdecay", "D_{0-}^{WHdec}", defaultnbins, 0, 1),
        Discriminant("D_0hplus_HadWHdecay", "D_{0h+}^{WHdec}", defaultnbins, 0, 1),
        Discriminant("D_L1_HadWHdecay", "D_{#Lambda1}^{WHdec}", defaultnbins, 0, 1),
#    ] + [
#        Discriminant(
#                     "D_g1{}_{}{}_HadWHdecay{}".format(i, gj, 4-i, prime),
#                     "D{}^[WHdec]_[g_[1]^[{}]g_[{}]^[{}]]".format("'" if prime else "", i, gj, 4-i).replace("[", "{").replace("]", "}"),
#                     defaultnbins,
#                     *minmax_g1jgik["WH", gj, i, prime]
#                    )
#            for prime in ("", "_prime")
#            for gj in ("g4", "g2", "g1prime2")
#            for i in range(1, 4)
    ] + [
        Discriminant("D_0minus_HadVH", "D_{0-}^{VH}", defaultnbins, 0, 1),
        Discriminant("D_CP_HadVH", "D_{CP}^{VH}", defaultnbins, -.4, .4),
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
    ] + [
        Discriminant("D_2jet_0plus", "D_{2jet}^{0+}", defaultnbins, 0, 1),
        Discriminant("D_2jet_0minus", "D_{2jet}^{0-}", defaultnbins, 0, 1),
        Discriminant("D_2jet_a2", "D_{2jet}^{0h+}", defaultnbins, 0, 1),
        Discriminant("D_2jet_L1", "D_{2jet}^{#Lambda1}", defaultnbins, 0, 1),
        Discriminant("D_HadZH_0plus", "D_{ZH}^{0+}", defaultnbins, 0, 1),
        Discriminant("D_HadZH_0minus", "D_{ZH}^{0-}", defaultnbins, 0, 1),
        Discriminant("D_HadZH_a2", "D_{ZH}^{0h+}", defaultnbins, 0, 1),
        Discriminant("D_HadZH_L1", "D_{ZH}^{#Lambda1}", defaultnbins, 0, 1),
        Discriminant("D_HadWH_0plus", "D_{WH}^{0+}", defaultnbins, 0, 1),
        Discriminant("D_HadWH_0minus", "D_{WH}^{0-}", defaultnbins, 0, 1),
        Discriminant("D_HadWH_a2", "D_{WH}^{0h+}", defaultnbins, 0, 1),
        Discriminant("D_HadWH_L1", "D_{WH}^{#Lambda1}", defaultnbins, 0, 1),
    ] + [
        Discriminant("ZZPt", "p_{T}^{ZZ}", defaultnbins, 0, 500),
    ]
}

del minmax_g1jgik

def discriminant(name):
    if isinstance(name, Discriminant): return name
    return discriminants[name]
