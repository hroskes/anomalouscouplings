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
        Discriminant("D_bkg_0plus", "D_{bkg}", defaultnbins, 0, 1),
        Discriminant("D_bkg_0plus_ResUp", "D_{bkg}^{ResUp}", defaultnbins, 0, 1),
        Discriminant("D_bkg_0plus_ResDown", "D_{bkg}^{ResDown}", defaultnbins, 0, 1),
        Discriminant("D_bkg_0plus_ScaleUp", "D_{bkg}^{ScaleUp}", defaultnbins, 0, 1),
        Discriminant("D_bkg_0plus_ScaleDown", "D_{bkg}^{ScaleDown}", defaultnbins, 0, 1),
        Discriminant("D_0minus_decay", "D_{0-}^{dec}", defaultnbins, 0, 1),
        Discriminant("D_CP_decay", "D_{CP}^{dec}", defaultnbins, -0.5, 0.5),
        Discriminant("D_g2_decay", "D_{0h+}^{dec}", defaultnbins, 0, 1),
        Discriminant("D_g1g2_decay", "D_{int}^{dec}", defaultnbins, 0, 1),
        Discriminant("D_g1prime2_decay", "D_{#Lambda1}^{dec}", defaultnbins, 0, 1),
        Discriminant("D_g1g1prime2_decay", "D_{#Lambda1int}^dec", defaultnbins, 0, 1),
        Discriminant("D_0minus_VBF", "D_{0-}^{VBF}", defaultnbins, 0, 1),
        Discriminant("D_CP_VBF", "D_{CP}^{VBF}", defaultnbins, -1, 1),
        Discriminant("D_g2_VBF", "D_{0h+}^{VBF}", defaultnbins, 0, 1),
        Discriminant("D_g1g2_VBF", "D_{int}^{VBF}", defaultnbins, -1, 1),
        Discriminant("D_g1prime2_VBF", "D_{#Lambda1}^{VBF}", defaultnbins, 0, 1),
        Discriminant("D_g1g1prime2_VBF", "D_{#Lambda1int}^{VBF}", defaultnbins, -1, 1),
        Discriminant("D_0minus_VBFdecay", "D_{0-}^{VBFdec}", defaultnbins, 0, 1),
        Discriminant("D_g2_VBFdecay", "D_{0h+}^{VBFdec}", defaultnbins, 0, 1),
        Discriminant("D_g1prime2_VBFdecay", "D_{#Lambda1}^{VBFdec}", defaultnbins, 0, 1),
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
        Discriminant("D_0minus_ZH_hadronic", "D_{0-}^{ZH}", defaultnbins, 0, 1),
        Discriminant("D_CP_ZH_hadronic", "D_{CP}^{ZH}", defaultnbins, -.4, .4),
        Discriminant("D_g2_ZH_hadronic", "D_{0h+}^{ZH}", defaultnbins, 0, 1),
        Discriminant("D_g1g2_ZH_hadronic", "D_{int}^{ZH}", defaultnbins, -1, .05),
        Discriminant("D_g1prime2_ZH_hadronic", "D_{#Lambda1}^{ZH}", defaultnbins, 0, 1),
        Discriminant("D_g1g1prime2_ZH_hadronic", "D_{#Lambda1int}^{ZH}", defaultnbins, 0, 1),
        Discriminant("D_0minus_ZHdecay_hadronic", "D_{0-}^{ZHdec}", defaultnbins, 0, 1),
        Discriminant("D_g2_ZHdecay_hadronic", "D_{0h+}^{ZHdec}", defaultnbins, 0, 1),
        Discriminant("D_g1prime2_ZHdecay_hadronic", "D_{#Lambda1}^{ZHdec}", defaultnbins, 0, 1),
#    ] + [
#        Discriminant(
#                     "D_g1{}_{}{}_ZHdecay_hadronic{}".format(i, gj, 4-i, prime),
#                     "D{}^[ZHdec]_[g_[1]^[{}]g_[{}]^[{}]]".format("'" if prime else "", i, gj, 4-i).replace("[", "{").replace("]", "}"),
#                     defaultnbins,
#                     *minmax_g1jgik["ZH", gj, i, prime]
#                    )
#            for prime in ("", "_prime")
#            for gj in ("g4", "g2", "g1prime2")
#            for i in range(1, 4)
    ] + [
        Discriminant("D_0minus_WH_hadronic", "D_{0-}^{WH}", defaultnbins, 0, 1),
        Discriminant("D_CP_WH_hadronic", "D_{CP}^{WH}", defaultnbins, -.1, .1),
        Discriminant("D_g2_WH_hadronic", "D_{0h+}^{WH}", defaultnbins, 0, 1),
        Discriminant("D_g1g2_WH_hadronic", "D_{int}^{WH}", defaultnbins, -1, .05),
        Discriminant("D_g1prime2_WH_hadronic", "D_{#Lambda1}^{WH}", defaultnbins, 0, 1),
        Discriminant("D_g1g1prime2_WH_hadronic", "D_{#Lambda1int}^{WH}", defaultnbins, 0, 1),
        Discriminant("D_0minus_WHdecay_hadronic", "D_{0-}^{WHdec}", defaultnbins, 0, 1),
        Discriminant("D_g2_WHdecay_hadronic", "D_{0h+}^{WHdec}", defaultnbins, 0, 1),
        Discriminant("D_g1prime2_WHdecay_hadronic", "D_{#Lambda1}^{WHdec}", defaultnbins, 0, 1),
#    ] + [
#        Discriminant(
#                     "D_g1{}_{}{}_WHdecay_hadronic{}".format(i, gj, 4-i, prime),
#                     "D{}^[WHdec]_[g_[1]^[{}]g_[{}]^[{}]]".format("'" if prime else "", i, gj, 4-i).replace("[", "{").replace("]", "}"),
#                     defaultnbins,
#                     *minmax_g1jgik["WH", gj, i, prime]
#                    )
#            for prime in ("", "_prime")
#            for gj in ("g4", "g2", "g1prime2")
#            for i in range(1, 4)
    ] + [
        Discriminant("D_0minus_VH_hadronic", "D_{0-}^{VH}", defaultnbins, 0, 1),
        Discriminant("D_CP_VH_hadronic", "D_{CP}^{VH}", defaultnbins, -.4, .4),
        Discriminant("D_g2_VH_hadronic", "D_{0h+}^{VH}", defaultnbins, 0, 1),
        Discriminant("D_g1g2_VH_hadronic", "D_{int}^{VH}", defaultnbins, -1, .05),
        Discriminant("D_g1prime2_VH_hadronic", "D_{#Lambda1}^{VH}", defaultnbins, 0, 1),
        Discriminant("D_g1g1prime2_VH_hadronic", "D_{#Lambda1int}^{VH}", defaultnbins, 0, 1),
        Discriminant("D_0minus_VHdecay_hadronic", "D_{0-}^{VHdec}", defaultnbins, 0, 1),
        Discriminant("D_g2_VHdecay_hadronic", "D_{0h+}^{VHdec}", defaultnbins, 0, 1),
        Discriminant("D_g1prime2_VHdecay_hadronic", "D_{#Lambda1}^{VHdec}", defaultnbins, 0, 1),
    ] + [
        Discriminant("ZZPt", "p_{T}^{ZZ}", defaultnbins, 0, 500)
    ]
}

del minmax_g1jgik

def discriminant(name):
    return discriminants[name]
