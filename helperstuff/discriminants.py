from collections import namedtuple

Discriminant = namedtuple("Discriminant", "name title bins min max")

minmax_g1jgik = {
    ("VBF", "g4", 1, ""): (-3, 3),
    ("VBF", "g4", 2, ""): (0, 5),
    ("VBF", "g4", 3, ""): (-4, 4),
    ("VBF", "g2", 1, ""): (-2, 3),
    ("VBF", "g2", 2, ""): (0, 5),
    ("VBF", "g2", 3, ""): (-4, 4),
    ("VBF", "g1prime2", 1, ""): (-1, 4),
    ("VBF", "g1prime2", 2, ""): (0, 1.5),
    ("VBF", "g1prime2", 3, ""): (-3, 3),

    ("VBF", "g4", 1, "_prime"): (-1, 1),
    ("VBF", "g4", 2, "_prime"): (0, 1),
    ("VBF", "g4", 3, "_prime"): (-.5, .5),
    ("VBF", "g2", 1, "_prime"): (-.5, 1),
    ("VBF", "g2", 2, "_prime"): (0, 1),
    ("VBF", "g2", 3, "_prime"): (-1, 1),
    ("VBF", "g1prime2", 1, "_prime"): (-.3, 1),
    ("VBF", "g1prime2", 2, "_prime"): (0, .5),
    ("VBF", "g1prime2", 3, "_prime"): (-.7, 1),

    ("ZHh", "g4", 1, ""): (-10, 10),
    ("ZHh", "g4", 2, ""): (0, 10),
    ("ZHh", "g4", 3, ""): (-10, 10),
    ("ZHh", "g2", 1, ""): (-10, 10),
    ("ZHh", "g2", 2, ""): (0, 10),
    ("ZHh", "g2", 3, ""): (-10, 10),
    ("ZHh", "g1prime2", 1, ""): (-10, 10),
    ("ZHh", "g1prime2", 2, ""): (0, 10),
    ("ZHh", "g1prime2", 3, ""): (-10, 10),

    ("ZHh", "g4", 1, "_prime"): (-1, 1),
    ("ZHh", "g4", 2, "_prime"): (0, 1),
    ("ZHh", "g4", 3, "_prime"): (-1, 1),
    ("ZHh", "g2", 1, "_prime"): (-1, 1),
    ("ZHh", "g2", 2, "_prime"): (0, 1),
    ("ZHh", "g2", 3, "_prime"): (-1, 1),
    ("ZHh", "g1prime2", 1, "_prime"): (-1, 1),
    ("ZHh", "g1prime2", 2, "_prime"): (0, 1),
    ("ZHh", "g1prime2", 3, "_prime"): (-1, 1),

    ("WHh", "g4", 1, ""): (-10, 10),
    ("WHh", "g4", 2, ""): (0, 10),
    ("WHh", "g4", 3, ""): (-10, 10),
    ("WHh", "g2", 1, ""): (-10, 10),
    ("WHh", "g2", 2, ""): (0, 10),
    ("WHh", "g2", 3, ""): (-10, 10),
    ("WHh", "g1prime2", 1, ""): (-10, 10),
    ("WHh", "g1prime2", 2, ""): (0, 10),
    ("WHh", "g1prime2", 3, ""): (-10, 10),

    ("WHh", "g4", 1, "_prime"): (-1, 1),
    ("WHh", "g4", 2, "_prime"): (0, 1),
    ("WHh", "g4", 3, "_prime"): (-1, 1),
    ("WHh", "g2", 1, "_prime"): (-1, 1),
    ("WHh", "g2", 2, "_prime"): (0, 1),
    ("WHh", "g2", 3, "_prime"): (-1, 1),
    ("WHh", "g1prime2", 1, "_prime"): (-1, 1),
    ("WHh", "g1prime2", 2, "_prime"): (0, 1),
    ("WHh", "g1prime2", 3, "_prime"): (-1, 1),
}

discriminants = {
    d.name: d for d in [
        Discriminant("D_bkg_0plus", "D_{bkg}", 20, 0, 1),
        Discriminant("D_bkg_0plus_ResUp", "D_{bkg}^{ResUp}", 20, 0, 1),
        Discriminant("D_bkg_0plus_ResDown", "D_{bkg}^{ResDown}", 20, 0, 1),
        Discriminant("D_bkg_0plus_ScaleUp", "D_{bkg}^{ScaleUp}", 20, 0, 1),
        Discriminant("D_bkg_0plus_ScaleDown", "D_{bkg}^{ScaleDown}", 20, 0, 1),
        Discriminant("D_0minus_decay", "D_{0-}^{dec}", 20, 0, 1),
        Discriminant("D_CP_decay", "D_{CP}^{dec}", 20, -0.5, 0.5),
        Discriminant("D_g2_decay", "D_{0h+}^{dec}", 20, 0, 1),
        Discriminant("D_g1g2_decay", "D_{int}^{dec}", 20, 0, 1),
        Discriminant("D_g1prime2_decay", "D_{#Lambda1}^{dec}", 20, 0, 1),
        Discriminant("D_g1g1prime2_decay", "D_{#Lambda1int}^dec", 20, 0, 1),
        Discriminant("D_0minus_VBF", "D_{0-}^{VBF}", 20, 0, 1),
        Discriminant("D_CP_VBF", "D_{CP}^{VBF}", 20, -1, 1),
        Discriminant("D_g2_VBF", "D_{0h+}^{VBF}", 20, 0, 1),
        Discriminant("D_g1g2_VBF", "D_{int}^{VBF}", 20, -1, 1),
        Discriminant("D_g1prime2_VBF", "D_{#Lambda1}^{VBF}", 20, 0, 1),
        Discriminant("D_g1g1prime2_VBF", "D_{#Lambda1int}^{VBF}", 20, -1, 1),
        Discriminant("D_0minus_VBFdecay", "D_{0-}^{VBFdec}", 20, 0, 1),
        Discriminant("D_g2_VBFdecay", "D_{0h+}^{VBFdec}", 20, 0, 1),
        Discriminant("D_g1prime2_VBFdecay", "D_{#Lambda1}^{VBFdec}", 20, 0, 1),
    ] + [
        Discriminant(
                     "D_g1{}_{}{}_VBFdecay{}".format(i, gj, 4-i, prime),
                     "D{}^[VBFdec]_[g_[1]^[{}]g_[{}]^[{}]]".format("'" if prime else "", i, gj, 4-i).replace("[", "{").replace("]", "}"),
                     20,
                     *minmax_g1jgik["VBF", gj, i, prime]
                    )
            for prime in ("", "_prime")
            for gj in ("g4", "g2", "g1prime2")
            for i in range(1, 4)
    ] + [
        Discriminant("D_0minus_ZH_hadronic", "D_{0-}^{ZHh}", 20, 0, 1),
        Discriminant("D_CP_ZH_hadronic", "D_{CP}^{ZHh}", 20, -1, 1),
        Discriminant("D_g2_ZH_hadronic", "D_{0h+}^{ZHh}", 20, 0, 1),
        Discriminant("D_g1g2_ZH_hadronic", "D_{int}^{ZHh}", 20, -1, 1),
        Discriminant("D_g1prime2_ZH_hadronic", "D_{#Lambda1}^{ZHh}", 20, 0, 1),
        Discriminant("D_g1g1prime2_ZH_hadronic", "D_{#Lambda1int}^{ZHh}", 20, -1, 1),
        Discriminant("D_0minus_ZHdecay_hadronic", "D_{0-}^{ZHhdec}", 20, 0, 1),
        Discriminant("D_g2_ZHdecay_hadronic", "D_{0h+}^{ZHhdec}", 20, 0, 1),
        Discriminant("D_g1prime2_ZHdecay_hadronic", "D_{#Lambda1}^{ZHhdec}", 20, 0, 1),
    ] + [
        Discriminant(
                     "D_g1{}_{}{}_ZHdecay_hadronic{}".format(i, gj, 4-i, prime),
                     "D{}^[ZHhdec]_[g_[1]^[{}]g_[{}]^[{}]]".format("'" if prime else "", i, gj, 4-i).replace("[", "{").replace("]", "}"),
                     20,
                     *minmax_g1jgik["ZHh", gj, i, prime]
                    )
            for prime in ("", "_prime")
            for gj in ("g4", "g2", "g1prime2")
            for i in range(1, 4)
    ] + [
        Discriminant("D_0minus_WH_hadronic", "D_{0-}^{WHh}", 20, 0, 1),
        Discriminant("D_CP_WH_hadronic", "D_{CP}^{WHh}", 20, -1, 1),
        Discriminant("D_g2_WH_hadronic", "D_{0h+}^{WHh}", 20, 0, 1),
        Discriminant("D_g1g2_WH_hadronic", "D_{int}^{WHh}", 20, -1, 1),
        Discriminant("D_g1prime2_WH_hadronic", "D_{#Lambda1}^{WHh}", 20, 0, 1),
        Discriminant("D_g1g1prime2_WH_hadronic", "D_{#Lambda1int}^{WHh}", 20, -1, 1),
        Discriminant("D_0minus_WHdecay_hadronic", "D_{0-}^{WHhdec}", 20, 0, 1),
        Discriminant("D_g2_WHdecay_hadronic", "D_{0h+}^{WHhdec}", 20, 0, 1),
        Discriminant("D_g1prime2_WHdecay_hadronic", "D_{#Lambda1}^{WHhdec}", 20, 0, 1),
    ] + [
        Discriminant(
                     "D_g1{}_{}{}_WHdecay_hadronic{}".format(i, gj, 4-i, prime),
                     "D{}^[WHhdec]_[g_[1]^[{}]g_[{}]^[{}]]".format("'" if prime else "", i, gj, 4-i).replace("[", "{").replace("]", "}"),
                     20,
                     *minmax_g1jgik["WHh", gj, i, prime]
                    )
            for prime in ("", "_prime")
            for gj in ("g4", "g2", "g1prime2")
            for i in range(1, 4)
    ]
}

del minmax_g1jgik

def discriminant(name):
    return discriminants[name]
