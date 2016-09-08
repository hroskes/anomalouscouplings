from collections import namedtuple

Discriminant = namedtuple("Discriminant", "name title bins min max")

minmax_g1jgik = {
    ("g4", 0, ""): (0, 1),
    ("g4", 4, ""): (0, 1),
    ("g2", 0, ""): (0, 1),
    ("g2", 4, ""): (0, 1),
    ("g1prime2", 0, ""): (0, 1),
    ("g1prime2", 4, ""): (0, 1),
    ("g4", 0, "_prime"): (0, 1),
    ("g4", 4, "_prime"): (0, 1),
    ("g2", 0, "_prime"): (0, 1),
    ("g2", 4, "_prime"): (0, 1),
    ("g1prime2", 0, "_prime"): (0, 1),
    ("g1prime2", 4, "_prime"): (0, 1),


    ("g4", 1, ""): (-3, 3),
    ("g4", 2, ""): (0, 5),
    ("g4", 3, ""): (-4, 4),
    ("g2", 1, ""): (-2, 3),
    ("g2", 2, ""): (0, 5),
    ("g2", 3, ""): (-4, 4),
    ("g1prime2", 1, ""): (-1, 4),
    ("g1prime2", 2, ""): (0, 1.5),
    ("g1prime2", 3, ""): (-3, 3),

    ("g4", 1, "_prime"): (-1, 1),
    ("g4", 2, "_prime"): (0, 1),
    ("g4", 3, "_prime"): (-.5, .5),
    ("g2", 1, "_prime"): (-.5, 1),
    ("g2", 2, "_prime"): (0, 1),
    ("g2", 3, "_prime"): (-1, 1),
    ("g1prime2", 1, "_prime"): (-.3, 1),
    ("g1prime2", 2, "_prime"): (0, .5),
    ("g1prime2", 3, "_prime"): (-.7, 1),
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
                     "D{}_[g_[1]^[{}]g_[{}]^[{}]]".format("'" if prime else "", i, gj, 4-i).replace("[", "{").replace("]", "}"),
                     20,
                     *minmax_g1jgik[gj, i, prime]
                    )
            for prime in ("", "_prime")
            for gj in ("g4", "g2", "g1prime2")
            for i in range(5)
    ]
}

del minmax_g1jgik

def discriminant(name):
    return discriminants[name]
