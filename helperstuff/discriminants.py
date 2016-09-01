from collections import namedtuple

Discriminant = namedtuple("Discriminant", "name title bins min max")

discriminants = {
    d.name: d for d in [
        Discriminant("D_bkg_0plus", "D_{bkg}", 50, 0, 1),
        Discriminant("D_bkg_0plus_ResUp", "D_{bkg}^{ResUp}", 50, 0, 1),
        Discriminant("D_bkg_0plus_ResDown", "D_{bkg}^{ResDown}", 50, 0, 1),
        Discriminant("D_bkg_0plus_ScaleUp", "D_{bkg}^{ScaleUp}", 50, 0, 1),
        Discriminant("D_bkg_0plus_ScaleDown", "D_{bkg}^{ScaleDown}", 50, 0, 1),
        Discriminant("D_0minus_decay", "D_{0-}^{dec}", 50, 0, 1),
        Discriminant("D_CP_decay", "D_{CP}^{dec}", 50, -0.5, 0.5),
        Discriminant("D_g2_decay", "D_{0h+}^{dec}", 50, 0, 1),
        Discriminant("D_g1g2_decay", "D_{int}^{dec}", 50, 0, 1),
        Discriminant("D_g1prime2_decay", "D_{#Lambda1}^{dec}", 50, 0, 1),
        Discriminant("D_g1g1prime2_decay", "D_{#Lambda1int}^dec", 50, 0, 1),
        Discriminant("D_0minus_VBF", "D_{0-}^{VBF}", 50, 0, 1),
        Discriminant("D_CP_VBF", "D_{CP}^{VBF}", 50, 0, 1),
        Discriminant("D_g2_VBF", "D_{0h+}^{VBF}", 50, 0, 1),
        Discriminant("D_g1g2_VBF", "D_{int}^{VBF}", 50, 0, 1),
        Discriminant("D_g1prime2_VBF", "D_{#Lambda1}^{VBF}", 50, 0, 1),
        Discriminant("D_g1g1prime2_VBF", "D_{#Lambda1int}^{VBF}", 50, 0, 1),
        Discriminant("D_0minus_VBFdecay", "D_{0-}^{VBFdec}", 50, 0, 1),
        Discriminant("D_g2_VBFdecay", "D_{0h+}^{VBFdec}", 50, 0, 1),
        Discriminant("D_g1prime2_VBFdecay", "D_{#Lambda1}^{VBFdec}", 50, 0, 1),
    ] + [
        Discriminant(
                     "D_g1{}_{}{}_VBFdecay{}".format(i, gj, 4-i, prime),
                     "D{}_[g_[1]^[{}]g_[{}]^[{}]]".format("'" if prime else "", i, gj, 4-i).replace("[", "{").replace("]", "}"),
                     50,
                     0 if i%2==0 else -1,
                     1
                    )
            for prime in ("", "_prime")
            for gj in ("g4", "g2", "g1prime2")
            for i in range(5)
    ]
}

def discriminant(name):
    return discriminants[name]
