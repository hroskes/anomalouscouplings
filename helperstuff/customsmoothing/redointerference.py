import ROOT

def redointerference(hsmooth, rawprojections, templatesandfactors, mirrorjsn, newf):
    hsmooth.Reset("M")
    for template, factor in templatesandfactors:
        hsmooth.Add(getattr(newf, template.templatename()), factor)

    if mirrorjsn:
        for key, value in mirrorjsn.iteritems():
            if key == "type":
                if value == "rescale": dorescale = True
                elif value == "mirror": dorescale = False
                else: raise ValueError("Unknown type {} in mirrorjsn".format(value))
            elif key == "axis":
                if value != 1: raise ValueError("Unknown axis {} in mirrorjsn".format(value))
            elif key == "antisymmetric":
                antisymmetric = value
            elif key == "factor":
                rescaleby = value
            else:
                raise ValueError("Unknown key {} in mirrorjsn".format(key))

        if dorescale:
            hsmooth.Scale(rescaleby)
        else:
            if antisymmetric: sign = -1
            else: sign = 1

            tmp = hsmooth.Clone("tmp"+hsmooth.GetName())
            for binx in range(1, tmp.GetNbinsX()+1):
                for biny in range(1, tmp.GetNbinsY()+1):
                    for binz in range(1, tmp.GetNbinsZ()+1):
                        hsmooth.SetBinContent(binx, biny, binz,
                                              tmp.GetBinContent(binx, biny, binz)
                                                + sign*tmp.GetBinContent(binx, tmp.GetNbinsY()-biny+1, binz)
                                             )

    return True #do make new control plots
