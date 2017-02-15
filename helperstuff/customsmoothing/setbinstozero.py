import ROOT

def SetBinContent(hist3d, binx, biny, binz, content):
    #wrapper to allow kwargs
    return hist3d.SetBinContent(binx, biny, binz, content)

def setbinstozero(hsmooth, rawprojections, nbinsonleft, nbinsonright):
    integral = hsmooth.Integral()
    for axis in (nbinsonleft.keys()+nbinsonright.keys()):
        nbins = rawprojections[axis].GetNbinsX()

        #http://stackoverflow.com/a/9359011/5228524
        nleft = nbinsonleft.get(axis, 0)
        nright = nbinsonright.get(axis, 0)
        binsto0 = [1+_ for _ in range(nleft)] + [nbins-_ for _ in range(nright)]

        thisaxis = "xyz"[axis]
        otheraxes = [axisname for axisname in "xyz" if axisname != thisaxis]
        assert len(otheraxes) == 2

        kwargs = {"binx": None, "biny": None, "binz": None, "content": 0}

        for kwargs["bin"+otheraxes[0]] in xrange(1, getattr(hsmooth, "GetNbins{}".format(otheraxes[0].upper()))()+1):
            for kwargs["bin"+otheraxes[1]] in xrange(1, getattr(hsmooth, "GetNbins{}".format(otheraxes[1].upper()))()+1):
                for kwargs["bin"+thisaxis] in binsto0:
                    hsmooth.SetBinContent(hsmooth, **kwargs)

    hsmooth.Scale(integral / hsmooth.Integral())

    return True #do make new control plots
