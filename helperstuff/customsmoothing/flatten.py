from collections import OrderedDict

import ROOT

def Integral(hist3d, binx1, binx2, biny1, biny2, binz1, binz2, option=""):
    #wrapper to allow kwargs
    return hist3d.Integral(binx1, binx2, biny1, biny2, binz1, binz2, option)
def SetBinContent(hist3d, binx, biny, binz, content):
    #wrapper to allow kwargs
    return hist3d.SetBinContent(binx, biny, binz, content)

def flatten(hsmooth, rawprojections, axes=[]):
    for axis in axes:
        proj = hsmooth.Projection(abs(axis))
        raw = rawprojections[axis]

        nbins = raw.GetNbinsX()
        assert nbins == proj.GetNbinsX()

        thisaxis = "xyz"[axis]
        otheraxes = [axisname for axisname in "xyz" if axisname != thisaxis]
        assert len(otheraxes) == 2

        otheraxisbins = [None, None]

        for thisaxisbins in xrange(1, nbins):
           for otheraxisbins[0] in xrange(1, getattr(hsmooth, "GetNbins{}".format(otheraxes[0].upper()))()+1):
              for otheraxisbins[1] in xrange(1, getattr(hsmooth, "GetNbins{}".format(otheraxes[1].upper()))()+1):
                integralkwargs = {}
                setbincontentkwargs = {}
                for otheraxis, otheraxisbin in zip(otheraxes, otheraxisbins):
                    integralkwargs["bin{}1".format(otheraxis)] = \
                    integralkwargs["bin{}2".format(otheraxis)] = \
                    setbincontentkwargs["bin{}".format(otheraxis)] = otheraxisbin

                integralkwargs["bin{}1".format(thisaxis)] = integralkwargs["bin{}2".format(thisaxis)] = -1

                sliceintegral = Integral(hsmooth, **integralkwargs)

                for binnumber, value in setcontent.iteritems():
                    setbincontentkwargs["bin{}".format(thisaxis)] = binnumber
                    setbincontentkwargs["content"] = sliceintegral / nbins

                    SetBinContent(hsmooth, **setbincontentkwargs)

    return True #do make new control plots
