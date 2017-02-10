from collections import OrderedDict
from itertools import izip

import ROOT

def Integral(hist3d, binx1, binx2, biny1, biny2, binz1, binz2, option=""):
    #wrapper to allow kwargs
    return hist3d.Integral(binx1, binx2, biny1, biny2, binz1, binz2, option)
def SetBinContent(hist3d, binx, biny, binz, content):
    #wrapper to allow kwargs
    return hist3d.SetBinContent(binx, biny, binz, content)

def reweightthingwithobviouspeak(hsmooth, rawprojections, axes=[], axesleft=[], axesright=[]):
    for axis in axesleft+axesright:
        if axis not in axes:
            raise ValueError("{} in axesleft or axesright but not in axes".format(axis))

    for axis in axes:
        leftpeak = axis in axesleft
        rightpeak = axis in axesright
        if not leftpeak and not rightpeak: raise ValueError("{} in axes but not in axesleft or axesright".format(axis))

        proj = hsmooth.Projection(abs(axis))
        raw = rawprojections[axis]

        content = raw.GetBinContent
        error = raw.GetBinError

        nbins = raw.GetNbinsX()
        assert nbins == proj.GetNbinsX()

        if leftpeak and (
                         content(1) < 1.5*content(2)
                         and (content(1) < content(2) or content(2) < 1.5*content(3))
                        ): raise ValueError("Histogram does not have an obvious left peak! {} {} {}".format(content(1), content(2), content(3)))
        if rightpeak and (
                          content(nbins) < 1.5*content(nbins-1)
                          and (content(nbins) < content(nbins-1) or content(nbins-1) < 1.5*content(nbins-2))
                         ): raise ValueError("Histogram does not have an obvious right peak! {} {} {}".format(content(nbins), content(nbins-1), content(nbins-2)))

        tailrange = [1, nbins]

        if leftpeak:
            tailrange[0] = 2
            while (content(tailrange[0])-error(tailrange[0]) > content(tailrange[0]+1)+error(tailrange[0]+1)
                   or content(tailrange[0]) > content(tailrange[0]+1) > content(tailrange[0]+2) > content(tailrange[0]+3)
                  ) and content(tailrange[0]+1) > content(1) / 4:
                tailrange[0] += 1

        if rightpeak:
            tailrange[1] = nbins-1
            while (content(tailrange[1])-error(tailrange[1]) > content(tailrange[1]-1)+error(tailrange[1]-1)
                   or content(tailrange[1]) > content(tailrange[1]-1) > content(tailrange[1]-2) > content(tailrange[1]-3)
                  ) and content(tailrange[1]-1) > content(nbins) / 4:
                tailrange[1] -= 1

        rawleftpeakintegral = rawrightpeakintegral = smoothleftpeakintegral = smoothrightpeakintegral = 0
        rawintegral = raw.Integral()
        smoothintegral = proj.Integral()
        if leftpeak:
            rawleftpeakintegral = raw.Integral(1, tailrange[0]-1)
            smoothleftpeakintegral = proj.Integral(1, tailrange[0]-1)
        rawtailintegral = raw.Integral(*tailrange)
        smoothtailintegral = proj.Integral(*tailrange)
        if rightpeak:
            rawrightpeakintegral = raw.Integral(tailrange[1]+1, nbins)
            smoothrightpeakintegral = raw.Integral(tailrange[1]+1, nbins)

        setcontent = OrderedDict()
        for i in range(1, tailrange[0]):
            setcontent[i] = raw.GetBinContent(i)
        for i in range(tailrange[0], tailrange[1]+1):
            setcontent[i] = (
                             proj.GetBinContent(i)
                                    * rawtailintegral / smoothtailintegral
                            )
        for i in range(tailrange[1]+1, nbins+1):
            setcontent[i] = raw.GetBinContent(i)

        thisaxis = "xyz"[axis]
        otheraxes = [axisname for axisname in "xyz" if axisname != thisaxis]
        assert len(otheraxes) == 2

        otheraxisbins = [None, None]

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
                sliceintegral_beforenormalization = sum(setcontent.values())
                ratio = sliceintegral / sliceintegral_beforenormalization

                for binnumber, value in setcontent.iteritems():
                    setbincontentkwargs["bin{}".format(thisaxis)] = binnumber
                    setbincontentkwargs["content"] = value*ratio

                    SetBinContent(hsmooth, **setbincontentkwargs)

    return True #do make new control plots
