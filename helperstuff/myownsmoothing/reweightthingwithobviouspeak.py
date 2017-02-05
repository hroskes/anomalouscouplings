import ROOT

from rootoverloads import histogramaxisnumbers

def Integral(hist3d, binx1, binx2, biny1, biny2, binz1, binz2, option=""):
    #wrapper to allow kwargs
    return hist3d.Integral(binx1, binx2, biny1, biny2, binz1, binz2, option)
def SetBinContent(hist3d, binx, biny, binz, content):
    #wrapper to allow kwargs
    return hist3d.SetBinContent(binx, biny, binz, content)

def reweightthingwithobviouspeak(hsmooth, rawprojections, axes=[], axesleft=[], axesright=[]):
    newh = h.Clone("h_reweightthingwithobviouspeak")

    for axis in axesleft, axesright:
        if axis not in axes:
            raise ValueError("{} in axesleft or axesright but not in axes".format(axis))

    for axis in axes:
        leftpeak = axis in axesleft
        rightpeak = axis in axesright
        if not leftpeak and not rightpeak: raise ValueError("{} in axes but not in axesleft or axesright".format(axis))

        proj = hsmooth.Projection(abs(axis))
        raw = rawprojections

        content = raw.GetBinContent
        error = raw.GetBinError

        nbins = h.GetNbinsX()

        if leftpeak and content(1) < 2*content(2): raise ValueError("Histogram does not have an obvious left peak! {} {}".format(content(1), content(2)))
        if rightpeak and content(nbins) < 2*content(nbins-1): raise ValueError("Histogram does not have an obvious right peak! {} {}".format(content(nbins), content(nbins-1)))

        tailrange = [1, nbins]

        if leftpeak:
            tailrange[0] = 2
            while content(tailrange[0])-error(tailrange[0]) > content(tailrange[0]+1)+error(tailrange[0]+1):
                tailrange[0] += 1

        if rightpeak:
            tailrange[1] = nbins-1
            while content(tailrange[1])-error(tailrange[1]) > content(tailrange[1]-1)+error(tailrange[1]-1):
                tailrange[1] -= 1

        leftpeakintegral = rightpeakintegral = 0
        if leftpeak:
            leftpeakintegral = raw.Integral(1, tailrange[0]-1)
        tailintegral = raw.Integral(*tailrange)
        smoothtailintegral = proj.Integral(*tailrange)
        if rightpeak:
            leftpeakintegral = raw.Integral(tailrange[1]+1, nbins)

        setcontent = OrderedDict()
        for i in range(1, tailrange[0]):
            setcontent[i] = raw.GetBinContent(i)
        for i in range(tailrange[0], tailrange[1]+1):
            setcontent[i] = proj.GetBinContent(i) * tailintegral / smoothtailintegral
        for i in range(tailrange[1]+1, nbins+1):
            setcontent[i] = raw.GetBinContent(i)

        thisaxis = "xyz".index(axis)
        otheraxes = [axisname for axisname in "xyz" if otheraxis != thisaxis]
        assert len(otheraxes) == 2

        for otheraxisbins in izip(xrange(getattr(h, "GetNbins{}".format(_.upper()))) for _ in otheraxes):
            integralkwargs = {}
            setbincontentkwargs = {}
            for otheraxis, otheraxisbin in zip(otheraxes, otheraxisbins):
                integralkwargs["bin{}1".format(otheraxis)] = \
                integralkwargs["bin{}2".format(otheraxis)] = \
                setbincontentkwargs["bin{}".format(otheraxis)] = otheraxisbin

                integralkwargs["bin{}1".format(thisaxis)] = integralkwargs["bin{}2".format(thisaxis)] = 0

                sliceintegral = Integral(h, **integralkwargs)
                sliceintegral_beforenormalization = sum(setcontent.values())
                ratio = sliceintegral / sliceintegral_beforenormalization

                for binnumber, value in setcontent.iteritems():
                     setbincontentkwargs["bin{}".format(thisaxis)] = binnumber
                     setbincontentkwargs["content"] = value*ratio

                     SetBinContent(newh, **setbincontentkwargs)

    return newh
