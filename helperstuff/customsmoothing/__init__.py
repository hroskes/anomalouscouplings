import ROOT

from rootoverloads import histogramaxisnumbers

from reweightthingwithobviouspeak import reweightthingwithobviouspeak

functions = {_.__name__: _ for _ in [reweightthingwithobviouspeak]}

def callfunction(hsmooth, rawprojections, **kwargs):
    name = kwargs["name"].lower()
    del kwargs["name"]
    return functions[name](hsmooth, rawprojections, **kwargs)

cache = []

def customsmoothing(hsmooth, rawprojections, **kwargs):
    print "customsmoothing({!r}, {!r}, **{!r})".format(hsmooth, rawprojections, kwargs)
    callfunction(hsmooth, rawprojections, **kwargs)
    cache.append(hsmooth)
    newcontrolplots = []
    for i, rawprojection in enumerate(rawprojections):
        smoothprojection = hsmooth.Projection(i)
        plotname = "control_{}_projAxis{}_afterCustomSmoothing".format(hsmooth.GetName(), i)
        rawname = plotname+"_raw"
        projname = plotname+"_proj"

        c = ROOT.TCanvas(plotname, plotname, 700, 700)
        rawprojection = rawprojection.Clone(rawname)
        smoothprojection.SetName(projname)

        smoothprojection.SetLineColor(ROOT.kRed)
        smoothprojection.SetLineWidth(2)
        rawprojection.SetLineColor(ROOT.kBlack)
        rawprojection.SetMarkerColor(ROOT.kBlack)
        rawprojection.SetMarkerStyle(20)

        rawprojection.Draw()
        smoothprojection.Draw("hist same")

        newcontrolplots.append(c)
        cache.append(c)
        cache.append(rawprojection)
        cache.append(smoothprojection)

    for _ in cache: print _.GetName()
    return newcontrolplots
