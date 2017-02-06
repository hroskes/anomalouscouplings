from rootoverloads import histogramaxisnumbers

from reweightthingwithobviouspeak import reweightthingwithobviouspeak

functions = {_.__name__: _ for _ in [reweightthingwithobviouspeak]}

def callfunction(hsmooth, rawprojections, **kwargs):
    name = kwargs["name"].lower()
    del kwargs["name"]
    return functions[name](*args, **kwargs)

cache = []

def customsmoothing(hsmooth, rawprojections, **kwargs):
    newh = callfunction(**kwargs)
    cache.append(newh)
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

        newcontrolplots.Append(c)
        cache.append(rawprojection)
        cache.append(smoothprojection)

    return newh, newcontrolplots
