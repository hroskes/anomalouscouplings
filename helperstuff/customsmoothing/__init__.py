import ROOT

from rootoverloads import histogramaxisnumbers

from reweightthingwithobviouspeak import reweightthingwithobviouspeak

def donothing(hsmooth, rawprojections, **kwargs):
   """do nothing"""

functions = {_.__name__: _ for _ in [reweightthingwithobviouspeak, donothing]}

def callfunction(hsmooth, rawprojections, **kwargs):
    if "name" in kwargs:
        name = kwargs["name"].lower()
        del kwargs["name"]
    else:
        if kwargs: raise ValueError("No function name given! {}".format(kwargs))
        name = "donothing"
    return functions[name](hsmooth, rawprojections, **kwargs)

def customsmoothing(hsmooth, rawprojections, templatedirectory, controlplotsdirectory, **kwargs):
    callfunction(hsmooth, rawprojections, **kwargs)
    newcontrolplots = []
    templatedirectory.cd()
    hsmooth.Write()
    for i, rawprojection in enumerate(rawprojections):
        smoothprojection = hsmooth.Projection(i)
        plotname = "control_{}_projAxis{}_afterCustomSmoothing".format(hsmooth.GetName(), i)
        rawname = plotname+"_raw"
        projname = plotname+"_proj"

        c = ROOT.TCanvas(plotname, plotname, 700, 700)
        rawprojection = rawprojection.Clone(rawname)
        smoothprojection.SetName(projname)

        rawprojection.Scale(smoothprojection.Integral() / rawprojection.Integral())

        smoothprojection.SetLineColor(ROOT.kRed)
        smoothprojection.SetLineWidth(2)
        rawprojection.SetLineColor(ROOT.kBlack)
        rawprojection.SetMarkerColor(ROOT.kBlack)
        rawprojection.SetMarkerStyle(20)

        rawprojection.Draw()
        smoothprojection.Draw("hist same")

        controlplotsdirectory.cd()
        c.Write()
