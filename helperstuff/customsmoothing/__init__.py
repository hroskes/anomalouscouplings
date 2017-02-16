import ROOT

from rootoverloads import histogramaxisnumbers

from redointerference import redointerference
from reweightthingwithobviouspeak import reweightthingwithobviouspeak
from setbinstozero import setbinstozero

def donothing(hsmooth, rawprojections, **kwargs):
   """do nothing"""
   return False #don't make new control plots

functions = {_.__name__: _ for _ in [donothing, redointerference, reweightthingwithobviouspeak, setbinstozero]}

def callfunction(hsmooth, rawprojections, **kwargs):
    if "name" in kwargs:
        name = kwargs["name"].lower()
        del kwargs["name"]
    else:
        if kwargs: raise ValueError("No function name given! {}".format(kwargs))
        name = "donothing"
    return functions[name](hsmooth, rawprojections, **kwargs)

def customsmoothing(hsmooth, rawprojections, templatedirectory, controlplotsdirectory, **kwargs):
    makenewcontrolplots = callfunction(hsmooth, rawprojections, **kwargs)
    templatedirectory.cd()
    hsmooth.Write()

    if not makenewcontrolplots: return

    newcontrolplots = []

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
