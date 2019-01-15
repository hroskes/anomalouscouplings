import ROOT

from rootoverloads import histogramaxisnumbers

from flatten import flatten
from justcopy import justcopy
from redointerference import redointerference
from reweightthingwithobviouspeak import reweightthingwithobviouspeak
from setbinstozero import setbinstozero
from seterrorforfloor import seterrorforfloor
from seterrortozero import seterrortozero
from useDbkgorthogonal import useDbkgorthogonal

def donothing(hsmooth, rawprojections, **kwargs):
  """do nothing"""
  return False #don't make new control plots

def multiplefunctions(hsmooth, rawprojections, listofkwargs):
  makenewcontrolplots = False
  for kwargs in listofkwargs:
    makenewcontrolplots = callfunction(hsmooth, rawprojections, **kwargs) or makenewcontrolplots  #don't reverse this!  it short circuits
  return makenewcontrolplots

functions = {_.__name__.lower(): _ for _ in [donothing, redointerference, reweightthingwithobviouspeak, setbinstozero, flatten, justcopy, useDbkgorthogonal, multiplefunctions, seterrorforfloor, seterrortozero]}

def callfunction(hsmooth, rawprojections, **kwargs):
    if "name" in kwargs:
        name = kwargs.pop("name").lower()
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

        if rawprojection.Integral() and smoothprojection.Integral():
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
