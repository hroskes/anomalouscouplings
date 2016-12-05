#!/usr/bin/env python

assert __name__ == "__main__"

from helperstuff import config
from helperstuff import constants
from helperstuff import style
from helperstuff.combinehelpers import gettemplate
from helperstuff.discriminants import discriminant
from helperstuff.enums import *
import helperstuff.rootoverloads.histogramfloor
from helperstuff.samples import *
import re
import ROOT
import os

#========================
#inputs
#weight, bins, min, max, category can be None
productionmode = "ggH"
disc           = "D_g2_decay"
weight         = samplewithfai(productionmode, "fa2", -.14, productionmodeforfai="ggH")
bins           = None
min            = None
max            = None
enrich         = False
masscut        = True
channel        = "2e2mu"
category       = "Untagged"
#========================

if isinstance(weight, SampleBase):
    weight = weight.MC_weight

if weight is not None:
    for name, value in constants.__dict__.iteritems():
        try:
            value = "{:g}".format(value)
        except ValueError:
            continue
        weight = re.sub(r"\b"+name+r"\b", value, weight)
        disc = re.sub(r"\b"+name+r"\b", value, disc)


hstack = ROOT.THStack()
legend = ROOT.TLegend(.6, .5, .9, .9)
cache = []

ROOT.gErrorIgnoreLevel = ROOT.kError

discname, title, discbins, discmin, discmax = discriminant(disc)
if bins is None:
    bins = discbins
if min is None:
    min = discmin
if max is None:
    max = discmax

c = ROOT.TCanvas()
hs = {}

def hypothesestouse():
    for hypothesis in hypotheses:
        if productionmode == "ggH" and hypothesis not in decayonlyhypotheses: continue
        if productionmode in ["VBF", "ZH", "WH"] and hypothesis not in prodonlyhypotheses: continue
        if productionmode in ["HJJ", "ttH"] and hypothesis not in hffhypotheses: continue
        yield hypothesis

for hypothesis in hypotheses:
    if hypothesis not in hypothesestouse():
        try:
            os.remove(os.path.join(config.plotsbasedir, "TEST", "{}.png".format(hypothesis)))
        except OSError:
            pass

for color, hypothesis in enumerate(hypothesestouse(), start=1):
    t = ROOT.TChain("candTree", "candTree")
    sample = Sample(productionmode, hypothesis, "160928")
    t.Add(sample.withdiscriminantsfile())
    hname = "h{}".format(hypothesis)

    weightname = weight if weight is not None else sample.weightname()

    weightfactors = [
                     weightname,
                     "{}>-998".format(discname),
                    ]

    if category is not None:
        weightfactors.append(" || ".join("(category=={})".format(_) for _ in Category(category).idnumbers))
    if enrich:
        weightfactors.append("D_bkg_0plus>0.5")
    if channel is not None:
        weightfactors.append("Z1Flav*Z2Flav=={}".format(Channel(channel).ZZFlav))
    if masscut:
        weightfactors.append("ZZMass > {}".format(config.m4lmin))
        weightfactors.append("ZZMass < {}".format(config.m4lmax))

    wt = "*".join("("+_+")" for _ in weightfactors)

    t.Draw("{}>>{}({},{},{})".format(discname, hname, bins, min, max), wt, "hist")
    h = hs[hypothesis] = getattr(ROOT, hname)
    h.GetXaxis().SetTitle(title)
    if isinstance(h, ROOT.TH1) and not isinstance(h, ROOT.TH2):
      h.SetBinContent(h.GetNbinsX(), h.GetBinContent(h.GetNbinsX()+1) + h.GetBinContent(h.GetNbinsX()))
      h.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1))
    h.Floor()
    hstack.Add(h)
    cache.append(h)
    legend.AddEntry(h, str(hypothesis), "l")
    h.SetLineColor(color)
    h.SetMarkerStyle(1)
    print "{:10} {:.3g}".format(hypothesis, h.Integral())
    try:
        os.makedirs(os.path.join(config.plotsbasedir, "TEST"))
    except OSError:
        pass
    c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "{}.png".format(hypothesis)))

hstack.Draw("histnostack")
hstack.GetXaxis().SetTitle(title)
c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "test.png"))

"""
hint = hs["0+"].Clone("hint")
hint.Add(hs["0-"])
hint.Scale(-.5)
hint.Add(hs["fa3prod0.5"])

hint.Draw("hist")
#c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "a1a3int.png"))
"""
