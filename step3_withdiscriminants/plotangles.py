#!/usr/bin/env python
from collections import namedtuple
from math import pi
import ROOT, rootoverloads
import style

t = ROOT.TChain("candTree")
t.Add("ggH0+170119.root")

Hypothesis = namedtuple("Hypothesis", "name weightname color")
Angle = namedtuple("Angle", "name title bins min max")

hypotheses = (
              Hypothesis("0^{+}", "MC_weight_ggH_g1", 1),
              None,
              Hypothesis("0^{-}", "MC_weight_ggH_g4", 2),
              Hypothesis("f_{a3}=0.5", "MC_weight_ggH_g1g4", 6),
              Hypothesis("0_{h}^{+}", "MC_weight_ggH_g2", 4),
              Hypothesis("f_{a2}=0.5", "MC_weight_ggH_g1g2", 7),
              Hypothesis("#Lambda_{1}", "MC_weight_ggH_g1prime2", ROOT.kGreen+3),
              Hypothesis("f_{#Lambda1}=0.5", "MC_weight_ggH_g1g1prime2", 3),
             )

angles = (
          Angle("Z1Mass", "m_{1}", 50, 40, 100),
          Angle("Z2Mass", "m_{2}", 50, 12, 60),
          Angle("costhetastar", "cos#theta*", 50, -1, 1),
          Angle("helcosthetaZ1", "cos#theta_{1}", 50, -1, 1),
          Angle("helcosthetaZ2", "cos#theta_{2}", 50, -1, 1),
          Angle("helphi", "#Phi", 50, -pi, pi),
          Angle("phistarZ1", "#Phi_{1}", 50, -pi, pi),
         )

h = {}
hstack = {}

addtolegend = True
legend = ROOT.TLegend(.5, .7, .9, .9)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetNColumns(2)

c1 = ROOT.TCanvas()

for angle in angles:
  anglename, title, bins, min, max = angle
  hstack[angle] = ROOT.THStack(anglename, "")

  for hypothesis in hypotheses:
    if hypothesis is None:
      if addtolegend:
        legend.AddEntry(0, "", "")
      continue

    name, weightname, color = hypothesis

    wt = "({})".format(weightname)
    draw = "{}>>{}({}, {}, {})".format(anglename, anglename+name, bins, min, max)
    t.Draw(draw, wt)
    hist = h[hypothesis,angle] = getattr(ROOT, anglename+name)

    hist.Scale(1/hist.Integral())
    hist.SetLineColor(color)

    hstack[angle].Add(hist)

    if addtolegend:
      legend.AddEntry(hist, name, "l")

  addtolegend = False

  hstack[angle].Draw("hist nostack")
  hstack[angle].GetXaxis().SetTitle(title)
  legend.Draw()

  for ext in "png eps root pdf".split():
    c1.SaveAs("~/www/TEST/oldtrees/{}.{}".format(anglename, ext))
