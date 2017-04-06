#!/usr/bin/env python

import array

import ROOT

import helperstuff.style

from helperstuff import config
from helperstuff.combinehelpers import getdatatree
from helperstuff.enums import channels
from projections import TemplateFromFile, TemplateSum

templates = sum(([
  TemplateFromFile("ggH", "Untagged", ch, config.productionforcombine, "", "fa3", "SM", "enrich"),
  TemplateFromFile("VBF", "Untagged", ch, config.productionforcombine, "", "fa3", "SM", "enrich"),
  TemplateFromFile("ZH", "Untagged", ch, config.productionforcombine, "", "fa3", "SM", "enrich"),
  TemplateFromFile("WH", "Untagged", ch, config.productionforcombine, "", "fa3", "SM", "enrich"),
  TemplateFromFile("ttH", "Untagged", ch, config.productionforcombine, "", "fa3", "SM", "enrich"),
  TemplateFromFile("ggZZ", "Untagged", ch, config.productionforcombine, "", "fa3", "enrich"),
  TemplateFromFile("qqZZ", "Untagged", ch, config.productionforcombine, "", "fa3", "enrich"),
  TemplateFromFile("ZX", "Untagged", ch, config.productionforcombine, "", "fa3", "enrich"),
  TemplateFromFile("VBF bkg", "Untagged", ch, config.productionforcombine, "", "fa3", "enrich"),
] for ch in channels), [])

t = TemplateSum("", *((t, 1) for t in templates))

c = ROOT.TCanvas()
t.Project3D("yx").Draw("colz")

discriminants = templates[0].discriminants

x, y = [], []
graphs = []

print t.h.Integral()

for style, ch in enumerate(("4e", "4mu", "2e2mu"), start=20):
  tree = getdatatree(ch, "Untagged", "fa3", config.productionforcombine)
  for entry in tree:
    if tree.D_bkg > 0.5:
      x.append(getattr(tree, discriminants[0].name))
      y.append(getattr(tree, discriminants[1].name))
  print len(x), len(y)
  g = ROOT.TGraph(len(x), array.array('d', x), array.array('d', y))
  g.SetMarkerStyle(style)
  graphs.append(g)
  g.Draw("P")

c.SaveAs("~/www/TEST/test.png")
c.SaveAs("~/www/TEST/test.pdf")
