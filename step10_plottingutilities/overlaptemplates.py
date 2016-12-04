#!/usr/bin/env python

assert __name__ == "__main__"

from helperstuff import config, style
import helperstuff.rootoverloads.histogramaxisnumbers
from helperstuff.templates import Template
import os
import ROOT

########################
#inputs
templates = [
             Template("fa3", "ggH", "0+", "2e2mu", "Untagged"),
             Template("fa3", "VBF", "0+", "2e2mu", "Untagged"),
             Template("fa3", "ZH", "0+", "2e2mu", "Untagged"),
            ]
axis = 2
########################

hstack = ROOT.THStack("hstack", "hstack")

for color, template in enumerate(templates, start=1):
    proj = template.gettemplate().Projection(axis)
    proj.Scale(1/proj.Integral())
    proj.SetLineColor(color)
    hstack.Add(proj)

c = ROOT.TCanvas()
hstack.Draw("hist nostack")
for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "test.{}".format(ext)))
