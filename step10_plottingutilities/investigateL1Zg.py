#!/usr/bin/env python

import os

import ROOT

from helperstuff import config
from helperstuff.samples import ReweightingSample, samplewithfai
from helperstuff.templates import TemplatesFile

from projections import ComponentTemplateSum, IntTemplateFromFile, TemplateFromFile

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

def investigateL1Zg():
  thesample = samplewithfai("ggH", "fL1Zg", .99)
  g1_mix, gi_mix = thesample.g1, thesample.ghzgs1prime2
  fainame = "f_{#Lambda1}^{Z#gamma, dec}"

  gi_ggHBSM = ReweightingSample("ggH", "L1Zg").ghzgs1prime2

  discriminants = TemplatesFile("2e2mu", "ggh", "fL1Zg", production, "Untagged").discriminants

  ggHSM     = TemplateFromFile(   0, "ggH", "Untagged", "fullrange", "rescalemixtures", production, "2e2mu", "fL1Zg", "0+_photoncut")
  ggHBSM    = TemplateFromFile(   0, "ggH", "Untagged", "fullrange", "rescalemixtures", production, "2e2mu", "fL1Zg", "L1Zg")
  ggHint    = IntTemplateFromFile(0, "ggH", "Untagged", "fullrange", "rescalemixtures", production, "2e2mu", "fL1Zg", "g11gi1")

  SM_name   = ComponentTemplateSum("ggH {}=0" .format(fainame), 1, ggHSM.Integral(), (ggHSM, 1))
  BSM_name  = ComponentTemplateSum("ggH {}=1" .format(fainame), 2, ggHSM.Integral(), (ggHBSM, 1))
  plus0p99 = ComponentTemplateSum("ggH {}=#plus0.99" .format(fainame), ROOT.kGreen+3, ggHSM.Integral(), (ggHSM, g1_mix**2), (ggHBSM, (gi_mix/gi_ggHBSM)**2), (ggHint,  g1_mix*gi_mix/gi_ggHBSM))
  minus0p99 = ComponentTemplateSum("ggH {}=#minus0.99" .format(fainame), 4, ggHSM.Integral(), (ggHSM, g1_mix**2), (ggHBSM, (-gi_mix/gi_ggHBSM)**2), (ggHint,  g1_mix*-gi_mix/gi_ggHBSM))

  templates = [ggHSM, ggHBSM, ggHint, SM_name, BSM_name, plus0p99, minus0p99]

  c1 = ROOT.TCanvas()
  legend = ROOT.TLegend(.65, .6, .9, .9)
  legend.SetBorderSize(0)
  legend.SetFillStyle(0)
  for template in templates:
    if template.color:
      if template.title:
        template.AddToLegend(legend)

  for i, discriminant in enumerate(discriminants):
    hstack = ROOT.THStack("{}".format(discriminant.name), discriminant.title)
    for template in templates:
      if template.color:
        hstack.Add(template.Projection(i), template.hstackoption)
    hstack.Draw("nostack")
    hstack.GetXaxis().SetTitle(discriminant.title)
    legend.Draw()

    saveasdir = os.path.join(config.plotsbasedir, "checkfL1Zg")

    try:
      os.makedirs(saveasdir)
    except OSError:
      pass
    for ext in "png eps root pdf".split():
      c1.SaveAs(os.path.join(saveasdir, "{}.{}".format(discriminant.name, ext)))

if __name__ == "__main__":
  investigateL1Zg()
