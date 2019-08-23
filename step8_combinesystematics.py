#!/usr/bin/env python

import argparse
if __name__ == "__main__":
  def __Analysis(*args, **kwargs):
    from helperstuff.enums import Analysis
    return Analysis(*args, **kwargs)
  parser = argparse.ArgumentParser()
  parser.add_argument("analysis", type=__Analysis, nargs="?")
  args = parser.parse_args()

import itertools, os

import numpy as np

from helperstuff import config

from helperstuff.combinehelpers import getrate
from helperstuff.enums import analyses, categories, channels, ProductionMode, productions, ShapeSystematic
from helperstuff.templates import IntTemplate, Template, TemplatesFile
from helperstuff.utilities import KeepWhileOpenFile, TFile

threshold = 100.

STXSuncertainties = "Mu", "Res", "Mig01", "Mig12", "VBF2j", "VBF3j", "PT60", "PT120", "qmtop"

def combinesystematics(channel, analysis, production, category, productionmode):
  templategroup = str(productionmode).lower()
  tfnominal = TemplatesFile(channel, analysis, production, category, templategroup)

  with TFile(tfnominal.templatesfile()) as fnominal:
    for systname in "ScaleUp", "ScaleDn", "ResUp", "ResDn", "JECUp", "JECDn", "THU_ggH_MuUp", "THU_ggH_ResUp", "THU_ggH_Mig01Up", "THU_ggH_Mig12Up", "THU_ggH_VBF2jUp", "THU_ggH_VBF3jUp", "THU_ggH_PT60Up", "THU_ggH_PT120Up", "THU_ggH_qmtopUp", "THU_ggH_MuDn", "THU_ggH_ResDn", "THU_ggH_Mig01Dn", "THU_ggH_Mig12Dn", "THU_ggH_VBF2jDn", "THU_ggH_VBF3jDn", "THU_ggH_PT60Dn", "THU_ggH_PT120Dn", "THU_ggH_qmtopDn":
      syst = ShapeSystematic(systname)
      if syst in ("ScaleUp", "ScaleDn", "ResUp", "ResDn") and not config.applym4lshapesystematics: continue
      if syst in ("JECUp", "JECDn") and not config.applyJECshapesystematics: continue
      if syst.isTHUggH and not config.applySTXSsystematics: continue

      if syst in ("JECUp", "JECDn") and category not in ("VBFtagged", "VHHadrtagged"): continue
      if syst.isTHUggH and not syst.applySTXStocategory(category): continue

      if syst.isTHUggH and productionmode != "ggH": continue

      tfsyst = TemplatesFile(channel, analysis, production, category, templategroup, syst)
      with TFile(tfsyst.templatesfile()) as fsyst:
        for hypothesis in syst.hypothesesforratio:
          numerator = getattr(fsyst, Template(tfsyst, productionmode, hypothesis).templatename())
          denominator = getattr(fnominal, Template(tfnominal, productionmode, hypothesis).templatename())

          ratio = numerator.Clone("ratio")
          ratio.Divide(denominator)

          for x, y, z in itertools.product(xrange(1, ratio.GetNbinsX()+1), xrange(1, ratio.GetNbinsY()+1), xrange(1, ratio.GetNbinsZ()+1)):
            if np.isclose(denominator.GetBinContent(x, y, z), 1e-10) or np.isclose(numerator.GetBinContent(x, y, z), 1e-10):
              ratio.SetBinContent(x, y, z, 1)
            if ratio.GetBinContent(x, y, z) > threshold or ratio.GetBinContent(x, y, z) < 1/threshold:
              if getrate(channel, analysis, production, category, productionmode, "fordata") * denominator.GetBinContent(x, y, z) / denominator.Integral() < 1e-3:
                ratio.SetBinContent(x, y, z, 1)
              else:
                raise ValueError("Huge or tiny ratio for ({}) / ({}) bin {} {} {}: {} / {} = {}".format(Template(tfsyst, productionmode, hypothesis), Template(tfnominal, productionmode, hypothesis), x, y, z, numerator.GetBinContent(x, y, z), denominator.GetBinContent(x, y, z), ratio.GetBinContent(x, y, z)))


          newsyst = ShapeSystematic(str(syst).replace("Up", hypothesis.combinename+"Up").replace("Dn", hypothesis.combinename+"Dn").replace("Down", hypothesis.combinename+"Down"))
          assert newsyst != syst, (syst, newsyst)
          newtfsyst = TemplatesFile(channel, analysis, production, category, templategroup, newsyst)

          newfilename = newtfsyst.templatesfile()

          with KeepWhileOpenFile(newfilename+".tmp") as kwof:
            if not kwof: continue
            if os.path.exists(newfilename): continue

            print production, analysis, channel, category, productionmode, newsyst
            cache = []
            with TFile(newfilename, "CREATE", deleteifbad=True) as newfsyst:
              for template in itertools.chain(tfnominal.templates(), tfnominal.inttemplates()):
                if isinstance(template, Template) and not template.hypothesis.ispure: continue
                newtemplate = getattr(fnominal, template.templatename()).Clone()
                newtemplate.Multiply(ratio)
                newtemplate.SetDirectory(newfsyst)
                cache.append(newtemplate)
            del cache

if __name__ == "__main__":
  for production in productions:
    for analysis in analyses:
      if analysis != args.analysis is not None: continue
      for channel in channels:
        for category in categories:
          for productionmode in ("ggH", "VBF", "VH", "ttH", "bbH"):
            productionmode = ProductionMode(productionmode)
            if productionmode != "ggH" and analysis.isdecayonly: continue
            if category != "Untagged" and analysis.isdecayonly: continue
            if category == "Boosted" and not analysis.useboosted: continue
            if category in ("VBF1jtagged", "VHLepttagged") and not analysis.usemorecategories: continue
            combinesystematics(channel, analysis, production, category, productionmode)
