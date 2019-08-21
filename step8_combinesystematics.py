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

from helperstuff.enums import analyses, categories, channels, ProductionMode, productions, ShapeSystematic
from helperstuff.templates import IntTemplate, Template, TemplatesFile
from helperstuff.utilities import KeepWhileOpenFile, TFile

def combinesystematics(channel, analysis, production, category, productionmode):
  templategroup = str(productionmode).lower()
  tfnominal = TemplatesFile(channel, analysis, production, category, templategroup)

  with TFile(tfnominal.templatesfile()) as fnominal:
    for syst in "ScaleUp", "ScaleDn", "ResUp", "ResDn", "JECUp", "JECDn":
      if syst in ("ScaleUp", "ScaleDn", "ResUp", "ResDn") and not config.applym4lshapesystematics: continue
      if syst in ("JECUp", "JECDn") and not config.applyJECshapesystematics: continue

      if syst in ("JECUp", "JECDn") and category not in ("VBFtagged", "VHHadrtagged"): continue

      tfsyst = TemplatesFile(channel, analysis, production, category, templategroup, syst)
      with TFile(tfsyst.templatesfile()) as fsyst:
        for hypothesis in ShapeSystematic(syst).hypothesesforratio:
          numerator = getattr(fsyst, Template(tfsyst, productionmode, hypothesis).templatename())
          denominator = getattr(fnominal, Template(tfnominal, productionmode, hypothesis).templatename())

          ratio = numerator.Clone("ratio")
          ratio.Divide(denominator)

          for x, y, z in itertools.product(xrange(1, ratio.GetNbinsX()+1), xrange(1, ratio.GetNbinsY()+1), xrange(1, ratio.GetNbinsZ()+1)):
            if np.isclose(denominator.GetBinContent(x, y, z), 1e-10):
              ratio.SetBinContent(x, y, z, 1)
            if ratio.GetBinContent(x, y, z) > 1000:
              raise ValueError("Huge ratio for ({}) / ({}) bin {} {} {}: {} / {} = {}".format(Template(tfsyst, productionmode, hypothesis), Template(tfnominal, productionmode, hypothesis), x, y, z, numerator.GetBinContent(x, y, z), denominator.GetBinContent(x, y, z), ratio.GetBinContent(x, y, z)))


          newsyst = ShapeSystematic(syst.replace("Up", hypothesis.combinename+"Up").replace("Dn", hypothesis.combinename+"Dn"))
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
