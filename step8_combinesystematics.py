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

from helperstuff import config

from helperstuff.enums import analyses, categories, channels, ProductionMode, productions, ShapeSystematic
from helperstuff.templates import IntTemplate, Template, TemplatesFile
from helperstuff.utilities import KeepWhileOpenFile, TFile

def combinesystematics(channel, analysis, production, category, productionmode):
  templategroup = str(productionmode).lower()
  tfnominal = TemplatesFile(channel, analysis, production, category, templategroup)

  with TFile(tfnominal.templatesfile()) as fnominal:
    for syst in "ScaleUp", "ScaleDn", "ResUp", "ResDn":
      tfsyst = TemplatesFile(channel, analysis, production, category, templategroup, syst)
      with TFile(tfsyst.templatesfile()) as fsyst:
        for hypothesis in ShapeSystematic(syst).hypothesesforratio:
          numerator = getattr(fsyst, Template(tfsyst, productionmode, hypothesis).templatename())
          denominator = getattr(fnominal, Template(tfnominal, productionmode, hypothesis).templatename())
          ratio = numerator.Clone("ratio")
          ratio.Divide(denominator)

          newsyst = ShapeSystematic(syst.replace("Up", hypothesis.combinename+"Up").replace("Dn", hypothesis.combinename+"Dn"))
          newtfsyst = TemplatesFile(channel, analysis, production, category, templategroup, newsyst)

          newfilename = newtfsyst.templatesfile()

          with KeepWhileOpenFile(newfilename+".tmp") as kwof:
            if not kwof: continue
            if os.path.exists(newfilename): continue

            print production, analysis, channel, category, productionmode, newsyst
            with TFile(newfilename, "CREATE", deleteifbad=True) as newfsyst:
              for template in itertools.chain(tfnominal.templates(), tfnominal.inttemplates()):
                if isinstance(template, Template) and not template.hypothesis.ispure: continue
                newtemplate = getattr(fnominal, template.templatename()).Clone()
                newtemplate.Multiply(ratio)
                newtemplate.SetDirectory(newfsyst)

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
