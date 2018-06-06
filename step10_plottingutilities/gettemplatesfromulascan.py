#!/usr/bin/env python

import os
from helperstuff import constants
from helperstuff.templates import IntTemplate, Template, templatesfiles
from helperstuff.utilities import KeepWhileOpenFile, TFile


def gettemplatefromulascan(template):
  folder = "/work-zfs/lhc/usarica/hep/SpinWidthPaper_2015/{year}Width/CMSSW_9_4_3/src/HiggsWidth_PostICHEP/Analysis/test/output/LHC_13TeV/Templates/{templatedate}/FinalTemplates/Stage1/{coupling}"

  kwargs = {}
  kwargses = [kwargs]

  if template.production == "180416_Ulascan":
    if template.analysis in ("fa2", "fL1"):
      kwargs["templatedate"] = 180520
    else:
      kwargs["templatedate"] = 180509
  if template.production == "180224_Ulascan":
    if template.analysis in ("fa2", "fL1"):
      kwargs["templatedate"] = 180520
    else:
      kwargs["templatedate"] = 180509

  if template.analysis == "fa3":
    kwargs["coupling"] = "a3"
    gi = constants.g4HZZ
  elif template.analysis == "fa2":
    kwargs["coupling"] = "a2"
    gi = constants.g2HZZ
  elif template.analysis == "fL1":
    kwargs["coupling"] = "L1"
    gi = constants.g1prime2HZZ

  basename = "HtoZZ{channel}_{category}_FinalTemplates_{productionmode}_{systematic}.root"
  templatename = "T_{productionmode}{_hypothesis}"

  kwargs["channel"] = template.channel
  kwargs["category"] = str(template.category).replace("VBFt", "JJVBFT").replace("VHHadrt", "HadVHT")
  if template.shapesystematic == "":
    kwargs["systematic"] = "Nominal"
  kwargs["year"] = template.production.year

  if template.productionmode == "ggH":
    kwargs["productionmode"] = "ggZZ"
    if isinstance(template, Template) and template.hypothesis == "0+":
      kwargs["_hypothesis"] = "_Sig"
      scalingpower = 0
    if isinstance(template, IntTemplate) and template.interferencetype == "g11gi1":
      kwargs["_hypothesis"] = "_Sig_ai1_1_Re"
      scalingpower = 1
    if isinstance(template, Template) and template.hypothesis in ("0-", "a2", "a3", "L1", "L1Zg"):
      kwargs["_hypothesis"] = "_Sig_ai1_2"
      scalingpower = 2

  if template.productionmode in ("VBF", "ZH", "WH"):
    kwargs["productionmode"] = template.productionmode
    if isinstance(template, Template) and template.hypothesis == "0+":
      kwargs["_hypothesis"] = "_Sig"
      scalingpower = 0
    if isinstance(template, IntTemplate) and template.interferencetype == "g13gi1":
      kwargs["_hypothesis"] = "_Sig_ai1_1_Re"
      scalingpower = 1
    if isinstance(template, IntTemplate) and template.interferencetype == "g12gi2":
      kwargs["_hypothesis"] = "_Sig_ai1_2_PosDef"
      scalingpower = 2
    if isinstance(template, IntTemplate) and template.interferencetype == "g11gi3":
      kwargs["_hypothesis"] = "_Sig_ai1_3_Re"
      scalingpower = 3
    if isinstance(template, Template) and template.hypothesis in ("0-", "a2", "a3", "L1", "L1Zg"):
      kwargs["_hypothesis"] = "_Sig_ai1_4"
      scalingpower = 4

  if template.productionmode == "ggZZ":
    kwargs["productionmode"] = "ggZZ"
    kwargs["_hypothesis"] = "_Bkg"
    scalingpower = 0

  if template.productionmode == "qqZZ":
    kwargs["productionmode"] = "bkg_qqzz"
    kwargs["_hypothesis"] = ""
    scalingpower = 0

  if template.productionmode == "ZX":
    kwargs["productionmode"] = "zjets"
    kwargs["_hypothesis"] = ""
    scalingpower = 0

  if template.productionmode == "VBF bkg":
    kwargs["productionmode"] = "VBF"
    kwargs["_hypothesis"] = "_Bkg"
    scalingpower = 0
    kwargses = [kwargs, kwargs.copy(), kwargs.copy()]
    kwargses[1]["productionmode"] = "ZH"
    kwargses[2]["productionmode"] = "WH"

  h = None

  for kwargs in kwargses:
    filename = os.path.join(folder, basename).format(**kwargs)
    with TFile(filename) as f:
      _ = getattr(f, templatename.format(**kwargs))
      _.Sumw2()
      if h is None:
        h = _
        h.SetDirectory(0)
      else:
        h.Add(_)

  print h.Integral(),
  if h.Integral(): h.Scale(h.Integral("width") / h.Integral())
  h.Scale(1 / gi ** scalingpower)
  if template.productionmode == "ZX":
    txtfilename = "/work-zfs/lhc/usarica/hep/SpinWidthPaper_2015/RunII_Combination/CMSSW_9_4_3/src/HZZ4l/CreateWidthDatacards/cards_180516_Onshell_{coupling}_13TeV_{year}/HCG/13TeV_{year}/hzz{channel}_{category}.txt".format(**kwargs)
    with open(txtfilename) as f:
      for line in f:
        if line.startswith("rate "):
          multipliers = line.split()[1:]
          assert all(float(_)==1 for _ in multipliers[:-1]), line
          print multipliers[-1],
          h.Scale(float(multipliers[-1]))
          break
  print h.Integral()
  h.SetName(template.templatename())

  return h

def gettemplatesfromulascan(tf):
  filename = tf.templatesfile()
  with KeepWhileOpenFile(filename+".tmp") as kwof:
    print tf
    if not kwof: return
    if os.path.exists(filename): return
    try:
      cache = []
      with TFile(filename, "CREATE", write=True) as f:
        for template in tf.templates() + tf.inttemplates():
          if isinstance(template, Template) and template.hypothesis and not template.hypothesis.ispure: continue
          h = gettemplatefromulascan(template)
          cache.append(h)
          h.SetDirectory(f)
    except:
      if os.path.exists(filename):
        os.remove(filename)
      raise

def getalltemplatesfromulascan():
  for tf in templatesfiles:
    if "Ulascan" not in str(tf.production): continue
    if tf.analysis not in ("fa3", "fa2", "fL1"): continue
    if tf.templategroup == "tth": continue
    gettemplatesfromulascan(tf)

if __name__ == "__main__":
  getalltemplatesfromulascan()
