#!/usr/bin/env python

from collections import namedtuple
import itertools
import logging
from math import sqrt
import os
import yaml

from helperstuff import config
from helperstuff.categorization import MultiCategorization
from helperstuff.enums import AlternateWeight, alternateweights, analyses, categories, Category, channels, flavors, ProductionMode, pythiasystematics
from helperstuff.samples import ReweightingSample, ReweightingSamplePlus, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import MultiplyCounter
from helperstuff.yields import count, totalrate, YieldSystematicValue, YieldValue

from categorysystematics import findsystematic

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

SampleCount = namedtuple("SampleCount", "productionmode samples")

def writeyields():
  tosamples_foryields = [
    SampleCount(ProductionMode("ggH"), {ReweightingSamplePlus("ggH", "0+", "MINLO")}),
    SampleCount(ProductionMode("VBF"), {ReweightingSamplePlus("VBF", "0+", "POWHEG")}),
    SampleCount(ProductionMode("ZH"), {ReweightingSamplePlus("ZH", "0+", "POWHEG")}),
    SampleCount(ProductionMode("WH"), {ReweightingSamplePlus("WplusH", "0+", "POWHEG"), ReweightingSamplePlus("WminusH", "0+", "POWHEG")}),
    SampleCount(ProductionMode("ttH"), {ReweightingSamplePlus("ttH", "0+", "Hff0+", "POWHEG")}),
    SampleCount(ProductionMode("qqZZ"), {ReweightingSample("qqZZ")}),
    SampleCount(ProductionMode("ggZZ"), {ReweightingSample("ggZZ", flavor) for flavor in flavors}),
    SampleCount(ProductionMode("VBF bkg"), {ReweightingSample("VBF bkg", flavor) for flavor in ("2e2mu", "4e", "4mu")}),
    SampleCount(ProductionMode("ZX"), {ReweightingSample("ZX")}),
  ]

  categorizations = [_ for _ in TreeWrapper.categorizations if isinstance(_, MultiCategorization)]

  for productionmode, samples in tosamples_foryields:
    print productionmode
    result = MultiplyCounter()
    samplegroups = [samples]
    if productionmode.issignal:
      samplegroups += [{ReweightingSamplePlus(s, systematic) for s in samples} for systematic in pythiasystematics]

    for usesamples in samplegroups:
      tmpresult = MultiplyCounter()
      for tosample in usesamples:
        usealternateweights = alternateweights
        if (productionmode in ("ggZZ", "VBF bkg", "ZX")
             or hasattr(tosample, "pythiasystematic") and tosample.pythiasystematic is not None):
          usealternateweights = [AlternateWeight("1")]

        tmpresult += count({Sample(tosample, production)}, {tosample}, categorizations, usealternateweights)

      for key in tmpresult:
        if len(key) == 3: continue
        elif 4 <= len(key) <= 5:
#          try:
            assert key[0] in usesamples
            tmpresult[key] /= sum(tmpresult[tosample,key[1],key[2]] for tosample in usesamples)
#          except ZeroDivisionError:
#            pass
        else: assert False

      result += tmpresult


    #if productionmode.issignal():
    #  result += count(
    #                  {ReweightingSample(productionmode, "0+")},
    #                  {Sample(productionmode, h, production) for h in productionmode.generatedhypotheses},
    #                  TreeWrapper.categorizations,
    #                 )

    for analysis in analyses:
      categorization = {_ for _ in categorizations if _.category_function_name == "category_"+analysis.categoryname}
      assert len(categorization) == 1, categorization
      categorization = categorization.pop()

      total = totalrate(productionmode, production, 1.0)
      for channel, category in itertools.product(channels, categories):
        YieldValue(channel, category, analysis, productionmode).value = total*sum(result[tosample, categorization, AlternateWeight("1"), category, channel] for tosample in samples)

      #same for all categories and channels
      #from yaml
      for systname in "pdf_Higgs_gg_yield", "pdf_Higgs_qq_yield", "pdf_Higgs_ttH_yield", "BRhiggs_hzz4l", "QCDscale_ggZH_yield", "QCDscale_ggVV_bonly_yield", "lumi_13TeV":
        for category, channel in itertools.product(categories, channels):
          syst = YieldSystematicValue(channel, category, analysis, productionmode, systname)
          values = syst.yieldsystematic.getfromyaml()
          if productionmode.yamlsystname in values:
            syst.value = values[productionmode.yamlsystname]
          else:
            syst.value = None

      #same for all categories
      #from yaml
      for systname in "CMS_eff_e", "CMS_eff_m":
        for category, channel in itertools.product(categories, channels):
          syst = YieldSystematicValue(channel, category, analysis, productionmode, systname)
          if channel == "4e" and systname == "CMS_eff_m" or channel == "4mu" and systname == "CMS_eff_e":
            values = {}
          else:
            values = syst.yieldsystematic.getfromyaml(channel=channel)
          if productionmode.yamlsystname in values:
            syst.value = values[productionmode.yamlsystname]
          else:
            syst.value = None

      #Z+X, for now
      for systname in "CMS_zz2e2mu_zjets", "CMS_zz4e_zjets", "CMS_zz4mu_zjets":
        for category, channel in itertools.product(categories, channels):
          syst = YieldSystematicValue(channel, category, analysis, productionmode, systname)
          if str(channel) in systname and productionmode == "ZX":
            syst.value = 0.4
          else:
            syst.value = None

      #same for all channels
      for category in categories:
        #variations on category definition: JEC and btagging
        nominal = sum(result[tosample, categorization, AlternateWeight("1"), category] for tosample in samples)
        nominaluntagged = sum(result[tosample, categorization, AlternateWeight("1"), Category("Untagged")] for tosample in samples)
        nominalyield = sum(result[tosample, categorization, AlternateWeight("1")] for tosample in samples)
        if productionmode == "ZX":
          JECUp = JECDn = btSFUp = btSFDn = 1
        else:
          JECUp = sum(result[tosample, findsystematic(categorizations, categorization, "JECUp", "Nominal"), AlternateWeight("1"), category] for tosample in samples) / nominal
          JECDn = sum(result[tosample, findsystematic(categorizations, categorization, "JECDn", "Nominal"), AlternateWeight("1"), category] for tosample in samples) / nominal
          btSFUp = sum(result[tosample, findsystematic(categorizations, categorization, "Nominal", "bTagSFUp"), AlternateWeight("1"), category] for tosample in samples) / nominal
          btSFDn = sum(result[tosample, findsystematic(categorizations, categorization, "Nominal", "bTagSFDn"), AlternateWeight("1"), category] for tosample in samples) / nominal

        for channel in channels:
          YieldSystematicValue(channel, category, analysis, productionmode, "JEC").value = (JECUp, JECDn)
          YieldSystematicValue(channel, category, analysis, productionmode, "BTag").value = (btSFUp, btSFDn)

        allQCDsystematicnames = [ProductionMode(_).QCDsystematicname for _ in ("ggH", "VBF", "ZH", "WH", "ttH", "qqZZ")]

        #QCD weight variations
        for systname in allQCDsystematicnames:
          if productionmode == "ggZZ":
            for channel in channels:
              YieldSystematicValue(channel, category, analysis, productionmode, systname+"_yield").value = YieldSystematicValue(channel, category, analysis, "ggH", systname+"_yield").value
              YieldSystematicValue(channel, category, analysis, productionmode, systname+"_cat").value = YieldSystematicValue(channel, category, analysis, "ggH", systname+"_cat").value
          elif productionmode == "VBF bkg":
            for channel in channels:
              YieldSystematicValue(channel, category, analysis, productionmode, systname+"_yield").value = YieldSystematicValue(channel, category, analysis, "VBF", systname+"_yield").value
              YieldSystematicValue(channel, category, analysis, productionmode, systname+"_cat").value = YieldSystematicValue(channel, category, analysis, "VBF", systname+"_cat").value
          elif systname == productionmode.QCDsystematicname:
            muYieldUp, muYieldDn = {}, {}
            muUp, muDn = {}, {}
            muUpUntagged, muDnUntagged = {}, {}

            muYieldUp["R"] = sum(result[tosample, categorization, AlternateWeight("muRUp")] for tosample in samples) / nominalyield
            muYieldUp["F"] = sum(result[tosample, categorization, AlternateWeight("muFUp")] for tosample in samples) / nominalyield
            muYieldDn["R"] = sum(result[tosample, categorization, AlternateWeight("muRDn")] for tosample in samples) / nominalyield
            muYieldDn["F"] = sum(result[tosample, categorization, AlternateWeight("muFDn")] for tosample in samples) / nominalyield
            muUp["R"] = sum(result[tosample, categorization, AlternateWeight("muRUp"), category] for tosample in samples) / nominal
            muUp["F"] = sum(result[tosample, categorization, AlternateWeight("muFUp"), category] for tosample in samples) / nominal
            muDn["R"] = sum(result[tosample, categorization, AlternateWeight("muRDn"), category] for tosample in samples) / nominal
            muDn["F"] = sum(result[tosample, categorization, AlternateWeight("muFDn"), category] for tosample in samples) / nominal

            #figure out which to combine with which
            muUpUntagged["R"] = sum(result[tosample, categorization, AlternateWeight("muRUp"), Category("Untagged")] for tosample in samples) / nominaluntagged
            muUpUntagged["F"] = sum(result[tosample, categorization, AlternateWeight("muFUp"), Category("Untagged")] for tosample in samples) / nominaluntagged
            muDnUntagged["R"] = sum(result[tosample, categorization, AlternateWeight("muRDn"), Category("Untagged")] for tosample in samples) / nominaluntagged
            muDnUntagged["F"] = sum(result[tosample, categorization, AlternateWeight("muFDn"), Category("Untagged")] for tosample in samples) / nominaluntagged

            if muYieldUp["R"] > 1 > muYieldDn["R"]:
              pass
            elif muYieldDn["R"] > 1 > muYieldUp["R"]:
              muYieldUp["R"], muYieldDn["R"] = muYieldDn["R"], muYieldUp["R"]
            else:
              assert False, (muYieldUp["R"], muYieldDn["R"])

            if muYieldUp["F"] > 1 > muYieldDn["F"]:
              pass
            elif muYieldDn["F"] > 1 > muYieldUp["F"]:
              muYieldUp["F"], muYieldDn["F"] = muYieldDn["F"], muYieldUp["F"]
            else:
              assert False, (muYieldUp["F"], muYieldDn["F"])

            if muUpUntagged["R"] > muDnUntagged["R"]:
              pass
            elif muDnUntagged["R"] > muUpUntagged["R"]:
              muUp["R"], muDn["R"] = muDn["R"], muUp["R"]
              muUpUntagged["R"], muDnUntagged["R"] = muDnUntagged["R"], muUpUntagged["R"]
            else:
              assert False, (muUpUntagged["R"], muDnUntagged["R"])

            if muUpUntagged["F"] > muDnUntagged["F"]:
              pass
            elif muDnUntagged["F"] > muUpUntagged["F"]:
              muUp["F"], muDn["F"] = muDn["F"], muUp["F"]
              muUpUntagged["F"], muDnUntagged["F"] = muDnUntagged["F"], muUpUntagged["F"]
            else:
              assert False, (muUpUntagged["F"], muDnUntagged["F"])

            if muUpUntagged["F"] > 1 and muDnUntagged["F"] > 1 or muUpUntagged["F"] < 1 and muDnUntagged["F"] < 1:
              logging.warning("muF up and down go in the same direction {} {}".format(muUpUntagged, muDnUntagged))
            if muUpUntagged["R"] > 1 and muDnUntagged["R"] > 1 or muUpUntagged["R"] < 1 and muDnUntagged["R"] < 1:
              logging.warning("muR up and down go in the same direction {} {}".format(muUpUntagged, muDnUntagged))

            QCDYieldUp = 1 + sqrt((muYieldUp["R"]-1)**2 + (muYieldUp["F"]-1)**2)
            QCDYieldDn = 1 - sqrt((muYieldDn["R"]-1)**2 + (muYieldDn["F"]-1)**2)

            if   muUp["R"] >= 1 and muUp["F"] >= 1 and muDn["R"] <= 1 and muDn["F"] <= 1: signup, signdn = +1, -1
            elif muUp["R"] <= 1 and muUp["F"] <= 1 and muDn["R"] >= 1 and muDn["F"] >= 1: signup, signdn = -1, +1
            else:
              logging.warning("muR and muF go in different directions in {} {}: {} {}".format(productionmode, category, muUp, muDn))
              if   sum(muUp.values()) > sum(muDn.values()): signup, signdn = +1, -1
              elif sum(muUp.values()) < sum(muDn.values()): signup, signdn = -1, +1
              else: assert False

            assert signup != signdn

            QCDUp = 1 + signup*sqrt((muUp["R"]-1)**2 + (muUp["F"]-1)**2)
            QCDDn = 1 + signdn*sqrt((muDn["R"]-1)**2 + (muDn["F"]-1)**2)

            for channel in channels:
              YieldSystematicValue(channel, category, analysis, productionmode, systname+"_cat").value = (QCDUp, QCDDn)
              YieldSystematicValue(channel, category, analysis, productionmode, systname+"_yield").value = (QCDYieldUp, QCDYieldDn)
          else:
            for channel in channels:
              YieldSystematicValue(channel, category, analysis, productionmode, systname+"_cat").value = None
              YieldSystematicValue(channel, category, analysis, productionmode, systname+"_yield").value = None

        #pythia scale and tune
        if productionmode in ("ggH", "VBF", "ZH", "WH", "ttH"):
          scaleup = sum(result[ReweightingSamplePlus(tosample, "ScaleUp"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
          scaledn = sum(result[ReweightingSamplePlus(tosample, "ScaleDn"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
          tuneup = sum(result[ReweightingSamplePlus(tosample, "TuneUp"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
          tunedn = sum(result[ReweightingSamplePlus(tosample, "TuneDn"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
        else:
          scaleup = scaledn = tuneup = tunedn = 1
        for channel in channels:
          YieldSystematicValue(channel, category, analysis, productionmode, "PythiaScale").value = (scaleup, scaledn)
          YieldSystematicValue(channel, category, analysis, productionmode, "PythiaTune").value = (tuneup, tunedn)

    YieldValue.writedict()
    YieldSystematicValue.writedict()

if __name__ == "__main__":
  writeyields()
