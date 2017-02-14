#!/usr/bin/env python

from collections import namedtuple
import itertools
import os
import yaml

from helperstuff import config
from helperstuff.enums import analyses, categories, channels, flavors, ProductionMode, pythiasystematics
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

  categorizations = TreeWrapper.categorizations

  for productionmode, samples in tosamples_foryields:
    print productionmode
    result = MultiplyCounter()
    samplegroups = [samples]
    if productionmode.issignal:
      samplegroups += [{ReweightingSamplePlus(s, systematic) for s in samples} for systematic in pythiasystematics]

    for usesamples in samplegroups:
      tmpresult = MultiplyCounter()
      for tosample in usesamples:
        ####################################
        if isinstance(tosample, ReweightingSamplePlus) and tosample.alternategenerator == "MINLO" and tosample.pythiasystematic is None:
          try:
            Sample(tosample, production)
          except ValueError:
            pass
          else:
            assert False, "Delete this section!"
          continue
        ####################################
        tmpresult += count({Sample(tosample, production)}, {tosample}, TreeWrapper.categorizations)

      for key in tmpresult:
        if isinstance(key, ReweightingSample): continue
        elif isinstance(key, tuple):
#          try:
            assert key[0] in usesamples
            tmpresult[key] /= sum(tmpresult[tosample] for tosample in usesamples)
#          except ZeroDivisionError:
#            pass
        else: assert False

      result += tmpresult


    ####################################
    if productionmode == "ggH":
      tosample = samples.copy().pop()
      try:
        Sample(tosample, production)
      except ValueError:
        pass
      else:
        assert False, "Delete this section!"
      for ctgrztn, category in itertools.product(categorizations, categories):
        result[tosample,ctgrztn,category] = .25 * sum(
            result[ReweightingSamplePlus(tosample, _), ctgrztn, category]
               for _ in pythiasystematics
        )
        for channel in channels:
          result[tosample,ctgrztn,category,channel] = .25 * sum(
              result[ReweightingSamplePlus(tosample, _), ctgrztn, category, channel]
                 for _ in pythiasystematics
          )
    ####################################

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
        YieldValue(channel, category, analysis, productionmode).value = total*sum(result[tosample, categorization, category, channel] for tosample in samples)

      #same for all categories and channels
      for systname in "pdf_Higgs_gg", "pdf_Higgs_qq", "pdf_Higgs_ttH", "pdf_Higgs_qq", "BRhiggs_hzz4l", "QCDscale_ggZH", "QCDscale_ggVV_bonly", "lumi_13TeV":
        for category, channel in itertools.product(categories, channels):
          syst = YieldSystematicValue(channel, category, analysis, productionmode, systname)
          values = syst.yieldsystematic.getfromyaml()
          if productionmode.yamlsystname in values:
            syst.value = values[productionmode.yamlsystname]
          else:
            syst.value = None

      #same for all categories
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

      #same for all channels
      for category in categories:
        nominal = sum(result[tosample, categorization, category] for tosample in samples)
        if productionmode == "ZX":
          JECUp = JECDn = btSFUp = btSFDn = 1
        else:
          JECUp = sum(result[tosample, findsystematic(categorizations, categorization, "JECUp", "Nominal"), category] for tosample in samples) / nominal
          JECDn = sum(result[tosample, findsystematic(categorizations, categorization, "JECDn", "Nominal"), category] for tosample in samples) / nominal
          btSFUp = sum(result[tosample, findsystematic(categorizations, categorization, "Nominal", "bTagSFUp"), category] for tosample in samples) / nominal
          btSFDn = sum(result[tosample, findsystematic(categorizations, categorization, "Nominal", "bTagSFDn"), category] for tosample in samples) / nominal
        for channel in channels:
          YieldSystematicValue(channel, category, analysis, productionmode, "JEC").value = (JECUp, JECDn)
          YieldSystematicValue(channel, category, analysis, productionmode, "BTag").value = (btSFUp, btSFDn)

    YieldValue.writeyieldsdict()
    YieldSystematicValue.writeyieldsystematicsdict()

if __name__ == "__main__":
  writeyields()
