#!/usr/bin/env python

from collections import namedtuple
import itertools
import logging
from math import sqrt
import os
import yaml

from helperstuff import config
from helperstuff.categorization import MultiCategorization, NoCategorization
from helperstuff.combinehelpers import gettemplate
from helperstuff.enums import AlternateWeight, analyses, categories, Category, Channel, channels, flavors, ProductionMode, pythiasystematics
from helperstuff.samples import ReweightingSample, ReweightingSamplePlus, ReweightingSampleWithFlavor, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import MultiplyCounter
from helperstuff.yields import count, totalrate, YieldSystematicValue, YieldValue

from categorysystematics import findsystematic

SampleCount = namedtuple("SampleCount", "productionmode samples")

def writeyields():
  if config.LHE: raise ValueError("For LHE you want writeyields_LHE()")

  for production in {_.productionforrate for _ in config.productionsforcombine}:
    tosamples_foryields = [
      SampleCount(ProductionMode("VBF"), {ReweightingSamplePlus("VBF", "0+", "POWHEG")}),
      SampleCount(ProductionMode("ZH"), {ReweightingSamplePlus("ZH", "0+", "POWHEG")}),
      SampleCount(ProductionMode("WH"), {ReweightingSamplePlus("WplusH", "0+", "POWHEG"), ReweightingSamplePlus("WminusH", "0+", "POWHEG")}),
      SampleCount(ProductionMode("ttH"), {ReweightingSamplePlus("ttH", "0+", "Hff0+", "POWHEG")}),
      SampleCount(ProductionMode("ggH"), {ReweightingSamplePlus("ggH", "0+", "POWHEG")}),
      SampleCount(ProductionMode("qqZZ"), {ReweightingSample("qqZZ"), ReweightingSamplePlus("qqZZ", "ext")}),
      SampleCount(ProductionMode("ggZZ"), {ReweightingSampleWithFlavor("ggZZ", flavor) for flavor in flavors}),
    ]
    if production.year == 2016:
      tosamples_foryields.append(SampleCount(ProductionMode("VBF bkg"), {ReweightingSampleWithFlavor("VBF bkg", flavor) for flavor in ("2e2mu", "4e", "4mu")}))
    if config.usedata:
      tosamples_foryields.append(SampleCount(ProductionMode("ZX"), {ReweightingSample("ZX")}))

    categorizations = [_ for _ in TreeWrapper.categorizations if isinstance(_, (MultiCategorization, NoCategorization))]

    result = MultiplyCounter()
    for productionmode, samples in tosamples_foryields:
      print productionmode
      samplegroups = [samples]
      if productionmode.issignal:
        samplegroups += [{ReweightingSamplePlus(s, systematic) for s in samples} for systematic in pythiasystematics
                            if systematic.hassample(production.year)]
        if production.year == 2017:
          samplegroups += [{ReweightingSamplePlus(s, "ext")} for s in samples]

      for usesamples in samplegroups:
        tmpresult = MultiplyCounter()
        for tosample in usesamples:
          if (productionmode in ("ggZZ", "VBF bkg", "ZX")
               or hasattr(tosample, "pythiasystematic") and tosample.pythiasystematic is not None):
            usealternateweights = [AlternateWeight("1")]
          elif productionmode.issignal and tosample.extension == "ext":
            usealternateweights = [AlternateWeight("PythiaScaleUp"), AlternateWeight("PythiaScaleDn")]
          else:
            usealternateweights = productionmode.alternateweights(production.year)

          tmpresult += count({Sample(tosample, production)}, {tosample}, categorizations, usealternateweights)

        result += tmpresult

      #if productionmode.issignal():
      #  result += count(
      #                  {ReweightingSample(productionmode, "0+")},
      #                  {Sample(productionmode, h, production) for h in productionmode.generatedhypotheses},
      #                  TreeWrapper.categorizations,
      #                 )

    result.freeze()


    for productionmode, samples in tosamples_foryields:
      print productionmode
      for analysis in analyses:
        if analysis.isdecayonly and productionmode in ("VBF", "ZH", "WH", "ttH"): continue
        categorization = {_ for _ in categorizations if _.category_function_name == "category_"+analysis.categoryname}
        assert len(categorization) == 1, categorization
        categorization = categorization.pop()

        total = totalrate(productionmode, production, 1.0)
        if analysis.isdecayonly and productionmode == "ggH": total += sum(totalrate(_, production, 1.0) for _ in ("VBF", "ZH", "WH", "ttH"))
        for channel, category in itertools.product(channels, categories):
          YieldValue(channel, category, analysis, productionmode, production).value = total * (
            sum(result[tosample, categorization, AlternateWeight("1"), category, channel] for tosample in samples)
          ) / (
            sum(result[tosample, categorization, AlternateWeight("1"), ca, ch] for tosample in samples for ca in categories for ch in channels)
          )

        #same for all categories and channels
        #from yaml
        for systname in "BRhiggs_hzz4l", "QCDscale_ggVV_bonly":
          for category, channel in itertools.product(categories, channels):
            if analysis.isdecayonly and category != "Untagged": continue
            syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
            syst.value = syst.yieldsystematic.valuefromyaml(productionmode, channel=channel)

        for systname in "lumi_13TeV_2016", "lumi_13TeV_2017":
          syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
          syst.value = syst.yieldsystematic.hardcodedvalue(production)

        #same for all categories
        #from yaml
        for systname in "CMS_eff_e", "CMS_eff_m":
          for category, channel in itertools.product(categories, channels):
            if analysis.isdecayonly and category != "Untagged": continue
            syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
            if channel == "4e" and systname == "CMS_eff_m" or channel == "4mu" and systname == "CMS_eff_e":
              values = {}
            else:
              values = syst.yieldsystematic.getfromyaml(channel=channel)
            if productionmode.yamlsystname in values:
              syst.value = values[productionmode.yamlsystname]
            else:
              syst.value = None

        #Z+X, for now
        for systname in "CMS_hzz4l_zz2e2mu_zjets", "CMS_hzz4l_zz4e_zjets", "CMS_hzz4l_zz4mu_zjets":
          for category, channel in itertools.product(categories, channels):
            if analysis.isdecayonly and category != "Untagged": continue
            syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
            if str(channel) in systname and productionmode == "ZX":
              syst.value = 1.4
            else:
              syst.value = None

        #same for all channels
        for category in categories:
          if analysis.isdecayonly and category != "Untagged": continue
          #variations on category definition: JEC and btagging
          nominal = sum(result[tosample, categorization, AlternateWeight("1"), category] for tosample in samples)
          nominaluntagged = sum(result[tosample, categorization, AlternateWeight("1"), Category("Untagged")] for tosample in samples)
          nominalyield = sum(result[tosample, categorization, AlternateWeight("1")] for tosample in samples)
          if productionmode == "ZX" or analysis.isdecayonly:
            JECUp = JECDn = btSFUp = btSFDn = 1
          else:
            JECUp = sum(result[tosample, findsystematic(categorizations, categorization, "JECUp", "Nominal"), AlternateWeight("1"), category] for tosample in samples) / nominal
            JECDn = sum(result[tosample, findsystematic(categorizations, categorization, "JECDn", "Nominal"), AlternateWeight("1"), category] for tosample in samples) / nominal
            btSFUp = sum(result[tosample, findsystematic(categorizations, categorization, "Nominal", "bTagSFUp"), AlternateWeight("1"), category] for tosample in samples) / nominal
            btSFDn = sum(result[tosample, findsystematic(categorizations, categorization, "Nominal", "bTagSFDn"), AlternateWeight("1"), category] for tosample in samples) / nominal

          for channel in channels:
            YieldSystematicValue(channel, category, analysis, productionmode, "JES", production).value = (JECUp, JECDn)
            YieldSystematicValue(channel, category, analysis, productionmode, "bTagSF", production).value = (btSFUp, btSFDn)

          #QCD muR weight variations

          for p in "ggH", "VBF", "ZH", "WH", "ttH", "qqZZ":
            p = ProductionMode(p)
            for systname, weight in (
                                     (p.QCDfacsystematicname, "muF"),
                                     (p.QCDrensystematicname, "muR"),
                                     (p.pdfvariationsystematicname, "PDF"),
                                     (p.pdfasmzsystematicname, "alphaS"), 
                                    ):
              if productionmode == "ggZZ":
                for channel in channels:
                  YieldSystematicValue(channel, category, analysis, productionmode, systname, production).value = YieldSystematicValue(channel, category, analysis, "ggH", systname, production).value
              elif productionmode == "VBF bkg" and analysis.isdecayonly:
                for channel in channels:
                  YieldSystematicValue(channel, category, analysis, productionmode, systname, production).value = 1
              elif productionmode == "VBF bkg":
                for channel in channels:
                  YieldSystematicValue(channel, category, analysis, productionmode, systname, production).value = YieldSystematicValue(channel, category, analysis, "VBF", systname, production).value
              elif systname in (p.QCDfacsystematicname, p.QCDrensystematicname, p.pdfvariationsystematicname, p.pdfasmzsystematicname):
                up = sum(result[tosample, categorization, AlternateWeight(weight+"Up"), category] for tosample in samples) / nominal)
                dn = sum(result[tosample, categorization, AlternateWeight(weight+"Dn"), category] for tosample in samples) / nominal)
                for channel in channels:
                  YieldSystematicValue(channel, category, analysis, productionmode, systname, production).value = (up, dn)
              else:
                for channel in channels:
                  YieldSystematicValue(channel, category, analysis, productionmode, systname, production).value = None

          #pythia scale and tune
          if productionmode in ("ggH", "VBF", "ZH", "WH", "ttH"):
            if production.year == 2016:
              scaleup = sum(result[ReweightingSamplePlus(tosample, "ScaleUp"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
              scaledn = sum(result[ReweightingSamplePlus(tosample, "ScaleDn"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
            elif production.year == 2017:
              scaleup = sum(result[tosample, categorization, AlternateWeight("PythiaScaleUp"), category] for tosample in samples) / nominal
              scaledn = sum(result[tosample, categorization, AlternateWeight("PythiaScaleUp"), category] for tosample in samples) / nominal
            tuneup = sum(result[ReweightingSamplePlus(tosample, "TuneUp"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
            tunedn = sum(result[ReweightingSamplePlus(tosample, "TuneDn"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
          else:
            scaleup = scaledn = tuneup = tunedn = 1
          for channel in channels:
            YieldSystematicValue(channel, category, analysis, productionmode, "PythiaScale", production).value = (scaleup, scaledn)
            YieldSystematicValue(channel, category, analysis, productionmode, "PythiaTune", production).value = (tuneup, tunedn)

          if productionmode == "qqZZ":
            EWcorrup = sum(result[tosample, categorization, AlternateWeight("EWcorrUp"), category] for tosample in samples) / nominal
            EWcorrdn = sum(result[tosample, categorization, AlternateWeight("EWcorrDn"), category] for tosample in samples) / nominal
          else:
            EWcorrup = EWcorrdn = 1
          for channel in channels:
            YieldSystematicValue(channel, category, analysis, productionmode, "EWcorr_VV", production).value = (EWcorrup, EWcorrdn)

      YieldValue.writedict()
      YieldSystematicValue.writedict()

def writeyields_LHE():
  if not config.LHE: raise ValueError("For non-LHE you want writeyields()")
  for analysis in analyses:
    if not analysis.isfL1fL1Zg: continue
    YieldValue("ggH", "2e2mu", "Untagged", analysis, production).value = 1.5258744890164022
    YieldValue("qqZZ", "2e2mu", "Untagged", analysis, production).value = 1.6250924214670268
  YieldValue.writedict()

if __name__ == "__main__":
  if config.LHE:
    writeyields_LHE()
  else:
    writeyields()
