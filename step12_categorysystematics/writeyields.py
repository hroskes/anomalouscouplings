#!/usr/bin/env python

from collections import namedtuple
import itertools
import logging
from math import sqrt
import os
import yaml

from helperstuff import config
from helperstuff.categorization import MultiCategorization
from helperstuff.enums import AlternateWeight, analyses, categories, Category, channels, flavors, ProductionMode, pythiasystematics
from helperstuff.samples import ReweightingSample, ReweightingSamplePlus, ReweightingSampleWithFlavor, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import MultiplyCounter
from helperstuff.yields import count, totalrate, YieldSystematicValue, YieldValue

from categorysystematics import findsystematic

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

SampleCount = namedtuple("SampleCount", "productionmode samples")

def writeyields():
  tosamples_foryields = [
    SampleCount(ProductionMode("ggH"), {ReweightingSamplePlus("ggH", "0+", "POWHEG")}),
    SampleCount(ProductionMode("VBF"), {ReweightingSamplePlus("VBF", "0+", "POWHEG")}),
    SampleCount(ProductionMode("ZH"), {ReweightingSamplePlus("ZH", "0+", "POWHEG")}),
    SampleCount(ProductionMode("WH"), {ReweightingSamplePlus("WplusH", "0+", "POWHEG"), ReweightingSamplePlus("WminusH", "0+", "POWHEG")}),
    SampleCount(ProductionMode("ttH"), {ReweightingSamplePlus("ttH", "0+", "Hff0+", "POWHEG")}),
    SampleCount(ProductionMode("qqZZ"), {ReweightingSample("qqZZ")}),
    SampleCount(ProductionMode("ggZZ"), {ReweightingSampleWithFlavor("ggZZ", flavor) for flavor in flavors}),
    SampleCount(ProductionMode("VBF bkg"), {ReweightingSampleWithFlavor("VBF bkg", flavor) for flavor in ("2e2mu", "4e", "4mu")}),
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
        if (productionmode in ("ggZZ", "VBF bkg", "ZX")
             or hasattr(tosample, "pythiasystematic") and tosample.pythiasystematic is not None):
          usealternateweights = [AlternateWeight("1")]
        else:
          usealternateweights = productionmode.alternateweights

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
      for systname in "pdf_Higgs_gg", "pdf_Higgs_qq", "pdf_Higgs_ttH", "pdf_qq", "BRhiggs_hzz4l", "QCDscale_ggVV_bonly", "lumi_13TeV", "QCDscale_ggH", "QCDscale_qqH", "QCDscale_VH", "QCDscale_ttH", "QCDscale_VV", "EWcorr_VV":
        for category, channel in itertools.product(categories, channels):
          syst = YieldSystematicValue(channel, category, analysis, productionmode, systname)
          syst.value = syst.yieldsystematic.valuefromyaml(productionmode, channel=channel)

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
          YieldSystematicValue(channel, category, analysis, productionmode, "JES").value = (JECUp, JECDn)
          YieldSystematicValue(channel, category, analysis, productionmode, "bTagSF").value = (btSFUp, btSFDn)

        allQCDsystematicnames = {ProductionMode(_).QCDsystematicname for _ in ("ggH", "VBF", "ZH", "WH", "ttH", "qqZZ")}
        allpdfsystematicnames = {ProductionMode(_).pdfsystematicname for _ in ("ggH", "VBF", "ZH", "WH", "ttH", "qqZZ")}

        #QCD weight variations
        for systname in allQCDsystematicnames | allpdfsystematicnames:
          if productionmode == "ggZZ":
            for channel in channels:
              YieldSystematicValue(channel, category, analysis, productionmode, systname).value = YieldSystematicValue(channel, category, analysis, "ggH", systname).value
          elif productionmode == "VBF bkg":
            for channel in channels:
              YieldSystematicValue(channel, category, analysis, productionmode, systname).value = YieldSystematicValue(channel, category, analysis, "VBF", systname).value
          elif systname == productionmode.QCDsystematicname or systname == productionmode.pdfsystematicname:
            if systname == productionmode.QCDsystematicname: first, second = "muR", "muF"
            if systname == productionmode.pdfsystematicname: first, second = "PDF", "alphaS"

            firstup, secondup, firstdn, seconddn = (AlternateWeight(_) for _ in (first+"Up", second+"Up", first+"Dn", second+"Dn"))

            up, dn = {}, {}
            upUntagged, dnUntagged = {}, {}

            up[1] = sum(result[tosample, categorization, firstup, category] for tosample in samples) / nominal
            up[2] = sum(result[tosample, categorization, secondup, category] for tosample in samples) / nominal
            dn[1] = sum(result[tosample, categorization, firstdn, category] for tosample in samples) / nominal
            dn[2] = sum(result[tosample, categorization, seconddn, category] for tosample in samples) / nominal

            #figure out which to combine with which
            upUntagged[1] = sum(result[tosample, categorization, firstup, Category("Untagged")] for tosample in samples) / nominaluntagged
            upUntagged[2] = sum(result[tosample, categorization, secondup, Category("Untagged")] for tosample in samples) / nominaluntagged
            dnUntagged[1] = sum(result[tosample, categorization, firstdn, Category("Untagged")] for tosample in samples) / nominaluntagged
            dnUntagged[2] = sum(result[tosample, categorization, seconddn, Category("Untagged")] for tosample in samples) / nominaluntagged

            if upUntagged[1] >= dnUntagged[1]:
              pass
            elif dnUntagged[1] > upUntagged[1]:
              up[1], dn[1] = dn[1], up[1]
              upUntagged[1], dnUntagged[1] = dnUntagged[1], upUntagged[1]
            else:
              assert False, (upUntagged[1], dnUntagged[1])

            if upUntagged[2] >= dnUntagged[2]:
              pass
            elif dnUntagged[2] > upUntagged[2]:
              up[2], dn[2] = dn[2], up[2]
              upUntagged[2], dnUntagged[2] = dnUntagged[2], upUntagged[2]
            else:
              assert False, (upUntagged[2], dnUntagged[2])

            if upUntagged[2] > 1 and dnUntagged[2] > 1 or upUntagged[2] < 1 and dnUntagged[2] < 1:
              logging.warning("{} up and down go in the same direction {} {}".format(second, upUntagged, dnUntagged))
            if upUntagged[1] > 1 and dnUntagged[1] > 1 or upUntagged[1] < 1 and dnUntagged[1] < 1:
              logging.warning("{} up and down go in the same direction {} {}".format(first, upUntagged, dnUntagged))

            if   up[1] >= 1 and up[2] >= 1 and dn[1] <= 1 and dn[2] <= 1: signup, signdn = +1, -1
            elif up[1] <= 1 and up[2] <= 1 and dn[1] >= 1 and dn[2] >= 1: signup, signdn = -1, +1
            else:
              logging.warning("{} and {} go in different directions in {} {}: {} {}".format(first, second, productionmode, category, up, dn))
              if   sum(up.values()) > sum(dn.values()): signup, signdn = +1, -1
              elif sum(up.values()) < sum(dn.values()): signup, signdn = -1, +1
              else: assert False

            assert signup != signdn

            QCDUp = 1 + signup*sqrt((up[1]-1)**2 + (up[2]-1)**2)
            QCDDn = 1 + signdn*sqrt((dn[1]-1)**2 + (dn[2]-1)**2)

            for channel in channels:
              YieldSystematicValue(channel, category, analysis, productionmode, systname).value = (QCDUp, QCDDn)
          else:
            for channel in channels:
              YieldSystematicValue(channel, category, analysis, productionmode, systname).value = None

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

        if productionmode == "qqZZ":
          EWcorrup = sum(result[tosample, categorization, AlternateWeight("EWcorrUp"), category] for tosample in samples) / nominal
          EWcorrdn = sum(result[tosample, categorization, AlternateWeight("EWcorrDn"), category] for tosample in samples) / nominal
        else:
          EWcorrup = EWcorrdn = 1
        for channel in channels:
          YieldSystematicValue(channel, category, analysis, productionmode, "EWcorr_VV_cat").value = (EWcorrup, EWcorrdn)

    YieldValue.writedict()
    YieldSystematicValue.writedict()

if __name__ == "__main__":
  writeyields()
