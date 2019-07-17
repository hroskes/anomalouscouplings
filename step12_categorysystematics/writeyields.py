#!/usr/bin/env python

if __name__ == "__main__":
  import argparse
  def __ProductionMode(*args, **kwargs):
    from helperstuff.enums import ProductionMode
    return ProductionMode(*args, **kwargs)
  def __Production(*args, **kwargs):
    from helperstuff.enums import Production
    return Production(*args, **kwargs)
  parser = argparse.ArgumentParser()
  parser.add_argument("--productionmode", action="append", type=__ProductionMode)
  parser.add_argument("--production", action="append", type=__Production)
  args = parser.parse_args()

from collections import namedtuple
import itertools
import logging
from math import exp, log, sqrt
import os
import pprint
import yaml

from TemplateBuilder.TemplateBuilder.moremath import weightedaverage

from helperstuff import config
from helperstuff.categorization import MultiCategorization, NoCategorization, SingleCategorizationgm4l
from helperstuff.combinehelpers import gettemplate
from helperstuff.enums import AlternateWeight, analyses, categories, Category, Channel, channels, flavors, ProductionMode, pythiasystematics
from helperstuff.samples import ReweightingSample, ReweightingSamplePlus, ReweightingSampleWithFlavor, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import deprecate, MultiplyCounter, sgn
from helperstuff.yields import count, YieldSystematicValue, YieldValue

from categorysystematics import findsystematic

SampleCount = namedtuple("SampleCount", "productionmode samples")

def writeyields(productionmodelist=None, productionlist=None):
  for production in sorted({_.productionforrate for _ in config.productionsforcombine if not _.LHE}):
    if productionlist and production not in productionlist: continue
    tosamples_foryields = [
      SampleCount(ProductionMode("VBF"), [[ReweightingSamplePlus("VBF", "0+", "POWHEG")]]),
      SampleCount(ProductionMode("VH"), [[ReweightingSamplePlus("ZH", "0+", "POWHEG")], [ReweightingSamplePlus("WplusH", "0+", "POWHEG")], [ReweightingSamplePlus("WminusH", "0+", "POWHEG")]]),
      SampleCount(ProductionMode("ttH"), [[ReweightingSamplePlus("ttH", "0+", "Hff0+", "POWHEG")]]),
      SampleCount(ProductionMode("bbH"), [[ReweightingSamplePlus("bbH", "0+")]]),
      SampleCount(ProductionMode("ggH"), [[ReweightingSamplePlus("ggH", "0+", "POWHEG")]]),
      SampleCount(ProductionMode("qqZZ"), [[ReweightingSample("qqZZ"), ReweightingSamplePlus("qqZZ", "ext")]]),
    ] + [
      SampleCount(ProductionMode("ggZZ"), [[ReweightingSampleWithFlavor("ggZZ", flavor)] for flavor in flavors]),
      SampleCount(ProductionMode("VBF bkg"), [[ReweightingSampleWithFlavor("VBF bkg", flavor)] for flavor in ("2e2mu", "4e", "4mu")])
    ][0:deprecate(1, 2019, 7, 20)] * (not production.GEN)

    if config.usedata and not production.GEN:
      tosamples_foryields.append(SampleCount(ProductionMode("ZX"), [[ReweightingSample("ZX")]]))

    categorizations = [_ for _ in TreeWrapper.categorizations if (isinstance(_, (MultiCategorization, NoCategorization)) or isinstance(_, SingleCategorizationgm4l) and _.hypothesis == "0+") and (not production.GEN or not _.issystematic)]

    result = MultiplyCounter()
    for productionmode, samplegroups in tosamples_foryields:
      if productionmodelist and productionmode not in productionmodelist: continue
      print productionmode
      for g in samplegroups:
        if productionmode.issignal and productionmode != "bbH" and not production.GEN:
          g += [
            ReweightingSamplePlus(s, systematic)
            for s in g
            for systematic in pythiasystematics
            if systematic.hassample(production.year)
          ]
          if production.year == 2017:
            g += [ReweightingSamplePlus(s, "ext") for s in g if s.pythiasystematic is None]
            if productionmode == "ggH":
              g += [ReweightingSamplePlus(s.reweightingsample, "MINLO") for s in g if s.pythiasystematic is None and s.extension is None]
        if productionmode == "qqZZ":
          if production.year == 2018:
            del g[:]
            g += [ReweightingSamplePlus("qqZZ", "ext1"), ReweightingSamplePlus("qqZZ", "ext2")]

      for usesamples in samplegroups:
        tmpresults = []
        for tosample in usesamples:
          if (productionmode in ("ggZZ", "VBF bkg", "ZX", "bbH")
               or hasattr(tosample, "pythiasystematic") and tosample.pythiasystematic is not None
               or hasattr(tosample, "alternategenerator") and tosample.alternategenerator == "MINLO"
               or production.GEN
             ):
            usealternateweights = [AlternateWeight("1")]
          elif productionmode.issignal and tosample.extension == "ext":
            usealternateweights = [AlternateWeight("PythiaScaleUp"), AlternateWeight("PythiaScaleDn")]
          else:
            usealternateweights = productionmode.alternateweights(production.year)

          tmpresults.append(count({Sample(tosample, production)}, {tosample}, categorizations, usealternateweights))

        try:
          result += MultiplyCounter({
            k:
            weightedaverage(
              tmpresult[k]
              for tmpresult in tmpresults
              if k in tmpresult
            )
            for k in frozenset.union(*(frozenset(tmpresult) for tmpresult in tmpresults))
          })
        except ZeroDivisionError:
          pprint.pprint([dict(_) for _ in tmpresults])
          raise

    result.freeze()

    for productionmode, samples in tosamples_foryields:
      if productionmodelist and productionmode not in productionmodelist: continue
      print productionmode
      for analysis in analyses:
        if analysis.isdecayonly and productionmode in ("VBF", "ZH", "WH", "VH", "ttH", "bbH"): continue
        categorization = {_ for _ in categorizations if _.category_function_name == "category_"+analysis.categoryname}
        assert len(categorization) == 1, categorization
        categorization = categorization.pop()

        for channel, category in itertools.product(channels, categories):
          if analysis.isdecayonly and category != "Untagged": continue
          yv = YieldValue(channel, category, analysis, productionmode, production)

          yvvalue = 0
          
          for g in samples:
            key = type(g[0])(*(getattr(g[0], needenum.enumname) for needenum in g[0].needenums if needenum.enumname != "extension")), categorization, AlternateWeight("1"), category, channel
            try:
              yvvalue += result[key].nominal_value
            except AttributeError:
              #result.key == 0, it's an int, so it doesn't have nominal_value
              if hasattr(g[0], "flavor") and str(channel) != g[0].flavor:
                  pass #e.g. ggZZ 2e2mu in the 4e channel: doesn't happen
              else:
                  #sanity check
                  pprint.pprint(dict(result))
                  print key
                  print result[key]
                  raise

          yv.value = yvvalue

        if production.GEN or deprecate(True, 2019, 7, 18): continue

        #same for all categories and channels
        for category, channel in itertools.product(categories, channels):
          if analysis.isdecayonly and category != "Untagged": continue

          syst = YieldSystematicValue(channel, category, analysis, productionmode, "BRhiggs_hzz4l", production)
          if productionmode.issignal:
            syst.value = 1.02
          else:
            syst.value = None

          syst = YieldSystematicValue(channel, category, analysis, productionmode, "CMS_eff_e", production)
          if productionmode == "ZX" or channel == "4mu":
            syst.value = None
          elif channel == "2e2mu":
            syst.value = 0.96, 1.039
          elif channel == "4e":
            syst.value = 0.914, 1.082

          syst = YieldSystematicValue(channel, category, analysis, productionmode, "CMS_eff_m", production)
          if productionmode == "ZX" or channel == "4e":
            syst.value = None
          elif channel == "2e2mu":
            syst.value = 1.025
          elif channel == "4mu":
            syst.value = 0.953, 1.046

          for systname in "lumi_13TeV_2016", "lumi_13TeV_2017":
            syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
            if productionmode == "ZX" or str(production.year) not in systname:
              syst.value = None
            elif production.year == 2016:
              syst.value = 1.026
            elif production.year == 2017:
              syst.value = 1.023

          for systname in "CMS_hzz4l_zz2e2mu_zjets", "CMS_hzz4l_zz4e_zjets", "CMS_hzz4l_zz4mu_zjets":
            syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
            if productionmode != "ZX" or str(channel) not in systname:
              syst.value = None
            elif production.year == 2016:
              syst.value = {
                (Channel("2e2mu"), Category(    "Untagged")): (0.746, 1.340),
                (Channel("2e2mu"), Category(   "VBFtagged")): (0.746, 1.341),
                (Channel("2e2mu"), Category("VHHadrtagged")): (0.746, 1.341),
                (Channel("4e"   ), Category(    "Untagged")): (0.758, 1.320),
                (Channel("4e"   ), Category(   "VBFtagged")): (0.758, 1.321),
                (Channel("4e"   ), Category("VHHadrtagged")): (0.758, 1.321),
                (Channel(  "4mu"), Category(    "Untagged")): (0.741, 1.350),
                (Channel(  "4mu"), Category(   "VBFtagged")): (0.740, 1.351),
                (Channel(  "4mu"), Category("VHHadrtagged")): (0.740, 1.352),
              }[channel, category]
            elif production.year == 2017:
              syst.value = {
                (Channel("2e2mu"), Category(    "Untagged")): (0.769, 1.300),
                (Channel("2e2mu"), Category(   "VBFtagged")): (0.768, 1.301),
                (Channel("2e2mu"), Category("VHHadrtagged")): (0.768, 1.301),
                (Channel("4e"   ), Category(    "Untagged")): (0.769, 1.300),
                (Channel("4e"   ), Category(   "VBFtagged")): (0.768, 1.301),
                (Channel("4e"   ), Category("VHHadrtagged")): (0.769, 1.301),
                (Channel(  "4mu"), Category(    "Untagged")): (0.769, 1.300),
                (Channel(  "4mu"), Category(   "VBFtagged")): (0.768, 1.302),
                (Channel(  "4mu"), Category("VHHadrtagged")): (0.768, 1.302),
              }[channel, category]
            else:
              assert False

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
            YieldSystematicValue(channel, category, analysis, productionmode, "CMS_scale_j_13TeV_{}".format(production.year), production).value = (JECDn, JECUp)
            YieldSystematicValue(channel, category, analysis, productionmode, "CMS_btag_comb_13TeV_{}".format(production.year), production).value = (btSFDn, btSFUp)
            for year in 2016, 2017, 2018:
              if year == production.year: continue
              YieldSystematicValue(channel, category, analysis, productionmode, "CMS_scale_j_13TeV_{}".format(year), production).value = None
              YieldSystematicValue(channel, category, analysis, productionmode, "CMS_btag_comb_13TeV_{}".format(year), production).value = None

          #QCD and PDF weight variations

          seen = set()
          for p in "ggH", "VBF", "ZH", "WH", "ttH", "qqZZ", "bbH":
            p = ProductionMode(p)
            for systname, weight in (
                                     (p.QCDfacsystematicname, "muF"),
                                     (p.QCDrensystematicname, "muR"),
                                     (p.pdfvariationsystematicname, "PDF"),
                                     (p.pdfasmzsystematicname, "alphaS"),
                                    ):
              if systname in seen: continue
              seen.add(systname)
              if productionmode == "ggZZ":
                for channel in channels:
                  syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
                  syst.value = YieldSystematicValue(channel, category, analysis, "ggH", systname, production).value
              elif productionmode == "VBF bkg" and analysis.isdecayonly:
                for channel in channels:
                  syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
                  syst.value = 1
              elif productionmode == "VBF bkg":
                for channel in channels:
                  syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
                  addsysts = [(YieldValue(channel, category, analysis, _, production).value,
                              YieldSystematicValue(channel, category, analysis, _, systname, production).value) for _ in "VBF", "ZH", "WH"]
                  addfirsts = [(_1, _2[0] if _2 is not None else 1) for _1, _2 in addsysts]
                  addseconds = [(_1, _2[1] if _2 is not None else 1) for _1, _2 in addsysts]
                  if addfirsts and addseconds:
                    first = sum(_1*_2 for _1, _2 in addfirsts) / sum(_1 for _1, _2 in addfirsts)
                    second = sum(_1*_2 for _1, _2 in addseconds) / sum(_1 for _1, _2 in addseconds)
                    #maxfirst = max(abs(log(_)) for _ in addfirsts)
                    #maxsecond = max(abs(log(_)) for _ in addseconds)
                    #addfirsts = [_ for _ in addfirsts if abs(log(_)) >= maxfirst/5]
                    #addseconds = [_ for _ in addseconds if abs(log(_)) >= maxsecond/5]
                    #assert len({sgn(log(_)) for _ in addfirsts}) <= 1, addfirsts
                    #assert len({sgn(log(_)) for _ in addseconds}) <= 1, addseconds
                    #first = exp(sqrt(sum(log(_)**2 for _ in addfirsts)))
                    #if addfirsts[0] < 1: first = 1/first
                    #second = exp(sqrt(sum(log(_)**2 for _ in addseconds)))
                    #if addseconds[0] < 1: second = 1/second
                    #print addfirsts, first
                  else:
                    first = second = 1
                  syst.value = first, second
              elif systname in (productionmode.QCDfacsystematicname, productionmode.QCDrensystematicname, productionmode.pdfvariationsystematicname, productionmode.pdfasmzsystematicname):
                if productionmode == "bbH":
                  dn, up = {
                    "QCDscale_ren_bbH": (1.128, 0.837),
                    "QCDscale_fac_bbH": (1.078, 0.96),
                    "pdf_asmz_Higgs_gg": (0.945, 1.075),
                    "pdf_variation_Higgs_gg": (1.113, 0.922),
                  }[systname]
                elif productionmode == "ttH" and systname == "pdf_asmz_Higgs_gg":
                  dn, up = 0.98, 1.02
                else:
                  up = sum(result[tosample, categorization, AlternateWeight(weight+"Up"), category] for tosample in samples) / nominal
                  dn = sum(result[tosample, categorization, AlternateWeight(weight+"Dn"), category] for tosample in samples) / nominal
                for channel in channels:
                  syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
                  syst.value = (dn, up)
              else:
                for channel in channels:
                  syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
                  syst.value = None

          #pythia scale and tune
          if productionmode in ("ggH", "VBF", "ZH", "WH", "ttH"):
            if production.year == 2016:
              scaleup = sum(result[ReweightingSamplePlus(tosample, "ScaleUp"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
              scaledn = sum(result[ReweightingSamplePlus(tosample, "ScaleDn"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
            elif production.year == 2017 or production.year == 2018:
              scaleup = sum(result[ReweightingSamplePlus(tosample, "ext"), categorization, AlternateWeight("PythiaScaleUp"), category] for tosample in samples) / nominal
              scaledn = sum(result[ReweightingSamplePlus(tosample, "ext"), categorization, AlternateWeight("PythiaScaleDn"), category] for tosample in samples) / nominal
            tuneup = sum(result[ReweightingSamplePlus(tosample, "TuneUp"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
            tunedn = sum(result[ReweightingSamplePlus(tosample, "TuneDn"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
          else:
            scaleup = scaledn = tuneup = tunedn = 1
          for channel in channels:
            YieldSystematicValue(channel, category, analysis, productionmode, "CMS_scale_pythia", production).value = (scaledn, scaleup)
            YieldSystematicValue(channel, category, analysis, productionmode, "CMS_tune_pythia", production).value = (tunedn, tuneup)

          if productionmode == "ggH":
            if production.year == 2016:
              minloup = minlodn = 1
            elif production.year == 2017:
              minloup = sum(result[ReweightingSamplePlus(tosample.reweightingsample, "MINLO"), categorization, AlternateWeight("1"), category] for tosample in samples) / nominal
              minlodn = 1 / minloup
          else:
            minloup = minlodn = 1
          for channel in channels:
            YieldSystematicValue(channel, category, analysis, productionmode, "QCDscale_ggH2in", production).value = (minlodn, minloup)

          if productionmode == "qqZZ":
            EWcorrup = sum(result[tosample, categorization, AlternateWeight("EWcorrUp"), category] for tosample in samples) / nominal
            EWcorrdn = sum(result[tosample, categorization, AlternateWeight("EWcorrDn"), category] for tosample in samples) / nominal
          else:
            EWcorrup = EWcorrdn = 1
          for channel in channels:
            YieldSystematicValue(channel, category, analysis, productionmode, "EWcorr_VV", production).value = (EWcorrdn, EWcorrup)

      YieldValue.writedict()
      YieldSystematicValue.writedict()

def writeyields_LHE():
  for production in sorted({_.productionforrate for _ in config.productionsforcombine if _.LHE}):
    for analysis in analyses:
      if not analysis.isfL1fL1Zg: continue
      YieldValue("ggH", "2e2mu", "Untagged", analysis, production).value = 1.5258744890164022
      YieldValue("qqZZ", "2e2mu", "Untagged", analysis, production).value = 1.6250924214670268
  YieldValue.writedict()

if __name__ == "__main__":
  writeyields_LHE()
  writeyields(productionmodelist=args.productionmode, productionlist=args.production)
