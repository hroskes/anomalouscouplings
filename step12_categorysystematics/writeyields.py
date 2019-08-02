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
from helperstuff.enums import AlternateWeight, analyses, categories, Category, Channel, channels, flavors, ProductionMode, PythiaSystematic, pythiasystematics
from helperstuff.samples import ReweightingSamplePlusWithFlavor as RSPWF, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import deprecate, MultiplyCounter, sgn
from helperstuff.yields import count, YieldSystematicValue, YieldValue

from categorysystematics import findsystematic

SampleCount = namedtuple("SampleCount", "productionmode samples")

def writeyields(productionmodelist=None, productionlist=None):
  for production in sorted({_.productionforrate for _ in config.productionsforcombine if not _.LHE}):
    if productionlist and production not in productionlist: continue
    print "Finding yields and category systematics for", production
    year = production.year

    tosamples_foryields = [
      SampleCount(ProductionMode("VBF"), [[RSPWF("VBF", "0+", "POWHEG")]]),
      SampleCount(ProductionMode("VH"), [[RSPWF("ZH", "0+", "POWHEG")], [RSPWF("WplusH", "0+", "POWHEG")], [RSPWF("WminusH", "0+", "POWHEG")]]),
      SampleCount(ProductionMode("ZH"), [[RSPWF("ZH", "0+", "POWHEG")]]),
      SampleCount(ProductionMode("WH"), [[RSPWF("WplusH", "0+", "POWHEG")], [RSPWF("WminusH", "0+", "POWHEG")]]),
      SampleCount(ProductionMode("ttH"), [[RSPWF("ttH", "0+", "Hff0+", "POWHEG")]]),
      SampleCount(ProductionMode("bbH"), [[RSPWF("bbH", "0+")]]),
      SampleCount(ProductionMode("ggH"), [[RSPWF("ggH", "0+", "POWHEG")]]),
      SampleCount(ProductionMode("qqZZ"), [[RSPWF("qqZZ"), RSPWF("qqZZ", "ext")]]),
    ] + [
      SampleCount(ProductionMode("ggZZ"), [[RSPWF("ggZZ", flavor)] for flavor in flavors]),
      SampleCount(ProductionMode("VBF bkg"), [[RSPWF("VBF bkg", flavor)] for flavor in ("2e2mu", "4e", "4mu")])
    ][0:deprecate(1, 2019, 8, 5)] * (not production.GEN)

    if config.usedata and not production.GEN:
      tosamples_foryields.append(SampleCount(ProductionMode("ZX"), [[RSPWF("ZX")]]))

    categorizations = [_ for _ in TreeWrapper.categorizations if (isinstance(_, (MultiCategorization, NoCategorization)) or isinstance(_, SingleCategorizationgm4l) and _.hypothesis == "0+") and (not production.GEN or not _.issystematic)]

    result = MultiplyCounter()

    print
    print "Iterating through the trees"
    print
    for productionmode, samplegroups in tosamples_foryields:
      if productionmodelist and productionmode not in productionmodelist: continue
      print productionmode
      for g in samplegroups:
        if productionmode.issignal and productionmode != "bbH" and not production.GEN:
          g += [
            RSPWF(s, systematic)
            for s in g
            for systematic in pythiasystematics
            if systematic.hassample(year)
          ]
          if year == 2017:
            g += [RSPWF(s, "ext") for s in g if s.pythiasystematic is None]
            if productionmode == "ggH":
              g += [RSPWF(s.reweightingsample, "MINLO") for s in g if s.pythiasystematic is None and s.extension is None]
        if productionmode == "qqZZ":
          if year == 2018:
            del g[:]
            g += [RSPWF("qqZZ", "ext1"), RSPWF("qqZZ", "ext2")]

      if productionmode in ("ZH", "WH") and (not productionmodelist or ProductionMode("VH") in productionmodelist): continue #it's already in the counter from VH

      for usesamples in samplegroups:
        tmpresults = []
        for tosample in usesamples:
          usealternateweights = Sample(tosample, production).alternateweights
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

    resultbypm = MultiplyCounter()

    print
    print "Collecting the results"
    print

    for productionmode, samples in tosamples_foryields:
      if productionmodelist and productionmode not in productionmodelist and not (productionmode in ("ZH", "WH") and "VH" in productionmodelist): continue
      if productionmode == "VH": continue
      print productionmode
      bkpresult = result.copy()
      for g in samples:
        for sample in g:
          sampleforkey = type(sample)(
            *(
              getattr(sample, needenum.enumname)
              for needenum in sample.needenums
              if needenum.enumname != "extension"
             )
          )
          for key in result.keys():
            if key[0] == sampleforkey:
              newkey = (
                (productionmode,)
                + tuple(
                  getattr(sampleforkey, needenum.enumname)
                  for needenum in sampleforkey.needenums
                  if getattr(sampleforkey, needenum.enumname) != getattr(g[0], needenum.enumname)
                  and getattr(sampleforkey, needenum.enumname) is not None
                )
                + key[1:]
              )
              value = result.pop(key)
              resultbypm[newkey] += value
              if productionmode in ("ZH", "WH"):
                newkey2 = (ProductionMode("VH"),) + newkey[1:]
                resultbypm[newkey2] += value

    if result:
      pprint.pprint(dict(result))
      assert False

    result = dict(resultbypm)
    del resultbypm

    print
    print "Calculating and writing to files"
    print

    for productionmode, samples in tosamples_foryields:
      if productionmodelist and productionmode not in productionmodelist: continue
      print productionmode
      for analysis in analyses:
        if analysis.isdecayonly and productionmode in ("VBF", "ZH", "WH", "VH", "ttH", "bbH"): continue
        categorization = {_ for _ in categorizations if _.category_function_name == "category_"+analysis.categoryname}
        assert len(categorization) == 1, categorization
        categorization = categorization.pop()

        for category in categories:
          if analysis.isdecayonly and category != "Untagged": continue
          sumcategories = [category]
          if not analysis.useboosted and not analysis.isdecayonly:
            if category == "Boosted": continue
            if category == "Untagged": sumcategories.append(Category("Boosted"))

          for channel in channels:
            yv = YieldValue(channel, category, analysis, productionmode, production)

            yvvalue = sum(result[productionmode, categorization, AlternateWeight("1"), cat, channel] for cat in sumcategories).nominal_value

            yv.value = yvvalue

            if production.GEN: continue

            syst = YieldSystematicValue(channel, category, analysis, productionmode, "hzz_br", production)
            if productionmode.isbkg:
              syst.value = None
            else:
              syst.value = 1.02

            syst = YieldSystematicValue(channel, category, analysis, productionmode, "CMS_eff_e", production)
            if productionmode == "ZX":
              syst.value = None
            else:
              syst.value = {
                (2016, "2e2mu"): (1.039, 0.960),
                (2016, "4e"):    (1.082, 0.914),
                (2016, "4mu"):   None,
                (2017, "2e2mu"): (1.058, 0.939),
                (2017, "4e"):    (1.125, 0.862),
                (2017, "4mu"):   None,
                (2018, "2e2mu"): 1.074,
                (2018, "4e"):    1.161,
                (2018, "4mu"):   None,
              }[year, str(channel)]

            syst = YieldSystematicValue(channel, category, analysis, productionmode, "CMS_eff_m", production)
            if productionmode == "ZX":
              syst.value = None
            else:
              syst.value = {
                (2016, "2e2mu"): (1.025, 0.975),
                (2016, "4e"):    None,
                (2016, "4mu"):   (1.046, 0.953),
                (2017, "2e2mu"): (1.030, 0.968),
                (2017, "4e"):    None,
                (2017, "4mu"):   (1.056, 0.937),
                (2018, "2e2mu"): (1.011, 0.992),
                (2018, "4e"):    None,
                (2018, "4mu"):   (1.016, 0.978),
              }[year, str(channel)]

            for systname in "lumi_13TeV_2016", "lumi_13TeV_2017", "lumi_13TeV_2018":
              syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
              if productionmode == "ZX" or str(year) not in systname:
                syst.value = None
              else:
                syst.value = {
                  2016: (1.026, 0.974),
                  2017: (1.023, 0.977),
                  2018: (1.025, 0.975),
                }[year]

            for yr in 2016, 2017, 2018:
              for chan in "2e2mu", "4e", "4mu":
                systname = "zjet_{}_{}".format(chan, yr)
                syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
                if productionmode != "ZX" or chan != channel or yr != year:
                  syst.value = None
                else:
                  syst.value = {
                    (2016, "2e2mu"): (1.152, 0.868),
                    (2016, "4e"):    (1.314, 0.728),
                    (2016, "4mu"):   (1.104, 0.899),
                    (2017, "2e2mu"): (1.330, 0.670),
                    (2017, "4e"):    (1.380, 0.640),
                    (2017, "4mu"):   (1.320, 0.680),
                    (2018, "2e2mu"): (1.300, 0.700),
                    (2018, "4e"):    (1.370, 0.630),
                    (2018, "4mu"):   (1.240, 0.760),
                  }[year, str(channel)]

          #variations on category definition: JEC and btagging
          nominal = sum(result[productionmode, categorization, AlternateWeight("1"), cat] for cat in sumcategories)
          if productionmode == "ZX" or analysis.isdecayonly:
            JECUp = JECDn = btSFUp = btSFDn = 1
          else:
            JECUp = (sum(result[productionmode, findsystematic(categorizations, categorization, "JECUp", "Nominal"), AlternateWeight("1"), cat] for cat in sumcategories) / nominal).nominal_value
            JECDn = (sum(result[productionmode, findsystematic(categorizations, categorization, "JECDn", "Nominal"), AlternateWeight("1"), cat] for cat in sumcategories) / nominal).nominal_value
            btSFUp = (sum(result[productionmode, findsystematic(categorizations, categorization, "Nominal", "bTagSFUp"), AlternateWeight("1"), cat] for cat in sumcategories) / nominal).nominal_value
            btSFDn = (sum(result[productionmode, findsystematic(categorizations, categorization, "Nominal", "bTagSFDn"), AlternateWeight("1"), cat] for cat in sumcategories) / nominal).nominal_value

          for channel in channels:
            YieldSystematicValue(channel, category, analysis, productionmode, "CMS_scale_j", production).value = (JECDn, JECUp)
            YieldSystematicValue(channel, category, analysis, productionmode, "CMS_btag_comb", production).value = (btSFDn, btSFUp)

          #QCD and PDF weight variations

          seen = set()
          for p in "ggH", "VBF", "ZH", "WH", "ttH", "qqZZ", "bbH":
            p = ProductionMode(p)
            for systname, weight in (
                                     (p.QCDmuFsystematicname, "muF"),
                                     (p.QCDmuRsystematicname, "muR"),
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
                  else:
                    first = second = 1
                  syst.value = first, second
              elif systname in (productionmode.QCDmuFsystematicname, productionmode.QCDmuRsystematicname, productionmode.pdfvariationsystematicname, productionmode.pdfasmzsystematicname):
                up = (sum(result[productionmode, categorization, AlternateWeight(weight+"Up"), cat] for cat in sumcategories) / nominal).nominal_value
                dn = (sum(result[productionmode, categorization, AlternateWeight(weight+"Dn"), cat] for cat in sumcategories) / nominal).nominal_value
                for channel in channels:
                  syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
                  syst.value = (dn, up)
              else:
                for channel in channels:
                  syst = YieldSystematicValue(channel, category, analysis, productionmode, systname, production)
                  syst.value = None

          #pythia scale and tune
          if productionmode in ("ggH", "VBF", "ZH", "WH", "ttH"):
            if year == 2016:
              scaleup = (sum(result[productionmode, PythiaSystematic("ScaleUp"), categorization, AlternateWeight("1"), cat] for cat in sumcategories) / nominal).nominal_value
              scaledn = (sum(result[productionmode, PythiaSystematic("ScaleDn"), categorization, AlternateWeight("1"), cat] for cat in sumcategories) / nominal).nominal_value
              if productionmode == "ttH": scaleup = scaledn = deprecate(1, 2019, 8, 10)
            elif year == 2017 or year == 2018:
              scaleup = (sum(result[productionmode, categorization, AlternateWeight("PythiaScaleUp"), cat] for cat in sumcategories) / nominal).nominal_value
              scaledn = (sum(result[productionmode, categorization, AlternateWeight("PythiaScaleDn"), cat] for cat in sumcategories) / nominal).nominal_value
            tuneup = (sum(result[productionmode, PythiaSystematic("TuneUp"), categorization, AlternateWeight("1"), cat] for cat in sumcategories) / nominal).nominal_value
            tunedn = (sum(result[productionmode, PythiaSystematic("TuneDn"), categorization, AlternateWeight("1"), cat] for cat in sumcategories) / nominal).nominal_value
          else:
            scaleup = scaledn = tuneup = tunedn = 1
          for channel in channels:
            YieldSystematicValue(channel, category, analysis, productionmode, "CMS_pythia_scale", production).value = (scaledn, scaleup)
            YieldSystematicValue(channel, category, analysis, productionmode, "CMS_pythia_tune", production).value = (tunedn, tuneup)

          if productionmode == "qqZZ":
            EWcorrup = (sum(result[productionmode, categorization, AlternateWeight("EWcorrUp"), cat] for cat in sumcategories) / nominal).nominal_value
            EWcorrdn = (sum(result[productionmode, categorization, AlternateWeight("EWcorrDn"), cat] for cat in sumcategories) / nominal).nominal_value
          else:
            EWcorrup = EWcorrdn = 1
          for channel in channels:
            YieldSystematicValue(channel, category, analysis, productionmode, "EWcorr_qqZZ", production).value = (EWcorrdn, EWcorrup)

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
