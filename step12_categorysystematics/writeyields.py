#!/usr/bin/env python

from collections import namedtuple
import itertools

from helperstuff import config
from helperstuff.enums import analyses, categories, channels, flavors, ProductionMode, pythiasystematics
from helperstuff.samples import ReweightingSample, ReweightingSamplePlus, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import MultiplyCounter
from helperstuff.yields import count, totalrate, YieldValue

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
    usesamples = samples.copy()
    if productionmode.issignal:
      usesamples |= {ReweightingSamplePlus(s, systematic) for systematic in pythiasystematics for s in samples}

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
      result += count({Sample(tosample, production)}, {tosample}, TreeWrapper.categorizations)

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
        if len(samples) == 1:
          tosample = list(samples)[0]
          YieldValue(channel, category, analysis, productionmode).value = total*result[tosample, categorization, category, channel]
        else:
          YieldValue(channel, category, analysis, productionmode).value = total*sum(result[tosample, categorization, category, channel]*Sample(tosample, production).xsec for tosample in samples) / sum(Sample(tosample, production).xsec for tosample in samples)

    YieldValue.writeyieldsdict()

if __name__ == "__main__":
  writeyields()
