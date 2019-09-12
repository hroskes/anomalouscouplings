#!/usr/bin/env python

if __name__ == "__main__":
  import argparse
  p = argparse.ArgumentParser()
  args = p.parse_args()

from abc import ABCMeta, abstractmethod, abstractproperty
from itertools import chain, product
import random
import sys

import numpy as np
import uncertainties

import ROOT

from helperstuff import config

from helperstuff.combinehelpers import getdatatree, getrate, Luminosity
from helperstuff.enums import Category, Channel, EnumItem, flavors, Hypothesis, JECSystematic, MultiEnum, MultiEnumABCMeta, MyEnum, ProductionMode
from helperstuff.samples import ReweightingSample, ReweightingSampleWithPdf, Sample
from helperstuff.utilities import KeyDefaultDict, MultiplyCounter, tfiles

from discriminantplots import getscaletos, gettrees, getweights, makehistogram, makehistogramnormalization, makehistogramrescaling, HypothesisLine, masscut

categories = list(Category(_) for _ in ("VBFtagged", "VHHadrtagged", "VHLepttagged", "VBF1jtagged", "Boosted", "Untagged"))
channels = list(Channel(_) for _ in ("4e", "4mu", "2e2mu"))
categorizationname = "category_0P_or_0M_or_a2_or_L1_or_L1Zg"

def categoryname(category):
  if category == "VBFtagged": return "VBF-2jet-tagged"
  if category == "VBF1jtagged": return "VBF-1jet-tagged"
  if category == "VHHadrtagged": return "$\V\PH$-hadronic-tagged"
  if category == "VHLepttagged": return "$\V\PH$-leptonic-tagged"
  if category == "Untagged": return "Untagged"
  if category == "Boosted": return "Boosted"
  assert False, category

class RowBaseBase(object):
  __metaclass__ = ABCMeta

  def __init__(self, *args, **kwargs):
    self.__categorydistribution = None
    super(RowBaseBase, self).__init__(*args, **kwargs)

  @property
  def categorydistribution(self):
    if self.__categorydistribution is None:
      self.__categorydistribution = self.getcategorydistribution()
    return self.__categorydistribution

  @abstractproperty
  def title(self): pass
  @abstractproperty
  def fmt(self): pass
  @abstractmethod
  def getlatex(self): pass
  @abstractmethod
  def getcategorydistribution(self): pass

class RowBase(RowBaseBase):
  def getlatex(self):
    result = "{}".format(self.title)
    for category in list(categories):
      result += " & "
      amount = self.categorydistribution[category]
      result += self.fmt.format(amount)
    return result

class Row(RowBase, MultiEnum):
  __metaclass__ = MultiEnumABCMeta
  enums = (ProductionMode, Hypothesis)
  enumname = "Row"
  def __init__(self, *args, **kwargs):
    self.__title = kwargs.pop("title", None)
    print self.__title
    super(Row, self).__init__(*args, **kwargs)
    if self.__title is None:
      self.__title = str(self.productionmode)

    kwargs = dict(
      xformula="1",
      binning=np.array([0., 2.]),
      linecolor=0, linestyle=0, linewidth=0, fillcolor=0, fillstyle=0, legendname="", legendlpf="",
      addonbottom=[],
      mirror=False,
    )
    if self.productionmode == "data":
      kwargs.update(
        trees=gettrees(("data",)),
        weightformulas=["1" for production in config.productionsforcombine],
        scaletos=None,
      )
    elif self.productionmode.isbkg:
      kwargs.update(
        trees=gettrees((str(self.productionmode),)),
        weightformulas=getweights((self.productionmode,), uselumi=self.productionmode != "ZX"),
        scaletos=None,
      )
    elif self.productionmode.issignal:
      hypothesisline = HypothesisLine(str(self.hypothesis), None, None, None, None, None)
      SMhypothesisline = HypothesisLine("0+", None, None, None, None, None)
      isL1Zg = self.hypothesis == "L1Zg"
      if self.productionmode == "ttH":
        otherargs = (("ttH", "Hff0+"),)
      elif self.productionmode == "VH" and isL1Zg:
        otherargs = (("ZH",),)
      else:
        otherargs = ((str(self.productionmode),),)

      if self.productionmode in ("ggH", "ttH", "bbH"):
        trees = gettrees(
          ("ggH", "0+"),
          ("bbH", "0+"),
          ("ttH", "0+", "Hff0+"),
        )
        SMweights = SMhypothesisline.ffHweights()
        weights = hypothesisline.ffHweights()
        SMhistogramnormalization = makehistogramnormalization(
          name="SMffHnormalization",
          trees=trees,
          weightformulas=SMweights,
          scaletos=getscaletos("ggH", "bbH", "ttH"),
        )
        histogramnormalization = makehistogramnormalization(
          name=str(self.hypothesis)+"ffHnormalization",
          trees=trees,
          weightformulas=weights,
          scaletos=getscaletos("ggH", "bbH", "ttH"),
          normalizationweightformulas=SMweights,
        )
      elif self.productionmode in ("VBF", "VH", "ZH", "WH"):
        SMargs = args = (
          ("VBF",),
          ("VH",),
        )
        if isL1Zg:
          args = (
            ("VBF",),
            ("ZH",),
          )
        SMweights = SMhypothesisline.getweights(*SMargs)
        SMtrees = SMhypothesisline.gettrees(*SMargs)
        weights = hypothesisline.getweights(*args)
        trees = hypothesisline.gettrees(*args)
        SMhistogramnormalization = makehistogramnormalization(
          name="SMVVHnormalization",
          trees=SMtrees,
          weightformulas=SMweights,
          scaletos=getscaletos(*[_[0] for _ in SMargs]),
        )
        histogramnormalization = makehistogramnormalization(
          name=str(self.hypothesis)+"VVHnormalization",
          trees=trees,
          weightformulas=weights,
          normalizationweightformulas=SMhypothesisline.getweights(*args),
          normalizationtrees=SMhypothesisline.gettrees(*args),
          scaletos=getscaletos(*[_[0] for _ in args], isL1Zg=isL1Zg),
        )
      rescaling = makehistogramrescaling(numerator=SMhistogramnormalization, denominator=histogramnormalization)


      kwargs.update(
        trees=hypothesisline.gettrees(*otherargs),
        weightformulas=hypothesisline.getweights(*otherargs),
        normalizationweightformulas=SMhypothesisline.getweights(*otherargs),
        normalizationtrees=SMhypothesisline.gettrees(*otherargs),
        scaletos=getscaletos(self.productionmode, isL1Zg=isL1Zg),
        rescalings=[rescaling for _ in hypothesisline.gettrees(*otherargs)]
      )
    else:
      assert False, self.productionmode

    self.__histograms = {
      category: makehistogram(
        name="dummy_{}_{}_{}".format(self.productionmode, self.hypothesis, category),
        cutformula=masscut + " && (" + " || ".join("{} == {}".format(categorizationname, c) for c in category.idnumbers) + ")",
        **kwargs
      ) for category in categories
    }

  def check(self, *args):
    dontcheck = []
    if self.productionmode == "data" or self.productionmode.isbkg:
      dontcheck.append(Hypothesis)
    super(Row, self).check(*args, dontcheck=dontcheck)

  @property
  def title(self):
    return self.__title
  @property
  def fmt(self):
    if self.productionmode == "data":
      return "{:d}"
    else:
      return "{:.1f}"

  def getcategorydistribution(self):
    result = MultiplyCounter()
    for category, histogram in self.__histograms.iteritems():
      histogram.makefinalhistogram()
      result[category] = histogram.histogram.Integral()
      if self.productionmode == "data":
        assert result[category] == int(result[category])
        result[category] = int(result[category])
    return result

class TotalRow(RowBase):
  def __init__(self, *rows, **kwargs):
    self.rows = rows
    self.__title = "Expected"
    for kw, kwarg in kwargs.iteritems():
      if kw == "title":
        self.__title = kwarg
      else:
        raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))
    super(TotalRow, self).__init__()
  @property
  def title(self):
    return self.__title
  def getcategorydistribution(self):
    return sum((row.categorydistribution for row in self.rows), MultiplyCounter())
  @property
  def fmt(self):
    return "{:.1f}"

class SlashRow(RowBaseBase):
  def __init__(self, *rows, **kwargs):
    self.rows = rows
    for kw, kwarg in kwargs.iteritems():
      if kw == "title":
        self.__title = kwarg
      else:
        raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))
    super(SlashRow, self).__init__()
  @property
  def title(self):
    return self.__title
  @property
  def fmt(self):
    return "{:.1f}"

  def getlatex(self):
    result = "{}".format(self.title)
    for category in categories:
      result += " & "
      parts = []
      for row in self.rows:
        total = 0
        total += row.categorydistribution[category]
        parts.append(self.fmt.format(total))
      assert len(parts) == 5
      result += "{} ({}/{}/{}/{})".format(*parts)
    return result

  def getcategorydistribution(self):
    return self.rows[0].categorydistribution

class Section(object):
  def __init__(self, title, *rows, **kwargs):
    self.c = ROOT.TCanvas()
    self.title = title
    self.rows = rows
    if not config.unblinddistributions:
      self.rows = tuple(row for row in self.rows if row.productionmode != "data")
  def getlatex(self):
    result = (r"\\"+"\n").join(_.getlatex() for _ in self.rows)
    return result
  def findrow(self, productionmode):
    def equalproductionmodes(a, row):
      try:
        ProductionMode(a)
        row.productionmode
      except (ValueError, AttributeError):
        return False
      return a == row.productionmode
    def equaltitles(a, row):
      try:
        ProductionMode(a)
      except ValueError:
        return a == row.title
      return False

    def isthisrow(a, row):
      return equalproductionmodes(a, row) or equaltitles(a, row)

    result = [row for row in self.rows if isthisrow(productionmode, row)]
    assert len(result) == 1, "{}\n{}".format(productionmode, [row.title for row in result])
    return result.pop()

def maketable():
  sections = [
    Section("SM",
      SlashRow(
        Row("VBF", "0+", title="VBF signal"),
        Row("VBF", "0-", title="VBF signal"),
        Row("VBF", "a2", title="VBF signal"),
        Row("VBF", "L1", title="VBF signal"),
        Row("VBF", "L1Zg", title="VBF signal"),
        title="VBF signal",
      ),
      SlashRow(
        Row("VH", "0+", title="VH signal"),
        Row("VH", "0-", title="VH signal"),
        Row("VH", "a2", title="VH signal"),
        Row("VH", "L1", title="VH signal"),
        Row("VH", "L1Zg", title="VH signal"),
        title=r"$\V\PH$ signal",
      ),
      SlashRow(
        Row("ggH", "0+", title="ggH signal"),
        Row("ggH", "0-", title="ggH signal"),
        Row("ggH", "a2", title="ggH signal"),
        Row("ggH", "L1", title="ggH signal"),
        Row("ggH", "L1Zg", title="ggH signal"),
        title=r"$\Pg\Pg\to\PH$ signal"
      ),
      SlashRow(
        Row("ttH", "0+", title="ttH signal"),
        Row("ttH", "0-", title="ttH signal"),
        Row("ttH", "a2", title="ttH signal"),
        Row("ttH", "L1", title="ttH signal"),
        Row("ttH", "L1Zg", title="ttH signal"),
        title=r"$\ttH$ signal"
      ),
      SlashRow(
        Row("bbH", "0+", title="bbH signal"),
        Row("bbH", "0-", title="bbH signal"),
        Row("bbH", "a2", title="bbH signal"),
        Row("bbH", "L1", title="bbH signal"),
        Row("bbH", "L1Zg", title="bbH signal"),
        title=r"$\bbH$ signal"
      ),
    ),
    Section("bkg",
      Row("ggZZ", title=r"\Pg\Pg\to4\ell bkg."),
      Row("qqZZ", title=r"$\qqbar\to4\ell$ bkg."),
      Row("ZX", title=r"$\Z\!+\!\X$ bkg."),
    ),
  ]
  sections.append(
    Section("Total",
      SlashRow(
        TotalRow(*chain((_.rows[0] for _ in sections[0].rows), sections[1].rows), title="Total expected"),
        TotalRow(*chain((_.rows[1] for _ in sections[0].rows), sections[1].rows), title="Total expected"),
        TotalRow(*chain((_.rows[2] for _ in sections[0].rows), sections[1].rows), title="Total expected"),
        TotalRow(*chain((_.rows[3] for _ in sections[0].rows), sections[1].rows), title="Total expected"),
        TotalRow(*chain((_.rows[4] for _ in sections[0].rows), sections[1].rows), title="Total expected"),
        title="Total expected",
      ),
      Row("data", title="Total observed"),
    )
  )

  result = ""
  result += r"\begin{{tabular}}{{{}}}".format("l" + "".join("c"*len(categories)) + "") + "\n"
  result += r"\hline" + "\n"
  result += " & " + " & ".join(list(r"{"+categoryname(_)+"}" for _ in categories)) + r" \\" + "\n"
  result += r"\hline" + "\n"
  joiner = r"\\"+"\n"+r"\hline"+"\n"

  result += joiner.join(section.getlatex() for section in sections) + "\n"
  result += r"\\\hline" + "\n"
  result += r"\end{tabular}" + "\n"

  print result

if __name__ == "__main__":
  maketable()
