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

categories = list(Category(_) for _ in ("Untagged", "Boosted", "VBF1jtagged", "VBFtagged", "VHLepttagged", "VHHadrtagged"))
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
    result += " & "
    amount = sum(self.categorydistribution[c] for c in list(categories))
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
        trees=gettrees(("data", "shift_pm4l")),
        weightformulas=["1" for production in config.productionsforcombine],
        scaletos=None,
      )
    elif self.productionmode.isbkg or self.productionmode == "tqH":
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
      return "{:.2f}"

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
    return "{:.2f}"

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
    return "{:.2f}"

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
    result += " & "
    parts = []
    for row in self.rows:
      total = 0
      total += sum(row.categorydistribution[c] for c in categories)
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
  signalrows = [
#    SlashRow(
      Row("ggH", "0+", title=r"$\ggH$ sig."),
#      Row("ggH", "0-", title=r"$\ggH$ signal"),
#      Row("ggH", "a2", title=r"$\ggH$ signal"),
#      Row("ggH", "L1", title=r"$\ggH$ signal"),
#      Row("ggH", "L1Zg", title=r"$\ggH$ signal"),
#      title=r"$\ggH$ signal"
#    ),
    SlashRow(
      Row("VBF", "0+", title="VBF sig."),
      Row("VBF", "0-", title="VBF sig."),
      Row("VBF", "a2", title="VBF sig."),
      Row("VBF", "L1", title="VBF sig."),
      Row("VBF", "L1Zg", title="VBF sig."),
      title="VBF sig.",
    ),
    SlashRow(
      Row("WH", "0+", title=r"$\WH$ sig."),
      Row("WH", "0-", title=r"$\WH$ sig."),
      Row("WH", "a2", title=r"$\WH$ sig."),
      Row("WH", "L1", title=r"$\WH$ sig."),
      Row("WH", "L1Zg", title=r"$\WH$ sig."),
      title=r"$\WH$ sig.",
    ),
    SlashRow(
      Row("ZH", "0+", title=r"$\ZH$ sig."),
      Row("ZH", "0-", title=r"$\ZH$ sig."),
      Row("ZH", "a2", title=r"$\ZH$ sig."),
      Row("ZH", "L1", title=r"$\ZH$ sig."),
      Row("ZH", "L1Zg", title=r"$\ZH$ sig."),
      title=r"$\ZH$ sig.",
    ),
#    SlashRow(
      Row("ttH", "0+", title=r"$\ttH$ sig."),
#      Row("ttH", "0-", title=r"$\ttH$ sig."),
#      Row("ttH", "a2", title=r"$\ttH$ sig."),
#      Row("ttH", "L1", title=r"$\ttH$ sig."),
#      Row("ttH", "L1Zg", title=r"$\ttH$ sig."),
#      title=r"$\ttH$ sig."
#    ),
#    SlashRow(
      Row("bbH", "0+", title=r"$\bbH$ sig."),
#      Row("bbH", "0-", title=r"$\bbH$ sig."),
#      Row("bbH", "a2", title=r"$\bbH$ sig."),
#      Row("bbH", "L1", title=r"$\bbH$ sig."),
#      Row("bbH", "L1Zg", title=r"$\bbH$ sig."),
#      title=r"$\bbH$ sig."
#    ),
  ]
  signalrows.append(
    SlashRow(
      TotalRow(*(_.rows[0] if isinstance(_, SlashRow) else _ for _ in signalrows), title="Signal expected"),
      TotalRow(*(_.rows[1] if isinstance(_, SlashRow) else _ for _ in signalrows), title="Signal expected"),
      TotalRow(*(_.rows[2] if isinstance(_, SlashRow) else _ for _ in signalrows), title="Signal expected"),
      TotalRow(*(_.rows[3] if isinstance(_, SlashRow) else _ for _ in signalrows), title="Signal expected"),
      TotalRow(*(_.rows[4] if isinstance(_, SlashRow) else _ for _ in signalrows), title="Signal expected"),
      title="Signal expected",
    )
  )
  sections = [
    Section("SM",
      *signalrows
    ),
    Section("bkg",
      Row("qqZZ", title=r"$\qqbar\to4\ell$ bkg."),
      Row("ggZZ", title=r"$\Pg\Pg\to4\ell$ bkg."),
      Row("EW", title="EW bkg."),
      Row("ZX", title=r"$\Z\!+\!\X$ bkg."),
    ),
  ]
  sections.append(
    Section("Total",
      SlashRow(
        TotalRow(signalrows[-1].rows[0], *sections[1].rows, title="Total expected"),
        TotalRow(signalrows[-1].rows[1], *sections[1].rows, title="Total expected"),
        TotalRow(signalrows[-1].rows[2], *sections[1].rows, title="Total expected"),
        TotalRow(signalrows[-1].rows[3], *sections[1].rows, title="Total expected"),
        TotalRow(signalrows[-1].rows[4], *sections[1].rows, title="Total expected"),
        title="Total expected",
      ),
      Row("data", title="Total observed"),
    )
  )

  result = ""
  result += r"\begin{{tabular}}{{{}}}".format("l" + "".join("c"*len(categories)) + "") + "\n"
  result += r"\hline" + "\n"
  result += " & " + " & ".join(list(r"{"+categoryname(_)+"}" for _ in categories)+["Total"]) + r" \\" + "\n"
  result += r"\hline" + "\n"
  joiner = r"\\"+"\n"+r"\hline"+"\n"

  result += joiner.join(section.getlatex() for section in sections) + "\n"
  result += r"\\\hline" + "\n"
  result += r"\end{tabular}" + "\n"

  print result

if __name__ == "__main__":
  maketable()
