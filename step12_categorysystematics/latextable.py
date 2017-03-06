#!/usr/bin/env python

from abc import ABCMeta, abstractmethod, abstractproperty
from itertools import product
import random
import sys

import ROOT

from helperstuff import config

from helperstuff.combinehelpers import getdatatree, getrate
from helperstuff.enums import analyses, Analysis, Category, Channel, flavors, HffHypothesis, Hypothesis, JECSystematic, MultiEnum, MultiEnumABCMeta, ProductionMode
from helperstuff.samples import ReweightingSample, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import KeyDefaultDict, MultiplyCounter, tfiles

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

categories = list(Category(_) for _ in ("VBFtagged", "VHHadrtagged", "Untagged"))
channels = list(Channel(_) for _ in ("4e", "4mu", "2e2mu"))

def categoryname(category):
  if category == "VBFtagged": return "VBF-jets"
  if category == "VHHadrtagged": return "$\V\PH$-jets"
  if category == "Untagged": return "untagged"

def gettree(productionmode):
  t = ROOT.TChain("candTree")
  for sample in productionmode.allsamples(production):
    t.Add(sample.withdiscriminantsfile())
  return t

trees = KeyDefaultDict(gettree)

class RowBase(object):
  __metaclass__ = ABCMeta

  def __init__(self, *args, **kwargs):
    self.__categorydistribution = None
    super(RowBase, self).__init__(*args, **kwargs)

  @property
  def categorydistribution(self):
    if self.__categorydistribution is None:
      self.__categorydistribution = self.getcategorydistribution()
    return self.__categorydistribution

  def scale(self, scaleby):
    self.categorydistribution
    self.__categorydistribution *= scaleby

  def getlatex(self, dochannels=True):
    result = " & {}".format(self.title)
    for category in categories:
      result += " & "
      total = 0
      for channel in channels:
        channel = Channel(channel)
        total += self.categorydistribution[channel, category]
        if dochannels:
          result += self.fmt.format(self.categorydistribution[channel, category])+"/"
      result += self.fmt.format(total)
    return result

  @abstractproperty
  def title(self): pass
  @abstractproperty
  def fmt(self): pass
  @abstractmethod
  def getcategorydistribution(self): pass

class Row(RowBase, MultiEnum):
  __metaclass__ = MultiEnumABCMeta
  enums = (ProductionMode, Hypothesis, Analysis, HffHypothesis)
  enumname = "Row"
  def __init__(self, *args, **kwargs):
    self.__title = None
    if "title" in kwargs:
      self.__title = kwargs["title"]
      del kwargs["title"]
    super(Row, self).__init__(*args, **kwargs)
    if self.__title is None:
      self.__title = str(self.productionmode)
  def check(self, *args):
    dontcheck = []
    if self.productionmode == "data" or self.productionmode.isbkg:
      dontcheck.append(Hypothesis)
    if self.productionmode not in ("ttH", "HJJ"):
      dontcheck.append(HffHypothesis)
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

  @property
  def tree(self): return trees[self.productionmode]
  @property
  def reweightingsample(self): return ReweightingSample(self.productionmode, self.hypothesis, self.hffhypothesis)

  @property
  def scalefactor(self):
    result = 1
    if self.productionmode in ("VBF", "ZH", "WH"):
      result *= (
                  (ReweightingSample(self.productionmode, self.hypothesis).xsec 
                      / sum(ReweightingSample(_, self.hypothesis).xsec for _ in ("VBF", "ZH", "WH")))
                 /
                  (ReweightingSample(self.productionmode, "SM"           ).xsec 
                      / sum(ReweightingSample(_, "SM"           ).xsec for _ in ("VBF", "ZH", "WH")))
                )
    if self.productionmode in ("VBF", "ZH", "WH", "ggH", "ttH"):
      result /= sum(
                    Sample.effectiveentries(
                                            reweightfrom=reweightfrom,
                                            reweightto=self.reweightingsample,
                                           )
                     for reweightfrom in self.productionmode.allsamples(production)
                   )
    if self.productionmode != "data":
      result *= config.expectedscanluminosity
    return result

  @property
  def rowchannels(self):
    return [RowChannel(self, channel) for channel in channels]

  def getcategorydistribution(self):
    result = MultiplyCounter()
    for rowchannel in self.rowchannels:
      result += rowchannel.getcategorydistribution()
    return result

  @property
  def weightname(self):
    if self.productionmode in ("ggZZ", "VBF bkg"): return ReweightingSample(self.productionmode, "2e2mu").weightname()
    if self.productionmode.isbkg: return ReweightingSample(self.productionmode).weightname()
    return self.reweightingsample.weightname()

class TotalRow(RowBase):
  def __init__(self, *rows):
    self.rows = rows
    super(TotalRow, self).__init__()
  @property
  def title(self):
    return "Expected"
  def getcategorydistribution(self):
    return sum((row.categorydistribution for row in self.rows), MultiplyCounter())
  @property
  def fmt(self):
    return "{:.1f}"

class RowChannel(Row, MultiEnum):
  enums = (Row, Channel)
  def getcategorydistribution(self):
    if self.productionmode == "data":
      return MultiplyCounter({
                              (self.channel, category): getdatatree(self.analysis, category, self.channel, production).GetEntries()
                                for category in categories
                            })
    if self.productionmode.isbkg or self.hypothesis == "SM":
      result = MultiplyCounter({(self.channel, category): getrate(self.channel, category, production, "forexpectedscan", self.analysis, self.productionmode) for category in categories})
    else:
      t = self.tree
      weightparts = [
                     "ZZMass>{}".format(config.m4lmin),
                     "ZZMass<{}".format(config.m4lmax),
                     "Z1Flav*Z2Flav=={}".format(self.ZZFlav),
                     self.weightname,
                    ]

      wt = " * ".join("("+_+")" for _ in weightparts)
      hname = "h{}".format(random.getrandbits(100))
      success = t.Draw("category_{}>>{}".format(self.analysis.categoryname, hname), wt)
      if success == 0: return MultiplyCounter()
      h = getattr(ROOT, hname)
      result = MultiplyCounter()
      for i in range(h.GetNbinsX()):
        result[self.channel, Category.fromid(i)] += h.GetBinContent(i+1)
      result *= self.scalefactor
    return result

  @property
  def ZZFlav(self):
    result = self.channel.ZZFlav
    if self.productionmode == "ZX": result *= -1
    return result

class Section(object):
  def __init__(self, title, *rows):
    self.c = ROOT.TCanvas()
    self.title = title
    self.rows = rows
    if not config.unblinddistributions:
      self.rows = tuple(row for row in self.rows if row.productionmode != "data")
  def getlatex(self, dochannels=True):
    result = r"\multirow{{{}}}{{*}}{{{}}}".format(len(self.rows), self.title)
    result += (r"\\\cline{{2-{}}}".format(len(categories)+2)+"\n").join(_.getlatex(dochannels=dochannels) for _ in self.rows)
    return result
  def findrow(self, productionmode):
    result = [row for row in self.rows if row.productionmode == productionmode]
    assert len(result) == 1
    return result.pop()

def scalerows(scalethese, tothese):
  wanttotal = sum(row.categorydistribution[channel, category] for row in tothese for channel in channels for category in categories)
  previoustotal = sum(row.categorydistribution[channel, category] for row in scalethese for channel in channels for category in categories)
  for row in scalethese:
    row.scale(wanttotal / previoustotal)

def total(rows):
  return sum(row.categorydistribution[channel, category] for row in rows for channel in channels for category in categories)

def maketable(analysis, dochannels=True):
  analysis = Analysis(analysis)
  sections = [
    Section("SM",
      Row(analysis, "VBF", analysis.purehypotheses[0]),
      Row(analysis, "ZH", analysis.purehypotheses[0], title=r"$\Z\PH$"),
      Row(analysis, "WH", analysis.purehypotheses[0], title=r"$\PW\PH$"),
      Row(analysis, "ggH", analysis.purehypotheses[0], title=r"$\Pg\Pg\to\PH$"),
      Row(analysis, "ttH", analysis.purehypotheses[0], "Hff0+", title=r"$\ttbar\PH$"),
    ),
    Section("${}=1$".format(analysis.title(latex=True)),
      Row(analysis, "VBF", analysis.purehypotheses[1]),
      Row(analysis, "ZH", analysis.purehypotheses[1], title=r"$\Z\PH$"),
      Row(analysis, "WH", analysis.purehypotheses[1], title=r"$\PW\PH$"),
      Row(analysis, "ggH", analysis.purehypotheses[1], title=r"$\Pg\Pg\to\PH$"),
      Row(analysis, "ttH", analysis.purehypotheses[1], "Hff0+", title=r"$\ttbar\PH$"),
    ),
    Section("bkg",
      Row(analysis, "qqZZ", title=r"$\qqbar\to4\ell$"),
      Row(analysis, "ggZZ", title=r"$\Pg\Pg\to4\ell$"),
      Row(analysis, "VBFbkg", title=r"VBF/$\V\V\V$"),
      Row(analysis, "ZX", title=r"$\Z+\X$"),
    )
  ]
  sections.append(
    Section("Total",
      TotalRow(*(sections[0].rows+sections[2].rows)),
      Row(analysis, "data", title="Observed"),
    )
  )

  scalerows([sections[1].findrow(p) for p in ("VBF", "ZH", "WH")], [sections[0].findrow(p) for p in ("VBF", "ZH", "WH")])
  scalerows([sections[1].findrow(p) for p in ("ggH", "ttH")], [sections[0].findrow(p) for p in ("ggH", "ttH")])

  print r"\begin{{tabular}}{{{}}}".format("|" + "|".join("c"*(len(categories)+2)) + "|")
  print r"\hline"
  print " & & " + " & ".join(categoryname(_) for _ in categories) + r"\\\hline\hline"
  print (r"\\\hline"+"\n").join(section.getlatex(dochannels=dochannels) for section in sections)
  print r"\\\hline"
  print r"\end{tabular}"

if __name__ == "__main__":
  maketable(sys.argv[1], dochannels=True)
