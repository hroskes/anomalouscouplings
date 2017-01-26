#!/usr/bin/env python

from itertools import product
import random

import ROOT

from helperstuff import config

from helperstuff.combinehelpers import getrate
from helperstuff.enums import analyses, Analysis, Category, Channel, channels, flavors, HffHypothesis, Hypothesis, JECSystematic, MultiEnum, ProductionMode
from helperstuff.samples import ReweightingSample, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import KeyDefaultDict, MultiplyCounter, tfiles

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

categories = list(Category(_) for _ in ("VBFtagged", "VHHadrtagged", "Untagged"))

def categoryname(category):
  if category == "VBFtagged": return "VBF jets"
  if category == "VHHadrtagged": return "\V\PH jets"
  if category == "Untagged": return "Inclusive"

def gettree(productionmode):
  t = ROOT.TChain("candTree")
  for sample in productionmode.allsamples(production):
    t.Add(sample.withdiscriminantsfile())
  return t

trees = KeyDefaultDict(gettree)

class Row(MultiEnum):
  enums = (ProductionMode, Hypothesis, Analysis, HffHypothesis)
  enumname = "Row"
  def __init__(self, *args, **kwargs):
    self.title = None
    if "title" in kwargs:
      self.title = kwargs["title"]
      del kwargs["title"]
    super(Row, self).__init__(*args, **kwargs)
    if self.title is None:
      self.title = str(self.productionmode)
  def check(self, *args):
    dontcheck = []
    if self.productionmode.isbkg:
      dontcheck.append(Hypothesis)
    if self.productionmode not in ("ttH", "HJJ"):
      dontcheck.append(HffHypothesis)
    super(Row, self).check(*args, dontcheck=dontcheck)

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

  def getlatex(self):
    categorydistribution = self.getcategorydistribution()
    result = " & {}".format(self.title)
    for category in categories:
      result += " & "
      total = 0
      for channel in "2e2mu", "4e", "4mu":
        channel = Channel(channel)
        total += categorydistribution[channel, category]
        result += "{:.2f}/".format(categorydistribution[channel, category])
      result += "{:.2f}".format(total)
    result += r"\\"
    return result

class RowChannel(Row, MultiEnum):
  enums = (Row, Channel)
  def getcategorydistribution(self):
    if self.productionmode.isbkg or (self.productionmode != "ttH" and self.hypothesis == "SM"):
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
      for i in range(6):
        if i+1 <= h.GetNbinsX():
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
  def getlatex(self):
    result = r"\multirow{{{}}}{{*}}{{{}}}".format(len(self.rows), self.title)
    result += "\n".join(_.getlatex() for _ in self.rows)
    return result

def maketable(analysis):
  analysis = Analysis(analysis)
  sections = [
    Section("SM",
      Row(analysis, "VBF", analysis.purehypotheses[0]),
      Row(analysis, "ZH", analysis.purehypotheses[0], title=r"\Z\PH"),
      Row(analysis, "WH", analysis.purehypotheses[0], title=r"\W\PH"),
      Row(analysis, "ggH", analysis.purehypotheses[0], title=r"\Pg\Pg\to\H"),
      Row(analysis, "ttH", analysis.purehypotheses[0], "Hff0+", title=r"\ttbar\PH"),
    ),
    Section("{}=1".format(analysis.title),
      Row(analysis, "VBF", analysis.purehypotheses[1]),
      Row(analysis, "ZH", analysis.purehypotheses[1], title=r"\Z\PH"),
      Row(analysis, "WH", analysis.purehypotheses[1], title=r"\W\PH"),
      Row(analysis, "ggH", analysis.purehypotheses[1], title=r"\Pg\Pg\to\H"),
      Row(analysis, "ttH", analysis.purehypotheses[1], "Hff0+", title=r"\ttbar\PH"),
    ),
    Section("bkg",
      Row(analysis, "qqZZ", title=r"\qqbar\to\Z\Z"),
      Row(analysis, "ggZZ", title=r"\Pg\Pg\to\Z\Z"),
      Row(analysis, "VBFbkg", title=r"VBF/\V\V\V bkg"),
      Row(analysis, "ZX", title=r"Z+X"),
    )
  ]
  print r"\begin{{tabular}}{{{}}}".format("|" + "|".join("c"*(len(categories)+2)) + "|")
  print " & & " + " & ".join(categoryname(_) for _ in categories)
  for section in sections:
    print section.getlatex()
  print r"\end{tabular}"

if __name__ == "__main__":
  maketable("fa3")
