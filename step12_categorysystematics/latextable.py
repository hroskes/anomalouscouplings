#!/usr/bin/env python

if __name__ == "__main__":
  import argparse
  p = argparse.ArgumentParser()
  g = p.add_mutually_exclusive_group(required=True)
  g.add_argument("--HIG17011", action="store_true")
  g.add_argument("--HIG17011PAS", action="store_true")
  g.add_argument("--HIG18002", action="store_true")
  p.add_argument("analysis", choices="fa3 fa2 fL1 fL1Zg".split())
  args = p.parse_args()

from abc import ABCMeta, abstractmethod, abstractproperty
from itertools import product
import random
import sys

import ROOT
import uncertainties

from helperstuff import config

from helperstuff.combinehelpers import getdatatree, getrate, getrate2015, Luminosity
from helperstuff.enums import analyses, Analysis, Category, Channel, EnumItem, flavors, HffHypothesis, Hypothesis, JECSystematic, MultiEnum, MultiEnumABCMeta, MyEnum, ProductionMode
from helperstuff.samples import ReweightingSample, ReweightingSampleWithPdf, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import KeyDefaultDict, MultiplyCounter, tfiles

class TableType(MyEnum):
  enumname = "tabletype"
  enumitems = (
    EnumItem("HIG17011PAS"),
    EnumItem("HIG17011"),
    EnumItem("HIG18002"),
  )

categories = list(Category(_) for _ in ("VBFtagged", "VHHadrtagged", "Untagged"))
channels = list(Channel(_) for _ in ("4e", "4mu", "2e2mu"))

def categoryname(category, tabletype):
  if tabletype in ("HIG17011PAS", "HIG17011"):
    if category == "VBFtagged": return "VBF-jets"
    if category == "VHHadrtagged": return "$\V\!\PH$-jets"
    if category == "Untagged": return "untagged"
  if tabletype == "HIG18002":
    if category == "VBFtagged": return "VBF-tagged"
    if category == "VHHadrtagged": return "$\V\PH$-tagged"
    if category == "Untagged": return "Untagged"

def gettree(productionmodeandproduction):
  productionmode, production = productionmodeandproduction
  t = ROOT.TChain("candTree")
  for sample in productionmode.allsamples(production):
    t.Add(sample.withdiscriminantsfile())
  return t

trees = KeyDefaultDict(gettree)

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

  def scale(self, scaleby):
    self.categorydistribution
    for channel in channels:
      for category in categories:
        assert (channel, category) in self.__categorydistribution
        self.__categorydistribution[channel, category] *= scaleby

  @abstractproperty
  def title(self): pass
  @abstractproperty
  def fmt(self): pass
  @abstractmethod
  def getlatex(self, tabletype): pass
  @abstractmethod
  def getcategorydistribution(self): pass

class RowBase(RowBaseBase):
  def getlatex(self, tabletype):
    result = ""
    if tabletype == "HIG17011PAS":
      result += " & "
    result = "{}".format(self.title)
    for category in list(categories)+[2015]:
      if tabletype == "HIG18002" and isinstance(category, int) and category == 2015: continue
      result += " & "
      total = 0
      for channel in channels:
        channel = Channel(channel)
        amount = self.categorydistribution[channel, category]
#        if tabletype == "HIG18002" and category == "Untagged": amount += self.categorydistribution[channel, 2015]
        total += amount
        if tabletype == "HIG17011PAS":
          result += self.fmt.format(amount)+"/"
      result += self.fmt.format(total)
      if tabletype == "HIG18002":
        result += " & xx"
    return result

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

  def tree(self, production): return trees[self.productionmode, production]
  @property
  def reweightingsample(self): return ReweightingSample(self.productionmode, self.hypothesis, self.hffhypothesis)

  def scalefactor(self, production):
    result = 1
    if self.productionmode in ("VBF", "ZH", "WH"):
      result *= (
                  (ReweightingSampleWithPdf(self.productionmode, self.hypothesis, production).xsec 
                      / sum(ReweightingSampleWithPdf(_, self.hypothesis, production).xsec for _ in ("VBF", "ZH", "WH")))
                 /
                  (ReweightingSample(self.productionmode, "SM"           ).xsec 
                      / sum(ReweightingSample(_, "SM"           ).xsec for _ in ("VBF", "ZH", "WH")))
                )
    if self.productionmode in ("VBF", "ZH", "WH", "ggH", "ttH", "bbH"):
      result /= sum(
                    Sample.effectiveentries(
                                            reweightfrom=reweightfrom,
                                            reweightto=self.reweightingsample,
                                           )
                     for reweightfrom in self.productionmode.allsamples(production)
                   )
    if self.productionmode != "data":
      result *= float(Luminosity(production, "fordata"))
    result = uncertainties.nominal_value(result)
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
    assert len(rows) == 2
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

  def getlatex(self, tabletype):
    assert tabletype != "HIG17011PAS"

    result = "{}".format(self.title)
    for category in categories+[2015]:
      if tabletype == "HIG18002" and isinstance(category, int) and category == 2015: continue
      result += " & "
      parts = []
      for row in self.rows:
        total = 0
        for channel in channels:
          channel = Channel(channel)
          total += row.categorydistribution[channel, category]
#          if tabletype == "HIG18002" and category == "Untagged": total += row.categorydistribution[channel, 2015]
        parts.append(self.fmt.format(total))
      assert len(parts) == 2
      result += "{} ({})".format(*parts)
      if tabletype == "HIG18002":
        result += " & xx (xx)"
    return result

  def getcategorydistribution(self):
    return self.rows[0].categorydistribution

class RowChannel(Row, MultiEnum):
  enums = (Row, Channel)
  def getcategorydistributionproduction(self, production):
    if self.productionmode == "data":
      result = MultiplyCounter({
                                (self.channel, category): getdatatree(self.analysis, category, self.channel, production).GetEntries()
                                  for category in categories
                              })
    elif self.productionmode.isbkg or self.hypothesis == "SM":
      result = MultiplyCounter({(self.channel, category): getrate(self.channel, category, production, "fordata", self.analysis, self.productionmode) for category in categories})
    else:
      t = self.tree(production)
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
      result *= self.scalefactor(production)
    return result

  def getcategorydistribution(self):
    result = sum((self.getcategorydistributionproduction(production) for production in config.productionsforcombine), MultiplyCounter())
    result[self.channel, 2015] = getrate2015(self.channel, self.productionmode)
    return result

  @property
  def ZZFlav(self):
    result = self.channel.ZZFlav
    if self.productionmode == "ZX": result *= -1
    return result

class Section(object):
  def __init__(self, title, *rows, **kwargs):
    self.c = ROOT.TCanvas()
    self.title = title
    self.rows = rows
    if not config.unblinddistributions:
      self.rows = tuple(row for row in self.rows if row.productionmode != "data")
  def getlatex(self, tabletype):
    if tabletype == "HIG17011" or tabletype == "HIG18002":
      result = (r"\\"+"\n").join(_.getlatex(tabletype=tabletype) for _ in self.rows)
    elif tabletype == "HIG17011PAS":
      result = r"\multirow{{{}}}{{*}}{{{}}}".format(len(self.rows), self.title)
      result += (r"\\\cline{{2-{}}}".format(len(categories)+2)+"\n").join(_.getlatex(tabletype=tabletype) for _ in self.rows)
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

def scalerows(scalethese, tothese):
  scalethese = list(scalethese)
  tothese = list(tothese)
  wanttotal = sum(row.categorydistribution[channel, category] for row in tothese for channel in channels for category in categories)
  previoustotal = sum(row.categorydistribution[channel, category] for row in scalethese for channel in channels for category in categories)
  for row in scalethese:
    row.scale(wanttotal / previoustotal)

  sum2015 = sum(row.categorydistribution[channel, 2015] for row in scalethese for channel in channels)

  if any(isinstance(_[1], int) and _[1] == 2015 for _ in sum((row.categorydistribution.keys() for row in scalethese+tothese), [])):
    for row in scalethese:
      for channel in channels:
        row.categorydistribution[channel, 2015] = sum2015 / wanttotal * sum(row.categorydistribution[channel, category] for category in categories)

def scaleslashrows(rows):
  for row in rows: assert len(row.rows) == 2
  scalerows([row.rows[1] for row in rows], [row.rows[0] for row in rows])

def total(rows):
  return sum(row.categorydistribution[channel, category] for row in rows for channel in channels for category in categories)

def maketable(analysis, tabletype):
  analysis = Analysis(analysis)
  if tabletype == "HIG17011":
    sections = [
      Section("SM",
        SlashRow(
          Row(analysis, "VBF", analysis.purehypotheses[0], title="VBF signal"),
          Row(analysis, "VBF", analysis.purehypotheses[1], title="VBF signal"),
          title="VBF signal",
        ),
        SlashRow(
          Row(analysis, "ZH", analysis.purehypotheses[0], title=r"$\Z\PH$"),
          Row(analysis, "ZH", analysis.purehypotheses[1], title=r"$\Z\PH$"),
          title=r"$\Z\!\PH$ signal",
        ),
        SlashRow(
          Row(analysis, "WH", analysis.purehypotheses[0], title=r"$\PW\PH$"),
          Row(analysis, "WH", analysis.purehypotheses[1], title=r"$\PW\PH$"),
          title=r"$\PW\!\PH$ signal",
        ),
        SlashRow(
          Row(analysis, "ggH", analysis.purehypotheses[0], title=r"$\Pg\Pg\to\PH$"),
          Row(analysis, "ggH", analysis.purehypotheses[1], title=r"$\Pg\Pg\to\PH$"),
          title=r"$\Pg\Pg\to\PH$ signal"
        ),
        SlashRow(
          Row(analysis, "ttH", analysis.purehypotheses[0], "Hff0+", title=r"$\ttbar\PH$"),
          Row(analysis, "ttH", analysis.purehypotheses[1], "Hff0+", title=r"$\ttbar\PH$"),
          title=r"$\ttbar\PH$ signal"
        ),
      ),
      Section("bkg",
        Row(analysis, "qqZZ", title=r"$\qqbar\to4\ell$ bkg"),
        Row(analysis, "ggZZ", title=r"$\Pg\Pg\to4\ell$ bkg"),
        Row(analysis, "VBFbkg", title=r"VBF/$\V\V\V$ bkg"),
        Row(analysis, "ZX", title=r"$\Z\!+\!\X$ bkg"),
      ),
    ]
    sections.append(
      Section("Total",
        TotalRow(*(sections[0].rows+sections[1].rows), title="Total expected"),
        Row(analysis, "data", title="Total observed"),
      )
    )
    scaleslashrows((sections[0].findrow("VBF signal"), sections[0].findrow(r"$\Z\!\PH$ signal"), sections[0].findrow(r"$\PW\!\PH$ signal")))
    scaleslashrows([sections[0].findrow(r"$\Pg\Pg\to\PH$ signal"), sections[0].findrow(r"$\ttbar\PH$ signal")])
  elif tabletype == "HIG18002":
    sections = [
      Section("SM",
        SlashRow(
          Row(analysis, "VBF", analysis.purehypotheses[0], title="VBF signal"),
          Row(analysis, "VBF", analysis.purehypotheses[1], title="VBF signal"),
          title="VBF signal",
        ),
        SlashRow(
          Row(analysis, "ZH", analysis.purehypotheses[0], title=r"$\Z\PH$"),
          Row(analysis, "ZH", analysis.purehypotheses[1], title=r"$\Z\PH$"),
          title=r"$\Z\PH$ signal",
        ),
        SlashRow(
          Row(analysis, "WH", analysis.purehypotheses[0], title=r"$\PW\PH$"),
          Row(analysis, "WH", analysis.purehypotheses[1], title=r"$\PW\PH$"),
          title=r"$\PW\PH$ signal",
        ),
        SlashRow(
          Row(analysis, "ggH", analysis.purehypotheses[0], title=r"$\Pg\Pg\to\PH$"),
          Row(analysis, "ggH", analysis.purehypotheses[1], title=r"$\Pg\Pg\to\PH$"),
          title=r"$\Pg\Pg\to\PH$ signal"
        ),
        SlashRow(
          Row(analysis, "ttH", analysis.purehypotheses[0], "Hff0+", title=r"$\ttH$"),
          Row(analysis, "ttH", analysis.purehypotheses[1], "Hff0+", title=r"$\ttH$"),
          title=r"$\ttH$ signal"
        ),
        SlashRow(
          Row(analysis, "bbH", analysis.purehypotheses[0], title=r"$\bbH$"),
          Row(analysis, "bbH", analysis.purehypotheses[1], title=r"$\bbH$"),
          title=r"$\bbH$ signal"
        ),
      ),
      Section("bkg",
        Row(analysis, "qqZZ", title=r"$\qqbar\to4\ell$ bkg"),
        Row(analysis, "ggZZ", title=r"$\Pg\Pg\to4\ell$ bkg"),
        Row(analysis, "VBFbkg", title=r"VBF/$\V\V\V$ bkg"),
        Row(analysis, "ZX", title=r"$\Z\!+\!\X$ bkg"),
      ),
    ]
    sections.append(
      Section("Total",
        TotalRow(*(sections[0].rows+sections[1].rows), title="Total expected"),
        Row(analysis, "data", title="Total observed"),
      )
    )
    scaleslashrows((sections[0].findrow("VBF signal"), sections[0].findrow(r"$\Z\PH$ signal"), sections[0].findrow(r"$\PW\PH$ signal")))
    scaleslashrows([sections[0].findrow(r"$\Pg\Pg\to\PH$ signal"), sections[0].findrow(r"$\ttH$ signal"), sections[0].findrow(r"$\bbH$ signal")])
  elif tabletype == "HIG17011PRL":
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

  if tabletype == "HIG17011PAS":
    scalerows([sections[1].findrow(p) for p in ("VBF", "ZH", "WH")], [sections[0].findrow(p) for p in ("VBF", "ZH", "WH")])
    scalerows([sections[1].findrow(p) for p in ("ggH", "ttH")], [sections[0].findrow(p) for p in ("ggH", "ttH")])

  if tabletype == "HIG17011":
    print r"\begin{{tabular}}{{{}}}".format("l" + "".join("c"*(len(categories)+1)) + "")
    print r"\hline\hline"
    print " & " + " & ".join(list(categoryname(_, tabletype) for _ in categories)+["~~2015~~"]) + r"\\\hline"
    joiner = r"\\"+"\n"
  elif tabletype == "HIG17011PAS":
    print r"\begin{{tabular}}{{{}}}".format("|" + "|".join("c"*(len(categories)+2)) + "|")
    print r"\hline"
    print " & & " + " & ".join(categoryname(_, tabletype) for _ in categories) + r"\\\hline\hline"
    joiner = r"\\\hline"+"\n"
  elif tabletype == "HIG18002":
    print r"\begin{{tabular}}{{{}}}".format("l" + "".join("c"*(2*len(categories))) + "")
    print r"\hline"
    print " & " + " & ".join(list(r"\multicolumn{2}{c}{"+categoryname(_, tabletype)+"}" for _ in categories)) + r" \\"
    print " & " + " & ".join(list(r"\onshell & \offshell" for _ in categories)) + r"  \\\hline"
    joiner = r"\\\hline"+"\n"

  print joiner.join(section.getlatex(tabletype=tabletype) for section in sections)

  if tabletype == "HIG17011":
    print r"\\\hline\hline"
  else:
    print r"\\\hline"
  print r"\end{tabular}"

if __name__ == "__main__":
  tabletype = None
  for _ in TableType.enumitems:
    if any(getattr(args, name) for name in _.names):
      assert tabletype is None
      tabletype = TableType(_)
  maketable(args.analysis, tabletype)
