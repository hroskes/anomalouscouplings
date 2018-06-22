#!/usr/bin/env python

from abc import ABCMeta, abstractmethod, abstractproperty
from itertools import product
import random
import sys

import ROOT

from helperstuff import config

from helperstuff.combinehelpers import getdatatree, getrate
from helperstuff.enums import analyses, Analysis, Category, Channel, EnumItem, flavors, HffHypothesis, Hypothesis, JECSystematic, MultiEnum, MultiEnumABCMeta, MyEnum, ProductionMode
from helperstuff.samples import ReweightingSample, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import KeyDefaultDict, MultiplyCounter, tfiles
from helperstuff.yields import YieldSystematicValue

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

categories = list(Category(_) for _ in ("VBFtagged", "VHHadrtagged", "Untagged"))
channels = list(Channel(_) for _ in ("4e", "4mu", "2e2mu"))

def categoryname(category):
  if category == "VBFtagged": return "VBF-jets"
  if category == "VHHadrtagged": return "$\V\PH$-jets"
  if category == "Untagged": return "untagged"

class Column(MyEnum):
  enumname = "column"
  enumitems = (
    EnumItem("QCD"),
    EnumItem("PDF"),
    EnumItem("EWcorr_VV"),
    EnumItem("PythiaScale"),
    EnumItem("PythiaTune"),
    EnumItem("CMS_scale_j_13TeV_2016"),
    EnumItem("CMS_btag_comb_13TeV_2016"),
  )

  @property
  def title(self):
    if self == "QCD": return "QCD scale"
    if self == "PDF": return "PDF set"
    if self == "EWcorr_VV": return "EW corr."
    return str(self)

columns = Column.items()

class Row(MultiEnum):
  enums = (ProductionMode, Analysis, Category)
  enumname = "row"

  def getlatex(self):
    result = " & {}".format(categoryname(self.category))
    for column in columns:
      result += " & {}".format(Cell(self, column).getlatex())
    return result

class Cell(MultiEnum):
  enums = (Row, Column)
  @property
  def yieldsystematic(self):
    if self.column == "QCD": return self.productionmode.QCDsystematicname
    if self.column == "PDF": return self.productionmode.pdfsystematicname
    return str(self.column)
  def getlatex(self):
    result = {YieldSystematicValue(self.row, self.yieldsystematic, channel).latexstr for channel in channels}
    assert len(result) == 1, (self, result)
    return result.pop()

class Section(MultiEnum):
  enums = (ProductionMode, Analysis)
  def __init__(self, *args, **kwargs):
    self.title = None
    if "title" in kwargs:
      self.title = kwargs["title"]
      del kwargs["title"]
    super(Section, self).__init__(*args, **kwargs)
    if self.title is None:
      self.title = str(self.productionmode)

  @property
  def rows(self): return [Row(self, category) for category in categories]

  def getlatex(self):
    result = r"\multirow{{{}}}{{*}}{{{}}}".format(len(self.rows), self.title)
    result += (r"\\\cline{{2-{}}}".format(len(columns)+2)+"\n").join(_.getlatex() for _ in self.rows)
    return result

def maketable(analysis, dochannels=True):
  sections = [
    Section("ggH", analysis, title=r"$\Pg\Pg\to\PH$"),
    Section("VBF", analysis),
    Section("ZH", analysis, title=r"$\Z\PH$"),
    Section("WH", analysis, title=r"$\PW\PH$"),
    Section("ttH", analysis, title=r"$\ttbar\PH$"),
    Section("qqZZ", analysis, title=r"$\qqbar\to4\ell$"),
  ]

  print r"\begin{{tabular}}{{{}}}".format("|" + "|".join("c"*(len(columns)+2)) + "|")
  print r"\hline"
  print " & & " + " & ".join(_.title for _ in columns) + r"\\\hline\hline"
  print (r"\\\hline"+"\n").join(section.getlatex() for section in sections)
  print r"\\\hline"
  print r"\end{tabular}"

if __name__ == "__main__":
  maketable(sys.argv[1])
