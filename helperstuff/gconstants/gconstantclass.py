from array import array
from itertools import izip
import os
import re

from ..enums import EnumItem, Hypothesis, MultiEnum, MyEnum, ProductionMode
from ..utilities import cache, cd, TFile

class Process(MyEnum):
  enumname = "process"
  enumitems = (
    EnumItem("HZZ2e2mu"),
    EnumItem("VBF"),
    EnumItem("ZH"),
    EnumItem("WH"),
    EnumItem("VH"),
  )
  def __init__(self, value):
    if isinstance(value, ProductionMode):
      value = str(value)
    return super(Process, self).__init__(value)

  @property
  def validhypotheses(self):
    if self == "WH": result = "a3", "a2", "L1"
    else: result = "a3", "a2", "L1", "L1Zg"
    if self == "HZZ2e2mu": result += ("g2Zg", "g4Zg", "g2gg", "g4gg")
    return [Hypothesis(_) for _ in result]


class GConstant(MultiEnum):
  commit = "3abe12fb451e3335efa3e4922acdc88028045bc5"

  enumname = "gconstant"
  enums = Process, Hypothesis

  def check(self, *args):
    if self.hypothesis not in self.process.validhypotheses:
      raise ValueError("Bad hypothesis {} for {}".format(self.hypothesis, self.process))

  @property
  def filename(self):
    if self.hypothesis in ("g2Zg", "g4Zg", "g2gg", "g4gg"): return None
    return "{}_{}.root".format(self.process, self.hypothesis)

  @property
  def coupling(self):
    for _ in "g2", "g4", "L1", "L1Zgs", "g2Zgs", "g4Zgs", "g2gsgs", "g4gsgs":
      if self.hypothesis == _:
        return _
    assert False, self

  @property
  def url(self):
    basename = "gConstant_"+str(self.process)+"_"+self.coupling+".root"
    return os.path.join("https://github.com/usarica/HiggsWidth_PostICHEP/raw", self.commit, "Analysis/data/", basename)

  @property
  def splinename(self):
    result = "sp_tgfinal_{self.process}_SM_over_tgfinal_{self.process}_{self.coupling}".format(**locals())
    if self.process == "VH":
      result = re.sub("(tgfinal_)VH(_[0-9a-zA-Z]*)", r"\1ZH\2_plus_\1WH\2", result)
    if self.hypothesis == "L1Zg":
      result = result.replace("SM", "SM_photoncut").replace("WH_SM_photoncut", "WH_SM").replace("_plus_tgfinal_WH_L1Zgs", "")
    return result

  @property
  @cache
  def spline(self):
    with cd(os.path.dirname(__file__)), TFile(self.filename) as f:
      try:
        result = getattr(f, self.splinename)
      except AttributeError:
        f.ls()
        raise
      return result

  def getvalue(self, m4l):
    if self.hypothesis == "g2Zg":
      #if self.process == "HZZ2e2mu": return 2.24e-3 ** 0.5  #table 1 in 14-018
      if self.process == "HZZ2e2mu": return 0.005226666666666666 ** 0.5  #from Savvas
      if self.process in ("VH", "ZH", "WH"): return 0.014778325123152709 ** 0.5 #from Savvas
      if self.process == "VBF": return 0.05263157894736842
      assert False, self.process
    if self.hypothesis == "g4Zg":
      #if self.process == "HZZ2e2mu": return 2.72e-3 ** 0.5  #table 1 in 14-018
      if self.process == "HZZ2e2mu": return 0.01088 ** 0.5  #from Savvas
      if self.process in ("VH", "ZH", "WH"): return 0.017857142857142856 ** .5 #from Savvas
      if self.process == "VBF": return 0.05263157894736842
      assert False, self.process
    if self.hypothesis == "g2gg":
      #if self.process == "HZZ2e2mu": return 2.82e-3 ** 0.5  #table 1 in 14-018
      if self.process == "HZZ2e2mu": return 0.005474117647058824 ** 0.5  #from Savvas
      assert False, self.process
    if self.hypothesis == "g4gg":
      #if self.process == "HZZ2e2mu": return 2.88e-3 ** 0.5  #table 1 in 14-018
      if self.process == "HZZ2e2mu": return 0.005590588235294118 ** 0.5  #from Savvas
      assert False, self.process
    return self.spline.Eval(m4l)

def gconstants():
  for process in Process.items():
    for coupling in process.validhypotheses:
      yield GConstant(process, coupling)
