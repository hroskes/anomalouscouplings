import itertools
import json
from numbers import Number
import os
import yaml

import ROOT

from combinehelpers import Luminosity
import config
from enums import Analysis, Category, categories, Channel, channels, EnumItem, MultiEnum, MyEnum, ProductionMode
from samples import ReweightingSample, Sample
from utilities import getnesteddictvalue, MultiplyCounter, setnesteddictvalue

class YieldSystematic(MyEnum):
    enumname = "yieldsystematic"
    enumitems = (
                 EnumItem("BTag"),
                 EnumItem("JEC"),
                 EnumItem("pdf_Higgs_gg"),
                 EnumItem("pdf_Higgs_qq"),
                 EnumItem("pdf_Higgs_ttH"),
                 EnumItem("pdf_Higgs_qq"),
                 EnumItem("BRhiggs_hzz4l"),
                )

class YieldValue(MultiEnum):
    enumname = "yieldvalue"
    enums = (Analysis, Category, Channel, ProductionMode)

    yieldsfile = os.path.join(config.repositorydir, "data", "yields.json")

    @classmethod
    def getyieldsdict(cls, trycache=True):
      import globals
      if globals.yieldsdict_cache is None or not trycache:
        try:
          with open(cls.yieldsfile) as f:
            jsonstring = f.read()
        except IOError:
          try:
            os.makedirs(os.path.dirname(cls.yieldsfile))
          except OSError:
            pass
          with open(cls.yieldsfile, "w") as f:
            f.write("{}\n")
            jsonstring = "{}"
        globals.yieldsdict_cache = json.loads(jsonstring)
      return globals.yieldsdict_cache

    @classmethod
    def writeyieldsdict(cls):
      dct = cls.getyieldsdict()
      jsonstring = json.dumps(dct, sort_keys=True, indent=4, separators=(',', ': '))
      with open(cls.yieldsfile, "w") as f:
        f.write(jsonstring)

    @property
    def keys(self):
        return (
                str(self.analysis),
                str(self.category),
                str(self.channel),
                str(self.productionmode),
               )

    @property
    def value(self):
        return getnesteddictvalue(self.getyieldsdict(), *self.keys, default=None)

    @value.setter
    def value(self, value):
        setnesteddictvalue(self.getyieldsdict(), *self.keys, value=value)
        assert self.value == value

    def __float__(self):
        return self.__value

class YieldSystematicValue(MultiEnum):
    enumname = "yieldsystematicvalue"
    enums = (YieldSystematic, Analysis, Category, Channel, ProductionMode)

    yieldsystematicsfile = os.path.join(config.repositorydir, "data", "categorysystematics.json")

    @classmethod
    def getyieldsystematicsdict(cls, trycache=True):
      import globals
      if globals.yieldsystematicsdict_cache is None or not trycache:
        try:
          with open(cls.yieldsystematicsfile) as f:
            jsonstring = f.read()
        except IOError:
          try:
            os.makedirs(os.path.dirname(cls.yieldsystematicsfile))
          except OSError:
            pass
          with open(cls.yieldsystematicsfile, "w") as f:
            f.write("{}\n")
            jsonstring = "{}"
        globals.yieldsystematicsdict_cache = json.loads(jsonstring)
      return globals.yieldsystematicsdict_cache

    @classmethod
    def writeyieldsystematicsdict(cls):
      dct = cls.getyieldsystematicsdict()
      jsonstring = json.dumps(dct, sort_keys=True, indent=4, separators=(',', ': '))
      with open(cls.yieldsystematicsfile, "w") as f:
        f.write(jsonstring)

    @property
    def keys(self):
        return (
                str(self.yieldsystematic),
                str(self.analysis),
                str(self.category),
                str(self.channel),
                str(self.productionmode),
               )

    @property
    def value(self):
        return getnesteddictvalue(self.getyieldsystematicsdict(), *self.keys, default=None)

    @value.setter
    def value(self, value):
        origvalue = value
        if isinstance(value, basestring):
          try:
            value = float(value)
          except ValueError:
            if value.count("/") == 1:
              try:
                value = [float(_) for _ in value.split("/")]
              except ValueError:
                raise ValueError("string value {!r} can't be parsed as a number or 2 numbers with / in between".format(origvalue))

        if hasattr(value, "__len__"):
          if len(value) != 2:
            raise ValueError("value '{!r}' is a list or similar, but has length {} instead of 2!".format(origvalue, len(value)))
          if not all(isinstance(_, Number) for _ in value):
            raise ValueError("value '{!r}' doesn't contain only numbers!".format(origvalue))

          if list(value) == [1, 1]: value = None
          elif 1 in list(value):
            value = list(value)
            value.remove(1)
            assert len(value) == 1
            value = value[0]

        elif isinstance(value, Number) or value is None:
          pass
        else:
          raise ValueError("value '{!r}' should be None, a number, or a list (tuple, etc.) of length 2".format(self, origvalue))

        if value == 1: value = None
        setnesteddictvalue(self.getyieldsystematicsdict(), *self.keys, value=value)
        assert self.value == value

    def __str__(self):
        if self.value is None or self.value == 1:
            return "-"
        if isinstance(self.value, Number):
            return str(self.value)
        if not (hasattr(self.value, "__len__") and len(self.value) == 2):
            raise ValueError("{!r} value '{!r}' should be None, a number, or a list (tuple, etc.) of length 2".format(self, self.value))
        return "{}/{}".format(self.value[0], self.value[1])

class __TotalRate(MultiEnum):
  enums = [ProductionMode, Luminosity]
  @property
  def yamlrate(self):
    lumi = None
    result = 0
    tags = [tag.replace("Ichep16", "").replace("tagged", "Tagged") for category in categories for tag in category.names if "Ichep16" in tag]
    for channel in channels:
      filename = os.path.join(config.repositorydir, "helperstuff", "Datacards13TeV_Moriond2017", "LegoCards", "configs", "inputs", "yields_per_tag_category_13TeV_{}.yaml".format(channel))
      with open(filename) as f:
        y = yaml.load(f)
      with open(filename) as f:
        for line in f:
          if "fb-1" in line:
            if lumi is None:
              lumi = float(line.split("=")[1].split("fb-1")[0])
            if lumi != float(line.split("=")[1].split("fb-1")[0]):
              raise ValueError("Different lumis in yaml files! {} {}".format(lumi, float(line.split("=")[1].split("fb-1")[0])))
            break
        else:
          raise IOError("No luminosity in {}".format(filename))

      p = self.productionmode.yamlratename

      for tag in tags:
        try:
          result += float(y[tag][p]) * float(self.luminosity) / lumi
        except ValueError:
          result += eval(y[tag][p].replace("@0", "125")) * float(self.luminosity) / lumi

    return result

  @property
  def treerate(self):
    assert self.productionmode == "VBF bkg"
    result = 0
    c = ROOT.TCanvas()
    t = ROOT.TChain("candTree")
    t.SetBranchStatus("*", 0)
    t.SetBranchStatus("MC_weight_*", 1)
    t.SetBranchStatus("ZZMass", 1)
    weightname = set()
    for flavor in "2e2mu", "4e", "4mu":
      t.Add(Sample(self.productionmode, flavor, self.production).withdiscriminantsfile())
      weightname.add(Sample(self.productionmode, flavor, self.production).weightname())
    assert len(weightname) == 1, weightname
    weightname = weightname.pop()
    t.Draw("1", "{}*(ZZMass>{} && ZZMass<{})".format(weightname, config.m4lmin, config.m4lmax))
    return c.FindObject("htemp").Integral()

  @property
  def rate(self):
    if self.productionmode == "VBF bkg": return self.treerate
    else: return self.yamlrate

def totalrate(*args):
  return __TotalRate(*args).rate

def count(fromsamples, tosamples, categorizations):
    t = ROOT.TChain("candTree")
    for fromsample in fromsamples:
        t.Add(fromsample.withdiscriminantsfile())
    length = t.GetEntries()
    result = MultiplyCounter()
    t.SetBranchStatus("*", 0)
    t.SetBranchStatus("category_*", 1)
    t.SetBranchStatus("MC_weight_*", 1)
    t.SetBranchStatus("Z*Flav", 1)
    t.SetBranchStatus("ZZMass", 1)

    c = ROOT.TCanvas()
    for tosample, categorization in itertools.product(tosamples, categorizations):
        if tosample == ReweightingSample("WH", "L1Zg"): continue
        t.Draw(categorization.category_function_name+":abs(Z1Flav*Z2Flav)", "{}*(ZZMass>{} && ZZMass<{})".format(tosample.weightname(), config.m4lmin, config.m4lmax), "LEGO")
        h = c.FindObject("htemp")
        for i in range(6):
            for channel in channels:
                toadd = h.GetBinContent(h.FindBin(channel.ZZFlav, i))
                result[tosample, categorization, Category.fromid(i), channel] += toadd
                result[tosample, categorization, Category.fromid(i)] += toadd

        if tosample not in result:
            result[tosample] = h.Integral()

    return result
