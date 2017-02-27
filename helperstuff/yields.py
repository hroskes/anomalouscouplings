import itertools
import json
from math import sqrt
from numbers import Number
import os
import yaml

import ROOT

from combinehelpers import Luminosity
import config
from enums import Analysis, Category, categories, Channel, channels, EnumItem, MultiEnum, MultiEnumABCMeta, MyEnum, ProductionMode
from samples import ReweightingSample, Sample
from utilities import JsonDict, MultiplyCounter

class YieldSystematic(MyEnum):
    enumname = "yieldsystematic"
    enumitems = (
                 EnumItem("BTag"),
                 EnumItem("JEC"),
                 EnumItem("QCDscale_ggH_yield"),
                 EnumItem("QCDscale_ggH_cat"),
                 EnumItem("QCDscale_qqH_yield"),
                 EnumItem("QCDscale_qqH_cat"),
                 EnumItem("QCDscale_VH_yield"),
                 EnumItem("QCDscale_VH_cat"),
                 EnumItem("QCDscale_ttH_yield"),
                 EnumItem("QCDscale_ttH_cat"),
                 EnumItem("QCDscale_VV_yield"),
                 EnumItem("QCDscale_VV_cat"),
                 #EnumItem("QCDscale_ggZH_yield"),
                 EnumItem("QCDscale_ggVV_bonly_yield"),
                 EnumItem("pdf_Higgs_gg_yield"),
                 EnumItem("pdf_Higgs_qq_yield"),
                 EnumItem("pdf_Higgs_ttH_yield"),
                 EnumItem("BRhiggs_hzz4l"),
                 EnumItem("PythiaScale"),
                 EnumItem("PythiaTune"),
                 EnumItem("lumi_13TeV"),
                 EnumItem("CMS_eff_e"),
                 EnumItem("CMS_eff_m"),
                 EnumItem("CMS_zz4e_zjets"),
                 EnumItem("CMS_zz4mu_zjets"),
                 EnumItem("CMS_zz2e2mu_zjets"),
                )
    def yamlfilename(self, channel=None):
      if self in ("pdf_Higgs_gg_yield", "pdf_Higgs_qq_yield", "pdf_Higgs_ttH_yield", "BRhiggs_hzz4l", "QCDscale_ggVV_bonly_yield"):
        return os.path.join(config.repositorydir, "helperstuff", "Datacards13TeV_Moriond2017", "STXSCards", "configs", "inputs", "systematics_theory_13TeV.yaml")
      if self in ("lumi_13TeV",):
        return os.path.join(config.repositorydir, "helperstuff", "Datacards13TeV_Moriond2017", "STXSCards", "configs", "inputs", "systematics_expt_13TeV.yaml")
      if self in ("CMS_eff_e", "CMS_eff_m"):
        if channel is None:
          raise ValueError("Need to give channel for {}".format(self))
        return os.path.join(config.repositorydir, "helperstuff", "Datacards13TeV_Moriond2017", "STXSCards", "configs", "inputs", "systematics_13TeV_{}.yaml".format(Channel(channel)))
      raise ValueError("{} is not from yaml!".format(self))

    def getfromyaml(self, channel=None):
      systname = str(self).replace("_yield", "")
      filename = self.yamlfilename(channel=channel)
      with open(filename) as f:
        y = yaml.load(f)
      if systname not in y: raise ValueError("{} not in {}!".format(systname, filename))
      if "Any" in y[systname]:
        values = y[systname]["Any"]
      elif "UnTagged" in y[systname] and all(y[systname][k] == y[systname]["UnTagged"] for k in y[systname]):
        values = y[systname]["UnTagged"]
      else:
        raise ValueError("Any not in {} in {}, and not all categories are the same!".format(systname, filename))
      return values

    def valuefromyaml(self, productionmode, channel=None):
      dct = self.getfromyaml(channel=channel)
      if productionmode in ("WH", "ZH"):
        if "{}_lep".format(productionmode) not in dct and "{}_had".format(productionmode) not in dct: return None
        lep = dct["{}_lep".format(productionmode)]
        had = dct["{}_had".format(productionmode)]
        if lep == had: return lep
        assert "/" in lep and "/" in had
        up = float(lep.split("/")[0])-1, float(had.split("/")[0])-1
        dn = float(lep.split("/")[1])-1, float(had.split("/")[1])-1
        rates = [totalrate(productionmode, _, 1) for _ in ("lep", "had")]
        assert all(_>0 for _ in up) and all(_<0 for _ in dn)
        up = 1 + sqrt((up[0]*rates[0])**2+(up[1]*rates[1])**2) / sum(rates)
        dn = 1 - sqrt((dn[0]*rates[0])**2+(dn[1]*rates[1])**2) / sum(rates)
        return up, dn
      else:
        return dct.get(productionmode.yamlsystname, None)

class YieldValue(MultiEnum, JsonDict):
    __metaclass__ = MultiEnumABCMeta
    enumname = "yieldvalue"
    enums = (Analysis, Category, Channel, ProductionMode)

    dictfile = os.path.join(config.repositorydir, "data", "yields.json")

    @property
    def keys(self):
        return (
                str(self.analysis),
                str(self.category),
                str(self.channel),
                str(self.productionmode),
               )

    def __float__(self):
        return self.__value

class YieldSystematicValue(MultiEnum, JsonDict):
    __metaclass__ = MultiEnumABCMeta
    enumname = "yieldsystematicvalue"
    enums = (YieldSystematic, Analysis, Category, Channel, ProductionMode)

    dictfile = os.path.join(config.repositorydir, "data", "categorysystematics.json")

    @property
    def keys(self):
        return (
                str(self.yieldsystematic),
                str(self.analysis),
                str(self.category),
                str(self.channel),
                str(self.productionmode),
               )

    def setvalue(self, value):
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

          if all(_ >= 1 for _ in list(value)):
            value = max(value)
          elif all(_ <= 1 for _ in list(value)):
            value = min(value)

        elif isinstance(value, Number) or value is None:
          pass
        else:
          raise ValueError("{!r} value {!r} should be None, a number, or a list (tuple, etc.) of length 2".format(self, origvalue))

        if value == 1: value = None

        super(YieldSystematicValue, self).setvalue(value)

    def __str__(self):
        if self.value is None or self.value == 1:
            return "-"
        if isinstance(self.value, Number):
            return str(self.value)
        if not (hasattr(self.value, "__len__") and len(self.value) == 2):
            raise ValueError("{!r} value '{!r}' should be None, a number, or a list (tuple, etc.) of length 2".format(self, self.value))
        return "{}/{}".format(self.value[0], self.value[1])

class VDecay(MyEnum):
  enumname = "vdecay"
  enumitems = (
               EnumItem("had"),
               EnumItem("lep"),
              )

class _TotalRate(MultiEnum):
  enums = [ProductionMode, Luminosity, VDecay]
  def check(self, *args):
    dontcheck = []
    if self.vdecay is not None and self.productionmode not in ("WH", "ZH"):
      raise ValueError("Can't have VDecay for {}\n{}".format(self.productionmode, args))
    if self.vdecay is None:
      dontcheck.append(VDecay)
    return super(_TotalRate, self).check(*args, dontcheck=dontcheck)

  @property
  def yamlrate(self):
    lumi = None
    result = 0
    tags = [tag.replace("Mor17", "").replace("tagged", "Tagged") for category in categories for tag in category.names if "Mor17" in tag]
    for channel in channels:
      filename = os.path.join(config.repositorydir, "helperstuff", "Datacards13TeV_Moriond2017", "STXSCards", "configs", "inputs", "yields_per_tag_category_13TeV_{}.yaml".format(channel))
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

      pnames = [p for p in self.productionmode.yamlratenames if self.vdecay is None or str(self.vdecay) in p]
      assert pnames

      for tag, p in itertools.product(tags, pnames):
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
  return _TotalRate(*args).rate

def count(fromsamples, tosamples, categorizations, alternateweights):
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
    for _ in alternateweights:
      if _ != "1":
        t.SetBranchStatus(_.weightname, 1)

    c = ROOT.TCanvas()
    for tosample, categorization, alternateweight in itertools.product(tosamples, categorizations, alternateweights):
        if tosample.productionmode == "WH" and tosample.hypothesis == "L1Zg": continue
        if alternateweight.issystematic and categorization.issystematic: continue
        t.Draw(categorization.category_function_name+":abs(Z1Flav*Z2Flav)", "{}*(ZZMass>{} && ZZMass<{})*{}".format(tosample.weightname(), config.m4lmin, config.m4lmax, alternateweight.weightname), "LEGO")
        h = c.FindObject("htemp")
        for i in range(h.GetNbinsY()):
            for channel in channels:
                toadd = h.GetBinContent(h.FindBin(channel.ZZFlav, i))
                result[tosample, categorization, alternateweight, Category.fromid(i), channel] += toadd
                result[tosample, categorization, alternateweight, Category.fromid(i)] += toadd
#        for k, v in result.iteritems(): print k, v
#        assert 0

        result[tosample,categorization,alternateweight] = h.Integral()

    return result
