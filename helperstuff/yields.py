from collections import Sequence
import itertools
import json
from math import sqrt
from numbers import Number
import os
import yaml

import ROOT

from combinehelpers import Luminosity
import config
from enums import Analysis, Category, categories, Channel, channels, EnumItem, MultiEnum, MultiEnumABCMeta, MyEnum, Production, ProductionMode
from samples import ReweightingSample, ReweightingSamplePlus, Sample
from utilities import cache, JsonDict, MultiplyCounter, TFile

class YieldSystematic(MyEnum):
    enumname = "yieldsystematic"
    enumitems = (
                 EnumItem("bTagSF"),
                 EnumItem("JES"),
                 EnumItem("QCDscale_fac_ggH"),
                 EnumItem("QCDscale_ren_ggH"),
                 EnumItem("QCDscale_fac_qqH"),
                 EnumItem("QCDscale_ren_qqH"),
                 EnumItem("QCDscale_fac_VH"),
                 EnumItem("QCDscale_ren_VH"),
                 EnumItem("QCDscale_fac_ttH"),
                 EnumItem("QCDscale_ren_ttH"),
                 EnumItem("QCDscale_fac_VV"),
                 EnumItem("QCDscale_ren_VV"),
                 EnumItem("EWcorr_VV"),
                 EnumItem("QCDscale_ggVV_bonly"),
                 EnumItem("pdf_Higgs_gg"),
                 EnumItem("pdf_Higgs_qq"),
                 EnumItem("pdf_Higgs_ttH"),
                 EnumItem("pdf_qq"),
                 EnumItem("BRhiggs_hzz4l"),
                 EnumItem("PythiaScale"),
                 EnumItem("PythiaTune"),
                 EnumItem("lumi_13TeV_2016"),
                 EnumItem("lumi_13TeV_2017"),
                 EnumItem("CMS_eff_e"),
                 EnumItem("CMS_eff_m"),
                 EnumItem("CMS_hzz4l_zz4e_zjets"),
                 EnumItem("CMS_hzz4l_zz4mu_zjets"),
                 EnumItem("CMS_hzz4l_zz2e2mu_zjets"),
                )

    def yamlfilename(self, channel=None):
      if self in ("pdf_Higgs_gg", "pdf_Higgs_qq", "pdf_Higgs_ttH", "pdf_qq", "BRhiggs_hzz4l", "QCDscale_ggVV_bonly", "QCDscale_ggH", "QCDscale_qqH", "QCDscale_VH", "QCDscale_ttH", "QCDscale_VV", "EWcorr_VV"):
        return os.path.join(config.repositorydir, "helperstuff", "Datacards13TeV_Moriond2017", "STXSCards", "configs", "inputs", "systematics_theory_13TeV.yaml")
      if self in ("CMS_eff_e", "CMS_eff_m"):
        if channel is None:
          raise ValueError("Need to give channel for {}".format(self))
        return os.path.join(config.repositorydir, "helperstuff", "Datacards13TeV_Moriond2017", "STXSCards", "configs", "inputs", "systematics_13TeV_{}.yaml".format(Channel(channel)))
      raise ValueError("{} is not from yaml!".format(self))

    def getfromyaml(self, channel=None):
      systname = str(self)
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
        assert False, (lep, had)
      else:
        return dct.get(productionmode.yamlsystname, None)

    def hardcodedvalue(self, production):
      if self == "lumi_13TeV_2016" and production.year == 2016: return 1.026
      if self == "lumi_13TeV_2017" and production.year == 2017: return 1.023
      return None

class YieldValue(MultiEnum, JsonDict):
    __metaclass__ = MultiEnumABCMeta
    enumname = "yieldvalue"
    enums = (Analysis, Category, Channel, ProductionMode, Production)

    dictfile = os.path.join(config.repositorydir, "data", "yields.json")

    @property
    def keys(self):
        return (
                str(self.production.productionforrate),
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
    enums = (YieldSystematic, Analysis, Category, Channel, ProductionMode, Production)

    dictfile = os.path.join(config.repositorydir, "data", "categorysystematics.json")

    @property
    def keys(self):
        return (
                str(self.production.productionforrate),
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
                value = tuple(float(_) for _ in value.split("/"))
              except ValueError:
                raise ValueError("string value {!r} can't be parsed as a number or 2 numbers with / in between".format(origvalue))

        if isinstance(value, Sequence):
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

    def getvalue(self):
        result = super(YieldSystematicValue, self).getvalue()
        if isinstance(result, list) and len(result) == 2:
          result = tuple(result)
        return result

    @property
    def latexstr(self):
        if self.value is None or self.value == 1:
            return "-"
        elif isinstance(self.value, Number):
            if abs(self.value-1) < .001: return "-"
            return "$\pm{:.1%}$".format(abs(self.value-1)).replace("%", r"\%")
        elif isinstance(self.value, Sequence) and len(self.value) == 2:
            result = [_-1 for _ in self.value]
            result = [float("{:.3f}".format(_)) for _ in result]
            if result[0] == result[1] == 0: return "-"
            if result[0] == 0: result.reverse()
            if result[1] == 0 or result[0] == abs(result[1]):
                return "$\pm{:.1%}$".format(abs(result[0])).replace("%", r"\%")
            if result[1] == abs(result[0]):
                return "$\mp{:.1%}$".format(abs(result[0])).replace("%", r"\%")
            return "${:+.1%}/{:+.1%}$".format(*(_-1 for _ in self.value)).replace("%", r"\%")

        raise ValueError("{!r} value '{!r}' should be None, a number, or a list (tuple, etc.) of length 2".format(self, self.value))

    def __str__(self):
        if self.value is None or self.value == 1:
            return "-"
        if isinstance(self.value, Number):
            return str(self.value)
        if not (isinstance(self.value, Sequence) and len(self.value) == 2):
            raise ValueError("{!r} value '{!r}' should be None, a number, or a list (tuple, etc.) of length 2".format(self, self.value))
        return "{}/{}".format(self.value[0], self.value[1])

class _TotalRate(MultiEnum):
  enums = [ProductionMode, Luminosity]

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

      pnames = [p for p in self.productionmode.yamlratenames]
      assert pnames

      for tag, p in itertools.product(tags, pnames):
        try:
          result += float(y[tag][p]) * float(self.luminosity) / lumi
        except ValueError:
          result += eval(y[tag][p].replace("@0", "125")) * float(self.luminosity) / lumi

    return result

  @property
  @cache
  def treerate(self):
    if self.productionmode == "WH":
      return sum(type(self)(self.luminosity, _).treerate for _ in ("WplusH", "WminusH"))

    result = 0
    c = ROOT.TCanvas()
    t = ROOT.TChain("candTree")
    t.SetBranchStatus("*", 0)
    t.SetBranchStatus("MC_weight_*", 1)
    t.SetBranchStatus("ZZMass", 1)
    t.SetBranchStatus("xsec", 1)
    t.SetBranchStatus("genxsec", 1)
    t.SetBranchStatus("genBR", 1)
    if self.productionmode in ("ggH", "ggZZ"): t.SetBranchStatus("KFactor_QCD_ggZZ_Nominal", 1)
    weightname = set()
    if self.productionmode == "VBF bkg":
      for flavor in "2e2mu", "4e", "4mu":
        t.Add(Sample(self.productionmode, flavor, self.production).withdiscriminantsfile())
        weightname.add(Sample(self.productionmode, flavor, self.production).weightname())
    elif self.productionmode in ("ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH"):
      if self.productionmode == "ggH":
        s = Sample(self.productionmode, "0+", "POWHEG", self.production)
        t.Add(Sample(s).withdiscriminantsfile())
        weightname.add("({}) * (genxsec * genBR * KFactor_QCD_ggZZ_Nominal / xsec)".format(s.weightname()))
      else:
        s = Sample(self.productionmode, "0+", "POWHEG", self.production, "Hff0+" if self.productionmode == "ttH" else None)
        t.Add(s.withdiscriminantsfile())
        weightname.add("({}) * (genxsec * genBR / xsec)".format(s.weightname()))
    else:
      assert False
    assert len(weightname) == 1, weightname
    weightname = weightname.pop()
    t.Draw("1", "{}*(ZZMass>{} && ZZMass<{})".format(weightname, config.m4lmin, config.m4lmax))
    return c.FindObject("htemp").Integral() * float(self.luminosity)

  @property
  def ratefromUlascan(self):
    from templates import Template
    rate = 0
    for ca, ch in itertools.product(channels, categories):
      t = Template(self.productionmode, str(self.production)+"_Ulascan", "0+" if self.productionmode.issignal else None, ca, ch, "fa2")
      with TFile(t.templatesfile.templatesfile()) as f:
        h = getattr(f, t.templatename())
        rate += h.Integral()
    if self.productionmode != "ZX": rate *= float(self.luminosity)
    return rate

  @property
  @cache
  def rate(self):
    if self.productionmode != "ttH": return self.ratefromUlascan
    elif self.productionmode == "VBF bkg": return self.treerate
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
    if any(_.productionmode in ("ggH", "ggZZ", "qqZZ") for _ in fromsamples): t.SetBranchStatus("KFactor_*", 1)
    t.SetBranchStatus("genxsec", 1)
    t.SetBranchStatus("genBR", 1)
    for _ in alternateweights:
      if _ in ("1", "EWcorrUp", "EWcorrDn", "PythiaScaleUp", "PythiaScaleDown"):
        pass
      else:
        t.SetBranchStatus(_.weightname, 1)
    if "PythiaScaleUp" in alternateweights or "PythiaScaleDown" in alternateweights:
      t.SetBranchStatus("PythiaWeight_*sr_muR0p25", 1)
      t.SetBranchStatus("PythiaWeight_*sr_muR4", 1)

    c = ROOT.TCanvas()
    for tosample, categorization, alternateweight in itertools.product(tosamples, categorizations, alternateweights):
        try:
            if tosample.productionmode == "WH" and tosample.hypothesis == "L1Zg": continue
            if alternateweight.issystematic and categorization.issystematic: continue
            t.Draw(categorization.category_function_name+":abs(Z1Flav*Z2Flav)", "{}*(ZZMass>{} && ZZMass<{})*{}".format(tosample.weightname(), config.m4lmin, config.m4lmax, alternateweight.weightname), "LEGO")
            h = c.FindObject("htemp")
            if tosample.productionmode == "ggH":
                t.GetEntry(0)
                assert all(_.productionmode == "ggH" for _ in fromsamples)
                h.Scale(t.genxsec * t.genBR * getattr(t, alternateweight.kfactorname) / tosample.SMxsec)
            for i in range(h.GetNbinsY()):
                for channel in channels:
                    toadd = h.GetBinContent(h.FindBin(channel.ZZFlav, i))
#                    print i, Category.fromid(i), channel, toadd
                    result[tosample, categorization, alternateweight, Category.fromid(i), channel] += toadd
                    result[tosample, categorization, alternateweight, Category.fromid(i)] += toadd
#            for k, v in result.iteritems(): print k, v
#            assert 0

            result[tosample,categorization,alternateweight] = h.Integral()
        except:
            if tosample == ReweightingSamplePlus("ZH", "0+", "POWHEG", "TuneUp") or tosample == ReweightingSamplePlus("ZH", "0+", "POWHEG", "TuneDn"):
              result[tosample,categorization,alternateweight] = result[ReweightingSamplePlus("ZH", "0+", "POWHEG"), categorization,alternateweight]
            else:
              print tosample, categorization, alternateweight
              raise

    return result

def check():
  filenames = set()
  for syst in YieldSystematic.items():
    for channel in channels:
      try:
        filenames.add(syst.yamlfilename(channel))
      except ValueError:
        pass
  y = {}
  for filename in filenames:
    with open(filename) as f:
      y.update(yaml.load(f))
  for key in y:
    if key == "CMS_zz4l_bkg_kdShape": continue
    if key.endswith("_cat") and key.replace("_cat", "") in y: continue
    try:
      YieldSystematic(key)
    except:
      raise ValueError("Don't have systematic for {}".format(key))

check()
del check
