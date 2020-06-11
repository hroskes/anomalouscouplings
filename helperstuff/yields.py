import array
from collections import Sequence
import itertools
import json
from math import sqrt
from numbers import Number
import os
import yaml

import numpy as np
import uncertainties

import ROOT

from combinehelpers import Luminosity
import config
from enums import AlternateWeight, Analysis, Category, categories, Channel, channels, EnumItem, MultiEnum, MultiEnumABCMeta, MyEnum, Production, ProductionMode
from samples import ReweightingSample, ReweightingSamplePlus, Sample
from utilities import cache, deprecate, JsonDict, MultiplyCounter, TFile

class YieldSystematic(MyEnum):
    enumname = "yieldsystematic"
    enumitems = (
                 EnumItem("CMS_btag_comb"),
                 EnumItem("CMS_scale_j", "CMS_scale_j_0PM", "JEC0PM"),
                 EnumItem("CMS_scale_j_0M", "JEC0M"),
                 EnumItem("CMS_scale_j_0PH", "JEC0PH"),
                 EnumItem("CMS_scale_j_0L1", "JEC0L1"),
                 EnumItem("CMS_scale_j_0L1Zg", "JEC0L1Zg"),
                 EnumItem("QCDscale_muF_ggH"),
                 EnumItem("QCDscale_muF_qqH"),
                 EnumItem("QCDscale_muF_VH"),
                 EnumItem("QCDscale_muF_ttH"),
                 EnumItem("QCDscale_muF_bbH"),
                 EnumItem("QCDscale_muF_VV"),
                 EnumItem("QCDscale_muR_ggH"),
                 EnumItem("QCDscale_muR_qqH"),
                 EnumItem("QCDscale_muR_VH"),
                 EnumItem("QCDscale_muR_ttH"),
                 EnumItem("QCDscale_muR_bbH"),
                 EnumItem("QCDscale_muR_VV"),
                 EnumItem("EWcorr_qqZZ"),
                 EnumItem("pdf_Higgs_gg"),
                 EnumItem("pdf_Higgs_qqbar"),
                 EnumItem("pdf_qqbar"),
                 EnumItem("pdf_As_Higgs_gg"),
                 EnumItem("pdf_As_Higgs_qqbar"),
                 EnumItem("pdf_As_qqbar"),
                 EnumItem("hzz_br"),
                 EnumItem("CMS_pythia_scale"),
                 EnumItem("CMS_pythia_tune"),
                 EnumItem("lumi_13TeV_2016"),
                 EnumItem("lumi_13TeV_2017"),
                 EnumItem("lumi_13TeV_2018"),
                 EnumItem("CMS_eff_e"),
                 EnumItem("CMS_eff_m"),
                 EnumItem("zjet_2e2mu", "zjet_2e2mu_2016"),
                 EnumItem("zjet_4e", "zjet_4e_2016"),
                 EnumItem("zjet_4mu", "zjet_4mu_2016"),
                 EnumItem("zjet_2e2mu_2017"),
                 EnumItem("zjet_4e_2017"),
                 EnumItem("zjet_4mu_2017"),
                 EnumItem("zjet_2e2mu_2018"),
                 EnumItem("zjet_4e_2018"),
                 EnumItem("zjet_4mu_2018"),
                 EnumItem("THU_ggH_Mu"),
                 EnumItem("THU_ggH_Res"),
                 EnumItem("THU_ggH_Mig01"),
                 EnumItem("THU_ggH_Mig12"),
                 EnumItem("THU_ggH_VBF2j"),
                 EnumItem("THU_ggH_VBF3j"),
                 EnumItem("THU_ggH_PT60"),
                 EnumItem("THU_ggH_PT120"),
                 EnumItem("THU_ggH_qmtop"),
                )

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
        if self.copyfromotheryieldsystematicvalue is not None: return

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

          value = list(value)
          for i, v in enumerate(value[:]):
            if np.isclose(v, 1, rtol=3e-4):
              value[i] = 1

          if all(_ == 1 for _ in value):
            value = 1
          elif np.isclose(*value):
            raise ValueError(str(value))

        elif isinstance(value, Number):
          if np.isclose(value, 1): value = 1
        elif value is None:
          pass
        else:
          raise ValueError("{!r} value {!r} should be None, a number, or a list (tuple, etc.) of length 2".format(self, origvalue))

        if value == 1: value = None

        super(YieldSystematicValue, self).setvalue(value)

    def getvalue(self):
        if self.copyfromotheryieldsystematicvalue is not None: return self.copyfromotheryieldsystematicvalue.getvalue()
        result = super(YieldSystematicValue, self).getvalue()
        if isinstance(result, list) and len(result) == 2:
          result = tuple(result)
        return result

    @property
    def downvalue(self):
      value = self.value
      if value is None: return 1
      if isinstance(value, Number): return 1./value
      assert len(value) == 2
      return value[0]

    @property
    def upvalue(self):
      value = self.value
      if value is None: return 1
      if isinstance(value, Number): return value
      assert len(value) == 2
      return value[1]

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

    @property
    def copyfromotheryieldsystematicvalue(self):
        kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in type(self).needenums}
        if self.production == "190821_2016" and self.yieldsystematic == "CMS_pythia_scale":
            kwargs["production"] = "190821_2017"
        elif self.yieldsystematic in ("JEC0PH", "JEC0M", "JEC0L1", "JEC0L1Zg"):
            kwargs["yieldsystematic"] = "JEC0PM"
        else:
            return None
        return type(self)(*kwargs.itervalues())

def count(fromsamples, tosamples, categorizations, alternateweights):
    assert len(tosamples) == 1
    GEN = {_.production.GEN for _ in fromsamples}
    assert len(GEN) == 1
    GEN = GEN.pop()
    productionmode = {_.productionmode for _ in fromsamples}
    assert len(productionmode) == 1
    productionmode = productionmode.pop()

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
    if productionmode in ("ggH", "ggZZ", "qqZZ"): t.SetBranchStatus("KFactor_*", 1)
    if productionmode != "ZX":
      t.SetBranchStatus("xsec", 1)
      t.SetBranchStatus("genxsec", 1)
      t.SetBranchStatus("genBR", 1)
      if not GEN: t.SetBranchStatus("LHEweight_*", 1)
    for _ in alternateweights:
      if _ in ("1", "EWcorrUp", "EWcorrDn", "PythiaScaleUp", "PythiaScaleDown"):
        pass
      else:
        t.SetBranchStatus("MC_weight_nominal", 1)
    if "PythiaScaleUp" in alternateweights or "PythiaScaleDown" in alternateweights:
      t.SetBranchStatus("PythiaWeight_*sr_muR0p25", 1)
      t.SetBranchStatus("PythiaWeight_*sr_muR4", 1)
    if any(_.isTHUggH for _ in alternateweights):
      t.SetBranchStatus("qcd_ggF_uncertSF", 1)

    c = ROOT.TCanvas()

    for tosample, categorization, alternateweight in itertools.product(tosamples, categorizations, alternateweights):
        if getattr(tosample, "extension", None) is not None:
            kwargs = {enum.enumname: getattr(tosample, enum.enumname) for enum in type(tosample).needenums}
            kwargs["extension"] = None
            tosample = type(tosample)(*kwargs.itervalues())
        try:
            if productionmode == "WH" and tosample.hypothesis == "L1Zg": continue
            if alternateweight.issystematic and categorization.issystematic: continue
            t.Draw(categorization.category_function_name+":abs(Z1Flav*Z2Flav)", "MC_weight_nominal*(ZZMass>{} && ZZMass<{})*{}".format(config.m4lmin, config.m4lmax, alternateweight.weightname), "LEGO")
            h = c.FindObject("htemp")
            t.GetEntry(0)
            if productionmode != "ZX":
                h.Scale(t.xsec / (t.genxsec * t.genBR))
            for i in range(h.GetNbinsY()):
                for channel in channels:
                    binnumber = h.FindBin(channel.ZZFlav, i)
                    toadd = uncertainties.ufloat(h.GetBinContent(binnumber), h.GetBinError(binnumber))
                    if not toadd.n and not toadd.s: continue
#                    print i, Category.fromid(i), channel, toadd
                    result[tosample, categorization, alternateweight, Category.fromid(i), channel] += toadd
                    result[tosample, categorization, alternateweight, Category.fromid(i)] += toadd
#            for k, v in result.iteritems(): print k, v
#            assert 0

            error = array.array("d", [0])
            result[tosample,categorization,alternateweight] = uncertainties.ufloat(h.IntegralAndError(1, h.GetNbinsX(), 1, h.GetNbinsY(), error), error[0])
        except:
            print tosample, categorization, alternateweight
            for fromsample in fromsamples:
              print fromsample.withdiscriminantsfile()
            raise

    if GEN: #scale to numbers from 18-002 https://arxiv.org/pdf/1901.00174.pdf
        if productionmode == "VBF":
            scaleto = 4.7 + 0.3 + 5.7
        elif productionmode == "WplusH":
            scaleto = (.3 + .7 + 2.1) * 8.400E-01 / (8.400E-01 + 5.328E-01)  #ratio from YR4
        elif productionmode == "WminusH":
            scaleto = (.3 + .7 + 2.1) * 5.328E-01 / (8.400E-01 + 5.328E-01)  #ratio from YR4
        elif productionmode == "ZH":
            scaleto = .2 + .5 + 1.5
        elif productionmode == "ggH":
            scaleto = 5.5 + 3.2 + 98.9 + .2 + .1 + 1.1 + .1 + .1 + 1.1 #include ttH and bbH
        elif productionmode == "qqZZ":
            scaleto = 1.6 + 1.5 + 120.3 + 5.2 + 3.0 + 46.3 + .8 + .3 + 12.7 + .2 + .1 + 1.5  #include Z+X, gg, VV backgrounds
        else:
            raise ValueError(productionmode)
        scaleto /= 77.5
        scaleby = scaleto / result[tosample,categorization,AlternateWeight("1")]
        print tosample, scaleby
        result *= scaleby

    return result
