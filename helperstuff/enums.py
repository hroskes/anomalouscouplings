import abc
from collections import OrderedDict
from itertools import combinations, permutations
import os
import re

import config
from utilities import cache, deprecate, generatortolist_condition, tfiles

class EnumItem(object):
    def __init__(self, name, *other):
        self.name = name
        self.names = tuple([name] + list(other))

    def __str__(self):
        return self.name
    def __hash__(self):
        return hash(self.names)

    def __eq__(self, other):
        if isinstance(other, (int, basestring)):
            return other in self.names
        if isinstance(other, type(self)):
            return str(self) == str(other)
        return NotImplemented
    def __ne__(self, other):
        return not self == other

class MyEnum(object):
    def __init__(self, value):
        if isinstance(value, (type(self), EnumItem)):
            value = str(value)
        for item in self.enumitems:
            if value in item.names and not isinstance(value, MyEnum):
                self.item = item
                break
        else:
            raise ValueError("{} is not a member of enum {}!  Valid choices:\n".format(value, type(self).__name__)
                               + "\n".join(" aka ".join(str(name) for name in item.names) for item in self.enumitems))

    def __str__(self):
        return str(self.item)
    def __repr__(self):
        return "{}({})".format(type(self).__name__, repr(self.item.name))

    def __eq__(self, other):
        if isinstance(other, MyEnum) and not isinstance(other, type(self)) and not isinstance(self, type(other)):
            return False
        try:
            other = type(self)(other)
            return self.item == other.item
        except ValueError:
            if other is None or (type(other) == str and other == ""):
                return False
            raise
    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.item)

    @classmethod
    def items(cls, condition=lambda item: True):
        return [cls(item) for item in cls.enumitems if condition(cls(item))]

    @property
    def names(self):
        return self.item.names

class Channel(MyEnum):
    enumname = "channel"
    enumitems = (
                 EnumItem("2e2mu", "2mu2e"),
                 EnumItem("4mu"),
                 EnumItem("4e"),
                )
    @property
    def ZZFlav(self):
        if self == "2e2mu": return 13*13*11*11
        if self == "4e":    return 11*11*11*11
        if self == "4mu":   return 13*13*13*13
        assert False

class Flavor(MyEnum):
    enumname = "flavor"
    enumitems = (
                 EnumItem("2e2mu", "2mu2e"),
                 EnumItem("4mu"),
                 EnumItem("4e"),
                 EnumItem("4tau"),
                 EnumItem("2mu2tau", "2tau2mu"),
                 EnumItem("2e2tau", "2tau2e"),
                )
    @property
    def hastaus(self):
        return self in ("2e2tau", "2mu2tau", "4tau")
    @property
    def countersbin(self):
        if self == "4mu": return 2
        if self == "4e": return 3
        if self == "2e2mu": return 4
        if self == "4tau": return 9
        assert False

def mixturepermutations_4d(enumitemsalreadythere):
  def inner():
    pures = "a3", "a2", "L1", "L1Zg"

    #pure hypotheses
    yield "0+",
    for _ in pures:
      yield _,

    for proddec in "dec", "prod", "proddec":
      kwargs = {"proddec": proddec, "sign": "-" if proddec == "proddec" else ""}
      #mixtures with SM
      for kwargs["ai"] in pures:
        yield "f{ai}{proddec}{sign}0.5".format(**kwargs),
      #mixtures with each other
      for kwargs["ai"], kwargs["aj"] in combinations(pures, 2):
        yield (
          "f{ai}{proddec}0.5f{aj}{proddec}{sign}0.5".format(**kwargs),
          "f{aj}{proddec}{sign}0.5f{ai}{proddec}0.5".format(**kwargs),
        )

      if proddec == "dec": yield None, "decay is done"

      #to get the terms with 3 unique couplings:
      #  there are 3 different new terms for each set of 3 couplings,
      #  since exactly one of the couplings is squared
      #  so we need 3 more mixtures
      #  prod, dec, proddec (with a minus sign somewhere for variety) will work
      for kwargs["ai"], kwargs["aj"] in combinations(pures, 2):
        yield (
          "f{ai}{proddec}0.33f{aj}{proddec}{sign}0.33".format(**kwargs),
          "f{aj}{proddec}{sign}0.33f{ai}{proddec}0.33".format(**kwargs),
        )

      for kwargs["ai"], kwargs["aj"], kwargs["ak"] in combinations(pures, 3):
        yield tuple(
          "".join(permutation).format(**kwargs)
            for permutation
            in permutations(("f{ai}{proddec}0.33", "f{aj}{proddec}0.33", "f{ak}{proddec}{sign}0.33"))
        )

    #for each set of 4 couplings, there's only one possible new term
    #proddec is the way to go, since it includes both
    #I'll keep it general using combinations, just in case that ever becomes necessary
    #it's not clear that a minus sign would help, so I'll leave it out
    kwargs = {}
    for kwargs["ai"], kwargs["aj"], kwargs["ak"] in combinations(pures, 3):
      yield tuple(
        "".join(permutation).format(**kwargs)
          for permutation
          in permutations(("f{ai}proddec0.25", "f{aj}proddec0.25", "f{ak}proddec0.25"))
      )
    for kwargs["ai"], kwargs["aj"], kwargs["ak"], kwargs["al"] in combinations(pures, 4):
      yield tuple(
        "".join(permutation).format(**kwargs)
          for permutation
          in permutations(("f{ai}proddec0.25", "f{aj}proddec0.25", "f{ak}proddec0.25", "f{al}proddec0.25"))
      )

  result = []
  allofthem = 0
  for enumitemnames in inner():
    if enumitemnames == (None, "decay is done"):
      assert allofthem == len(list(combinations(range(6), 2))), (allofthem, len(list(combinations(range(6), 2))))
      continue

    allofthem += 1
    enumitemnames_modified = [_ for _ in enumitemnames]
    enumitemnames_modified += [re.sub("(?<!prod)dec", "", _) for _ in enumitemnames_modified]
    if any(_ in enumitemsalreadythere for _ in enumitemnames_modified):
      continue
    result.append(EnumItem(*enumitemnames_modified))

  assert allofthem == len(list(combinations(range(8), 4))), (allofthem, len(list(combinations(range(8), 4))))
  return tuple(result)

class Hypothesis(MyEnum):
    enumname = "hypothesis"
    enumitems = (
                 EnumItem("0+", "SM", "scalar", "0PM", "a1", "g1"),
                 EnumItem("a2", "0h+", "0PH", "g2"),
                 EnumItem("0-", "a3", "PS", "pseudoscalar", "0M", "g4"),
                 EnumItem("L1", "Lambda1", "0L1"),
                 EnumItem("L1Zg", "0L1Zg", "L1Zgs"),
                 EnumItem("fa20.5", "fa2dec0.5", "fa2+0.5", "fa2dec+0.5"),
                 EnumItem("fa30.5", "fa3dec0.5", "fa3+0.5", "fa3dec+0.5"),
                 EnumItem("fL10.5", "fL1dec0.5", "fL1+0.5", "fL1dec+0.5"),
                 EnumItem("fL1Zg0.5", "fL1Zgdec0.5", "fL1Zg+0.5", "fL1Zgdec+0.5"),
                 EnumItem("fa2prod0.5", "fa2prod+0.5"),
                 EnumItem("fa3prod0.5", "fa3prod+0.5"),
                 EnumItem("fL1prod0.5", "fL1prod+0.5"),
                 EnumItem("fL1Zgprod0.5", "fL1Zgprod+0.5"),
                 EnumItem("fa2proddec0.5", "fa2proddec+0.5"),
                 EnumItem("fa3proddec0.5", "fa3proddec+0.5"),
                 EnumItem("fL1proddec0.5", "fL1proddec+0.5"),
                 EnumItem("fL1Zgproddec0.5", "fL1Zgproddec+0.5"),
                 EnumItem("fa2dec-0.5", "fa2-0.5"),
                 EnumItem("fa3dec-0.5", "fa3-0.5"),
                 EnumItem("fL1dec-0.5", "fL1-0.5"),
                 EnumItem("fL1Zgdec-0.5", "fL1Zg-0.5"),
                 EnumItem("fa2prod-0.5"),
                 EnumItem("fa3prod-0.5"),
                 EnumItem("fL1prod-0.5"),
                 EnumItem("fL1Zgprod-0.5"),
                 EnumItem("fa2proddec-0.5"),
                 EnumItem("fa3proddec-0.5"),
                 EnumItem("fL1proddec-0.5"),
                 EnumItem("fL1Zgproddec-0.5"),
                 EnumItem("fa2dec-0.9", "fa2-0.9"),
                )
    enumitems = enumitems + mixturepermutations_4d(enumitems)

    @property
    def ispure(self):
        return self in ("0+", "0-", "0h+", "L1", "L1Zg")
    @property
    def couplingname(self):
        if self == "0+": return "g1"
        if self == "a2": return "g2"
        if self == "0-": return "g4"
        if self == "L1": return "g1prime2"
        if self == "L1Zg": return "ghzgs1prime2"
        assert False, self
    @property
    def combinename(self):
        for _ in "0PM", "0PH", "0M", "0L1", "0L1Zg":
            if self == _:
                return _
        assert False, self

class HffHypothesis(MyEnum):
    enumname = "hffhypothesis"
    enumitems = (
                 EnumItem("Hff0+"),
                 EnumItem("Hff0-"),
                 EnumItem("fCP0.5"),
                )

class ProductionMode(MyEnum):
    enumname = "productionmode"
    enumitems = (
                 EnumItem("ggH"),
                 EnumItem("VBF", "qqH"),
                 EnumItem("HJJ", "H+jj"),
                 EnumItem("ZH"),
                 EnumItem("WH"),
                 EnumItem("ttH"),
                 EnumItem("qqZZ", "bkg_qqzz"),
                 EnumItem("ggZZ", "bkg_ggzz"),
                 EnumItem("VBFbkg", "VBF bkg", "bkg_vbf"),
                 EnumItem("ZX", "bkg_zjets"),
                 EnumItem("data"),
                 EnumItem("WplusH"),
                 EnumItem("WminusH"),
                 EnumItem("bbH"),
                 EnumItem("tqH"),
                )
    @property
    def combinename(self):
        for name in "ggH", "qqH", "ZH", "WH", "ttH", "bbH", "bkg_qqzz", "bkg_ggzz", "bkg_vbf", "bkg_zjets":
            if self == name:
                return name
        assert False
    @property
    def yamlratenames(self):
        if self == "VBF":
            return ["qqH"]
        elif self == "ZX":
            return ["zjets"]
        elif self in ["ZH", "WH"]:
            return ["{}_{}".format(self, dec) for dec in ("lep", "had")]
        return [str(self)]
    @property
    def yamlratenames2015(self):
        if self == "VBF":
            return ["qqH"]
        elif self == "ZX":
            return ["zjets"]
        return [str(self)]
    @property
    def yamlsystname(self):
        if self in ("VBF", "VBF bkg"):
            return "qqH"
        elif self == "ZX":
            return "zjets"
        return str(self)
    @property
    def isbkg(self):
        if self in ("ggH", "VBF", "ZH", "WH", "ttH", "HJJ", "WplusH", "WminusH", "bbH", "tqH"):
            return False
        elif self in ("ggZZ", "qqZZ", "VBF bkg", "ZX"):
            return True
        assert False
    @property
    def issignal(self):
        return not self.isbkg
    @property
    @cache
    def validhypotheses(self):
        if self == "ggH":
            return Hypothesis.items(lambda x: x in decayonlyhypotheses)
        if self in ("ttH", "HJJ", "bbH"):
            return Hypothesis.items(lambda x: x in decayonlyhypotheses)
        if self == "VBF":
            return Hypothesis.items(lambda x: x in proddechypotheses)
        if self == "ZH":
            return Hypothesis.items(lambda x: x in proddechypotheses)
        if self == "WH":
            return Hypothesis.items(lambda x: x in proddechypotheses)
        if self in ("WplusH", "WminusH", "tqH"):
            return Hypothesis.items(lambda x: x == "0+")
        assert False
    def generatedhypotheses(self, production):
        if not production.LHE:
            if self == "ggH":
                return Hypothesis.items(lambda x: x in ("0+", "0-", "a2", "L1", "fa30.5", "fa20.5", "fL10.5"))
            if self in ("VBF", "ZH", "WH"):
                return Hypothesis.items(lambda x: x in ("0+", "0-", "a2", "L1", "fa3prod0.5", "fa2prod0.5", "fL1prod0.5"))
            if self in ("WplusH", "WminusH", "ttH", "HJJ", "bbH", "tqH"):
                return Hypothesis.items(lambda x: x == "0+")
        else:
            if self == "ggH":
                return Hypothesis.items(lambda x: x in ("0+", "L1", "L1Zg", "fL10.5", "fL1Zg0.5", "fL10.5fL1Zg0.5"))
        assert False, self
    @generatortolist_condition(lambda x: tfiles[x.withdiscriminantsfile()].candTree.GetEntries())
    def allsamples(self, production):
        from samples import Sample
        from utilities import tfiles
        if self == "VBF bkg":
            for flavor in "2e2mu", "4e", "4mu":
                yield Sample(self, flavor, production)
        elif self == "ggZZ":
            for flavor in flavors:
                yield Sample(self, flavor, production)
        elif self.isbkg:
            yield Sample(self, production)
        else:
            hff = None
            if self in ("HJJ", "ttH"): hff = "Hff0+"
            for h in self.generatedhypotheses(production):
                yield Sample(self, h, hff, production)

    @property
    def QCDfacsystematicname(self):
      if self == "ggH": return "QCDscale_fac_ggH"
      if self == "qqH": return "QCDscale_fac_qqH"
      if self in ("ZH", "WH"): return "QCDscale_fac_VH"
      if self == "ttH": return "QCDscale_fac_ttH"
      if self == "bbH": return "QCDscale_fac_bbH"
      if self == "qqZZ": return "QCDscale_fac_VV"
      return None

    @property
    def QCDrensystematicname(self):
      if self == "ggH": return "QCDscale_ren_ggH"
      if self == "qqH": return "QCDscale_ren_qqH"
      if self in ("ZH", "WH"): return "QCDscale_ren_VH"
      if self == "ttH": return "QCDscale_ren_ttH"
      if self == "bbH": return "QCDscale_ren_bbH"
      if self == "qqZZ": return "QCDscale_ren_VV"
      return None

    @property
    def pdfvariationsystematicname(self):
      if self in ("ggH", "ttH", "bbH"): return "pdf_variation_Higgs_gg"
      if self in ("qqH", "ZH", "WH"): return "pdf_variation_Higgs_qqbar"
      if self == "qqZZ": return "pdf_variation_qqbar"
      return None

    @property
    def pdfasmzsystematicname(self):
      if self in ("ggH", "ttH", "bbH"): return "pdf_asmz_Higgs_gg"
      if self in ("qqH", "ZH", "WH"): return "pdf_asmz_Higgs_qqbar"
      if self == "qqZZ": return "pdf_asmz_qqbar"
      return None

    def workspaceshapesystematics(self, category):
      result = []
      if self in ("ggH", "qqH", "ZH", "WH", "ttH", "bbH"):
        if (
            config.applym4lshapesystematicsUntagged and category == "Untagged"
            or config.applym4lshapesystematicsVBFVHtagged and category != "Untagged"
            or config.applym4lshapesystematicsggH and self == "ggH"
            or config.applym4lshapesystematicsggHUntagged and self == "ggH" and category == "Untagged"
            or config.applym4lshapesystematicsdiagonal and (
                                                            self == "ggH" and category == "Untagged"
                                                            or self == "VBF" and category == "VBFtagged"
                                                            or self in ("ZH", "WH") and category == "VHHadrtagged"
                                                           )
           ):
          if config.combinem4lshapesystematics:
            result += ["ScaleRes"]
          else:
            result += ["Scale", "Res"]
      if self == "ggH" and category in ("VBFtagged", "VHHadrtagged") and config.applyMINLOsystematics:
        result += ["MINLO"]
      if self == "ggH" and category in ("VBFtagged", "VHHadrtagged") and config.applyJECshapesystematics:
        result += ["CMS_scale_j_13TeV_2016", "CMS_scale_j_13TeV_2017"]
      return [WorkspaceShapeSystematic(_) for _ in result]

    def alternateweights(self, year):
      def yearcondition(systematic):
        if systematic == "PythiaScaleUp" or systematic == "PythiaScaleDown" and year == 2016: return False
        return True

      if self in ("ggH", "qqH", "ZH", "WH", "ttH"):
        return AlternateWeight.items(lambda x: x!="EWcorrUp" and x!="EWcorrDn" and yearcondition(x))
      if self == "qqZZ":
        return AlternateWeight.items(yearcondition)
      assert False

class WorkspaceShapeSystematic(MyEnum):
    enumname = "workspaceshapesystematic"
    enumitems = (
                 EnumItem("CMS_res_", "Res"),
                 EnumItem("CMS_scale_", "Scale"),
                 EnumItem("CMS_scaleres_", "ScaleRes"),
                 EnumItem("QCDscale_ggH2in", "MINLO"),
                 EnumItem("CMS_scale_j_13TeV_2016"),
                 EnumItem("CMS_scale_j_13TeV_2017"),
                )
    @property
    def isperchannel(self):
        if self in ("Res", "Scale", "ScaleRes"): return True
        if self in ("MINLO", "CMS_scale_j_13TeV_2016", "CMS_scale_j_13TeV_2017"): return False
        assert False, self

    @property
    def years(self):
        if self == "CMS_scale_j_13TeV_2016": return 2016,
        if self == "CMS_scale_j_13TeV_2017": return 2017,
        if self in ("Res", "Scale", "ScaleRes", "MINLO"): return 2016, 2017
        assert False, self

    @property
    def nickname(self):
      for _ in "Res", "Scale", "ScaleRes", "MINLO":
        if self == _:
          return _
      if self in ("CMS_scale_j_13TeV_2016", "CMS_scale_j_13TeV_2017"): return "JEC"
      return str(self)

class SystematicDirection(MyEnum):
    enumname = "systematicdirection"
    enumitems = (
                 EnumItem("Up"),
                 EnumItem("Down", "Dn"),
                )

class ShapeSystematic(MyEnum):
    enumname = "shapesystematic"
    enumitems = (
                 EnumItem(""),
                 EnumItem("ResUp"),
                 EnumItem("ResDown"),
                 EnumItem("ScaleUp"),
                 EnumItem("ScaleDown"),
                 EnumItem("ScaleResUp", "ResScaleUp"),
                 EnumItem("ScaleResDown", "ResScaleDown"),
                 EnumItem("JECUp"),
                 EnumItem("JECDn", "JECDown"),
                 EnumItem("ZXUp"),
                 EnumItem("ZXDown", "ZXDn"),
                 EnumItem("MINLO_SM"),
                 EnumItem("MINLOUp"),
                 EnumItem("MINLODn", "MINLODown"),
                )
    @property
    def Dbkgappendname(self):
        for _ in ("ScaleUp", "ScaleDown", "ResUp", "ResDown"):
            if self == _:
                return "_" + _
        return ""
    @property
    def categoryappendname(self):
        for _ in ("JECUp", "JECDn"):
            if self == _:
                return "_" + _
        if self in ("", "ResUp", "ResDown", "ScaleUp", "ScaleDown", "MINLO_SM", "ZXUp", "ZXDn"):
            return ""
        assert False, self
    def appliesto(self, templategroup):
        if self == "":
            return True
        if self in ("JECUp", "JECDn"):
            return templategroup in ("ggh", "zh", "wh")
        if self in ("ResUp", "ResDown", "ScaleUp", "ScaleDown", "ScaleResUp", "ScaleResDown"):
            return templategroup in ("ggh", "vbf", "zh", "wh", "tth", "bbh")
        if self in ("ZXUp", "ZXDown"):
            return templategroup == "bkg"
        if self in ("MINLO_SM", "MINLOUp", "MINLODn"):
            return templategroup == "ggh"
        assert False

class JECSystematic(MyEnum):
    enumname = "jecsystematic"
    enumitems = (
                 EnumItem("JECNominal", "Nominal"),
                 EnumItem("JECUp", "Up"),
                 EnumItem("JECDn", "Down", "JecDown"),
                )
    @property
    def appendname(self):
        if self == "JECNominal": return ""
        return "_" + str(self)
    @property
    def njetsappendname(self):
        return self.appendname.replace("JEC", "jec")

class BTagSystematic(MyEnum):
    enumname = "btagsystematic"
    enumitems = (
                 EnumItem("bTagSFNominal", "Nominal"),
                 EnumItem("bTagSFUp", "Up"),
                 EnumItem("bTagSFDn", "Dn"),
                )
    @property
    def appendname(self):
        if self == "bTagSFNominal": return ""
        return "_" + str(self)
    @property
    def njetsappendname(self):
        return "_"+str(self).replace("Nominal", "")

class TemplateGroup(MyEnum):
    enumname = "templategroup"
    enumitems = (
                 EnumItem("ggh"),
                 EnumItem("vbf"),
                 EnumItem("zh"),
                 EnumItem("wh"),
                 EnumItem("tth"),
                 EnumItem("bbh"),
                 EnumItem("background", "bkg"),
                 EnumItem("DATA"),
                )

class Analysis(MyEnum):
    enumname = "analysis"
    enumitems = (
                 EnumItem("fa3"),
                 EnumItem("fa2"),
                 EnumItem("fL1"),
                 EnumItem("fL1Zg"),
                 EnumItem("fL1fL1Zg_DL1_DL1L1Zgint"),
                 EnumItem("fL1fL1Zg", "fL1fL1Zg_DL1_DL1Zgint"),
                 EnumItem("fL1fL1Zg_DeR_DeLeR"),
                 EnumItem("fL1fL1Zg_DeR_DeLint"),
                 EnumItem("fL1fL1Zg_m1_m2"),
                 EnumItem("fL1fL1Zg_m1_phi"),
                 EnumItem("fL1fL1Zg_m2_phi"),
                 EnumItem("fa3_STXS"),
                 EnumItem("fa3_multiparameter_nodbkg"),
                 EnumItem("fa3_only6bins"),
                 EnumItem("fa3fa2fL1fL1Zg"),
                )
    def title(self, latex=False, superscript=None):
        if self.dimensions > 1: return self.fais[0].title(latex=latex, superscript=superscript)

        if self.fais == ("fa3",):
            result = "f_{{a3}}"
        elif self == "fa2":
            result = "f_{{a2}}"
        elif self == "fL1":
            result = "f_{{{Lambda}1}}"
        elif self == "fL1Zg":
            result = "f_{{{Lambda}1}}^{{{Z}{gamma}}}"
        else:
            assert False

        if latex:
            repmap = {"Lambda": r"\Lambda", "Z": r"\Z", "gamma": r"\gamma"}
        else:
            repmap = {"Lambda": "#Lambda", "Z": "Z", "gamma": r"#gamma"}
        result = result.format(**repmap)

        if superscript is not None:
            if "^" in result:
                result = result.rstrip("}") + ", " + str(superscript) + "}"
            else:
                result += "^{" + str(superscript) + "}"

        return result

    @property
    def phi(self):
        if self.dimensions > 1: return self.fais[0].phi
        if self.fais == ("fa3",):
            return "#phi_{a3}"
        if self == "fa2":
            return "#phi_{a2}"
        if self == "fL1":
            return "#phi_{#Lambda1}"
        if self == "fL1Zg":
            return "#phi_{#Lambda1}^{Z#gamma}"
        assert False
    @property
    def phi_lower(self):
        return self.phi.replace("{", "{#lower[-0.25]{").replace("}", "}}")
    @property
    def couplingname(self):
        if self.fais == ("fa3",): return "g4"
        if self == "fa2": return "g2"
        if self == "fL1": return "g1prime2"
        if self == "fL1Zg": return "ghzgs1prime2"
        assert False, self
    @property
    def couplingnames(self):
        if self.isfL1fL1Zg: return "g1prime2", "ghzgs1prime2"
        if self == "fa3fa2fL1fL1Zg": return "g4", "g2", "g1prime2", "ghzgs1prime2"
        if self.dimensions == 1: return self.couplingname,
        assert False, self
    @property
    def couplingtitle(self):
        if self.fais == ("fa3",): return "a_{3}"
        if self == "fa2": return "a_{2}"
        if self == "fL1": return "#Lambda_{1}"
        if self == "fL1Zg": return "#Lambda_{1}^{Z#gamma}"
    @property
    def latexfai(self):
        if self.fais == ("fa3",):
            return r"\fcospAC{3}"
        if self == "fa2":
            return r"\fcospAC{2}"
        if self == "fL1":
            return r"\fcospLC{1}"
        if self == "fL1Zg":
            return r"\fcospLZGs"
        assert False, self
    @property
    def purehypotheses(self):
        if self.fais == ("fa3",):
            return Hypothesis("0+"), Hypothesis("0-")
        if self == "fa2":
            return Hypothesis("0+"), Hypothesis("a2")
        if self == "fL1":
            return Hypothesis("0+"), Hypothesis("L1")
        if self == "fL1Zg":
            return Hypothesis("0+"), Hypothesis("L1Zg")
        if self.isfL1fL1Zg:
            return Hypothesis("0+"), Hypothesis("L1"), Hypothesis("L1Zg")
        if self == "fa3fa2fL1fL1Zg":
            return Hypothesis("0+"), Hypothesis("a3"), Hypothesis("a2"), Hypothesis("L1"), Hypothesis("L1Zg")
        assert False, self
    @property
    def mixdecayhypothesis(self):
        if self.fais == ("fa3",):
            return Hypothesis("fa3dec0.5")
        if self == "fa2":
            return Hypothesis("fa2dec0.5")
        if self == "fL1":
            return Hypothesis("fL1dec0.5")
        if self == "fL1Zg":
            return Hypothesis("fL1Zgdec0.5")
        assert False
    @property
    def mixprodhypothesis(self):
        if self.fais == ("fa3",):
            return Hypothesis("fa3prod0.5")
        if self == "fa2":
            return Hypothesis("fa2prod0.5")
        if self == "fL1":
            return Hypothesis("fL1prod0.5")
        if self == "fL1Zg":
            return Hypothesis("fL1Zgprod0.5")
        assert False
    @property
    def categoryname(self):
        if self.fais == ("fa3",): return "0P_or_0M"
        if self == "fa2": return "0P_or_a2"
        if self == "fL1": return "0P_or_L1"
        if self == "fL1Zg": return "0P_or_L1Zg"
        if self.isdecayonly: return "nocategorization"
        if self == "fa3fa2fL1fL1Zg": return "0P_or_0M_or_a2_or_L1_or_L1Zg"
        assert False, self
    @property
    def dimensions(self):
        return len(self.fais)
    @property
    def isdecayonly(self):
        if self.isfL1fL1Zg: return True
        return False
    @property
    def doGEN(self):
        if self == "fa3_multiparameter_nodbkg": return True
        if self == "fa3_only6bins": return True
        if self in ("fa2", "fa3", "fa3_STXS", "fL1", "fL1Zg"): return False
        if self.isfL1fL1Zg: return False
        if self == "fa3fa2fL1fL1Zg": return False
        assert False, self
    @property
    def doCMS(self):
        if self in ("fa2", "fa3", "fL1", "fL1Zg"): return True
        if self.isfL1fL1Zg: return False
        if self in ("fa3_STXS", "fa3_multiparameter_nodbkg", "fa3fa2fL1fL1Zg"): return False
        assert False, self
    @property
    def fais(self):
        if self in ("fa3_STXS", "fa3_multiparameter_nodbkg", "fa3_only6bins"): return "fa3",
        if self in ("fa3", "fa2", "fL1", "fL1Zg"): return self,
        if self.isfL1fL1Zg: return Analysis("fL1"), Analysis("fL1Zg")
        if self == "fa3fa2fL1fL1Zg": return Analysis("fa3"), Analysis("fa2"), Analysis("fL1"), Analysis("fL1Zg")
        assert False, self
    @property
    def isfL1fL1Zg(self):
        return "fL1fL1Zg" in str(self) and self != "fa3fa2fL1fL1Zg"

    @property
    def usehistogramsforcombine(self):
        if self == "fa3_multiparameter_nodbkg": return True
        if self == "fa3_only6bins": return True
        if self == "fa3fa2fL1fL1Zg": return True
        return False

class Production(MyEnum):
    enumname = "production"
    enumitems = (
                 EnumItem("170203"),
                 EnumItem("170222"),
                 EnumItem("170712"),
                 EnumItem("170825"),
                 EnumItem("180121"),
                 EnumItem("180224"),
                 EnumItem("180224_10bins"),
                 EnumItem("180224_newdiscriminants"),
                 EnumItem("180416"),
                 EnumItem("180530", "180531_2016"),
                 EnumItem("180530_Ulascan"),
                 EnumItem("180531", "180531_2017"),
                 EnumItem("180531_Ulascan"),
                 EnumItem("180721", "180721_2016"),
                 EnumItem("180721_Ulascan"),
                 EnumItem("180722", "180721_2017"),
                 EnumItem("180722_Ulascan"),
                 EnumItem("LHE_170509"),
                 EnumItem("GEN_Meng"),
                 EnumItem("GEN_181119"),
                )
    def __cmp__(self, other):
        return cmp(str(self), str(type(self)(other)))
    def CJLSTdir(self):
        if "_" in str(self) and "LHE" not in str(self) and "GEN" not in str(self): return type(self)(str(self).split("_")[0]).CJLSTdir()
        if self == "170203":
            if config.host == "lxplus":
                return "root://lxcms03//data3/Higgs/170203"
            elif config.host == "MARCC":
                return "/work-zfs/lhc/CJLSTtrees/170203"
        if self == "170222":
            if config.host == "lxplus":
                return "root://lxcms03//data3/Higgs/170222"
            elif config.host == "MARCC":
                return "/work-zfs/lhc/CJLSTtrees/170222"
        if self == "170712":
            if config.host == "lxplus":
                return "root://lxcms03//data3/Higgs/170623"
            elif config.host == "MARCC":
                return "/work-zfs/lhc/CJLSTtrees/170623"
        if self == "170825":
            if config.host == "MARCC":
                return "/work-zfs/lhc/CJLSTtrees/170825_Heshy_fL1fL1Zg"
        if self == "180121":
            if config.host == "lxplus":
                return "root://lxcms03//data3/Higgs/180121"
            elif config.host == "MARCC":
                return "/work-zfs/lhc/CJLSTtrees/180121"
        if self in ("180224", "180224_10bins", "180224_newdiscriminants"):
            if config.host == "lxplus":
                return "/eos/user/u/usarica/CJLST/4l/180224"
            elif config.host == "MARCC":
                return "/work-zfs/lhc/CJLSTtrees/180224"
        if self == "180416":
            if config.host == "lxplus":
                return "root://lxcms03//data3/Higgs/180416"
            elif config.host == "MARCC":
                return "/work-zfs/lhc/CJLSTtrees/180416"
        if self in ("180531_2016", "180530_Ulascan"):
            if config.host == "lxplus":
                return "root://lxcms03//data3/Higgs/180531_2016"
            elif config.host == "MARCC":
                return "/work-zfs/lhc/CJLSTtrees/180530"
        if self in ("180531", "180531_Ulascan"):
            if config.host == "lxplus":
                return "root://lxcms03//data3/Higgs/180531"
            elif config.host == "MARCC":
                return "/work-zfs/lhc/CJLSTtrees/180531"
        if self in ("180721_2016", "180721_Ulascan"):
            if config.host == "lxplus":
                return "root://lxcms03//data3/Higgs/180721_2016"
            elif config.host == "MARCC":
                return "/work-zfs/lhc/CJLSTtrees/180721"
        if self in ("180721_2017", "180722_Ulascan"):
            if config.host == "lxplus":
                return "root://lxcms03//data3/Higgs/180721_2017"
            elif config.host == "MARCC":
                return "/work-zfs/lhc/CJLSTtrees/180722"
        if self == "GEN_Meng":
            if config.host == "lxplus":
                return "/eos/cms/store/user/xiaomeng/writeupTrees/"
            elif config.host == "MARCC":
                assert False
        if self == "GEN_181119":
            if config.host == "lxplus":
                assert False
            elif config.host == "MARCC":
                return "/work-zfs/lhc/GENtrees/181119_2017MC"
        assert False
    def CJLSTdir_anomalous(self):
        return self.CJLSTdir()
    def CJLSTdir_data(self):
        if self == "170712":
            if config.host == "lxplus":
                return "root://lxcms03//data3/Higgs/170712_Data2017"
            elif config.host == "MARCC":
                return "/work-zfs/lhc/CJLSTtrees/170712_Data2017"
        return self.CJLSTdir()
    def CJLSTdir_anomalous_VBF(self):
        return self.CJLSTdir()
    def CJLSTdir_anomalous_VH(self):
        return self.CJLSTdir()
    def CJLSTdir_MINLO(self):
        return self.CJLSTdir()
    @property
    def dataluminosity(self):
        if self in ("170203", "170222", "170825", "180121", "180224", "180531_2016", "180721_2016"): return 35.921875594646
        if self in ("180416", "180531", "180721_2017"): return 41.529343499127
        if self.LHE or self.GEN: return 300

        if self in ("180224_10bins", "180224_newdiscriminants") or "_Ulascan" in str(self):
            return type(self)(str(self).split("_")[0]).dataluminosity

        assert False
    def __int__(self):
        return int(str(self))
    @property
    def year(self):
        if "_" in str(self) and "LHE" not in str(self) and "GEN" not in str(self): return type(self)(str(self).split("_")[0]).year
        if self <= "180224" or self == "180531_2016" or self == "180721_2016":
            return 2016
        if self == "180416" or self == "180531" or self == "180721_2017" or self == "GEN_Meng" or self == "GEN_181119":
            return 2017
        assert False
    @property
    def productionforsmoothingparameters(self):
        if self == "170222": return type(self)("170203")
        if self == "170825": return type(self)("170203")
        if self == "180121": return type(self)("170203")
        if self == "180224": return type(self)("170203")

        if self == "180224_10bins": return type(self)("180416")
        if self == "180224_newdiscriminants": return type(self)("180416")

        if self == "180531": return type(self)("180416")
        if self == "180531_2016": return type(self)("180416")  #10 bins

        if self == "180721_2016": return type(self)("180416")
        if self == "180721_2017": return type(self)("180416")
        return self
    @property
    def productionforrate(self):
        if self in ("180224_10bins", "180224_newdiscriminants") or "_Ulascan" in str(self):
            return type(self)(str(self).split("_")[0])
        return self
    @property
    def LHE(self):
        return "LHE" in str(self)
    @property
    def GEN(self):
        return "GEN" in str(self)
    @property
    def pdf(self):
      if self.year == 2016: return "NNPDF30_lo_as_0130"
      if self.year == 2017: return "NNPDF31_lo_as_0130"
      assert False, self

class Category(MyEnum):
    """
    For now just 3 categories, VBF2j, VH hadronic, and dump everything else into untagged
    """
    enumname = "category"
    enumitems = (
                 EnumItem("Untagged", "UntaggedMor17", "VBF1jTaggedMor17", "VHLeptTaggedMor17", "ttHTaggedMor17", "VHMETTaggedMor17", "UntaggedMor18", "VBF1jTaggedMor18", "VHLeptTaggedMor18", "ttHHadrTaggedMor18", "ttHLeptTaggedMor18", "VHMETTaggedMor18", ""),
                 EnumItem("VHHadrtagged", "VHHadrTaggedMor17", "VHHadrTaggedMor18"),
                 EnumItem("VBFtagged", "VBF2jTaggedMor17", "VBF2jTaggedMor18"),
                )
    @property
    def idnumbers(self):
        """
        returns a list of ints corresponding to the C++ enums corresponding to this category
        (defined in Category.h)
        """
        import CJLSTscripts
        return {getattr(CJLSTscripts, name) for name in self.item.names if "Mor" in name}

    @property
    def yamlname(self):
        return str(self).replace("tagged", "Tagged")

    def __contains__(self, other):
        self.checkidnumbers()
        return other in self.idnumbers

    @property
    def mainsignals(self):
        if self == "Untagged": return [ProductionMode("ggH")]
        if self == "VBFtagged": return [ProductionMode("VBF")]
        if self == "VHHadrtagged": return [ProductionMode("ZH"), ProductionMode("WH")]

    @property
    def nameonplot(self):
        if self == "Untagged": return "Untagged"
        if self == "VBFtagged": return "VBF-tagged"
        if self == "VHHadrtagged": return "VH-tagged"

    @classmethod
    def fromid(cls, number):
        for category in cls.items():
            if number in category:
                return category
        raise ValueError("Invalid id {}".format(number))

    @classmethod
    @cache
    def checkidnumbers(cls):
        for cat1, cat2 in combinations(cls.enumitems, 2):
            assert not cls(cat1).idnumbers & cls(cat2).idnumbers

class AlternateGenerator(MyEnum):
    enumname = "alternategenerator"
    enumitems = (
                 EnumItem("POWHEG"),
                 EnumItem("MINLO"),
                 EnumItem("NNLOPS"),
                )

class Extension(MyEnum):
    enumname = "extension"
    enumitems = (
                 EnumItem("ext"),
                )

class PythiaSystematic(MyEnum):
    enumname = "pythiasystematic"
    enumitems = (
                 EnumItem("ScaleUp"),
                 EnumItem("ScaleDn", "ScaleDown"),
                 EnumItem("TuneUp"),
                 EnumItem("TuneDn", "TuneDown"),
                )

    @property
    def appendname(self):
        return "_" + str(self).lower().replace("dn", "down")
    def hassample(self, year):
        if self in ("TuneUp", "TuneDn") and year in (2016, 2017): return True
        if self in ("ScaleUp", "ScaleDn") and year == 2016: return True
        if self in ("ScaleUp", "ScaleDn") and year == 2017: return False
        assert False, (self, year)

class AlternateWeight(MyEnum):
    enumname = "alternateweight"
    enumitems = (
                 EnumItem("1"),
                 EnumItem("muRUp"),
                 EnumItem("muRDn", "muRDown"),
                 EnumItem("muFUp"),
                 EnumItem("muFDn", "muFDown"),
                 EnumItem("PDFUp"),
                 EnumItem("PDFDn", "PDFDown"),
                 EnumItem("alphaSUp"),
                 EnumItem("alphaSDn", "alphaSDown"),
                 EnumItem("EWcorrUp"),
                 EnumItem("EWcorrDn", "EWcorrDown"),
                 EnumItem("PythiaScaleUp"),
                 EnumItem("PythiaScaleDn", "PythiaScaleDown"),
                )
    @property
    def issystematic(self): return self != "1"
    @property
    def weightname(self):
      if self == "1": return "1"
      if self == "muRUp": return "LHEweight_QCDscale_muR2_muF1"
      if self == "muRDn": return "LHEweight_QCDscale_muR0p5_muF1"
      if self == "muFUp": return "LHEweight_QCDscale_muR1_muF2"
      if self == "muFDn": return "LHEweight_QCDscale_muR1_muF0p5"
      if self == "PDFUp": return "LHEweight_PDFVariation_Up"
      if self == "PDFDn": return "LHEweight_PDFVariation_Dn"
      if self == "alphaSUp": return "LHEweight_AsMZ_Up"
      if self == "alphaSDn": return "LHEweight_AsMZ_Dn"
      if self == "EWcorrUp": return "(1 + KFactor_EW_qqZZ_unc/KFactor_EW_qqZZ)"
      if self == "EWcorrDn": return "(1 - KFactor_EW_qqZZ_unc/KFactor_EW_qqZZ)"
      if self == "PythiaScaleUp": return "(PythiaWeight_isr_muR4 * PythiaWeight_fsr_muR4)"
      if self == "PythiaScaleDn": return "(PythiaWeight_isr_muR0p25 * PythiaWeight_fsr_muR0p25)"
      assert False
    @property
    def kfactorname(self):
      if self in ("1", "PythiaScaleUp", "PythiaScaleDn"): return "KFactor_QCD_ggZZ_Nominal"
      if self == "muRUp": return "KFactor_QCD_ggZZ_QCDScaleUp"
      if self == "muRDn": return "KFactor_QCD_ggZZ_QCDScaleDn"
      if self == "muFUp": return "KFactor_QCD_ggZZ_PDFScaleUp"
      if self == "muFDn": return "KFactor_QCD_ggZZ_PDFScaleDn"
      if self == "PDFUp": return "KFactor_QCD_ggZZ_PDFReplicaUp"
      if self == "PDFDn": return "KFactor_QCD_ggZZ_PDFReplicaDn"
      if self == "alphaSUp": return "KFactor_QCD_ggZZ_AsUp"
      if self == "alphaSDn": return "KFactor_QCD_ggZZ_AsDn"
      assert False

channels = Channel.items()
JECsystematics = JECSystematic.items()
btagsystematics = BTagSystematic.items()
pythiasystematics = PythiaSystematic.items()
flavors = Flavor.items()
hypotheses = Hypothesis.items()
decayonlyhypotheses = Hypothesis.items(lambda x: all("prod" not in name and "0.33" not in name and "0.25" not in name for name in x.names))
prodonlyhypotheses = Hypothesis.items(lambda x: all("dec" not in name for name in x.names))
proddechypotheses = Hypothesis.items(lambda x: True)
purehypotheses = Hypothesis.items(lambda x: x.ispure)
hffhypotheses = HffHypothesis.items()
productionmodes = ProductionMode.items()
analyses = Analysis.items(lambda x: x.doGEN)
config.productionsforcombine = type(config.productionsforcombine)(Production(production) for production in config.productionsforcombine)
if len(config.productionsforcombine) == 1:
    config.productionforcombine = Production(config.productionforcombine)
productions = Production.items(lambda x: x in config.productionsforcombine)
categories = Category.items()

_ = [""]
if config.applym4lshapesystematicsUntagged or config.applym4lshapesystematicsVBFVHtagged:
    _ += ["ResUp", "ResDown", "ScaleUp", "ScaleDown"]
    if config.combinem4lshapesystematics:
        _ += ["ScaleResUp", "ScaleResDown"]
if config.applyZXshapesystematicsUntagged or config.applyZXshapesystematicsVBFVHtagged:
    _ += ["ZXUp", "ZXDown"]
if config.applyJECshapesystematics:
    _ += ["JECUp", "JECDown"]
if config.applyMINLOsystematics:
    _ += ["MINLO_SM", "MINLOUp", "MINLODn"]
shapesystematics = ShapeSystematic.items(lambda x: x in _)
treeshapesystematics = ShapeSystematic.items(lambda x: x in _ and x in ("", "ResUp", "ResDown", "ScaleUp", "ScaleDown", "JECUp", "JECDown", "MINLO_SM") + ("ZXUp", "ZXDown")*config.usenewZXsystematics)
del _

class MetaclassForMultiEnums(type):
    def __new__(cls, clsname, bases, dct):
        try:
            enums = dct["enums"]
        except KeyError:
            for base in bases:
                try:
                    enums = base.enums
                    break
                except AttributeError:
                    pass
            else:
                raise TypeError("MultiEnum class {} does not define enums!".format(clsname))
            dct["enums"] = enums
        dct["needenums"] = needenums = list(enums)
        dct["subenums"] = subenums = OrderedDict()
        while True:
            for enum in needenums[:]:
                if issubclass(enum, MultiEnum):
                    #replace it with its subenums
                    needenums.remove(enum)
                    needenums += [a for a in enum.enums if a not in needenums]
                    #break out of for in order to continue in the while, to allow for more deeply nested MultiEnums
                    subenums[enum] = enum.enums
                    break
                elif issubclass(enum, MyEnum):
                    if enum not in needenums:
                        needenums.append(enum)
                else:
                    raise TypeError("{} is not a MyEnum or a MultiEnum! (in class {})".format(enum, clsname))
            else:
                break
        for _ in "enums", "needenums":
            dct[_] = tuple(dct[_])
        return super(MetaclassForMultiEnums, cls).__new__(cls, clsname, bases, dct)

class MultiEnumABCMeta(MetaclassForMultiEnums, abc.ABCMeta):
    """
    needed to resolve conflict
    http://code.activestate.com/recipes/204197-solving-the-metaclass-conflict/
    except don't need all their fancy stuff
    the only function in MetaClassForMultiEnums is __new__ and it calls super
    """

class MultiEnum(object):
    __metaclass__ = MetaclassForMultiEnums
    enums = []
    def __init__(self, *args):
        for enum in self.enums:
            setattr(self, enum.enumname, None)
        enumsdict = {enum: None for enum in self.needenums}

        argscopy = list(args)
        while True:
            for i, arg in enumerate(argscopy[:]):
                if isinstance(arg, MultiEnum):
                    argscopy[i:i+1] = [getattr(arg, enum.enumname) for enum in arg.enums]
                    break
            else:
                break
        while None in argscopy:
            argscopy.remove(None)

        for enum in self.needenums:
            for i, arg in enumerate(argscopy[:]):
                try:
                    tmp = enum(arg)
                except ValueError:
                    pass
                else:
                    if enumsdict[enum] is not None: raise TypeError("Multiple arguments provided that fit {} ({} and {})".format(enum, enumsdict[enum], arg))
                    enumsdict[enum] = tmp
                    del argscopy[i]
                    break
        if argscopy:
            raise ValueError("Extra arguments {}\n{}".format(argscopy, args))

        self.applysynonyms(enumsdict)

        for enum in self.enums:
            if issubclass(enum, MyEnum):
                setattr(self, enum.enumname, enumsdict[enum])
            if issubclass(enum, MultiEnum):
                setattr(self, enum.enumname, enum(*(v for k, v in enumsdict.iteritems() if k in enum.needenums and v is not None)))

        for subenum, subsubenums in self.subenums.iteritems():
            for subsubenum in subsubenums:
                setattr(self, subsubenum.enumname, getattr(getattr(self, subenum.enumname), subsubenum.enumname))

        self.check(*args)
        for enum in self.enums:
            if isinstance(enum, MultiEnum):
                enum.check(*args)
        self.items = tuple(getattr(self, enum.enumname) for enum in self.enums)

    def __eq__(self, other):
        if type(self) != type(other): raise ValueError("Comparing two different types of MultiEnums!")
        return self.items == other.items
    def __ne__(self, other):
        return not self == other
    def __hash__(self):
        return hash(self.items)

    def __str__(self):
        return " ".join(str(item) for item in self.items if item is not None)
    def __repr__(self):
        return "{}({})".format(type(self).__name__, ", ".join(repr(_.item.name if isinstance(_, MyEnum) else _) for _ in self.items if _ is not None))

    def applysynonyms(self, enumsdict):
        pass

    def check(self, *args, **kwargs):
        dontcheck = ()
        for kw, kwarg in kwargs.iteritems():
            if kw == "dontcheck":
                dontcheck = tuple(kwarg)
        for enum in self.enums:
            if enum not in dontcheck and getattr(self, enum.enumname) is None:
                raise ValueError("No option provided for {}\n{}".format(enum.enumname, args))
