import abc
from collections import OrderedDict
import os

import config
import constants
from utilities import generatortolist_condition, tfiles

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
    def moriondcardfile(self):
        return os.path.join(config.repositorydir, "helperstuff/moriondcards/hzz4lcard_{}_0.txt".format(self))

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

class Hypothesis(MyEnum):
    enumname = "hypothesis"
    enumitems = (
                 EnumItem("0+", "SM", "scalar"),
                 EnumItem("0+_photoncut", "SM_photoncut"),
                 EnumItem("a2", "0h+"),
                 EnumItem("0-", "a3", "PS", "pseudoscalar"),
                 EnumItem("L1", "Lambda1"),
                 EnumItem("L1Zg"),
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
    @property
    def ispure(self):
        return self in ("0+", "0+_photoncut", "0-", "0h+", "L1", "L1Zg")
    @property
    def couplingname(self):
        if self in ("0+", "0+_photoncut"): return "g1"
        if self == "a2": return "g2"
        if self == "0-": return "g4"
        if self == "L1": return "g1prime2"
        if self == "L1Zg": return "ghzgs1prime2"
        assert False
    @property
    def photoncut(self):
        if self in ("0+_photoncut", "L1Zg"): return True
        if self in ("0+", "0-", "a2", "L1"): return False
        for b in "prod", "dec", "proddec":
            for c in "+-":
                for a in "fa2", "fa3", "fL1":
                    if self == "{}{}{}0.5".format(a, b, c): return False
                if self == "{}{}{}0.5".format("fL1Zg", b, c): return True
        if self == "fa2dec-0.9": return False
        assert False

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
                )
    @property
    def combinename(self):
        for name in "ggH", "qqH", "ZH", "WH", "ttH", "bkg_qqzz", "bkg_ggzz", "bkg_vbf", "bkg_zjets":
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
    def yamlsystname(self):
        if self in ("VBF", "VBF bkg"):
            return "qqH"
        elif self == "ZX":
            return "zjets"
        return str(self)
    @property
    def isbkg(self):
        if self in ("ggH", "VBF", "ZH", "WH", "ttH", "HJJ", "WplusH", "WminusH"):
            return False
        elif self in ("ggZZ", "qqZZ", "VBF bkg", "ZX"):
            return True
        assert False
    @property
    def issignal(self):
        return not self.isbkg
    @property
    def validhypotheses(self):
        if self in ("ggH", "ttH", "HJJ"):
            return Hypothesis.items(lambda x: x in decayonlyhypotheses)
        if self == "VBF":
            return Hypothesis.items(lambda x: x in proddechypotheses)
        if self == "ZH":
            return Hypothesis.items(lambda x: x in proddechypotheses)
        if self == "WH":
            return Hypothesis.items(lambda x: x in proddechypotheses)
        if self in ("WplusH", "WminusH"):
            return Hypothesis.items(lambda x: x == "0+")
        assert False
    @property
    def generatedhypotheses(self):
        if self == "ggH":
            return Hypothesis.items(lambda x: x in ("0+", "0-", "a2", "L1", "fa30.5", "fa20.5", "fL10.5"))
        if self in ("VBF", "ZH", "WH"):
            return Hypothesis.items(lambda x: x in ("0+", "0-", "a2", "L1", "fa3prod0.5", "fa2prod0.5", "fL1prod0.5"))
        if self in ("WplusH", "WminusH", "ttH", "HJJ"):
            return Hypothesis.items(lambda x: x == "0+")
        assert False
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
            for h in self.generatedhypotheses:
                yield Sample(self, h, hff, production)

    @property
    def QCDsystematicname(self):
      if self == "ggH": return "QCDscale_ggH_cat"
      if self == "qqH": return "QCDscale_qqH_cat"
      if self in ("ZH", "WH"): return "QCDscale_VH_cat"
      if self == "ttH": return "QCDscale_ttH_cat"
      if self == "qqZZ": return "QCDscale_VV_cat"
      return None

    @property
    def pdfsystematicname(self):
      if self == "ggH": return "pdf_Higgs_gg_cat"
      if self in ("qqH", "ZH", "WH"): return "pdf_Higgs_qq_cat"
      if self == "ttH": return "pdf_Higgs_ttH_cat"
      if self == "qqZZ": return "pdf_qq_cat"
      return None

    def workspaceshapesystematics(self, category):
      result = []
      if self in ("ggH", "qqH", "ZH", "WH", "ttH"):
        if config.applym4lshapesystematicsUntagged and category == "Untagged" or config.applym4lshapesystematicsVBFVHtagged and category != "Untagged":
          if config.combinem4lshapesystematics:
            result += ["ScaleRes"]
          else:
            result += ["Scale", "Res"]
      if self == "ggH" and category in ("VBFtagged", "VHHadrtagged") and config.applyMINLOsystematics:
        result += ["MINLO"]
      return [WorkspaceShapeSystematic(_) for _ in result]

    @property
    def alternateweights(self):
      if self in ("ggH", "qqH", "ZH", "WH", "ttH"):
         return AlternateWeight.items(lambda x: x!="EWcorrUp" and x!="EWcorrDn")
      if self == "qqZZ":
         return AlternateWeight.items()
      assert False

class WorkspaceShapeSystematic(MyEnum):
    enumname = "workspaceshapesystematic"
    enumitems = (
                 EnumItem("Res"),
                 EnumItem("Scale"),
                 EnumItem("ScaleRes"),
                 EnumItem("MINLO"),
                )
    @property
    def isperchannel(self):
        if self in ("Res", "Scale", "ScaleRes"): return True
        return False

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
                 EnumItem("JECDown", "JECDn"),
                 EnumItem("ZXUp"),
                 EnumItem("ZXDown", "ZXDn"),
                 EnumItem("MINLO_SM"),
                 EnumItem("MINLOUp"),
                 EnumItem("MINLODn", "MINLODown"),
                )
    def appendname(self):
        if self in ("ScaleUp", "ScaleDown", "ResUp", "ResDown"): return "_" + str(self)
        return ""
    def D_bkg(self, title=False):
        from discriminants import discriminant
        return discriminant("D_bkg"+self.appendname())
    def appliesto(self, templategroup):
        if self == "":
            return True
        if self in ("ResUp", "ResDown", "ScaleUp", "ScaleDown", "ScaleResUp", "ScaleResDown"):
            return templategroup in ("ggh", "vbf", "zh", "wh", "tth")
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
                )
    def title(self, latex=False, superscript=None):
        if self == "fa3":
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
        if self == "fa3":
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
        if self == "fa3": return "g4"
        if self == "fa2": return "g2"
        if self == "fL1": return "g1prime2"
        if self == "fL1Zg": return "ghzgs1prime2"
    @property
    def purehypotheses(self):
        if self == "fa3":
            return Hypothesis("0+"), Hypothesis("0-")
        if self == "fa2":
            return Hypothesis("0+"), Hypothesis("a2")
        if self == "fL1":
            return Hypothesis("0+"), Hypothesis("L1")
        if self == "fL1Zg":
            return Hypothesis("0+_photoncut"), Hypothesis("L1Zg")
    @property
    def mixdecayhypothesis(self):
        if self == "fa3":
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
        if self == "fa3":
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
        if self == "fa3": return "0P_or_0M"
        if self == "fa2": return "0P_or_a2"
        if self == "fL1": return "0P_or_L1"
        if self == "fL1Zg": return "0P_or_L1Zg"
        assert False
    @property
    def photoncut(self):
        if self == "fL1Zg": return True
        if self in ("fa2", "fa3", "fL1"): return False
        assert False

class Production(MyEnum):
    enumname = "production"
    enumitems = (
                 EnumItem("170203"),
                 EnumItem("170222"),
                )
    def __cmp__(self, other):
        return cmp(str(self), str(type(self)(other)))
    def CJLSTdir(self):
        if self == "170203":
            if config.host == "lxplus":
                return "root://lxcms03//data3/Higgs/170203"
        if self == "170222":
            if config.host == "lxplus":
                return "root://lxcms03//data3/Higgs/170222"
        assert False
    def CJLSTdir_anomalous(self):
        return self.CJLSTdir()
    def CJLSTdir_data(self):
        return self.CJLSTdir()
    def CJLSTdir_anomalous_VBF(self):
        return self.CJLSTdir()
    def CJLSTdir_anomalous_VH(self):
        return self.CJLSTdir()
    def CJLSTdir_MINLO(self):
        return self.CJLSTdir()
    @property
    def dataluminosity(self):
        if self in ("170203", "170222"): return 35.8671
        assert False
    def __int__(self):
        return int(str(self))
    @property
    def year(self):
        if "170222" <= self:
            return 2016
        assert False
    @property
    def productionforsmoothingparameters(self):
        if self == "170222": return type(self)("170203")
        return self

class Category(MyEnum):
    """
    For now just 3 categories, VBF2j, VH hadronic, and dump everything else into untagged
    """
    enumname = "category"
    enumitems = (
                 EnumItem("Untagged", "UntaggedMor17", "VBF1jTaggedMor17", "VHLeptTaggedMor17", "ttHTaggedMor17", "VHMETTaggedMor17"),
                 EnumItem("VHHadrtagged", "VHHadrTaggedMor17"),
                 EnumItem("VBFtagged", "VBF2jTaggedMor17"),
                )
    @property
    def idnumbers(self):
        """
        returns a list of ints corresponding to the C++ enums corresponding to this category
        (defined in Category.h)
        """
        import CJLSTscripts
        return [getattr(CJLSTscripts, name) for name in self.item.names if "Mor17" in name]

    @property
    def yamlname(self):
        return str(self).replace("tagged", "Tagged")

    def __contains__(self, other):
        return other in self.idnumbers

    @property
    def mainsignals(self):
        if self == "Untagged": return [ProductionMode("ggH")]
        if self == "VBFtagged": return [ProductionMode("VBF")]
        if self == "VHHadrtagged": return [ProductionMode("ZH"), ProductionMode("WH")]

    @classmethod
    def fromid(cls, number):
        for category in cls.items():
            if number in category:
                return category
        raise ValueError("Invalid id {}".format(number))

class AlternateGenerator(MyEnum):
    enumname = "alternategenerator"
    enumitems = (
                 EnumItem("POWHEG"),
                 EnumItem("MINLO"),
                 EnumItem("NNLOPS"),
                )

class PythiaSystematic(MyEnum):
    enumname = "pythiasystematic"
    enumitems = (
                 EnumItem("ScaleUp"),
                 EnumItem("ScaleDn"),
                 EnumItem("TuneUp"),
                 EnumItem("TuneDn"),
                )

    @property
    def appendname(self):
        return "_" + str(self).lower().replace("dn", "down")

class AlternateWeight(MyEnum):
    enumname = "alternateweight"
    enumitems = (
                 EnumItem("1"),
                 EnumItem("muRUp"),
                 EnumItem("muRDn"),
                 EnumItem("muFUp"),
                 EnumItem("muFDn"),
                 EnumItem("PDFUp"),
                 EnumItem("PDFDn"),
                 EnumItem("alphaSUp"),
                 EnumItem("alphaSDn"),
                 EnumItem("EWcorrUp"),
                 EnumItem("EWcorrDn"),
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
      assert False

channels = Channel.items()
JECsystematics = JECSystematic.items()
btagsystematics = BTagSystematic.items()
pythiasystematics = PythiaSystematic.items()
flavors = Flavor.items()
hypotheses = Hypothesis.items()
decayonlyhypotheses = Hypothesis.items(lambda x: x in ("0+", "0+_photoncut", "a2", "0-", "L1", "L1Zg", "fa2dec0.5", "fa3dec0.5", "fL1dec0.5", "fL1Zgdec0.5", "fa2dec-0.5", "fa3dec-0.5", "fL1dec-0.5", "fL1Zgdec-0.5", "fa2dec-0.9"))
prodonlyhypotheses = Hypothesis.items(lambda x: x in ("0+", "0+_photoncut", "a2", "0-", "L1", "L1Zg", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5", "fL1Zgprod0.5", "fa2prod-0.5", "fa3prod-0.5", "fL1prod-0.5", "fL1Zgprod-0.5"))
proddechypotheses = Hypothesis.items()
purehypotheses = Hypothesis.items(lambda x: x.ispure)
hffhypotheses = HffHypothesis.items()
productionmodes = ProductionMode.items()
analyses = Analysis.items()
config.productionsforcombine = type(config.productionsforcombine)(Production(production) for production in config.productionsforcombine)
productions = Production.items(lambda x: x in ("170203", "170222"))
#productions = Production.items(lambda x: x in config.productionsforcombine)
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
treeshapesystematics = ShapeSystematic.items(lambda x: x in _ and x in ("", "ResUp", "ResDown", "ScaleUp", "ScaleDown", "JECUp", "JECDown"))
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
