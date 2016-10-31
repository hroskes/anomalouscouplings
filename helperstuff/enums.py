import abc
from collections import OrderedDict
import config
import constants
import os

class EnumItem(object):
    def __init__(self, name, *other):
        self.name = name
        self.names = tuple([name] + list(other))

    def __str__(self):
        return self.name
    def __hash__(self):
        return hash(self.names)

    def __eq__(self, other):
        if type(other) == int or type(other) == str:
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

class Hypothesis(MyEnum):
    enumname = "hypothesis"
    enumitems = (
                 EnumItem("0+", "SM", "scalar"),
                 EnumItem("a2", "0h+"),
                 EnumItem("0-", "a3", "PS", "pseudoscalar"),
                 EnumItem("L1", "Lambda1"),
                 EnumItem("fa20.5", "fa2dec0.5"),
                 EnumItem("fa30.5", "fa3dec0.5"),
                 EnumItem("fL10.5", "fL1dec0.5"),
                 EnumItem("fa2prod0.5"),
                 EnumItem("fa3prod0.5"),
                 EnumItem("fL1prod0.5"),
                 EnumItem("fa2proddec-0.5"),
                 EnumItem("fa3proddec-0.5"),
                 EnumItem("fL1proddec-0.5"),
                )
    @property
    def ispure(self):
        return self in ("0+", "0-", "0h+", "L1")

class ProductionMode(MyEnum):
    enumname = "productionmode"
    enumitems = (
                 EnumItem("ggH"),
                 EnumItem("VBF", "qqH"),
                 EnumItem("H+jj", "HJJ"),
                 EnumItem("ZH"),
                 EnumItem("WH"),
                 EnumItem("ttH"),
                 EnumItem("qqZZ", "bkg_qqzz"),
                 EnumItem("ggZZ", "bkg_ggzz"),
                 EnumItem("VBFbkg", "VBF bkg", "bkg_vbf"),
                 EnumItem("ZX", "bkg_zjets"),
                 EnumItem("data"),
                )
    @property
    def combinename(self):
        import combinehelpers
        for name in combinehelpers.datacardprocessline.split():
            if self == name:
                return name
        assert False
    @property
    def yamlratenames(self):
        if self == "ggH":
            return ["ggH", "ttH"]
        elif self == "VBF":
            return ["qqH"]
        elif self == "ZX":
            return ["zjets"]
        return [str(self)]
    @property
    def yamlsystematicsname(self):
        if self == "ggH":
            return "ggH"
        elif self == "VBF":
            return "qqH"
        elif self == "ZX":
            return "zjets"
        elif self == "VBF bkg":
            return "qqZZ"
        return str(self)

class Systematic(MyEnum):
    enumname = "systematic"
    enumitems = (
                 EnumItem(""),
                 EnumItem("ResUp"),
                 EnumItem("ResDown"),
                 EnumItem("ScaleUp"),
                 EnumItem("ScaleDown"),
                 EnumItem("ScaleResUp", "ResScaleUp"),
                 EnumItem("ScaleResDown", "ResScaleDown"),
                 EnumItem("ZXUp"),
                 EnumItem("ZXDown"),
                )
    def appendname(self):
        if self == "": return ""
        return "_" + str(self)
    def D_bkg_0plus(self, title=False):
        from discriminants import discriminant
        return discriminant("D_bkg_0plus"+self.appendname())
    def appliesto(self, templategroup):
        if templategroup in ("ggh", "vbf", "zh", "wh"):
            return self in ("", "ResUp", "ResDown", "ScaleUp", "ScaleDown", "ScaleResUp", "ScaleResDown")
        elif templategroup == "bkg":
            return self in ("", "ZXUp", "ZXDown")
        elif templategroup == "DATA":
            return self in ("", )
        assert False


class TemplateGroup(MyEnum):
    enumname = "templategroup"
    enumitems = (
                 EnumItem("ggh"),
                 EnumItem("vbf"),
                 EnumItem("zh"),
                 EnumItem("wh"),
                 EnumItem("background", "bkg"),
                 EnumItem("DATA"),
                )

class Analysis(MyEnum):
    enumname = "analysis"
    enumitems = (
                 EnumItem("fa3"),
                 EnumItem("fa2"),
                 EnumItem("fL1"),
                )
    @property
    def title(self):
        if self == "fa3":
            return "f_{a3}"
        if self == "fa2":
            return "f_{a2}"
        if self == "fL1":
            return "f_{#Lambda1}"
        assert False
    @property
    def phi(self):
        if self == "fa3":
            return "#phi_{a3}"
        if self == "fa2":
            return "#phi_{a2}"
        if self == "fL1":
            return "#phi_{#Lambda1}"
        assert False
    @property
    def phi_lower(self):
        return self.phi.replace("{", "{#lower[-0.25]{").replace("}", "}}")
    def purediscriminant(self, title=False):
        if not title:
            if self == "fa3":
                return "D_0minus_decay"
            if self == "fa2":
                return "D_g2_decay"
            if self == "fL1":
                return "D_g1prime2_decay"
        else:
            if self == "fa3":
                return "D_{0-}"
            if self == "fa2":
                return "D_{0h+}"
            if self == "fL1":
                return "D_{#Lambda1}"
        assert False
    def mixdiscriminant(self, title=False):
        if not title:
            if self == "fa3":
                return "D_CP_decay"
            if self == "fa2":
                return "D_g1g2_decay"
            if self == "fL1":
                return "D_g2_decay"
        else:
            if self == "fa3":
                return "D_{CP}"
            if self == "fa2":
                return "D_{int}"
            if self == "fL1":
                return "D_{0h+}"
        assert False
    def interfxsec(self):
        if self == "fa3":
            return constants.JHUXSggH2L2la1a3
        elif self == "fa2":
            return constants.JHUXSggH2L2la1a2
        elif self == "fL1":
            return constants.JHUXSggH2L2la1L1
        assert False
    @property
    def couplingname(self):
        if self == "fa3": return "g4"
        if self == "fa2": return "g2"
        if self == "fL1": return "g1prime2"
    @property
    def purehypotheses(self):
        if self == "fa3":
            return "0+", "0-"
        if self == "fa2":
            return "0+", "a2"
        if self == "fL1":
            return "0+", "L1"
    @property
    def mixdecayhypothesis(self):
        if self == "fa3":
            return "fa3dec0.5"
        if self == "fa2":
            return "fa2dec0.5"
        if self == "fL1":
            return "fL1dec0.5"

class Production(MyEnum):
    enumname = "production"
    enumitems = (
                 EnumItem("160225"),
                 EnumItem("160624"),
                 EnumItem("160714"),
                 EnumItem("160720"),
                 EnumItem("160725", "160726"),
                 EnumItem("160729"),
                 EnumItem("160901"),
                 EnumItem("160909"),
                 EnumItem("160919"),
                 EnumItem("160928"),
                )
    def __cmp__(self, other):
        return cmp(str(self), str(type(self)(other)))
    def CJLSTdir(self):
        if self == "160225":
            return "root://lxcms03//data3/Higgs/160225/"
        if self == "160624":
            return "root://lxcms03//data3/Higgs/160624/"
        if self == "160714":
            return "root://lxcms03//data3/Higgs/160714/"
        if self == "160720":
            return "root://lxcms03//data3/Higgs/160720/"
        if self == "160725" or self == "160729":
            return "root://lxcms03//data3/Higgs/160726/"
        if self == "160901":
            return "/afs/cern.ch/work/h/hroskes/reweighting_CJLST/CMSSW_8_0_8/src/ZZAnalysis/AnalysisStep/test/prod/VBFanomalous_prodMEs_fix/AAAOK"
        if self == "160909":
            return "/afs/cern.ch/work/h/hroskes/reweighting_CJLST/CMSSW_8_0_8/src/ZZAnalysis/AnalysisStep/test/prod/VHanomalous/AAAOK"
        if self == "160919":
            return "/afs/cern.ch/work/h/hroskes/reweighting_CJLST/CMSSW_8_0_8/src/ZZAnalysis/AnalysisStep/test/prod/ggHVBFVHanomalous_bkp/AAAOK"
        if self == "160928":
            return "/afs/cern.ch/work/h/hroskes/reweighting_CJLST/CMSSW_8_0_8/src/ZZAnalysis/AnalysisStep/test/prod/ggHVBFVHanomalous/AAAOK"
        assert False
    def CJLSTdir_anomalous(self):
        if self < "160624":
            return type(self)("160624").CJLSTdir_anomalous()
        if self == "160714":
            return "/afs/cern.ch/work/h/hroskes/public/CJLST/CMSSW_8_0_8/src/ZZAnalysis/AnalysisStep/test/prod/anomalous/PT13TeV"
        if self == "160720":
            return "root://lxcms03//data3/Higgs/160718/"
        return self.CJLSTdir()
    def CJLSTdir_data(self):
        if self == "160714":
            return "root://lxcms03//data3/Higgs/160716/"
        if self == "160725":
            return "root://lxcms03//data3/Higgs/160725/"
        if self == "160729":
            return "root://lxcms03//data3/Higgs/160729_complete/"
        return self.CJLSTdir()
    def CJLSTdir_anomalous_VBF(self):
        if self < "160901":
            return type(self)("160901").CJLSTdir()
        return self.CJLSTdir()
    def CJLSTdir_anomalous_VH(self):
        if self < "160909":
            return type(self)("160909").CJLSTdir()
        return self.CJLSTdir()
    @property
    def useMELAv2(self):
        if self in ("160225", "160624"):
            return False
        return True
    @property
    def dataluminosity(self):
        if self == "160225": return 2.8
        if self == "160714": return 7.65
        if self == "160720": return 9.2
        if "160725" <= self: return 12.9
        assert False
    def __int__(self):
        return int(str(self))
    @property
    def year(self):
        if self == "160225":
            return 2015
        if "160624" <= self:
            return 2016
        assert False

class BlindStatus(MyEnum):
    enumname = "blindstatus"
    enumitems = (
                 EnumItem("unblind"),
                 EnumItem("blind"),
                )

class Category(MyEnum):
    """
    For now just 3 categories, VBF2j, VH hadronic, and dump everything else into untagged
    """
    enumname = "category"
    enumitems = (
                 EnumItem("Untagged", "UntaggedIchep16", "VBF1jTaggedIchep16", "VHLeptTaggedIchep16", "ttHTaggedIchep16"),
                 EnumItem("VHHadrtagged", "VHHadrTaggedIchep16"),
                 EnumItem("VBFtagged", "VBF2jTaggedIchep16"),
                )
    @property
    def idnumbers(self):
        """
        returns a list of ints corresponding to the C++ enums corresponding to this category
        (defined in Category.h)
        """
        import CJLSTscripts
        return [getattr(CJLSTscripts, name) for name in self.item.names if "Ichep16" in name]

    @property
    def yamlname(self):
        return str(self).replace("tagged", "Tagged")

    def __contains__(self, other):
        return other in self.idnumbers

class WhichProdDiscriminants(MyEnum):
    """
    D_bkg and D_(0minus/0hplus/L1)_VBFdecay, but which third one?
    """
    enumname = "whichproddiscriminants"
    enumitems = (
                 EnumItem("D_int_decay"),
                 EnumItem("D_int_prod"),
                 EnumItem("D_g11gi3"),
                 EnumItem("D_g12gi2"),
                 EnumItem("D_g13gi1"),
                 EnumItem("D_g11gi3_prime"),
                 EnumItem("D_g12gi2_prime"),
                 EnumItem("D_g13gi1_prime"),
                )

channels = Channel.items()
if config.applyshapesystematics:
    systematics = Systematic.items()
    treesystematics = Systematic.items(lambda x: x in ("", "ResUp", "ResDown", "ScaleUp", "ScaleDown"))
else:
    systematics = treesystematics = Systematic.items(lambda x: x == "")
flavors = Flavor.items()
hypotheses = Hypothesis.items()
decayonlyhypotheses = Hypothesis.items(lambda x: x in ("0+", "a2", "0-", "L1", "fa20.5", "fa30.5", "fL10.5"))
prodonlyhypotheses = Hypothesis.items(lambda x: x in ("0+", "a2", "0-", "L1", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5"))
proddechypotheses = Hypothesis.items(lambda x: x in ("0+", "a2", "0-", "L1", "fa2dec0.5", "fa3dec0.5", "fL1dec0.5", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5", "fa2proddec-0.5", "fa3proddec-0.5", "fL1proddec-0.5"))
productionmodes = ProductionMode.items()
analyses = Analysis.items()
#productions = Production.items(lambda x: x in ("160225", "160729"))
config.productionsforcombine = type(config.productionsforcombine)(Production(production) for production in config.productionsforcombine)
productions = Production.items(lambda x: x in config.productionsforcombine)
blindstatuses = BlindStatus.items()
categories = Category.items()
whichproddiscriminants = WhichProdDiscriminants.items(lambda x: x == "D_int_prod")
#whichproddiscriminants = WhichProdDiscriminants.items()

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
        return self.items == other.items
    def __ne__(self, other):
        return not self == other
    def __hash__(self):
        return hash(self.items)

    def __str__(self):
        return " ".join(str(item) for item in self.items if item is not None)

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
