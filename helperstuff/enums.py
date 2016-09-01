from collections import OrderedDict
import config
import constants
from filemanager import tfiles
from itertools import product
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

class Hypothesis(MyEnum):
    enumname = "hypothesis"
    enumitems = (
                 EnumItem("0+", "SM", "scalar"),
                 EnumItem("a2", "0h+"),
                 EnumItem("0-", "PS", "pseudoscalar"),
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

class ProductionMode(MyEnum):
    enumname = "productionmode"
    enumitems = (
                 EnumItem("ggH"),
                 EnumItem("VBF"),
                 EnumItem("H+jj", "HJJ"),
                 EnumItem("ZH"),
                 EnumItem("WminusH"),
                 EnumItem("WplusH"),
                 EnumItem("ttH"),
                 EnumItem("qqZZ"),
                 EnumItem("ggZZ"),
                 EnumItem("ZX"),
                 EnumItem("data"),
                )

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
        if templategroup in ("ggh", "vbf"):
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
    def signaltemplates(self, *args):
        return [Template(sample, self, *args) for sample in self.signalsamples()]
    def interfxsec(self):
        if self == "fa3":
            return constants.JHUXSggH2L2la1a3
        elif self == "fa2":
            return constants.JHUXSggH2L2la1a2
        elif self == "fL1":
            return constants.JHUXSggH2L2la1L1
        assert False

class Production(MyEnum):
    enumname = "production"
    enumitems = (
                 EnumItem("160225"),
                 EnumItem("160624"),
                 EnumItem("160714"),
                 EnumItem("160720"),
                 EnumItem("160725", "160726"),
                 EnumItem("160729"),
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
        if self <= "160729":
            return "/afs/cern.ch/work/h/hroskes/reweighting_CJLST/CMSSW_8_0_8/src/ZZAnalysis/AnalysisStep/test/prod/VBFanomalous_prodMEs/PT13TeV"
        return self.CJLSTdir()
    @property
    def useMELAv2(self):
        if self in ("160225", "160624"):
            return False
        return True
    @property
    def release(self):
        if self == "160225":
            return self.Release("76X")
        elif self in ("160624", "160714", "160720", "160725", "160729"):
            return self.Release("80X")
        assert False
    @property
    def dataluminosity(self):
        if self == "160225": return 2.8
        if self == "160714": return 7.65
        if self == "160720": return 9.2
        if self in ("160725", "160729"): return 12.9
        assert False
    def __int__(self):
        return int(str(self))
    @property
    def year(self):
        if self == "160225":
            return 2015
        if self in ("160624", "160714", "160720", "160725", "160729"):
            return 2016
        assert False

    #put this in here to avoid me getting really confused
    class Release(MyEnum):
        enumitems = (
                     EnumItem("76X"),
                     EnumItem("80X"),
                    )
        def __int__(self):
            return int(str(self).replace("X", ""))

class BlindStatus(MyEnum):
    enumname = "blindstatus"
    enumitems = (
                 EnumItem("unblind"),
                 EnumItem("blind"),
                )

class AnalysisType(MyEnum):
    enumname = "analysistype"
    enumitems = (
                 EnumItem("decayonly"),
                 EnumItem("prod+dec"),
                )

class Category(MyEnum):
    """
    For now just 2 categories, VBF2j and dump everything else into untagged
    """
    enumname = "category"
    enumitems = (
                 EnumItem("UntaggedIchep16", "VBF1jTaggedIchep16", "VHLeptTaggedIchep16", "VHHadrTaggedIchep16", "ttHTaggedIchep16"),
                 EnumItem("VBF2jTaggedIchep16"),
                )
    @property
    def idnumbers(self):
        """
        returns a list of ints corresponding to the C++ enums corresponding to this category
        (defined in Category.h)
        """
        import CJLSTscripts
        return [getattr(CJLSTscripts, name) for name in self.item.names]

class WhichProdDiscriminants(MyEnum):
    """
    D_bkg and D_(0minus/0hplus/L1)_VBFdecay, but which third one?
    """
    enumname = "whichproddiscriminants"
    enumitems = (
                 EnumItem("D_int_decay"),
                 EnumItem("D_int_VBF"),
                 EnumItem("g11gi3"),
                 EnumItem("g12gi2"),
                 EnumItem("g13gi1"),
                 EnumItem("g11gi3_prime"),
                 EnumItem("g12gi2_prime"),
                 EnumItem("g13gi1_prime"),
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
productions = Production.items(lambda x: x in ["160729"])
config.productionsforcombine = type(config.productionsforcombine)(Production(production) for production in config.productionsforcombine)
blindstatuses = BlindStatus.items()
categories = Category.items()
whichproddiscriminants = WhichProdDiscriminants.items()

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
        dct["needenums"] = needenums = enums[:]
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

class TemplatesFile(MultiEnum):
    enumname = "templatesfile"
    enums = [Channel, Systematic, TemplateGroup, Analysis, Production, BlindStatus, AnalysisType, WhichProdDiscriminants, Category]

    def check(self, *args):
        dontcheck = []

        if self.systematic is None:
            self.systematic = Systematic("")

        if self.category is None:
            self.category = Category("UntaggedIchep16")

        if self.templategroup != "DATA":
            if self.blindstatus is not None:
                raise ValueError("Can't blind MC!\n{}".format(args))
            dontcheck.append(BlindStatus)

        if self.analysistype == "decayonly":
            if self.whichproddiscriminants is not None:
                raise ValueError("Don't provide whichproddiscriminants for decayonly\n{}".format(args))
            if self.category != "UntaggedIchep16":
                raise ValueError("Only one category for decayonly analysis\n{}".format(args))
            dontcheck.append(WhichProdDiscriminants)

        if self.category == "UntaggedIchep16":
            if self.whichproddiscriminants is not None:
                raise ValueError("Don't provide whichproddiscriminants for untagged category\n{}".format(args))
            dontcheck.append(WhichProdDiscriminants)

        super(TemplatesFile, self).check(*args, dontcheck=dontcheck)

        if not self.systematic.appliesto(self.templategroup):
            raise ValueError("Systematic {} does not apply to {}\n{}".format(self.systematic, self.templategroup, args))

    def jsonfile(self):
        folder = os.path.join(config.repositorydir, "step5_json")

        nameparts = ["templates", self.templategroup, self.analysis, self.whichproddiscriminants, self.channel, self.categorynamepart, self.systematic, self.production, self.blindnamepart]

        nameparts = [str(x) for x in nameparts if x]
        result = os.path.join(folder, "_".join(x for x in nameparts if x) + ".json")

        return result

    def templatesfile(self):
        folder = os.path.join(config.repositorydir, "step7_templates")

        nameparts = ["templates", self.templategroup, self.analysis, self.channel, self.categorynamepart, self.systematic, self.production, self.blindnamepart]

        nameparts = [str(x) for x in nameparts]
        result = os.path.join(folder, "_".join(x for x in nameparts if x) + ".root")

        return result

    def signalsamples(self):
        from samples import ReweightingSample

        if self.templategroup == "ggh":
            if self.analysis == "fa3":
                reweightingsamples = [ReweightingSample("ggH", "0+"), ReweightingSample("ggH", "0-"), ReweightingSample("ggH", "fa30.5")]
            if self.analysis == "fa2":
                reweightingsamples = [ReweightingSample("ggH", "0+"), ReweightingSample("ggH", "a2"), ReweightingSample("ggH", "fa20.5")]
            if self.analysis == "fL1":
                reweightingsamples = [ReweightingSample("ggH", "0+"), ReweightingSample("ggH", "L1"), ReweightingSample("ggH", "fL10.5")]

        elif self.templategroup == "vbf":
            if self.analysis == "fa3":
                reweightingsamples = [ReweightingSample("VBF", "0+"), ReweightingSample("VBF", "0-"), ReweightingSample("VBF", "fa3prod0.5"), ReweightingSample("VBF", "fa3dec0.5"), ReweightingSample("VBF", "fa3proddec-0.5")]
            if self.analysis == "fa2":
                reweightingsamples = [ReweightingSample("VBF", "0+"), ReweightingSample("VBF", "a2"), ReweightingSample("VBF", "fa2prod0.5"), ReweightingSample("VBF", "fa2dec0.5"), ReweightingSample("VBF", "fa2proddec-0.5")]
            if self.analysis == "fL1":
                reweightingsamples = [ReweightingSample("VBF", "0+"), ReweightingSample("VBF", "L1"), ReweightingSample("VBF", "fL1prod0.5"), ReweightingSample("VBF", "fL1dec0.5"), ReweightingSample("VBF", "fL1proddec-0.5")]

        return reweightingsamples

    def templates(self):
        if self.templategroup in ["ggh", "vbf"]:
            return [Template(self, sample) for sample in self.signalsamples()]
        elif self.templategroup == "bkg":
            return [Template(self, productionmode) for productionmode in ("qqZZ", "ggZZ", "ZX")]
        elif self.templategroup == "DATA":
            return [Template(self, "data")]

    @property
    def bkgdiscriminant(self):
        return self.systematic.D_bkg_0plus()

    @property
    def purediscriminant(self):
        from discriminants import discriminant

        if self.category == "UntaggedIchep16":
            if self.analysis == "fa3":
                return discriminant("D_0minus_decay")
            if self.analysis == "fa2":
                return discriminant("D_g2_decay")
            if self.analysis == "fL1":
                return discriminant("D_g1prime2_decay")
        
        if self.category == "VBF2jTaggedIchep16":
            if self.analysis == "fa3":
                return discriminant("D_0minus_VBFdecay")
            if self.analysis == "fa2":
                return discriminant("D_g2_VBFdecay")
            if self.analysis == "fL1":
                return discriminant("D_g1prime2_VBFdecay")
        

    @property
    def mixdiscriminant(self):
        from discriminants import discriminant

        if self.category == "UntaggedIchep16" or (self.category == "VBF2jTaggedIchep16" and self.whichproddiscriminants == "D_int_decay"):
            if self.analysis == "fa3":
                return discriminant("D_CP_decay")
            if self.analysis == "fa2":
                return discriminant("D_g1g2_decay")
            if self.analysis == "fL1":
                return discriminant("D_g2_decay")
        
        if self.category == "VBF2jTaggedIchep16" and self.whichproddiscriminants == "D_int_VBF":
            if self.analysis == "fa3":
                return discriminant("D_CP_VBF")
            if self.analysis == "fa2":
                return discriminant("D_g1g2_VBF")
            if self.analysis == "fL1":
                return discriminant("D_g2_VBF")

        for i, prime in product(range(1, 4), ("", "_prime")):
            if self.category == "VBF2jTaggedIchep16" and self.whichproddiscriminants == "g1{}gi{}{}".format(i, 4-i, prime):
                if self.analysis == "fa3":
                    return discriminant("D_g1{}_g4{}_VBFdecay{}".format(i, 4-i, prime))
                if self.analysis == "fa2":
                    return discriminant("D_g1{}_g2{}_VBFdecay{}".format(i, 4-i, prime))
                if self.analysis == "fL1":
                    return discriminant("D_g1{}_g1prime2{}_VBFdecay{}".format(i, 4-i, prime))

        assert False

    @property
    def discriminants(self):
        return [self.purediscriminant, self.mixdiscriminant, self.bkgdiscriminant]

    @property
    def blind(self):
        assert self.templategroup == "DATA"
        return self.blindstatus == "blind"

    @property
    def blindnamepart(self):
        if self.templategroup != "DATA": return ""
        if self.blind: return "blind"
        return ""

    @property
    def categorynamepart(self):
        if self.category == "UntaggedIchep16":
            return ""
        if self.category == "VBF2jTaggedIchep16":
            return "VBFtag"

templatesfiles = []
def tmp():
    for systematic in treesystematics:
        for channel in channels:
            for production in productions:
                for analysis in analyses:
                    for category in categories:
                        for w in whichproddiscriminants:
                            if category == "UntaggedIchep16":
                                if w == whichproddiscriminants[0]:
                                    w = None
                                else:
                                    continue
                            templatesfiles.append(TemplatesFile(channel, systematic, "ggh", analysis, production, category, "prod+dec", w))
                            templatesfiles.append(TemplatesFile(channel, systematic, "vbf", analysis, production, category, "prod+dec", w))
                            if systematic == "":
                                templatesfiles.append(TemplatesFile(channel, "bkg", analysis, production, category, "prod+dec", w))
                            for blindstatus in blindstatuses:
                                if systematic == "" and (blindstatus == "blind" or config.unblinddistributions):
                                    templatesfiles.append(TemplatesFile(channel, "DATA", analysis, production, blindstatus, category, "prod+dec", w))
tmp()
del tmp

class Template(MultiEnum):
    enums = [TemplatesFile, ProductionMode, Hypothesis]
    enumname = "template"

    def applysynonyms(self, enumsdict):
        if enumsdict[TemplateGroup] is None:
            if enumsdict[ProductionMode] == "ggH":
                enumsdict[TemplateGroup] = "ggh"
            elif enumsdict[ProductionMode] == "VBF":
                enumsdict[TemplateGroup] = "vbf"
            elif enumsdict[ProductionMode] in ("qqZZ", "ggZZ", "ZX"):
                enumsdict[TemplateGroup] = "bkg"
            elif enumsdict[ProductionMode] == "data":
                enumsdict[TemplateGroup] = "DATA"
            else:
                assert False

    def check(self, *args):
        from samples import ReweightingSample
        if self.productionmode is None:
            raise ValueError("No option provided for productionmode\n{}".format(args))
        elif self.productionmode in ("ggH", "VBF"):
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if ReweightingSample(self.productionmode, self.hypothesis) not in self.templatesfile.signalsamples():
                raise ValueError("{} {} is not used to make templates for {} {}!\n{}".format(self.productionmode, self.hypothesis, self.templategroup, self.analysis, args))
            if self.templategroup != str(self.productionmode).lower():
                raise ValueError("{} is not {}!\n{}".format(self.productionmode, self.templategroup, args))
        elif self.productionmode in ("ggZZ", "qqZZ", "ZX"):
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.templategroup != "bkg":
                raise ValueError("{} is not {}!\n{}".format(self.hypothesis, self.templategroup, args))
        elif self.productionmode == "data":
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.templategroup != "DATA":
                raise ValueError("{} is not {}!\n{}".format(self.hypothesis, self.templategroup, args))
        else:
            raise ValueError("No templates for {}\n{}".format(self.productionmode, args))

    def templatefile(self, run1=False):
        return self.templatesfile.templatesfile(run1)

    def templatename(self, final=True):
        if self.productionmode == "ggH":
            if self.hypothesis == "0+":
                name = "template0PlusAdapSmooth"
            elif self.hypothesis == "0-":
                name = "template0MinusAdapSmooth"
            elif self.hypothesis == "a2":
                name = "template0HPlusAdapSmooth"
            elif self.hypothesis == "L1":
                name = "template0L1AdapSmooth"
            elif self.hypothesis in ("fa20.5", "fa30.5", "fL10.5"):
                if final:
                    name = "templateIntAdapSmooth"
                else:
                    name = "templateMixAdapSmooth"
        elif self.productionmode == "VBF":
            if self.hypothesis == "0+":
                name = "template0PlusAdapSmooth"
            elif self.hypothesis == "0-":
                name = "template0MinusAdapSmooth"
            elif self.hypothesis == "a2":
                name = "template0HPlusAdapSmooth"
            elif self.hypothesis == "L1":
                name = "template0L1AdapSmooth"
            elif self.hypothesis in ("fa2dec0.5", "fa3dec0.5", "fL1dec0.5"):
                name = "templateMixDecayAdapSmooth"
            elif self.hypothesis in ("fa2prod0.5", "fa3prod0.5", "fL1prod0.5"):
                name = "templateMixProdAdapSmooth"
            elif self.hypothesis in ("fa2proddec-0.5", "fa3proddec-0.5", "fL1proddec-0.5"):
                name = "templateMixProdDecPiAdapSmooth"
        elif self.productionmode in ("ggZZ", "qqZZ", "ZX"):
            name = "template{}AdapSmooth".format(self.productionmode)
        elif self.productionmode == "data":
            name = "datatemplate"

        name

        if self.analysis == "fa3" and final and self.productionmode != "data":
            name += "Mirror"

        return name

    def title(self):
        if self.productionmode == "ggH":
            return "{} {}".format(self.productionmode, self.hypothesis)
        if self.productionmode == "ggZZ" or self.productionmode == "qqZZ":
            return str(self.productionmode).replace("ZZ", "#rightarrowZZ")
        if self.productionmode == "ZX":
            return "Z+X"
        if self.productionmode == "data":
            return "data"
        assert False

    def gettemplate(self):
        f = tfiles[self.templatefile()]
        try:
            return getattr(f, self.templatename())
        except AttributeError:
            raise IOError("No template {} in {}".format(self.templatename(), self.templatefile()))

    def weightname(self):
        from samples import ReweightingSample
        if self.productionmode == "ggZZ":
            return ReweightingSample(self.productionmode, "2e2mu").weightname()
        if self.hypothesis is not None:
            return ReweightingSample(self.productionmode, self.hypothesis).weightname()
        if self.productionmode == "data":
            return None
        return ReweightingSample(self.productionmode).weightname()

    def reweightfrom(self):
        from samples import Sample
        if self.productionmode == "ggH":
            if self.analysis in ("fa2", "fa3"):
                result=[
                        Sample(self.production, "ggH", "0+"),
                        Sample(self.production, "ggH", "a2"),
                        Sample(self.production, "ggH", "0-"),
                        Sample(self.production, "ggH", "L1"),
                        Sample(self.production, "ggH", "fa20.5"),
                        Sample(self.production, "ggH", "fa30.5"),
                        #Sample(self.production, "ggH", "fL10.5"),   #NOT fL1 for now
                       ]
            if self.analysis == "fL1":
                if self.hypothesis in ("0+", "L1"):
                    result=[
                            Sample(self.production, "ggH", "0+"),
                            Sample(self.production, "ggH", "a2"),
                            Sample(self.production, "ggH", "0-"),
                            Sample(self.production, "ggH", "L1"),
                            Sample(self.production, "ggH", "fa20.5"),
                            Sample(self.production, "ggH", "fa30.5"),
                            #Sample(self.production, "ggH", "fL10.5"),   #NOT fL1 for now
                           ]
                elif self.hypothesis == "fL10.5":
                    result=[
                            #Sample(self.production, "ggH", "0+"),
                            Sample(self.production, "ggH", "a2"),
                            #Sample(self.production, "ggH", "0-"),
                            Sample(self.production, "ggH", "L1"),
                            Sample(self.production, "ggH", "fa20.5"),
                            Sample(self.production, "ggH", "fa30.5"),
                            Sample(self.production, "ggH", "fL10.5"),
                           ]
        if self.productionmode == "VBF" and self.analysistype == "prod+dec":
            if self.hypothesis == "0+":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("0+", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       ]
            if self.hypothesis == "a2":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("0-", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5") #0- and L1 are borderline
                       ]
            if self.hypothesis == "0-":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("0-", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       ]
            if self.hypothesis == "L1":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("L1", "fL1prod0.5")   #????? what happened?
                       ]
            if self.hypothesis == "fa2dec0.5":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("0-", "a2", "L1", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5") #0+, 0-, and L1 are borderline
                       ]
            if self.hypothesis == "fa3dec0.5":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("0-", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       ]
            if self.hypothesis == "fL1dec0.5":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("fL1prod0.5",)
                       ]
            if self.hypothesis == "fa2prod0.5":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("0+", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       ]
            if self.hypothesis == "fa3prod0.5":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("0+", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")  #0+ and a2 are borderline
                       ]
            if self.hypothesis == "fL1prod0.5":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("0+", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5") #????
                       ]
            if self.hypothesis == "fa2proddec-0.5":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       ]
            if self.hypothesis == "fa3proddec-0.5":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("0-", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       ]
            if self.hypothesis == "fL1proddec-0.5":
                return [
                        Sample(self.production, "VBF", hypothesis)
                            for hypothesis in ("a2", "fa3prod0.5", "fL1prod0.5") #????
                       ]
        if self.productionmode in ("qqZZ", "ZX"):
            result = [Sample(self.production, self.productionmode)]
        if self.productionmode == "ggZZ":
            result = [Sample(self.production, self.productionmode, flavor) for flavor in flavors]
        if self.productionmode == "data":
            result = [Sample(self.production, self.productionmode, self.blindstatus)]
        result = [sample for sample in result if tfiles[sample.withdiscriminantsfile()].candTree.GetEntries() != 0]
        assert result
        return result

    @property
    def scalefactor(self):
        from samples import ReweightingSample
        if self.templategroup in ("bkg", "DATA"): return 1
        if self.productionmode in ("VBF", "ggH"):
            result = ReweightingSample(self.productionmode, self.hypothesis).xsec / ReweightingSample(self.productionmode, "SM").xsec
        result /= len(self.reweightfrom())
        return result

    def domirror(self, final=True):
        return self.analysis == "fa3" and not (not final and self.hypothesis in ("fa30.5", "fa3prod0.5", "fa3proddec-0.5")) and self.productionmode != "data"

    @property
    def discriminants(self):
        return self.templatesfile.discriminants

    @property
    def binning(self):
        result = sum(([d.bins, d.min, d.max] for d in self.discriminants), [])
        for i in 1, 2, 4, 5, 7, 8:
            result[i] = float(result[i])
        return result

    def puretemplatestosubtract(self):
        if self.templategroup == "ggh":
            if self.analysis == "fa2" and self.hypothesis == "fa20.5":
                return (Template(self.templatesfile, "ggH", "0+"), Template(self.templatesfile, "ggH", "a2"))
            if self.analysis == "fa3" and self.hypothesis == "fa30.5":
                return (Template(self.templatesfile, "ggH", "0+"), Template(self.templatesfile, "ggH", "0-"))
            if self.analysis == "fL1" and self.hypothesis == "fL10.5":
                return (Template(self.templatesfile, "ggH", "0+"), Template(self.templatesfile, "ggH", "L1"))
        assert False

    def smoothentriesperbin(self):
        if self.analysistype == "decayonly":
            if self.analysis == "fL1":
                if self.channel in ["2e2mu", "4e"]:
                    if self.productionmode == "ZX":
                        return 20
                    if self.productionmode == "ggZZ":
                        return 100
                    if self.productionmode == "qqZZ":
                        return 70
                if self.channel == "4mu":
                    if self.productionmode == "ZX":
                        return 9
                    if self.productionmode == "ggZZ":
                        return 100
                    if self.productionmode == "qqZZ":
                        return 100
            if self.analysis == "fa2":
                if self.channel == "2e2mu":
                    if self.productionmode == "ZX":
                        return 45
                    if self.productionmode == "ggZZ":
                        return 20
                    if self.productionmode == "qqZZ":
                        return 30
                if self.channel == "4e":
                    if self.productionmode == "ZX":
                        return 6   #too much and the peak in axis 1 at .9 goes away, less (or reweight) and the bump at .4 comes back
                    if self.productionmode == "ggZZ":
                        return 100
                    if self.productionmode == "qqZZ":
                        if self.production == "160225":
                            return 25
                        return 15  #similar to Z+X
                if self.channel == "4mu":
                    if self.productionmode == "ZX":
                        return 9
                    if self.productionmode == "ggZZ":
                        return 150
                    if self.productionmode == "qqZZ":
                        if self.production == "160225":
                            return 50
                        return 100
            if self.analysis == "fa3":
                if self.channel == "2e2mu":
                    if self.productionmode == "ZX":
                        return 20
                    if self.productionmode == "ggZZ":
                        return 500
                    if self.productionmode == "qqZZ":
                        return 60
                if self.channel == "4e":
                    if self.productionmode == "ZX":
                        return 10
                    if self.productionmode == "ggZZ":
                        return 100
                    if self.productionmode == "qqZZ":
                        return 50
                if self.channel == "4mu":
                    if self.productionmode == "ZX":
                        return 10
                    if self.productionmode == "ggZZ":
                        return 100
                    if self.productionmode == "qqZZ":
                        return 100
            if self.productionmode == "ZX":
                return 5
            elif self.templategroup == "bkg":
                return 20
        return 50

    def reweightaxes(self):
        if self.analysistype == "decayonly":
            if self.channel == "2e2mu" and self.productionmode == "ggZZ" and self.analysis == "fa3":
                return [1, 2]
            if self.channel == "2e2mu" and self.productionmode == "ggZZ" and self.analysis == "fa2":
                return [1, 2]
            if self.channel == "2e2mu" and self.productionmode == "ggZZ" and self.analysis == "fL1":
                return [0, 2]
            if self.channel == "2e2mu" and self.productionmode == "qqZZ" and self.analysis == "fa2":
                return [2]
            if self.channel == "2e2mu" and self.productionmode == "qqZZ" and self.analysis == "fL1":
                return [0, 2]
            if self.channel == "2e2mu" and self.productionmode == "qqZZ" and self.analysis == "fa3":
                return [1, 2]
            if self.channel == "2e2mu" and self.productionmode == "ZX"   and self.analysis == "fL1":
                return [0, 2]
            if self.channel == "2e2mu" and self.productionmode == "ZX"   and self.analysis == "fa3":
                return [2]
            if self.channel == "2e2mu" and self.productionmode == "ZX"   and self.analysis == "fa2":
                return [2]
            if self.channel == "4e"    and self.productionmode == "ggZZ" and self.analysis == "fa3":
                return [1, 2]
            if self.channel == "4e"    and self.productionmode == "qqZZ" and self.analysis == "fa2":
                if self.production == "160225":
                    return [2]
                return [0, 2]
            if self.channel == "4e"    and self.productionmode == "qqZZ" and self.analysis == "fa3":
                return [1, 2]
            if self.channel == "4e"    and self.productionmode == "ZX"   and self.analysis == "fL1":
                return [0, 2]
            if self.channel == "4e"    and self.productionmode == "ZX"   and self.analysis == "fa3":
                return [1]
            if self.channel == "4e"    and self.productionmode == "ZX"   and self.analysis == "fa2":
                return []
            if self.channel == "4mu"   and self.productionmode == "qqZZ" and self.analysis in ["fa3", "fa2"]:
                return [1, 2]
            if self.channel == "4mu"   and self.productionmode == "qqZZ" and self.analysis == "fL1":
                if self.production != "160225":
                    return [0, 2]
            if self.channel == "4mu"   and self.productionmode == "ZX":
                return []
        return [0, 1, 2]

    def getjson(self):
        jsn = {
               "templates": [
                 {
                   "name": self.templatename(final=False),
                   "files": [os.path.basename(sample.withdiscriminantsfile()) for sample in self.reweightfrom()],
                   "tree": "candTree",
                   "variables": [d.name for d in self.discriminants],
                   "weight": self.weightname(),
                   "selection": self.selection,
                   "assertion": "D_0minus_decay >= 0. && D_0minus_decay <= 1.",
                   "binning": {
                     "type": "fixed",
                     "bins": self.binning,
                   },
                   "conserveSumOfWeights": True,
                   "postprocessing": [
                     {"type": "smooth", "kernel": "adaptive", "entriesperbin": self.smoothentriesperbin()},
                     {"type": "reweight", "axes": self.reweightaxes()},
                     {"type": "rescale","factor": self.scalefactor},
                   ],
                   "filloverflows": True,
                  },
                ],
              }

        if self.domirror(final=False):
            mirrorjsn = {
                          "templates":[
                            {
                              "name":self.templatename(final=True),
                              "templatesum":[
                                {"name":self.templatename(final=False),"factor":1.0},
                              ],
                              "postprocessing":[
                                {"type":"mirror", "axis":1},
                                {"type":"floor"},
                              ],
                            },
                          ],
                        }
            jsn["templates"] += mirrorjsn["templates"]
        else:
            jsn["templates"][0]["postprocessing"].append({"type": "floor"})

        if self.hypothesis in ["fa30.5", "fa20.5", "fL10.5"] and self.templategroup == "ggh":
            intjsn = {
                       "templates":[
                         {
                           "name": self.templatename(final=True),
                           "templatesum":[
                             {"name":self.templatename(final=False),"factor":1.0},
                           ] + [
                             {"name":t.templatename(final=False),"factor":-1.0} for t in self.puretemplatestosubtract()
                           ],
                           "postprocessing":[
                           ],
                         },
                       ],
                     }
            if self.domirror(final=True):
                intjsn["templates"][0]["postprocessing"].append({"type":"mirror", "antisymmetric":True, "axis":1})
            jsn["templates"] += intjsn["templates"]

        if self.productionmode == "data":
            del jsn["templates"][0]["postprocessing"]
            del jsn["templates"][0]["weight"]

        return jsn

    @property
    def blind(self):
        if self.templategroup != "DATA": assert False
        return self.templatesfile.blind

    @property
    def selection(self):
        result = "ZZMass>{} && ZZMass<{} && Z1Flav*Z2Flav == {}".format(config.m4lmin, config.m4lmax, self.ZZFlav)
        return result

    @property
    def ZZFlav(self):
        result = self.channel.ZZFlav
        if self.productionmode == "ZX": result *= -1
        return result

class DataTree(MultiEnum):
    enums = [Channel, Production]
    enumname = "datatree"
    @property
    def originaltreefile(self):
        from samples import Sample
        return Sample("data", self.production, "unblind").withdiscriminantsfile()
    @property
    def treefile(self):
        return os.path.join(config.repositorydir, "step7_templates", "data_{}_{}.root".format(self.production, self.channel))
    def passescut(self, t):
        return abs(t.Z1Flav * t.Z2Flav) == self.channel.ZZFlav and config.m4lmin < t.ZZMass < config.m4lmax and config.unblindscans

datatrees = []
for channel in channels:
    for production in productions:
        datatrees.append(DataTree(channel, production))



class SubtractProduction(MyEnum):
    enumname = "subtractproduction"
    enumitems = (
                 EnumItem("subtract160720"),
                )
    subtracttree = None
    def passescut(self, t):
        if self.subtracttree is None:
            from samples import Sample
            self.subtracttree = tfiles[Sample("data", "unblind", str(self).replace("subtract", "")).withdiscriminantsfile()].candTree
        run, event, lumi = t.RunNumber, t.EventNumber, t.LumiNumber
        for t2 in self.subtracttree:
            if (run, event, lumi) == (t2.RunNumber, t2.EventNumber, t2.LumiNumber):
                return False
        return True

class SubtractDataTree(DataTree, MultiEnum):
    enums = DataTree.enums + (SubtractProduction,)

    @property
    def treefile(self):
        return os.path.join(config.repositorydir, "step7_templates", "data_{}_{}_{}.root".format(self.production, self.channel, self.subtractproduction))
    def check(self, *args, **kwargs):
        return super(SubtractDataTree, self).check(*args, **kwargs)
    def passescut(self, t):
        return super(SubtractDataTree, self).passescut(t) and self.subtractproduction.passescut(t)

for channel in channels:
    if "160729" in productions:
        datatrees.append(SubtractDataTree("160729", "subtract160720", channel))
