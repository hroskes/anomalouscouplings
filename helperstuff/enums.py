from collections import OrderedDict
import config
import constants
from filemanager import tfiles
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
                 EnumItem("fa20.5"),
                 EnumItem("fa30.5"),
                 EnumItem("fL10.5", "flambda10.5"),
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
        if not title:
            return "D_bkg_0plus"+self.appendname()
        else:
            return "D_{bkg}^{"+self.appendname().replace("_", "")+"}"
    def appliesto(self, signalorbkg):
        if signalorbkg == "signal":
            return self in ("", "ResUp", "ResDown", "ScaleUp", "ScaleDown", "ScaleResUp", "ScaleResDown")
        elif signalorbkg == "bkg":
            return self in ("", "ZXUp", "ZXDown")
        elif signalorbkg == "DATA":
            return self in ("", )
        assert False


class SignalOrBkg(MyEnum):
    enumname = "signalorbkg"
    enumitems = (
                 EnumItem("signal", "sig"),
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
    def mixdiscriminantmin(self):
        if self == "fa3":
            return -.5
        elif self == "fa2":
            return 0.
        else:
            assert False
    def mixdiscriminantmax(self):
        if self == "fa3":
            return .5
        elif self == "fa2":
            return 1.
        else:
            assert False
    def signalsamples(self):
        from samples import ReweightingSample
        if self == "fa3":
            return [ReweightingSample("ggH", "0+"), ReweightingSample("ggH", "0-"), ReweightingSample("ggH", "fa30.5")]
        elif self == "fa2":
            return [ReweightingSample("ggH", "0+"), ReweightingSample("ggH", "a2"), ReweightingSample("ggH", "fa20.5")]
        elif self == "fL1":
            return [ReweightingSample("ggH", "0+"), ReweightingSample("ggH", "L1"), ReweightingSample("ggH", "fL10.5")]
        else:
            assert False
    def signaltemplates(self, *args):
        return [Template(sample, self, *args) for sample in self.signalsamples()]
    def interfxsec(self):
        if self == "fa3":
            return constants.JHUXS2L2la1a3 - 2*constants.JHUXS2L2la1
        elif self == "fa2":
            return constants.JHUXS2L2la1a2 - 2*constants.JHUXS2L2la1
        elif self == "fL1":
            return constants.JHUXS2L2la1L1 - 2*constants.JHUXS2L2la1
        assert False

class Production(MyEnum):
    enumname = "production"
    enumitems = (
                 EnumItem("160225"),
                 EnumItem("160624"),
                 EnumItem("160714"),
                 EnumItem("160720"),
                 EnumItem("160725", "160726"),
                )
    def CJLSTdir(self):
        if self == "160225":
            return "root://lxcms03//data3/Higgs/160225/"
        if self == "160624":
            return "root://lxcms03//data3/Higgs/160624/"
        if self == "160714":
            return "root://lxcms03//data3/Higgs/160714/"
        if self == "160720":
            return "root://lxcms03//data3/Higgs/160720/"
        if self == "160725":
            return "root://lxcms03//data3/Higgs/160726/"
        assert False
    def CJLSTdir_anomalous(self):
        if self == "160225":
            return "/afs/cern.ch/work/h/hroskes/reweighting_CJLST/CMSSW_7_6_3_patch2/src/ZZAnalysis/AnalysisStep/test/prod/AnomalousCouplingsReweighting/PT13TeV"
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
        return self.CJLSTdir()
    @property
    def useMELAv2(self):
        if self in ("160225", "160624"):
            return False
        if self in ("160714", "160720", "160725"):
            return True
        assert False
    @property
    def release(self):
        if self == "160225":
            return self.Release("76X")
        elif self in ("160624", "160714", "160720", "160725"):
            return self.Release("80X")
        assert False
    @property
    def dataluminosity(self):
        if self == "160225": return 2.8
        if self == "160714": return 7.65
        if self == "160720": return 9.2
        if self == "160725": return 12.9
        assert False
    def __int__(self):
        return int(str(self))
    @property
    def year(self):
        if self == "160225":
            return 2015
        if self in ("160624", "160714", "160720", "160725"):
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

channels = Channel.items()
systematics = Systematic.items()
treesystematics = Systematic.items(lambda x: x in ("", "ResUp", "ResDown", "ScaleUp", "ScaleDown"))
flavors = Flavor.items()
hypotheses = Hypothesis.items()
productionmodes = ProductionMode.items()
analyses = Analysis.items()
productions = Production.items(lambda x: x in ("160225", "160725"))
config.productionsforcombine = type(config.productionsforcombine)(Production(production) for production in config.productionsforcombine)
blindstatuses = BlindStatus.items()

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
    enums = [Channel, Systematic, SignalOrBkg, Analysis, Production, BlindStatus]

    def check(self, *args):
        if self.systematic is None:
            self.systematic = Systematic("")
        if self.blindstatus is None and self.signalorbkg == "DATA":
            raise ValueError("No option provided for blind!\n{}".format(args))
        if self.blindstatus is not None and self.signalorbkg != "DATA":
            raise ValueError("Can't blind MC!\n{}".format(args))
        super(TemplatesFile, self).check(*args, dontcheck=[BlindStatus])
        if not self.systematic.appliesto(self.signalorbkg):
            raise ValueError("Systematic {} does not apply to {}\n{}".format(self.systematic, self.signalorbkg, args))

    def jsonfile(self):
        folder = os.path.join(config.repositorydir, "step5_json")

        if self.signalorbkg == "bkg":
            result = os.path.join(folder, "templates_{}_{}_bkg_{}".format(self.analysis, self.channel, self.production))
        elif self.signalorbkg == "DATA":
            result = os.path.join(folder, "datatemplates_{}_{}_{}".format(self.analysis, self.channel, self.production))
        elif self.signalorbkg == "signal":
            result = os.path.join(folder, "templates_{}_{}{}_{}".format(self.analysis, self.channel, self.systematic.appendname(), self.production))
        else:
            assert False

        if self.signalorbkg == "DATA" and self.blind:
            result += "_blind"
        result += ".json"
        return result

    def templatesfile(self, run1=False):
        if run1:
            assert self.analysis == "fa3"
            return "/afs/cern.ch/work/x/xiaomeng/public/forChris/{}_fa3Adap_new{}.root".format(self.channel, "_bkg" if self.signalorbkg == "bkg" else self.systematic.appendname())

        folder = os.path.join(config.repositorydir, "step7_templates")

        if self.signalorbkg == "bkg":
            result = os.path.join(folder, "{}_bkg{}_{}Adap_{}".format(self.channel, self.systematic.appendname(), self.analysis, self.production))
        elif self.signalorbkg == "DATA":
            result = os.path.join(folder, "data_{}_{}_{}".format(self.channel, self.analysis, self.production))
        elif self.signalorbkg == "signal":
            result = os.path.join(folder, "{}{}_{}Adap_{}".format(self.channel, self.systematic.appendname(), self.analysis, self.production))
        else:
            assert False

        if self.signalorbkg == "DATA" and self.blind:
            result += "_blind"
        result += ".root"
        return result

    def signalsamples(self):
        from samples import Sample
        return [Sample(reweightingsample, self.production) for reweightingsample in self.analysis.signalsamples()]

    def templates(self):
        if self.signalorbkg == "signal":
            return [Template(self, sample.productionmode, sample.hypothesis) for sample in self.signalsamples()]
        elif self.signalorbkg == "bkg":
            return [Template(self, productionmode) for productionmode in ("qqZZ", "ggZZ", "ZX")]
        elif self.signalorbkg == "DATA":
            return [Template(self, "data")]

    @property
    def blind(self):
        if self.signalorbkg != "DATA": assert False
        return self.blindstatus == "blind"

templatesfiles = []
def tmp():
    for systematic in treesystematics:
        for channel in channels:
            for production in productions:
                for analysis in analyses:
                    templatesfiles.append(TemplatesFile(channel, systematic, "signal", analysis, production))
                    if systematic == "":
                        templatesfiles.append(TemplatesFile(channel, "bkg", analysis, production))
                    for blindstatus in blindstatuses:
                        if systematic == "" and (blindstatus == "blind" or config.unblinddistributions):
                            templatesfiles.append(TemplatesFile(channel, "DATA", analysis, production, blindstatus))
tmp()
del tmp

class Template(MultiEnum):
    enums = [TemplatesFile, ProductionMode, Hypothesis]
    enumname = "template"

    def applysynonyms(self, enumsdict):
        if enumsdict[SignalOrBkg] is None:
            if enumsdict[ProductionMode] == "ggH":
                enumsdict[SignalOrBkg] = "signal"
            elif enumsdict[ProductionMode] in ("qqZZ", "ggZZ", "ZX"):
                enumsdict[SignalOrBkg] = "bkg"
            elif enumsdict[ProductionMode] == "data":
                enumsdict[SignalOrBkg] = "DATA"
            else:
                assert False

    def check(self, *args):
        from samples import ReweightingSample
        if self.productionmode is None:
            raise ValueError("No option provided for productionmode\n{}".format(args))
        elif self.productionmode == "ggH":
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if ReweightingSample(self.productionmode, self.hypothesis) not in self.analysis.signalsamples():
                raise ValueError("Hypothesis {} is not used in analysis {}!\n{}".format(self.hypothesis, self.analysis, args))
            if self.signalorbkg != "signal":
                raise ValueError("{} is not {}!\n{}".format(self.productionmode, self.signalorbkg, args))
        elif self.productionmode in ("ggZZ", "qqZZ", "ZX"):
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.signalorbkg != "bkg":
                raise ValueError("{} is not {}!\n{}".format(self.hypothesis, self.signalorbkg, args))
        elif self.productionmode == "data":
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.signalorbkg != "DATA":
                raise ValueError("{} is not {}!\n{}".format(self.hypothesis, self.signalorbkg, args))
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
            else:
                assert False
        elif self.productionmode in ("ggZZ", "qqZZ", "ZX"):
            name = "template{}AdapSmooth".format(self.productionmode)
        elif self.productionmode == "data":
            name = "datatemplate"
        else:
            assert False

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
        if self.productionmode in ("qqZZ", "ZX"):
            result = [Sample(self.production, self.productionmode)]
        if self.productionmode == "ggZZ":
            result = [Sample(self.production, self.productionmode, flavor) for flavor in flavors]
        if self.productionmode == "data":
            result = [Sample(self.production, self.productionmode, self.blindstatus)]
        result = [sample for sample in result if tfiles[sample.withdiscriminantsfile()].candTree.GetEntries() != 0]
        assert result
        return result

    def scalefactor(self):
        if self.signalorbkg in ("bkg", "DATA"): return 1
        if self.hypothesis in ("0+", "0-", "a2", "L1"):
            result = 1.0
        elif self.hypothesis == "fa30.5":
            result = constants.JHUXS2L2la1a3 / constants.JHUXS2L2la1
        elif self.hypothesis == "fa20.5":
            result = constants.JHUXS2L2la1a2 / constants.JHUXS2L2la1
        elif self.hypothesis == "fL10.5":
            result = constants.JHUXS2L2la1L1 / constants.JHUXS2L2la1
        else:
            assert False
        result /= len(self.reweightfrom())
        return result

    def domirror(self, final=True):
        return self.analysis == "fa3" and not (not final and self.hypothesis == "fa30.5") and self.productionmode != "data"

    def discriminants(self):
        return (
                self.analysis.purediscriminant(),
                self.analysis.mixdiscriminant(),
                self.systematic.D_bkg_0plus(),
               )
    def binning(self):
        if self.analysis == "fa2":
            result = [50, 0, 1, 50, 0, 1, 50, 0, 1]
        elif self.analysis == "fa3":
            result = [50, 0, 1, 50, -0.5, 0.5, 50, 0, 1]
        elif self.analysis == "fL1":
            result = [50, 0, 1, 50, 0, 1, 50, 0, 1]
        else:
            assert False
        for i in 1, 2, 4, 5, 7, 8:
            result[i] = float(result[i])
        return result

    def puretemplatestosubtract(self):
        if self.analysis == "fa2" and self.hypothesis == "fa20.5":
            return (Template(self.templatesfile, "ggH", "0+"), Template(self.templatesfile, "ggH", "a2"))
        if self.analysis == "fa3" and self.hypothesis == "fa30.5":
            return (Template(self.templatesfile, "ggH", "0+"), Template(self.templatesfile, "ggH", "0-"))
        if self.analysis == "fL1" and self.hypothesis == "fL10.5":
            return (Template(self.templatesfile, "ggH", "0+"), Template(self.templatesfile, "ggH", "L1"))
        assert False

    def smoothentriesperbin(self):
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
        elif self.signalorbkg == "bkg":
            return 20
        return 50

    def reweightaxes(self):
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
                   "variables": self.discriminants(),
                   "weight": self.weightname(),
                   "selection": self.selection,
                   "assertion": "D_0minus_decay >= 0. && D_0minus_decay <= 1.",
                   "binning": {
                     "type": "fixed",
                     "bins": self.binning(),
                   },
                   "conserveSumOfWeights": True,
                   "postprocessing": [
                     {"type": "smooth", "kernel": "adaptive", "entriesperbin": self.smoothentriesperbin()},
                     {"type": "reweight", "axes": self.reweightaxes()},
                     {"type": "rescale","factor": self.scalefactor()},
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

        if self.hypothesis in ["fa30.5", "fa20.5", "fL10.5"]:
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
        if self.signalorbkg != "DATA": assert False
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
    if "160725" in productions:
        datatrees.append(SubtractDataTree("160725", "subtract160720", channel))
