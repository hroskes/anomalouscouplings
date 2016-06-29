import config
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
            if value in item.names:
                self.item = item
                break
        else:
            raise ValueError("%s is not a member of enum "%value + type(self).__name__ + "!  Valid choices:\n"
                               + "\n".join(" aka ".join(str(name) for name in item.names) for item in self.enumitems))

    def __str__(self):
        return str(self.item)

    def __eq__(self, other):
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

class Channel(MyEnum):
    enumname = "channel"
    enumitems = (
                 EnumItem("2e2mu", "2mu2e"),
                 EnumItem("4mu"),
                 EnumItem("4e"),
                )
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
    def appliesto(self, signalorbkg):
        if signalorbkg == "signal":
            return self in ("", "ResUp", "ResDown", "ScaleUp", "ScaleDown", "ScaleResUp", "ScaleResDown")
        elif signalorbkg == "bkg":
            return self in ("", "ZXUp", "ZXDown")
        assert False


class SignalOrBkg(MyEnum):
    enumname = "signalorbkg"
    enumitems = (
                 EnumItem("signal", "sig"),
                 EnumItem("background", "bkg"),
                )

class Analysis(MyEnum):
    enumname = "analysis"
    enumitems = (
                 EnumItem("fa3"),
                 EnumItem("fa2"),
                )
    def purediscriminant(self):
        if self == "fa3":
            return "D_0minus_decay"
        elif self == "fa2":
            return "D_g2_decay"
        elif self == "fL1":
            return "D_g1prime2_decay"
        else:
            assert False
    def mixdiscriminant(self):
        if self == "fa3":
            return "D_CP_decay"
        elif self == "fa2":
            return "D_g1g2_decay"
        elif self == "fL1":
            return "D_g1g1prime2_decay"
        else:
            assert False
    def signalsamples(self):
        from samples import Sample
        if self == "fa3":
            return (Sample("ggH", "0+"), Sample("ggH", "0-"), Sample("ggH", "fa30.5"))
        elif self == "fa2":
            return (Sample("ggH", "0+"), Sample("ggH", "a2"), Sample("ggH", "fa20.5"))
        elif self == "fL1":
            return (Sample("ggH", "0+"), Sample("ggH", "L1"), Sample("ggH", "fL10.5"))
        else:
            assert False
    def domirror(self):
        if self == "fa3": return True
        if self in ["fa2", "fL1"]: return False
        assert False

channels = [Channel(item) for item in Channel.enumitems]
systematics = [Systematic(item) for item in Systematic.enumitems]
treesystematics = [Systematic(item) for item in ("", "ResUp", "ResDown", "ScaleUp", "ScaleDown")]
flavors = [Flavor(item) for item in Flavor.enumitems]
hypotheses = [Hypothesis(item) for item in Hypothesis.enumitems]
productionmodes = [ProductionMode(item) for item in ProductionMode.enumitems]
analyses = [Analysis(item) for item in Analysis.enumitems]

class MultiEnum(object):
    multienums = []
    def __init__(self, *args):
        for enum in self.enums:
            setattr(self, enum.enumname, None)
        for arg in args:
            for enum in self.enums:
                try:
                    setattr(self, enum.enumname, enum(arg))
                    break
                except ValueError:
                    pass
            else:
                raise ValueError("{} is not a valid choice for any of: {}".format(arg, ", ".join(e.__name__ for e in self.enums)))
        for multienum in self.multienums:
            setattr(self, multienum.multienumname, multienum(*(getattr(self, enum.enumname) for enum in self.enums if enum in multienum.enums)))

        self.check(*args)
        self.items = tuple(getattr(self, enum.enumname) for enum in self.enums)

    def __eq__(self, other):
        return self.items == other.items
    def __ne__(self, other):
        return not self == other
    def __hash__(self):
        return hash(self.items)

    def __str__(self):
        return " ".join(str(item) for item in self.items if item is not None)

    def check(self, *args):
        for enum in self.enums:
            if getattr(self, enum.enumname) is None:
                raise ValueError("No option provided for {}\n{}".format(enum.enumname, args))



class TemplatesFile(MultiEnum):
    enums = [Channel, Systematic, SignalOrBkg, Analysis]

    def check(self, *args):
        if self.systematic is None:
            self.systematic = Systematic("")
        if not self.systematic.appliesto(self.signalorbkg):
            raise ValueError("Systematic {} does not apply to {}\n{}".format(self.systematic, self.signalorbkg, args))
        super(TemplatesFile, self).check(*args)

    def jsonfile(self):
        return os.path.join(config.repositorydir, "step5_json/templates_{}_{}{}.json".format(self.analysis, self.channel, "_bkg" if self.signalorbkg == "bkg" else self.systematic.appendname()))

    def templatesfile(self, run1=False):
        if not run1:
            return os.path.join(config.repositorydir, "step7_templates/{}{}{}_{}Adap_new.root".format(self.channel, "_bkg" if self.signalorbkg == "bkg" else "", self.systematic.appendname(), self.analysis))
        else:
            assert self.analysis == "fa3"
            return "/afs/cern.ch/work/x/xiaomeng/public/forChris/{}_fa3Adap_new{}.root".format(self.channel, "_bkg" if self.signalorbkg == "bkg" else self.systematic.appendname())
