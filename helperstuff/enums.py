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
    def jsonfile(self, systematic, isbkg):
        assert (not isbkg or systematic == "")
        return os.path.join(config.repositorydir, "step5_json/templates_{}{}.json".format(self, "_bkg" if isbkg else systematic.appendname()))
    def templatesfile(self, systematic, isbkg, run1=False):
        assert (not isbkg or systematic == "")
        if not run1:
            return os.path.join(config.repositorydir, "step7_templates/{}{}_fa3Adap_new.root".format(self, "_bkg" if isbkg else systematic.appendname()))
        else:
            return "/afs/cern.ch/work/x/xiaomeng/public/forChris/{}_fa3Adap_new{}.root".format(self, "_bkg" if isbkg else systematic.appendname())

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
                 EnumItem("WH"),
                 EnumItem("qqZZ"),
                 EnumItem("ggZZ"),
                 EnumItem("ZX"),
                 EnumItem("data"),
                )

class Systematic(MyEnum):
    enumitems = (
                 EnumItem(""),
                 EnumItem("ResUp"),
                 EnumItem("ResDown"),
                 EnumItem("ScaleUp"),
                 EnumItem("ScaleDown"),
                )
    def appendname(self):
        if self == "": return ""
        return "_" + str(self)

channels = [Channel(item) for item in Channel.enumitems]
systematics = [Systematic(item) for item in Systematic.enumitems]
flavors = [Flavor(item) for item in Flavor.enumitems]
hypotheses = [Hypothesis(item) for item in Hypothesis.enumitems]
productionmodes = [ProductionMode(item) for item in ProductionMode.enumitems]


class MultiEnum(object):
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
