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
    def D_bkg_0plus(self):
        return "D_bkg_0plus"+self.appendname()
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
                 EnumItem("fL1"),
                )
    def title(self):
        if self == "fa3":
            return "f_{a3}"
        if self == "fa2":
            return "f_{a2}"
        if self == "fL1":
            return "f_{#Lambda1}"
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
                return "D_{a2}"
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
                return "D_g1g1prime2_decay"
        else:
            if self == "fa3":
                return "D_{CP}"
            if self == "fa2":
                return "D_{a1a2}"
            if self == "fL1":
                return "D_{a1#Lambda1}"
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
        from samples import Sample
        if self == "fa3":
            return [Sample("ggH", "0+"), Sample("ggH", "0-"), Sample("ggH", "fa30.5")]
        elif self == "fa2":
            return [Sample("ggH", "0+"), Sample("ggH", "a2"), Sample("ggH", "fa20.5")]
        elif self == "fL1":
            return [Sample("ggH", "0+"), Sample("ggH", "L1"), Sample("ggH", "fL10.5")]
        else:
            assert False
    def signaltemplates(self, channel):
        return [Template(sample, self, channel) for sample in self.signalsamples()]

channels = [Channel(item) for item in Channel.enumitems]
systematics = [Systematic(item) for item in Systematic.enumitems]
treesystematics = [Systematic(item) for item in ("", "ResUp", "ResDown", "ScaleUp", "ScaleDown")]
flavors = [Flavor(item) for item in Flavor.enumitems]
hypotheses = [Hypothesis(item) for item in Hypothesis.enumitems]
productionmodes = [ProductionMode(item) for item in ProductionMode.enumitems]
analyses = [Analysis(item) for item in Analysis.enumitems]

class MetaclassForMultiEnums(type):
    def __new__(cls, clsname, bases, dct):
        enums = dct["enums"]
        dct["needenums"] = needenums = enums[:]
        while True:
            for enum in needenums[:]:
                if issubclass(enum, MultiEnum):
                    #replace it with its subenums
                    needenums.remove(enum)
                    needenums += [a for a in enum.enums if a not in needenums]
                    #break out of for in order to continue in the while, to allow for more deeply nested MultiEnums
                    break
                elif issubclass(enum, MyEnum):
                    if enum not in needenums:
                        needenums.append(enum)
                else:
                    raise TypeError("{} is not a MyEnum or a MultiEnum! (in class {})".format(enum, clsname))
            else:
                break
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

    def check(self, *args):
        for enum in self.enums:
            if getattr(self, enum.enumname) is None:
                raise ValueError("No option provided for {}\n{}".format(enum.enumname, args))



class TemplatesFile(MultiEnum):
    enumname = "templatesfile"
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

    def templates(self):
        if self.signalorbkg == "signal":
            return [Template(self, sample.productionmode, sample.hypothesis) for sample in self.analysis.signalsamples()]
        elif self.signalorbkg == "bkg":
            return [Template(self, productionmode) for productionmode in ("qqZZ", "ggZZ", "ZX")]

class Template(MultiEnum):
    enums = [TemplatesFile, ProductionMode, Hypothesis]

    def applysynonyms(self, enumsdict):
        if enumsdict[SignalOrBkg] is None:
            if enumsdict[ProductionMode] == "ggH":
                enumsdict[SignalOrBkg] = "signal"
            if enumsdict[ProductionMode] in ("qqZZ", "ggZZ", "ZX"):
                enumsdict[SignalOrBkg] = "bkg"

    def check(self, *args):
        from samples import Sample
        if self.productionmode is None:
            raise ValueError("No option provided for productionmode\n{}".format(args))
        elif self.productionmode == "ggH":
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for ggH productionmode\n{}".format(args))
            if Sample(self.productionmode, self.hypothesis) not in self.templatesfile.analysis.signalsamples():
                raise ValueError("Hypothesis {} is not used in analysis {}!\n{}".format(self.hypothesis, self.templatesfile.analysis, args))
            if self.templatesfile.signalorbkg == "bkg":
                raise ValueError("ggH is not bkg!\n{}".format(args))
        elif self.productionmode in ("ggZZ", "qqZZ", "ZX"):
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.templatesfile.signalorbkg == "sig":
                raise ValueError("{} is not signal!\n{}".format(self.hypothesis, args))
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
                name = "template0HPlusAdapSmooth"
            elif self.hypothesis in ("fa20.5", "fa30.5", "fL10.5"):
                if final:
                    name = "templateIntAdapSmooth"
                else:
                    name = "templateMixAdapSmooth"
            else:
                assert False
        elif self.productionmode in ("ggZZ", "qqZZ", "ZX"):
            name = "template{}AdapSmooth".format(self.productionmode)
        else:
            assert False

        if self.templatesfile.analysis == "fa3" and final:
            name += "Mirror"

        return name

    def title(self):
        if self.productionmode == "ggH":
            return "{} {}".format(self.productionmode, self.hypothesis)
        if self.productionmode == "ggZZ" or self.productionmode == "qqZZ":
            return str(self.productionmode).replace("ZZ", "#rightarrowZZ")
        if self.productionmode == "ZX":
            return "Z+X"

    def gettemplate(self):
        f = tfiles[self.templatefile()]
        try:
            return getattr(f, self.templatename())
        except AttributeError:
            raise IOError("No template {} in {}".format(self.templatename(), self.templatefile()))

    def weightname(self):
        from samples import Sample
        if self.productionmode == "ggZZ":
            return Sample(self.productionmode, "2e2mu").weightname()
        if self.hypothesis is not None:
            return Sample(self.productionmode, self.hypothesis).weightname()
        return Sample(self.productionmode).weightname()

    def reweightfrom(self):
        from samples import Sample
        if self.productionmode == "ggH":
            if self.templatesfile.analysis in ("fa2", "fa3"):
                return [
                        Sample("ggH", "0+"),
                        Sample("ggH", "a2"),
                        Sample("ggH", "0-"),
                        Sample("ggH", "L1"),
                        Sample("ggH", "fa20.5"),
                        Sample("ggH", "fa30.5"),
                        #Sample("ggH", "fL10.5"),   #NOT fL1 for now
                       ]
            if self.templatesfile.analysis == "fL1":
                if self.hypothesis in ("0+", "L1"):
                    return [
                            Sample("ggH", "0+"),
                            Sample("ggH", "a2"),
                            Sample("ggH", "0-"),
                            Sample("ggH", "L1"),
                            Sample("ggH", "fa20.5"),
                            Sample("ggH", "fa30.5"),
                            #Sample("ggH", "fL10.5"),   #NOT fL1 for now
                           ]
                elif self.hypothesis == "fL10.5":
                    return [
                            #Sample("ggH", "0+"),
                            Sample("ggH", "a2"),
                            #Sample("ggH", "0-"),
                            Sample("ggH", "L1"),
                            Sample("ggH", "fa20.5"),
                            Sample("ggH", "fa30.5"),
                            Sample("ggH", "fL10.5"),
                           ]
        if self.productionmode in ("qqZZ", "ZX"):
            return [Sample(self.productionmode)]
        if self.productionmode == "ggZZ":
            return [Sample(self.productionmode, flavor) for flavor in flavors]
        assert False

    def scalefactor(self):
        if self.templatesfile.signalorbkg == "bkg": return 1
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
        return self.templatesfile.analysis == "fa3" and not (not final and self.hypothesis == "fa30.5")

    def discriminants(self):
        return [
                self.templatesfile.analysis.purediscriminant(),
                self.templatesfile.analysis.mixdiscriminant(),
                self.templatesfile.systematic.D_bkg_0plus(),
               ]
    def binning(self):
        if self.templatesfile.analysis == "fa2":
            result = [25, 0, 1, 25, 0, 1, 25, 0, 1]
        elif self.templatesfile.analysis == "fa3":
            result = [25, 0, 1, 25, -0.5, 0.5, 25, 0, 1]
        elif self.templatesfile.analysis == "fL1":
            result = [25, 0, 1, 25, -0.5, 1, 25, 0, 1]
        else:
            assert False
        for i in 1, 2, 4, 5, 7, 8:
            result[i] = float(result[i])
        return result

    def puretemplatestosubtract(self):
        if self.templatesfile.analysis == "fa2" and self.hypothesis == "fa20.5":
            return (Template(self.templatesfile, "ggH", "0+"), Template(self.templatesfile, "ggH", "a2"))
        if self.templatesfile.analysis == "fa3" and self.hypothesis == "fa30.5":
            return (Template(self.templatesfile, "ggH", "0+"), Template(self.templatesfile, "ggH", "0-"))
        if self.templatesfile.analysis == "fL1" and self.hypothesis == "fL10.5":
            return (Template(self.templatesfile, "ggH", "0+"), Template(self.templatesfile, "ggH", "L1"))
        assert False

    def smoothentriesperbin(self):
        if self.productionmode == "ZX":
            return 5
        elif self.templatesfile.signalorbkg == "bkg":
            return 20
        else:
            return 50

    def getjson(self):
        jsn = {
               "templates": [
                 {
                   "name": self.templatename(final=False),
                   "files": [os.path.basename(sample.withdiscriminantsfile()) for sample in self.reweightfrom()],
                   "tree": "candTree",
                   "variables": self.discriminants(),
                   "weight": self.weightname(),
                   "selection": "ZZMass>105 && ZZMass<140 && Z1Flav*Z2Flav == {}".format(self.templatesfile.channel.ZZFlav()),
                   "assertion": "D_0minus_decay >= 0. && D_0minus_decay <= 1.",
                   "binning": {
                     "type": "fixed",
                     "bins": self.binning(),
                   },
                   "conserveSumOfWeights": True,
                   "postprocessing": [
                     {"type": "smooth", "kernel": "adaptive", "entriesperbin": self.smoothentriesperbin()},
                     {"type": "reweight", "axes": [0,1,2]},
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

        return jsn

if __name__ == "__main__":
    print Template("ggH", "4mu", "fa30.5", "fa3").gettemplate().Integral()
    print Template("ggH", "4mu", "fa20.5", "fa2").gettemplate().Integral()
