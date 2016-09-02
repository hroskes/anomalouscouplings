import abc
import config
from enums import Analysis, analyses, AnalysisType, BlindStatus, blindstatuses, Channel, channels, Category, categories, EnumItem, flavors, Hypothesis, MultiEnum, MyEnum, Production, ProductionMode, productions, Systematic, TemplateGroup, treesystematics, WhichProdDiscriminants, whichproddiscriminants
from filemanager import tfiles
from itertools import product
import numpy
import os
from samples import ReweightingSample, Sample

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

        nameparts = ["templates", self.templategroup, self.analysis, self.whichproddiscriminants, self.channel, self.categorynamepart, self.systematic, self.production, self.blindnamepart]

        nameparts = [str(x) for x in nameparts if x]
        result = os.path.join(folder, "_".join(x for x in nameparts if x) + ".root")

        return result

    def signalsamples(self):
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
            if config.usedata:
                return [Template(self, productionmode) for productionmode in ("qqZZ", "ggZZ", "ZX")]
            else:
                return [Template(self, productionmode) for productionmode in ("qqZZ", "ggZZ")]
        elif self.templategroup == "DATA":
            return [Template(self, "data")]
        assert False

    def inttemplates(self):
        if self.templategroup == "ggh":
            return [IntTemplate(self, "ggH", "g11gi1")]
        elif self.templategroup == "vbf":
            return [IntTemplate(self, "VBF", "g1{}gi{}".format(i, 4-i)) for i in (1, 2, 3)]
        elif self.templategroup in ("bkg", "DATA"):
            return []
        assert False

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
        assert False

    @property
    def invertedmatrix(self):
        assert self.templategroup == "vbf"

        if not hasattr(self, "__invertedmatrix"):

            ganomalous = self.analysis.couplingname
            matrix = numpy.matrix(
                                  [
                                   [
                                    template.sample.g1**(4-i) * getattr(template.sample, ganomalous)**i
                                        for i in range(5)
                                   ]
                                       for template in self.templates()
                                  ]
                                 )
            """
            matrix looks something like this:
            1,    0,      0,        0,      0
            0,    0,      0,        0,      g4^4
            g1^4, g1^3g4, g1^2g4^2, g1g4^3, g4^4
            g1^4, g1^3g4, g1^2g4^2, g1g4^3, g4^4
            g1^4, g1^3g4, g1^2g4^2, g1g4^3, g4^4

            invert the matrix, then multiply by the vector of templates.  This should give back
               templates for each respective term (g1^4g4^0, g1^3g4^1, ...)
            In the PDF, the templates need to be multiplied by (g1^i)(g4^(4-i))
            """
            self.__invertedmatrix = invertedmatrix = matrix.I
            #make sure the two pure templates can be used as is
            #these assertions should be equivalent to asserting that SMtemplate.g1 == anomaloustemplate.ganomalous == 1
            # and that the two pure samples are in the right places on the list
            assert invertedmatrix[0,0] == 1 and all(invertedmatrix[0,i] == 0 for i in range(1,5))
            assert invertedmatrix[4,1] == 1 and invertedmatrix[4,0] == 0 and all(invertedmatrix[4,i] == 0 for i in range(2,5))

        return self.__invertedmatrix

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
                                if systematic == "" and (blindstatus == "blind" or config.unblinddistributions) and config.usedata:
                                    templatesfiles.append(TemplatesFile(channel, "DATA", analysis, production, blindstatus, category, "prod+dec", w))
tmp()
del tmp

class TemplateBase(object):
    __metaclass__ = abc.ABCMeta
    def templatefile(self):
        return self.templatesfile.templatesfile()

    @abc.abstractmethod
    def templatename(self):
        pass

    def gettemplate(self):
        f = tfiles[self.templatefile()]
        try:
            return getattr(f, self.templatename())
        except AttributeError:
            raise IOError("No template {} in {}".format(self.templatename(), self.templatefile()))

    @property
    def discriminants(self):
        return self.templatesfile.discriminants

    @property
    def binning(self):
        result = sum(([d.bins, d.min, d.max] for d in self.discriminants), [])
        for i in 1, 2, 4, 5, 7, 8:
            result[i] = float(result[i])
        return result

    @abc.abstractmethod
    def getjson(self):
        pass

class MultiEnumABCMeta(MultiEnum.__metaclass__, TemplateBase.__metaclass__):
    """
    needed to resolve conflict
    http://code.activestate.com/recipes/204197-solving-the-metaclass-conflict/
    except don't need all their fancy stuff
    the only function in MetaClassForMultiEnums is __new__ and it calls super
    """

class Template(TemplateBase, MultiEnum):
    __metaclass__ = MultiEnumABCMeta
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

    def weightname(self):
        if self.productionmode == "ggZZ":
            return ReweightingSample(self.productionmode, "2e2mu").weightname()
        if self.hypothesis is not None:
            return ReweightingSample(self.productionmode, self.hypothesis).weightname()
        if self.productionmode == "data":
            return None
        return ReweightingSample(self.productionmode).weightname()

    def reweightfrom(self):
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
        if self.templategroup in ("bkg", "DATA"): return 1
        if self.productionmode in ("VBF", "ggH"):
            result = ReweightingSample(self.productionmode, self.hypothesis).xsec / ReweightingSample(self.productionmode, "SM").xsec
        result /= len(self.reweightfrom())
        return result

    def domirror(self):
        return self.analysis == "fa3" and self.hypothesis not in ("fa30.5", "fa3prod0.5", "fa3proddec-0.5") and self.productionmode != "data"

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
#                     {"type": "smooth", "kernel": "adaptive", "entriesperbin": self.smoothentriesperbin()},
#                     {"type": "reweight", "axes": self.reweightaxes()},
                     {"type": "rescale","factor": self.scalefactor},
                   ],
                   "filloverflows": True,
                  },
                ],
              }

        if self.domirror():
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
        result += " && (" + " || ".join("category == {}".format(c) for c in self.category.idnumbers) + ")"
        return result

    @property
    def ZZFlav(self):
        result = self.channel.ZZFlav
        if self.productionmode == "ZX": result *= -1
        return result

    @property
    def sample(self):
        return ReweightingSample(self.hypothesis, self.productionmode)

class IntTemplate(TemplateBase, MultiEnum):
    __metaclass__ = MultiEnumABCMeta
    class InterferenceType(MyEnum):
        enumname = "interferencetype"
        enumitems = (
                     EnumItem("g11gi1"),
                     EnumItem("g11gi3"),
                     EnumItem("g12gi2"),
                     EnumItem("g13gi1"),
                    )
    enums = [TemplatesFile, ProductionMode, InterferenceType]

    def check(self, *args):
        if self.productionmode == "ggH":
            if self.interferencetype != "g11gi1":
                raise ValueError("Invalid interferencetype {} for productionmode {}!\n{}".format(self.interferencetype, self.productionmode, args))
        elif self.productionmode == "VBF":
            if self.interferencetype not in ("g11gi3", "g12gi2", "g13gi1"):
                raise ValueError("Invalid interferencetype {} for productionmode {}!\n{}".format(self.interferencetype, self.productionmode, args))
        else:
            raise ValueError("Invalid productionmode {}!\n{}".format(self.productionmode, args))
        super(IntTemplate, self).check(*args)

    def templatename(self):
        if self.productionmode == "ggH":
            if self.interferencetype == "g11gi1":
                result = "templateIntAdapSmooth"
        if self.productionmode == "VBF":
            if self.interferencetype == "g11gi3":
                result = "templateg11{}3AdapSmooth".format(self.analysis.couplingname)
            if self.interferencetype == "g12gi2":
                result = "templateg12{}2AdapSmooth".format(self.analysis.couplingname)
            if self.interferencetype == "g13gi1":
                result = "templateg13{}1AdapSmooth".format(self.analysis.couplingname)
        if self.domirror:
            result += "Mirror"
        return result

    @property
    def domirror(self):
        if self.analysis != "fa3": return False
        if self.interferencetype in ("g11gi1", "g11gi3", "g13gi1"):
            return True
        elif self.interferencetype == "g12gi2":
            return False
        assert False

    def getjson(self):
        if self.interferencetype == "g11gi1":
            puretemplates = [Template(self.templatesfile, self.productionmode, h) for h in self.analysis.purehypotheses]
            mixtemplate = Template(self.templatesfile, self.productionmode, self.analysis.mixdecayhypothesis)
            intjsn = {
                       "templates":[
                         {
                           "name": self.templatename(),
                           "templatesum":[
                             {"name":mixtemplate.templatename(final=False),"factor":1.0},
                           ] + [
                             {"name":t.templatename(final=False),"factor":-1.0} for t in puretemplates
                           ],
                           "postprocessing":[],
                         },
                       ],
                     }

        if self.interferencetype in ("g11gi3", "g12gi2", "g13gi1"):
            g1exp, giexp = (int(i) for i in str(self.interferencetype).replace("g1", "").replace("gi", ""))
            invertedmatrix = self.templatesfile.invertedmatrix
            vectoroftemplates = self.templatesfile.templates()  #ok technically it's a list not a vector
            rowofinvertedmatrix = giexp  #first row is labeled 0
            templatesum = []
            for j, template in enumerate(vectoroftemplates):
                templatesum.append({
                                    "name": template.templatename(final=False),
                                    "factor": invertedmatrix[rowofinvertedmatrix,j],
                                   })
            intjsn = {
                       "templates":[
                         {
                           "name": self.templatename(),
                           "templatesum":templatesum,
                           "postprocessing":[],
                         },
                       ],
                     }

        if self.domirror:
            intjsn["templates"][0]["postprocessing"].append({"type":"mirror", "antisymmetric":True, "axis":1})
        return intjsn


class DataTree(MultiEnum):
    enums = [Channel, Production]
    enumname = "datatree"
    @property
    def originaltreefile(self):
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
