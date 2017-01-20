import abc
import config
from enums import Analysis, analyses, BlindStatus, blindstatuses, Channel, channels, Category, categories, EnumItem, flavors, Hypothesis, MultiEnum, MultiEnumABCMeta, MyEnum, prodonlyhypotheses, Production, ProductionMode, productions, ShapeSystematic, TemplateGroup, treeshapesystematics, WhichProdDiscriminants, whichproddiscriminants
from itertools import product
import json
import numpy
import os
from samples import ReweightingSample, Sample, SampleBasis
from utilities import cache, getnesteddictvalue, jsonloads, tfiles

class TemplatesFile(MultiEnum):
    enumname = "templatesfile"
    enums = [Channel, ShapeSystematic, TemplateGroup, Analysis, Production, BlindStatus, WhichProdDiscriminants, Category]

    def applysynonyms(self, enumsdict):
        if enumsdict[Production] is None and len(config.productionsforcombine) == 1:
            enumsdict[Production] = config.productionsforcombine[0]
        if enumsdict[WhichProdDiscriminants] is None and len(whichproddiscriminants) == 1:
            enumsdict[WhichProdDiscriminants] = whichproddiscriminants[0]
        super(TemplatesFile, self).applysynonyms(enumsdict)

    def check(self, *args):
        dontcheck = []

        if self.shapesystematic is None:
            self.shapesystematic = ShapeSystematic("")

        if self.category is None:
            self.category = Category("UntaggedIchep16")

        if self.category == "UntaggedIchep16":
            self.whichproddiscriminants = None

        if self.templategroup != "DATA":
            if self.blindstatus is not None:
                raise ValueError("Can't blind MC!\n{}".format(args))
            dontcheck.append(BlindStatus)

        if self.category == "UntaggedIchep16":
            if len(whichproddiscriminants) == 1: self.whichproddiscriminants = None
            if self.whichproddiscriminants is not None:
                raise ValueError("Don't provide whichproddiscriminants for untagged category\n{}".format(args))
            dontcheck.append(WhichProdDiscriminants)

        if self.analysis == "fL1":
            if len(whichproddiscriminants) == 1: self.whichproddiscriminants = None
            if self.whichproddiscriminants is not None:
                raise ValueError("Don't provide whichproddiscriminants for fL1\n{}".format(args))
            dontcheck.append(WhichProdDiscriminants)

        super(TemplatesFile, self).check(*args, dontcheck=dontcheck)

        if not self.shapesystematic.appliesto(self.templategroup):
            raise ValueError("ShapeSystematic {} does not apply to {}\n{}".format(self.shapesystematic, self.templategroup, args))

    def jsonfile(self, iteration=None):
        folder = os.path.join(config.repositorydir, "step5_json")
        if iteration is not None:
            folder = os.path.join(folder, "bkp_iter{}".format(iteration))

        nameparts = ["templates", self.templategroup, self.analysis, self.whichproddiscriminants, self.channel, self.categorynamepart, self.shapesystematic, self.production, self.blindnamepart]

        nameparts = [str(x) for x in nameparts if x]
        result = os.path.join(folder, "_".join(x for x in nameparts if x) + ".json")

        return result

    def templatesfile(self, iteration=None):
        folder = os.path.join(config.repositorydir, "step7_templates")
        if iteration is not None:
            folder = os.path.join(folder, "bkp_iter{}".format(iteration))
            if not os.path.exists(folder):
                raise IOError("No folder {}".format(folder))

        nameparts = ["templates", self.templategroup, self.analysis, self.whichproddiscriminants, self.channel, self.categorynamepart, self.shapesystematic if config.applyshapesystematics else "", self.production, self.blindnamepart]

        nameparts = [str(x) for x in nameparts if x]
        result = os.path.join(folder, "_".join(x for x in nameparts if x) + ".root")

        return result

    def controlplotsdir(self, *args, **kwargs):
        relpath = os.path.relpath(self.templatesfile(*args, **kwargs), os.path.join(config.repositorydir, "step7_templates"))
        assert ".." not in relpath
        return os.path.join(config.plotsbasedir, "templateprojections", "controlplots", relpath.replace(".root", "").replace("bkp_", ""))

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
                reweightingsamples = [ReweightingSample("VBF", "0+"), ReweightingSample("VBF", "L1"), ReweightingSample("VBF", "fL1prod0.5"), ReweightingSample("VBF", "fL1dec0.5"), ReweightingSample("VBF", "fL1proddec0.5")]

        elif self.templategroup == "zh":
            if self.analysis == "fa3":
                reweightingsamples = [ReweightingSample("ZH", "0+"), ReweightingSample("ZH", "0-"), ReweightingSample("ZH", "fa3prod0.5"), ReweightingSample("ZH", "fa3dec0.5"), ReweightingSample("ZH", "fa3proddec-0.5")]
            if self.analysis == "fa2":
                reweightingsamples = [ReweightingSample("ZH", "0+"), ReweightingSample("ZH", "a2"), ReweightingSample("ZH", "fa2prod0.5"), ReweightingSample("ZH", "fa2dec0.5"), ReweightingSample("ZH", "fa2proddec-0.5")]
            if self.analysis == "fL1":
                reweightingsamples = [ReweightingSample("ZH", "0+"), ReweightingSample("ZH", "L1"), ReweightingSample("ZH", "fL1prod0.5"), ReweightingSample("ZH", "fL1dec0.5"), ReweightingSample("ZH", "fL1proddec0.5")]

        elif self.templategroup == "wh":
            if self.analysis == "fa3":
                reweightingsamples = [ReweightingSample("WH", "0+"), ReweightingSample("WH", "0-"), ReweightingSample("WH", "fa3prod0.5"), ReweightingSample("WH", "fa3dec0.5"), ReweightingSample("WH", "fa3proddec-0.5")]
            if self.analysis == "fa2":
                reweightingsamples = [ReweightingSample("WH", "0+"), ReweightingSample("WH", "a2"), ReweightingSample("WH", "fa2prod0.5"), ReweightingSample("WH", "fa2dec0.5"), ReweightingSample("WH", "fa2proddec-0.5")]
            if self.analysis == "fL1":
                reweightingsamples = [ReweightingSample("WH", "0+"), ReweightingSample("WH", "L1"), ReweightingSample("WH", "fL1prod0.5"), ReweightingSample("WH", "fL1dec0.5"), ReweightingSample("WH", "fL1proddec0.5")]

        return reweightingsamples

    def templates(self):
        if self.templategroup in ["ggh", "vbf", "zh", "wh"]:
            return [Template(self, sample) for sample in self.signalsamples()]
        elif self.templategroup == "bkg":
            result = ["qqZZ", "ggZZ"]
            if config.usedata:
                result.append("ZX")
            result.append("VBF bkg")
            return [Template(self, productionmode) for productionmode in result]
        elif self.templategroup == "DATA":
            return [Template(self, "data")]
        assert False

    def inttemplates(self):
        if self.templategroup == "ggh":
            return [IntTemplate(self, "ggH", "g11gi1")]
        elif self.templategroup == "vbf":
            return [IntTemplate(self, "VBF", "g1{}gi{}".format(i, 4-i)) for i in (1, 2, 3)]
        elif self.templategroup == "zh":
            return [IntTemplate(self, "ZH", "g1{}gi{}".format(i, 4-i)) for i in (1, 2, 3)]
        elif self.templategroup == "wh":
            return [IntTemplate(self, "WH", "g1{}gi{}".format(i, 4-i)) for i in (1, 2, 3)]
        elif self.templategroup in ("bkg", "DATA"):
            return []
        assert False

    @property
    def bkgdiscriminant(self):
        return self.shapesystematic.D_bkg_0plus()

    @property
    def purediscriminant(self):
        from discriminants import discriminant

        if self.category == "UntaggedIchep16":
            if self.analysis == "fa3":
                return discriminant("D_0minus_decay")
            if self.analysis == "fa2":
                return discriminant("D_0hplus_decay")
            if self.analysis == "fL1":
                return discriminant("D_L1_decay")

        if self.category == "VBF2jTaggedIchep16":
            if self.analysis == "fa3":
                return discriminant("D_0minus_VBFdecay")
            if self.analysis == "fa2":
                return discriminant("D_0hplus_VBFdecay")
            if self.analysis == "fL1":
                return discriminant("D_L1_VBFdecay")

        if self.category == "VHHadrTaggedIchep16":
            if self.analysis == "fa3":
                return discriminant("D_0minus_HadVHdecay")
            if self.analysis == "fa2":
                return discriminant("D_0hplus_HadVHdecay")
            if self.analysis == "fL1":
                return discriminant("D_L1_HadVHdecay")

        assert False

    @property
    def mixdiscriminant(self):
        from discriminants import discriminant

        if self.category == "UntaggedIchep16" or self.whichproddiscriminants == "D_int_decay":
            if self.analysis == "fa3":
                return discriminant("D_CP_decay")
            if self.analysis == "fa2":
                return discriminant("D_int_decay")
            if self.analysis == "fL1":
                return discriminant("D_0hplus_decay")

        if self.category == "VBF2jTaggedIchep16" and self.whichproddiscriminants == "D_int_prod":
            if self.analysis == "fa3":
                return discriminant("D_CP_VBF")
            if self.analysis == "fa2":
                return discriminant("D_int_VBF")
            if self.analysis == "fL1":
                return discriminant("D_0hplus_VBF")

        if self.category == "VHHadrTaggedIchep16" and self.whichproddiscriminants == "D_int_prod":
            if self.analysis == "fa3":
                return discriminant("D_CP_HadVH")
            if self.analysis == "fa2":
                return discriminant("D_int_HadVH")
            if self.analysis == "fL1":
                return discriminant("D_0hplus_HadVH")

        if self.analysis == "fL1":
            if self.category == "VBF2jTaggedIchep16":
                return discriminant("D_0hplus_VBFdecay")
            if self.category == "VHHadrTaggedIchep16":
                return discriminant("D_0hplus_HadVHdecay")

        for i, prime in product(range(1, 4), ("", "_prime")):
            if self.category == "VBF2jTaggedIchep16" and self.whichproddiscriminants == "D_g1{}gi{}{}".format(i, 4-i, prime):
                if self.analysis == "fa3":
                    return discriminant("D_g1{}_g4{}_VBFdecay{}".format(i, 4-i, prime))
                if self.analysis == "fa2":
                    return discriminant("D_g1{}_g2{}_VBFdecay{}".format(i, 4-i, prime))
                if self.analysis == "fL1":
                    return discriminant("D_g1{}_g1prime2{}_VBFdecay{}".format(i, 4-i, prime))

            if self.category == "VHHadrTaggedIchep16" and self.whichproddiscriminants == "D_g1{}gi{}{}".format(i, 4-i, prime):
                if self.analysis == "fa3":
                    return discriminant("D_g1{}_g4{}_HadZHdecay{}".format(i, 4-i, prime))
                if self.analysis == "fa2":
                    return discriminant("D_g1{}_g2{}_HadZHdecay{}".format(i, 4-i, prime))
                if self.analysis == "fL1":
                    return discriminant("D_g1{}_g1prime2{}_HadZHdecay{}".format(i, 4-i, prime))

        assert False

    @property
    def discriminants(self):
        return (self.purediscriminant, self.mixdiscriminant, self.bkgdiscriminant)

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
        if self.category == "VHHadrTaggedIchep16":
            return "VHhadrtag"
        assert False

    @property
    @cache
    def invertedmatrix(self):
        ganomalous = self.analysis.couplingname
        productionmode = str(self.templategroup).upper()
        basis = SampleBasis([template.hypothesis for template in self.templates()], productionmode, self.analysis)
        invertedmatrix = basis.invertedmatrix
        if self.templategroup in ("vbf", "zh", "wh"):
            """
            basis.matrix looks something like this:
            1,    0,      0,        0,      0
            0,    0,      0,        0,      g4^4
            g1^4, g1^3g4, g1^2g4^2, g1g4^3, g4^4
            g1^4, g1^3g4, g1^2g4^2, g1g4^3, g4^4
            g1^4, g1^3g4, g1^2g4^2, g1g4^3, g4^4

            multiply inverted matrix by the vector of templates.  This should give back
               templates for each respective term (g1^4g4^0, g1^3g4^1, ...)
            In the PDF, the templates need to be multiplied by (g1^i)(g4^(4-i))
            """
            #make sure the two pure templates can be used as is
            #these assertions should be equivalent to asserting that SMtemplate.g1 == anomaloustemplate.ganomalous == 1
            # and that the two pure samples are in the right places on the list
            assert invertedmatrix[0,0] == 1 and all(invertedmatrix[0,i] == 0 for i in range(1,5))
            assert invertedmatrix[4,1] == 1 and invertedmatrix[4,0] == 0 and all(invertedmatrix[4,i] == 0 for i in range(2,5))

        if self.templategroup == "ggh":
            assert invertedmatrix[0,0] == 1 and all(invertedmatrix[0,i] == 0 for i in range(1,3))
            assert invertedmatrix[2,0] == 0 and invertedmatrix[2,2] == 0 #can't assert invertedmatrix[2,1] == 1 because different convention :(

        return invertedmatrix

def listfromiterator(function):
    return list(function())

@listfromiterator
def templatesfiles():
    for shapesystematic in treeshapesystematics:
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
                            elif analysis == "fL1":
                                if w == whichproddiscriminants[0]:
                                    w = None
                                else:
                                    continue
                            else:
                                if w != "D_int_prod":
                                    continue
                            yield TemplatesFile(channel, shapesystematic, "ggh", analysis, production, category, w)
                            yield TemplatesFile(channel, shapesystematic, "vbf", analysis, production, category, w)
                            yield TemplatesFile(channel, shapesystematic, "zh", analysis, production, category, w)
                            yield TemplatesFile(channel, shapesystematic, "wh", analysis, production, category, w)
                            if shapesystematic == "":
                                yield TemplatesFile(channel, "bkg", analysis, production, category, w)
                            for blindstatus in blindstatuses:
                                if shapesystematic == "" and (blindstatus == "blind" or config.unblinddistributions) and config.usedata:
                                    yield TemplatesFile(channel, "DATA", analysis, production, blindstatus, category, w)

class TemplateBase(object):
    __metaclass__ = abc.ABCMeta

    def applysynonyms(self, enumsdict):
        if enumsdict[TemplateGroup] is None:
            if enumsdict[ProductionMode] == "ggH":
                enumsdict[TemplateGroup] = "ggh"
            elif enumsdict[ProductionMode] == "VBF":
                enumsdict[TemplateGroup] = "vbf"
            elif enumsdict[ProductionMode] == "WH":
                enumsdict[TemplateGroup] = "wh"
            elif enumsdict[ProductionMode] == "ZH":
                enumsdict[TemplateGroup] = "zh"
            elif enumsdict[ProductionMode] in ("qqZZ", "ggZZ", "VBF bkg", "ZX"):
                enumsdict[TemplateGroup] = "bkg"
            elif enumsdict[ProductionMode] == "data":
                enumsdict[TemplateGroup] = "DATA"
            elif enumsdict[ProductionMode] is None:
                pass
            else:
                assert False
        super(TemplateBase, self).applysynonyms(enumsdict)

    def templatefile(self, *args, **kwargs):
        return self.templatesfile.templatesfile(*args, **kwargs)

    @abc.abstractmethod
    def templatename(self):
        pass

    def gettemplate(self, *args, **kwargs):
        f = tfiles[self.templatefile(**kwargs)]
        try:
            return getattr(f, self.templatename())
        except AttributeError:
            raise IOError("No template {} in {}".format(self.templatename(), self.templatefile(**kwargs)))

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

smoothingparametersfile = os.path.join(config.repositorydir, "step5_json", "smoothingparameters", "smoothingparameters.json")

class Template(TemplateBase, MultiEnum):
    __metaclass__ = MultiEnumABCMeta
    enums = [TemplatesFile, ProductionMode, Hypothesis]
    enumname = "template"

    def check(self, *args):
        if self.productionmode is None:
            raise ValueError("No option provided for productionmode\n{}".format(args))
        elif self.productionmode in ("ggH", "VBF", "ZH", "WH"):
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if ReweightingSample(self.productionmode, self.hypothesis) not in self.templatesfile.signalsamples():
                raise ValueError("{} {} is not used to make templates for {} {}!\n{}".format(self.productionmode, self.hypothesis, self.templategroup, self.analysis, args))
            if self.templategroup != str(self.productionmode).lower():
                raise ValueError("{} is not {}!\n{}".format(self.productionmode, self.templategroup, args))
        elif self.productionmode in ("ggZZ", "qqZZ", "VBF bkg", "ZX"):
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
        elif self.productionmode in ("VBF", "ZH", "WH"):
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
            elif self.hypothesis in ("fa2proddec-0.5", "fa3proddec-0.5"):
                name = "templateMixProdDecPiAdapSmooth"
            elif self.hypothesis in ("fa2dec-0.5",):
                name = "templateMixDecayPiAdapSmooth"
            elif self.hypothesis in ("fa2prod-0.5",):
                name = "templateMixProdPiAdapSmooth"
            elif self.hypothesis in ("fa2proddec0.5", "fL1proddec0.5"):
                name = "templateMixProdDecAdapSmooth"
        elif self.productionmode == "ZX" and not config.usedata:
            name = "templateqqZZAdapSmooth"
        elif self.productionmode in ("ggZZ", "qqZZ", "VBF bkg", "ZX"):
            name = "template{}AdapSmooth".format(self.productionmode)
        elif self.productionmode == "data":
            name = "datatemplate"

        name

        if self.domirror and final:
            name += "Mirror"

        return name

    def title(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH"):
            return "{} {}".format(self.productionmode, self.hypothesis)
        if self.productionmode == "ggZZ" or self.productionmode == "qqZZ":
            return str(self.productionmode).replace("ZZ", "#rightarrowZZ")
        if self.productionmode == "ZX":
            return "Z+X"
        if self.productionmode == "VBF bkg":
            return "VBF bkg"
        if self.productionmode == "data":
            return "data"
        assert False

    def weightname(self):
        if self.productionmode in ("ggZZ", "VBF bkg"):
            return ReweightingSample(self.productionmode, "2e2mu").weightname()
        if self.hypothesis is not None:
            return ReweightingSample(self.productionmode, self.hypothesis).weightname()
        if self.productionmode == "data":
            return None
        return ReweightingSample(self.productionmode).weightname()

    def reweightfrom(self):
        if self.productionmode == "ggH":
            if self.analysis in ("fa2", "fa3"):
                result={
                        Sample(self.production, "ggH", "0+"),
                        Sample(self.production, "ggH", "a2"),
                        Sample(self.production, "ggH", "0-"),
                        Sample(self.production, "ggH", "L1"),
                        Sample(self.production, "ggH", "fa20.5"),
                        Sample(self.production, "ggH", "fa30.5"),
                        #Sample(self.production, "ggH", "fL10.5"),   #NOT fL1 for now
                       }
            if self.analysis == "fL1":
                if self.hypothesis in ("0+", "L1"):
                    result={
                            Sample(self.production, "ggH", "0+"),
                            Sample(self.production, "ggH", "a2"),
                            Sample(self.production, "ggH", "0-"),
                            Sample(self.production, "ggH", "L1"),
                            Sample(self.production, "ggH", "fa20.5"),
                            Sample(self.production, "ggH", "fa30.5"),
                            #Sample(self.production, "ggH", "fL10.5"),   #NOT fL1 for now
                           }
                elif self.hypothesis == "fL10.5":
                    result={
                            #Sample(self.production, "ggH", "0+"),
                            Sample(self.production, "ggH", "a2"),
                            #Sample(self.production, "ggH", "0-"),
                            Sample(self.production, "ggH", "L1"),
                            Sample(self.production, "ggH", "fa20.5"),
                            Sample(self.production, "ggH", "fa30.5"),
                            Sample(self.production, "ggH", "fL10.5"),
                           }
        if self.productionmode == "VBF":
            if self.hypothesis == "0+":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       }
            if self.hypothesis == "a2":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0-", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5") #0- and L1 are borderline
                       }
            if self.hypothesis == "0-":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0-", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       }
            if self.hypothesis == "L1":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("L1", "fL1prod0.5", "0-", "fa3prod0.5")
                       }
            if self.hypothesis == "fa2dec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0-", "a2", "L1", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5") #0+, 0-, and L1 are borderline
                       }
            if self.hypothesis == "fa3dec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0-", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       }
            if self.hypothesis == "fL1dec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "0-", "a2", "L1", "fa2prod0.5", "fL1prod0.5", "fa3prod0.5")
                       }
            if self.hypothesis == "fa2prod0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       }
            if self.hypothesis == "fa3prod0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")  #0+ and a2 are borderline
                       }
            if self.hypothesis == "fL1prod0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5") #????
                       }
            if self.hypothesis == "fa2proddec-0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       }
            if self.hypothesis == "fa3proddec-0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0-", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       }
            if self.hypothesis == "fL1proddec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "0-", "a2", "L1", "fa2prod0.5", "fL1prod0.5", "fa3prod0.5")
                       }
            if self.hypothesis == "fa2dec-0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "a2", "fa2prod0.5", "fL1prod0.5")
                       }
            if self.hypothesis == "fa2prod-0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "fa3prod0.5", "fa2prod0.5", "fL1prod0.5")
                       }
            if self.hypothesis == "fa2proddec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "0-", "a2", "L1", "fa2prod0.5", "fL1prod0.5", "fa3prod0.5")
                       }
        if self.productionmode == "ZH":
            if self.hypothesis == "0+":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "0-", "a2", "fa2prod0.5", "fa3prod0.5")
                       }
            if self.hypothesis == "a2":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "0-", "a2", "fa2prod0.5", "fa3prod0.5")
                       }
            if self.hypothesis == "0-":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "0-", "a2", "fa2prod0.5", "fa3prod0.5")
                       }
            if self.hypothesis == "L1":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in prodonlyhypotheses
                       }
            if self.hypothesis == "fa2dec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "0-", "a2", "fa2prod0.5", "fa3prod0.5")
                       }
            if self.hypothesis == "fa3dec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "0-", "a2", "fa2prod0.5", "fa3prod0.5")
                       }
            if self.hypothesis == "fL1dec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("a2", "fa2prod0.5", "fL1prod0.5")  #not great
                       }
            if self.hypothesis == "fa2prod0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0-", "a2", "fa2prod0.5", "fa3prod0.5")
                       }
            if self.hypothesis == "fa3prod0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "0-", "a2", "fa2prod0.5", "fa3prod0.5")
                       }
            if self.hypothesis == "fL1prod0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in prodonlyhypotheses
                       }
            if self.hypothesis == "fa2proddec-0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "0-", "a2", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       }
            if self.hypothesis == "fa3proddec-0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "0-", "a2", "fa2prod0.5", "fa3prod0.5")
                       }
            if self.hypothesis == "fL1proddec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0-", "a2", "L1", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                       }
            if self.hypothesis == "fa2dec-0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ()
                       }
            if self.hypothesis == "fa2prod-0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ()
                       }
            if self.hypothesis == "fa2proddec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ()
                       }
        if self.productionmode == "WH":
            if self.hypothesis == "0+":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "a2", "fa2prod0.5", "0-", "fa3prod0.5")
                       }
            if self.hypothesis == "a2":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("a2", "fa2prod0.5", "L1", "a3", "fa3prod0.5")
                       }
            if self.hypothesis == "0-":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0-", "fa3prod0.5", "a2", "fa2prod0.5", "L1")
                       }
            if self.hypothesis == "L1":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "a2", "fa2prod0.5", "L1", "fL1prod0.5", "0-", "fa3prod0.5")
                       }
            if self.hypothesis == "fa2dec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("a2", "fa2prod0.5", "L1", "0-", "fa3prod0.5")
                       }
            if self.hypothesis == "fa3dec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("a2", "fa3prod0.5", "0-", "fa3prod0.5", "L1")
                       }
            if self.hypothesis == "fL1dec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0+", "a2", "fa2prod0.5", "0-", "fa3prod0.5", "L1", "fL1prod0.5",)
                       }
            if self.hypothesis == "fa2prod0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("a3", "fa3prod0.5", "a2", "fa2prod0.5", "L1", "fL1prod0.5")
                       }
            if self.hypothesis == "fa3prod0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0-", "fa3prod0.5", "a2", "fa2prod0.5", "L1")
                       }
            if self.hypothesis == "fL1prod0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("fa2prod0.5", "a3", "fa3prod0.5", "fL1prod0.5", "a2", "L1")
                       }
            if self.hypothesis == "fa2proddec-0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0-", "fa3prod0.5", "a2", "fa2prod0.5", "L1")
                       }
            if self.hypothesis == "fa3proddec-0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0-", "fa3prod0.5", "a2", "fa2prod0.5")
                       }
            if self.hypothesis == "fL1proddec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ("0-", "fa3prod0.5", "a2", "L1", "fL1prod0.5")
                       }
            if self.hypothesis == "fa2dec-0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ()
                       }
            if self.hypothesis == "fa2prod-0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ()
                       }
            if self.hypothesis == "fa2proddec0.5":
                result={
                        Sample(self.production, self.productionmode, hypothesis)
                            for hypothesis in ()
                       }
        if self.productionmode in ("qqZZ", "ZX"):
            result = {Sample(self.production, self.productionmode)}
        if self.productionmode == "ggZZ":
            result = {Sample(self.production, self.productionmode, flavor) for flavor in flavors}
        if self.productionmode == "VBF bkg":
            result = {Sample(self.production, self.productionmode, flavor) for flavor in flavors if not flavor.hastaus}
        if self.productionmode == "data":
            result = {Sample(self.production, self.productionmode, self.blindstatus)}
        result = {sample for sample in result if tfiles[sample.withdiscriminantsfile()].candTree.GetEntries() != 0}
        assert result
        return result

    @property
    def scalefactor(self):
        if self.templategroup in ("bkg", "DATA"): return 1
        if self.productionmode in ("VBF", "ggH", "ZH", "WH"):
            result = ReweightingSample(self.productionmode, self.hypothesis).xsec / ReweightingSample(self.productionmode, "SM").xsec
        result /= len(self.reweightfrom())
        return result

    @property
    def domirror(self):
        if self.analysis != "fa3": return False
        if self.productionmode == "data": return False
        if self.whichproddiscriminants in ("D_g12gi2", "D_g12gi2_prime"): return False

        if self.hypothesis in ("fa30.5", "fa3prod0.5", "fa3proddec-0.5"): return False
        if self.hypothesis in ("0+", "0-"): return True
        if self.productionmode in ("qqZZ", "ggZZ", "VBF bkg", "ZX"): return True
        assert False

    @staticmethod
    def getsmoothingparametersdict(trycache=True):
      import globals
      if globals.smoothingparametersdict_cache is None or not trycache:
        try:
          with open(smoothingparametersfile) as f:
            jsonstring = f.read()
        except IOError:
          try:
            os.makedirs(os.path.dirname(smoothingparametersfile))
          except OSError:
            pass
          with open(smoothingparametersfile, "w") as f:
            f.write("{}\n")
            jsonstring = "{}"
        globals.smoothingparametersdict_cache = json.loads(jsonstring)
      return globals.smoothingparametersdict_cache

    @property
    def smoothingparameters(self):
      keys = (
              str(self.productionmode),
              str(self.category),
              str(self.channel),
              str(self.analysis),
              str(self.whichproddiscriminants),
              str(self.hypothesis),
              str(self.shapesystematic),
              str(self.production),
             )
      return getnesteddictvalue(self.getsmoothingparametersdict(), *keys, default=[None, None, None])

    @property
    def smoothentriesperbin(self):
      return self.smoothingparameters[0]

    @property
    def reweightaxes(self):
      return self.smoothingparameters[1]

    @property
    def reweightrebin(self):
      return self.smoothingparameters[2]

    @property
    def postprocessingjson(self):
      result = []
      if self.smoothentriesperbin:
        result.append({"type": "smooth", "kernel": "adaptive", "entriesperbin": self.smoothentriesperbin})
        if self.reweightaxes:
          reweight = {"type": "reweight", "axes": self.reweightaxes}
          if self.reweightrebin and any(self.reweightrebin):
            reweight.update({"rebinning": self.reweightrebin})
          result.append(reweight)
      result.append({"type": "rescale", "factor": self.scalefactor})
      return result

    def getjson(self):
        jsn = {
               "templates": [
                 {
                   "name": self.templatename(final=False),
                   "files": sorted([os.path.basename(sample.withdiscriminantsfile()) for sample in self.reweightfrom()]),
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
                   "postprocessing": self.postprocessingjson,
                   "filloverflows": True,
                  },
                ],
              }

        if self.domirror:
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
    enumname = "inttemplate"
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
        elif self.productionmode in ("VBF", "ZH", "WH"):
            if self.interferencetype not in ("g11gi3", "g12gi2", "g13gi1"):
                raise ValueError("Invalid interferencetype {} for productionmode {}!\n{}".format(self.interferencetype, self.productionmode, args))
        else:
            raise ValueError("Invalid productionmode {}!\n{}".format(self.productionmode, args))
        super(IntTemplate, self).check(*args)

    def templatename(self):
        if self.productionmode == "ggH":
            if self.interferencetype == "g11gi1":
                result = "templateIntAdapSmooth"
        if self.productionmode in ("VBF", "ZH", "WH"):
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
    def mirrorjsn(self):
        if self.analysis != "fa3": return None
        if self.whichproddiscriminants in ("D_g12gi2", "D_g12gi2_prime"): return None

        #cross talk - production discriminants for the wrong category don't make sense
        if self.category in ("VBFtagged", "VHHadrtagged") and self.productionmode == "ggH":
            if self.whichproddiscriminants == "D_int_prod":
                #ggH has no production information, so mirror symmetric
                return {"type":"mirror", "antisymmetric":False, "axis":1}
            elif self.whichproddiscriminants == "D_int_decay":
                pass #continue to the end of the function
            else:
                #mix of production and decay information - can't mirror either way
                return None

        if (   self.category == "VBFtagged" and self.productionmode in ("ZH", "WH")
            or self.category == "VHHadrtagged" and self.productionmode == "VBF"
           ):
            if self.whichproddiscriminants == "D_int_decay":
                pass
            elif self.interferencetype in ("g11gi3", "g13gi1"):
                return None #Even for D_int_prod can't mirror, there is still CP violation even
                            # though it's not the way this ME expects
            #but for g12gi2 can mirror, since there's no CP violation

        #cross talk to the untagged category is ok, since the decay is the same

        if self.interferencetype in ("g11gi1", "g11gi3", "g13gi1"):
            return {"type":"mirror", "antisymmetric":True, "axis":1}
        elif self.interferencetype == "g12gi2":
            return {"type":"mirror", "antisymmetric":False, "axis":1}
        assert False

    @property
    def domirror(self):
        return bool(self.mirrorjsn)

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
            intjsn["templates"][0]["postprocessing"].append(self.mirrorjsn)
        return intjsn


class DataTree(MultiEnum):
    enums = [Channel, Production, Category]
    enumname = "datatree"
    @property
    def originaltreefile(self):
        return Sample("data", self.production, "unblind").withdiscriminantsfile()
    @property
    def treefile(self):
        return os.path.join(config.repositorydir, "step7_templates", "data_{}_{}_{}.root".format(self.production, self.channel, self.category))
    def passescut(self, t):
        return abs(t.Z1Flav * t.Z2Flav) == self.channel.ZZFlav and config.m4lmin < t.ZZMass < config.m4lmax and config.unblindscans and t.category in self.category

datatrees = []
for channel in channels:
    for production in productions:
        for category in categories:
            datatrees.append(DataTree(channel, production, category))



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
        return os.path.join(config.repositorydir, "step7_templates", "data_{}_{}_{}_{}.root".format(self.production, self.channel, self.category, self.subtractproduction))
    def check(self, *args, **kwargs):
        return super(SubtractDataTree, self).check(*args, **kwargs)
    def passescut(self, t):
        return super(SubtractDataTree, self).passescut(t) and self.subtractproduction.passescut(t)

#for channel in channels:
#    if "160729" in productions:
#        datatrees.append(SubtractDataTree("160729", "subtract160720", channel))
