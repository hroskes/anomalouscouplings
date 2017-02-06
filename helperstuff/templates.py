import abc
from itertools import product
import json
from math import isnan
import os
import ROOT

import numpy

import config
import customsmoothing
from enums import Analysis, analyses, Channel, channels, Category, categories, EnumItem, flavors, HffHypothesis, Hypothesis, MultiEnum, MultiEnumABCMeta, MyEnum, prodonlyhypotheses, Production, ProductionMode, productions, ShapeSystematic, TemplateGroup, treeshapesystematics
from samples import ReweightingSample, Sample, SampleBasis
from utilities import cache, getnesteddictvalue, is_almost_integer, jsonloads, tfiles

class TemplatesFile(MultiEnum):
    enumname = "templatesfile"
    enums = [Channel, ShapeSystematic, TemplateGroup, Analysis, Production, Category]

    def applysynonyms(self, enumsdict):
        if enumsdict[Production] is None and len(config.productionsforcombine) == 1:
            enumsdict[Production] = config.productionsforcombine[0]
        super(TemplatesFile, self).applysynonyms(enumsdict)

    def check(self, *args):
        dontcheck = []

        if self.shapesystematic is None:
            self.shapesystematic = ShapeSystematic("")

        if self.category is None:
            self.category = Category("UntaggedIchep16")

        super(TemplatesFile, self).check(*args, dontcheck=dontcheck)

        if not self.shapesystematic.appliesto(self.templategroup):
            raise ValueError("ShapeSystematic {} does not apply to {}\n{}".format(self.shapesystematic, self.templategroup, args))

    def jsonfile(self, iteration=None):
        folder = os.path.join(config.repositorydir, "step5_json")
        if iteration is not None:
            folder = os.path.join(folder, "bkp_iter{}".format(iteration))

        nameparts = ["templates", self.templategroup, self.analysis, self.channel, self.categorynamepart, self.shapesystematic, self.production]

        nameparts = [str(x) for x in nameparts if x]
        result = os.path.join(folder, "_".join(x for x in nameparts if x) + ".json")

        return result

    def templatesfile(self, iteration=None, firststep=False):
        folder = os.path.join(config.repositorydir, "step7_templates")
        if iteration is not None:
            folder = os.path.join(folder, "bkp_iter{}".format(iteration))
            if not os.path.exists(folder):
                raise IOError("No folder {}".format(folder))

        nameparts = ["templates", self.templategroup, self.analysis, self.channel, self.categorynamepart, self.shapesystematic if config.applyshapesystematics else "", self.production]
        if firststep: nameparts.append("firststep")

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
            if self.analysis == "fL1Zg":
                reweightingsamples = [ReweightingSample("ggH", "0+_photoncut"), ReweightingSample("ggH", "L1Zg"), ReweightingSample("ggH", "fL1Zg-0.5")]

        elif self.templategroup == "vbf":
            if self.analysis == "fa3":
                reweightingsamples = [ReweightingSample("VBF", "0+"), ReweightingSample("VBF", "0-"), ReweightingSample("VBF", "fa3prod0.5"), ReweightingSample("VBF", "fa3dec0.5"), ReweightingSample("VBF", "fa3proddec-0.5")]
            if self.analysis == "fa2":
                reweightingsamples = [ReweightingSample("VBF", "0+"), ReweightingSample("VBF", "a2"), ReweightingSample("VBF", "fa2prod0.5"), ReweightingSample("VBF", "fa2dec-0.5"), ReweightingSample("VBF", "fa2proddec-0.5")]
            if self.analysis == "fL1":
                reweightingsamples = [ReweightingSample("VBF", "0+"), ReweightingSample("VBF", "L1"), ReweightingSample("VBF", "fL1prod0.5"), ReweightingSample("VBF", "fL1dec0.5"), ReweightingSample("VBF", "fL1proddec0.5")]
            if self.analysis == "fL1Zg":
                reweightingsamples = [ReweightingSample("VBF", "0+_photoncut"), ReweightingSample("VBF", "L1Zg"), ReweightingSample("VBF", "fL1Zgprod0.5"), ReweightingSample("VBF", "fL1Zgdec0.5"), ReweightingSample("VBF", "fL1Zgproddec-0.5")]

        elif self.templategroup == "zh":
            if self.analysis == "fa3":
                reweightingsamples = [ReweightingSample("ZH", "0+"), ReweightingSample("ZH", "0-"), ReweightingSample("ZH", "fa3prod0.5"), ReweightingSample("ZH", "fa3dec0.5"), ReweightingSample("ZH", "fa3proddec-0.5")]
            if self.analysis == "fa2":
                reweightingsamples = [ReweightingSample("ZH", "0+"), ReweightingSample("ZH", "a2"), ReweightingSample("ZH", "fa2prod0.5"), ReweightingSample("ZH", "fa2dec-0.5"), ReweightingSample("ZH", "fa2proddec-0.5")]
            if self.analysis == "fL1":
                reweightingsamples = [ReweightingSample("ZH", "0+"), ReweightingSample("ZH", "L1"), ReweightingSample("ZH", "fL1prod0.5"), ReweightingSample("ZH", "fL1dec0.5"), ReweightingSample("ZH", "fL1proddec0.5")]
            if self.analysis == "fL1Zg":
                reweightingsamples = [ReweightingSample("ZH", "0+_photoncut"), ReweightingSample("ZH", "L1Zg"), ReweightingSample("ZH", "fL1Zgprod0.5"), ReweightingSample("ZH", "fL1Zgdec0.5"), ReweightingSample("ZH", "fL1Zgproddec-0.5")]

        elif self.templategroup == "wh":
            if self.analysis == "fa3":
                reweightingsamples = [ReweightingSample("WH", "0+"), ReweightingSample("WH", "0-"), ReweightingSample("WH", "fa3prod0.5"), ReweightingSample("WH", "fa3dec0.5"), ReweightingSample("WH", "fa3proddec-0.5")]
            if self.analysis == "fa2":
                reweightingsamples = [ReweightingSample("WH", "0+"), ReweightingSample("WH", "a2"), ReweightingSample("WH", "fa2prod0.5"), ReweightingSample("WH", "fa2dec-0.5"), ReweightingSample("WH", "fa2proddec-0.5")]
            if self.analysis == "fL1":
                reweightingsamples = [ReweightingSample("WH", "0+"), ReweightingSample("WH", "L1"), ReweightingSample("WH", "fL1prod0.5"), ReweightingSample("WH", "fL1dec0.5"), ReweightingSample("WH", "fL1proddec0.5")]
            if self.analysis == "fL1Zg":
                reweightingsamples = [ReweightingSample("WH", "0+_photoncut"), ReweightingSample("WH", "L1Zg"), ReweightingSample("WH", "fL1Zgprod0.5"), ReweightingSample("WH", "fL1Zgdec0.5"), ReweightingSample("WH", "fL1Zgproddec-0.5")]

        elif self.templategroup == "tth":
            if self.analysis == "fa3":
                reweightingsamples = [ReweightingSample("ttH", "0+", "Hff0+"), ReweightingSample("ttH", "0-", "Hff0+"), ReweightingSample("ttH", "fa30.5", "Hff0+")]
            if self.analysis == "fa2":
                reweightingsamples = [ReweightingSample("ttH", "0+", "Hff0+"), ReweightingSample("ttH", "a2", "Hff0+"), ReweightingSample("ttH", "fa20.5", "Hff0+")]
            if self.analysis == "fL1":
                reweightingsamples = [ReweightingSample("ttH", "0+", "Hff0+"), ReweightingSample("ttH", "L1", "Hff0+"), ReweightingSample("ttH", "fL10.5", "Hff0+")]
            if self.analysis == "fL1Zg":
                reweightingsamples = [ReweightingSample("ttH", "0+_photoncut", "Hff0+"), ReweightingSample("ttH", "L1Zg", "Hff0+"), ReweightingSample("ttH", "fL1Zg0.5", "Hff0+")]

        return reweightingsamples

    def templates(self):
        if self.templategroup in ["ggh", "vbf", "zh", "wh", "tth"]:
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
        elif self.templategroup == "tth":
            return [IntTemplate(self, "ttH", "g11gi1")]
        elif self.templategroup in ("bkg", "DATA"):
            return []
        assert False

    @property
    def bkgdiscriminant(self):
        return self.shapesystematic.D_bkg()

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
            if self.analysis == "fL1Zg":
                return discriminant("D_L1Zg_decay")

        if self.category == "VBF2jTaggedIchep16":
            if self.analysis == "fa3":
                return discriminant("D_0minus_VBFdecay")
            if self.analysis == "fa2":
                return discriminant("D_0hplus_VBFdecay")
            if self.analysis == "fL1":
                return discriminant("D_L1_VBFdecay")
            if self.analysis == "fL1Zg":
                return discriminant("D_L1Zg_VBFdecay")

        if self.category == "VHHadrTaggedIchep16":
            if self.analysis == "fa3":
                return discriminant("D_0minus_HadVHdecay")
            if self.analysis == "fa2":
                return discriminant("D_0hplus_HadVHdecay")
            if self.analysis == "fL1":
                return discriminant("D_L1_HadVHdecay")
            if self.analysis == "fL1Zg":
                return discriminant("D_L1Zg_HadVHdecay")

        assert False

    @property
    def mixdiscriminant(self):
        from discriminants import discriminant

        if self.category == "UntaggedIchep16":
            if self.analysis == "fa3":
                return discriminant("D_CP_decay")
            if self.analysis == "fa2":
                return discriminant("D_int_decay")
            if self.analysis == "fL1":
                return discriminant("D_0hplus_decay")
            if self.analysis == "fL1Zg":
                return discriminant("D_0hplus_decay")

        if self.category == "VBF2jTaggedIchep16":
            if self.analysis == "fa3":
                return discriminant("D_CP_VBF")
            if self.analysis == "fa2":
                return discriminant("D_int_VBF")
            if self.analysis == "fL1":
                return discriminant("D_0hplus_VBFdecay")
            if self.analysis == "fL1Zg":
                return discriminant("D_0hplus_VBFdecay")

        if self.category == "VHHadrTaggedIchep16":
            if self.analysis == "fa3":
                return discriminant("D_CP_HadVH")
            if self.analysis == "fa2":
                return discriminant("D_int_HadVH")
            if self.analysis == "fL1":
                return discriminant("D_0hplus_HadVHdecay")
            if self.analysis == "fL1Zg":
                return discriminant("D_0hplus_HadVHdecay")

        assert False

    @property
    def discriminants(self):
        return (self.purediscriminant, self.mixdiscriminant, self.bkgdiscriminant)

    @property
    def categorynamepart(self):
        if self.category == "UntaggedIchep16":
            return "Untagged"
        if self.category == "VBF2jTaggedIchep16":
            return "VBFtag"
        if self.category == "VHHadrTaggedIchep16":
            return "VHhadrtag"
        assert False

    @property
    @cache
    def invertedmatrix(self):
        ganomalous = self.analysis.couplingname
        productionmode = str(self.templategroup).upper().replace("GGH", "ggH").replace("TTH", "ttH")
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

        if self.templategroup in ("ggh", "tth"):
            assert invertedmatrix[0,0] == 1 and all(invertedmatrix[0,i] == 0 for i in range(1,3))
            assert invertedmatrix[2,0] == 0 and invertedmatrix[2,2] == 0 #can't assert invertedmatrix[2,1] == 1 because different convention :(

        return invertedmatrix

    def getjson(self):
        return {
                "inputDirectory": "step3_withdiscriminants/",
                "outputFile": self.templatesfile(firststep=self.hascustomsmoothing),
                "templates": sum((_.getjson() for _ in self.templates()+self.inttemplates()), []),
               }

    @property
    def hascustomsmoothing(self):
        result = any(_.hascustomsmoothing for _ in self.templates() + self.inttemplates())
        if self.inttemplates():
            assert not result, self
        return result

    def docustomsmoothing(self):
        if not self.hascustomsmoothing: return
        newf = ROOT.TFile(self.templatesfile(), "RECREATE")
        oldf = ROOT.TFile(self.templatesfile(firststep=True))
        newf.cd()

        thetemplates = {template: {getattr(oldf, template.templatename(final=False)), getattr(oldf, template.templatename(final=True))} for template in self.templates()+self.inttemplates()}

        controlplotsdir = newf.mkdir("controlPlots")
        controlplotsdir.cd()
        for key in oldf.controlPlots.GetListOfKeys():
            key.ReadObj().Write()

        for template in self.templates()+self.inttemplates():
            newcontrolplots = template.docustomsmoothing(newf, controlplotsdir)

def listfromiterator(function):
    return list(function())

@listfromiterator
def templatesfiles():
    for shapesystematic in treeshapesystematics:
        for channel in channels:
            for production in productions:
                for analysis in analyses:
                    for category in categories:
                        yield TemplatesFile(channel, shapesystematic, "ggh", analysis, production, category)
                        yield TemplatesFile(channel, shapesystematic, "vbf", analysis, production, category)
                        yield TemplatesFile(channel, shapesystematic, "zh", analysis, production, category)
                        yield TemplatesFile(channel, shapesystematic, "wh", analysis, production, category)
                        yield TemplatesFile(channel, shapesystematic, "tth", analysis, production, category)
                        if shapesystematic == "":
                            yield TemplatesFile(channel, "bkg", analysis, production, category)
                        if shapesystematic == "" and config.usedata:
                            yield TemplatesFile(channel, "DATA", analysis, production, category)

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
            elif enumsdict[ProductionMode] == "ttH":
                enumsdict[TemplateGroup] = "tth"
            elif enumsdict[ProductionMode] in ("qqZZ", "ggZZ", "VBF bkg", "ZX"):
                enumsdict[TemplateGroup] = "bkg"
            elif enumsdict[ProductionMode] == "data":
                enumsdict[TemplateGroup] = "DATA"
            elif enumsdict[ProductionMode] is None:
                pass
            else:
                assert False

        if enumsdict[HffHypothesis] is None:
            if enumsdict[ProductionMode] == "ttH":
                enumsdict[HffHypothesis] = "Hff0+"

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
    enums = [TemplatesFile, ProductionMode, Hypothesis, HffHypothesis]
    enumname = "template"

    def check(self, *args):
        if self.productionmode is None:
            raise ValueError("No option provided for productionmode\n{}".format(args))

        elif self.productionmode in ("ggH", "VBF", "ZH", "WH", "ttH"):
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if ReweightingSample(self.productionmode, self.hypothesis, self.hffhypothesis) not in self.templatesfile.signalsamples():
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

        if self.productionmode == "ttH":
            if self.hffhypothesis != "Hff0+":
                raise ValueError("{} is not used to make templates {}!\n{}".format(self.hffhypothesis, args))
        else:
            if self.hffhypothesis is not None:
               raise ValueError("HffHypothesis {} provided for {}!\n{}".format(self.hffhypothesis, self.productionmode, args))

    def templatename(self, final=True):
        if self.productionmode in ("ggH", "ttH"):
            if self.hypothesis in ("0+", "0+_photoncut"):
                name = "template0PlusAdapSmooth"
            elif self.hypothesis == "0-":
                name = "template0MinusAdapSmooth"
            elif self.hypothesis == "a2":
                name = "template0HPlusAdapSmooth"
            elif self.hypothesis == "L1":
                name = "template0L1AdapSmooth"
            elif self.hypothesis == "L1Zg":
                name = "template0L1ZgAdapSmooth"
            elif self.hypothesis in ("fa20.5", "fa30.5", "fL10.5", "fL1Zg0.5"):
                name = "templateMixAdapSmooth"
            elif self.hypothesis in ("fa2-0.5", "fa3-0.5", "fL1-0.5", "fL1Zg-0.5"):
                name = "templateMixPiAdapSmooth"
        elif self.productionmode in ("VBF", "ZH", "WH"):
            if self.hypothesis in ("0+", "0+_photoncut"):
                name = "template0PlusAdapSmooth"
            elif self.hypothesis == "0-":
                name = "template0MinusAdapSmooth"
            elif self.hypothesis == "a2":
                name = "template0HPlusAdapSmooth"
            elif self.hypothesis == "L1":
                name = "template0L1AdapSmooth"
            elif self.hypothesis == "L1Zg":
                name = "template0L1ZgAdapSmooth"
            elif self.hypothesis in ("fa2dec0.5", "fa3dec0.5", "fL1dec0.5", "fL1Zgdec0.5"):
                name = "templateMixDecayAdapSmooth"
            elif self.hypothesis in ("fa2prod0.5", "fa3prod0.5", "fL1prod0.5", "fL1Zgprod0.5"):
                name = "templateMixProdAdapSmooth"
            elif self.hypothesis in ("fa2proddec-0.5", "fa3proddec-0.5", "fL1proddec-0.5", "fL1Zgproddec-0.5"):
                name = "templateMixProdDecPiAdapSmooth"
            elif self.hypothesis in ("fa2dec-0.5", "fa3dec-0.5", "fL1dec-0.5", "fL1Zgdec-0.5"):
                name = "templateMixDecayPiAdapSmooth"
            elif self.hypothesis in ("fa2prod-0.5", "fa3prod-0.5", "fL1prod-0.5", "fL1Zgprod-0.5"):
                name = "templateMixProdPiAdapSmooth"
            elif self.hypothesis in ("fa2proddec0.5", "fa3proddec0.5", "fL1proddec0.5", "fL1Zgproddec0.5"):
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
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "ttH"):
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
        if self.productionmode == "ttH":
            return ReweightingSample(self.productionmode, self.hypothesis, self.hffhypothesis).weightname()
        if self.hypothesis is not None:
            return ReweightingSample(self.productionmode, self.hypothesis).weightname()
        if self.productionmode == "data":
            return None
        return ReweightingSample(self.productionmode).weightname()

    @property
    def categoryname(self):
        result = "category_"
        result += self.analysis.categoryname

        from treewrapper import TreeWrapper
        for categorization in TreeWrapper.categorizations:
            if categorization.category_function_name == result: return result
        assert False, "{} does not exist in TreeWrapper".format(result)

    def reweightfrom(self):
        if self.productionmode == "ggH":
            result={
                    Sample(self.production, self.productionmode, hypothesis)
                        for hypothesis in ("0+", "0-", "a2", "L1", "fa2dec0.5", "fa3dec0.5", "fL1dec0.5")
                   }
        if self.productionmode in ["VBF", "ZH", "WH"]:
            result={
                    Sample(self.production, self.productionmode, hypothesis)
                        for hypothesis in ("0+", "0-", "a2", "L1", "fa2prod0.5", "fa3prod0.5", "fL1prod0.5")
                   }
        if self.productionmode == "ttH":
            result={
                    Sample(self.production, self.productionmode, "0+", self.hffhypothesis)
                   }
        if self.productionmode in ("qqZZ", "ZX"):
            result = {Sample(self.production, self.productionmode)}
        if self.productionmode == "ggZZ":
            result = {Sample(self.production, self.productionmode, flavor) for flavor in flavors}
        if self.productionmode == "VBF bkg":
            result = {Sample(self.production, self.productionmode, flavor) for flavor in flavors if not flavor.hastaus}
        if self.productionmode == "data":
            result = {Sample(self.production, self.productionmode)}
        result = {sample for sample in result if tfiles[sample.withdiscriminantsfile()].candTree.GetEntries() != 0}
        assert result
        return result

    @property
    def scalefactor(self):
        if self.templategroup in ("bkg", "DATA"): return 1
        if self.productionmode in ("VBF", "ggH", "ZH", "WH"):
            result = ReweightingSample(self.productionmode, self.hypothesis).xsec / ReweightingSample(self.productionmode, "SM").xsec
        if self.productionmode == "ttH":
            result = ReweightingSample(self.productionmode, self.hypothesis, self.hffhypothesis).xsec / ReweightingSample(self.productionmode, "SM", "Hff0+").xsec
        result /= sum(
                      Sample.effectiveentries(
                                              reweightfrom=reweightfrom,
                                              reweightto=ReweightingSample(self.productionmode, self.hypothesis, self.hffhypothesis)
                                             )
                       for reweightfrom in self.reweightfrom()
                     )
        if isnan(result): result = 0
        return result

    @property
    def domirror(self):
        if self.analysis != "fa3": return False
        if self.productionmode == "data": return False

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

    @classmethod
    def writesmoothingparametersdict(cls):
      dct = cls.getsmoothingparametersdict()
      jsonstring = json.dumps(dct, sort_keys=True, indent=4, separators=(',', ': '))
      with open(smoothingparametersfile, "w") as f:
        f.write(jsonstring)

    @property
    def smoothingparameters(self):
      keys = (
              str(self.productionmode),
              str(self.category),
              str(self.channel),
              str(self.analysis),
              str(self.hypothesis),
              str(self.shapesystematic),
              str(self.production),
             )
      result = getnesteddictvalue(self.getsmoothingparametersdict(), *keys, default=[None, None, None])
      if len(result) == 3:
          result = [result, None]
      if result[1] is None:
          result[1] = {}
      return result

    @property
    def smoothentriesperbin(self):
      return self.smoothingparameters[0][0]

    @property
    def reweightaxes(self):
      return self.smoothingparameters[0][1]

    @property
    def reweightrebin(self):
      result = self.smoothingparameters[0][2]
      #validation
      if result is not None:
        if len(result) != 3:
          raise ValueError("len(reweightrebin) for {!r} != 3!\n{}".format(self, reweightrebin))
        for axis, (_, disc) in enumerate(zip(result, self.discriminants)):
          if _ is not None:
            if _ != sorted(_):
              raise ValueError("reweightrebin for {!r} axis {} is not sorted!\n{}".format(self, axis, _))
            if len(_) != len(set(_)):
              raise ValueError("reweightrebin for {!r} axis {} has duplicates!\n{}".format(self, axis, _))
            if min(_) != disc.min:
              raise ValueError("first entry {} of reweightrebin for {!r} axis {} is not the same as the discriminant minimum {}!\n{}".format(min(_), self, axis, disc.min, _))
            if max(_) != disc.max:
              raise ValueError("last entry {} of reweightrebin for {!r} axis {} is not the same as the discriminant maximum {}!\n{}".format(max(_), self, axis, disc.max, _))
            for bin_ in _:
              shouldbeint = (bin_-disc.min) * disc.bins / (disc.max-disc.min)
              if not is_almost_integer(shouldbeint):
                raise ValueError("({bin}-{min}) * {bins} / ({max}-{min}) = {!r} in reweightrebin for {!r} axis {} is not an integer\n{}".format(shouldbeint, self, axis, _, bin=bin_, min=disc.min, max=disc.max, bins=disc.bins))
      return result

    @property
    def customsmoothingkwargs(self):
      return self.smoothingparameters[1]

    @property
    def hascustomsmoothing(self):
        return bool(self.customsmoothingkwargs)

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
        jsn = [
               {
                "name": self.templatename(final=False),
                "files": sorted([os.path.basename(sample.withdiscriminantsfile()) for sample in self.reweightfrom()]),
                "tree": "candTree",
                "variables": [d.name for d in self.discriminants],
                "weight": self.weightname(),
                "selection": self.selection,
                "binning": {
                  "type": "fixed",
                  "bins": self.binning,
                },
                "conserveSumOfWeights": True,
                "postprocessing": self.postprocessingjson,
                "filloverflows": True,
               },
              ]

        if self.domirror:
            mirrorjsn = [
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
                        ]
            jsn += mirrorjsn
        else:
            jsn[0]["postprocessing"].append({"type": "floor"})

        if self.productionmode == "data":
            del jsn[0]["postprocessing"]
            del jsn[0]["weight"]

        return jsn

    @property
    def selection(self):
        result = "ZZMass>{} && ZZMass<{} && Z1Flav*Z2Flav == {}".format(config.m4lmin, config.m4lmax, self.ZZFlav)
        result += " && (" + " || ".join("{} == {}".format(self.categoryname, c) for c in self.category.idnumbers) + ")"
        return result

    @property
    def ZZFlav(self):
        result = self.channel.ZZFlav
        if self.productionmode == "ZX": result *= -1
        return result

    @property
    def sample(self):
        return ReweightingSample(self.hypothesis, self.productionmode)

    def docustomsmoothing(self, newf, controlplotsdir):
        f = ROOT.TFile(self.templatesfile.templatesfile(firststep=True))
        h = getattr(f, self.templatename(final=False))
        rawprojections = [f.Get("controlPlots/control_{}_projAxis{}_afterFill".format(self.templatename(final=False), i)).GetListOfPrimitives()[0]
                              for i in range(3)]
        try:
            customsmoothing.customsmoothing(h, rawprojections, newf, controlplotsdir, **self.customsmoothingkwargs)
        except:
            print "Error while smoothing {}:".format(self.templatename(final=False))
            raise
        if self.domirror:
            assert "axes" not in self.customsmoothingkwargs or 1 not in self.customsmoothingkwargs["axes"]
            hmirror = getattr(f, self.templatename(final=True))
            try:
                customsmoothing.customsmoothing(hmirror, rawprojections, newf, controlplotsdir, **self.customsmoothingkwargs)
            except:
                print "Error while custom smoothing {}:".format(self.templatename(final=True))
                raise

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
    enums = [TemplatesFile, ProductionMode, InterferenceType, HffHypothesis]

    def check(self, *args):

        dontcheck = []

        if self.productionmode in ("ggH", "ttH"):
            if self.interferencetype != "g11gi1":
                raise ValueError("Invalid interferencetype {} for productionmode {}!\n{}".format(self.interferencetype, self.productionmode, args))
        elif self.productionmode in ("VBF", "ZH", "WH"):
            if self.interferencetype not in ("g11gi3", "g12gi2", "g13gi1"):
                raise ValueError("Invalid interferencetype {} for productionmode {}!\n{}".format(self.interferencetype, self.productionmode, args))
        else:
            raise ValueError("Invalid productionmode {}!\n{}".format(self.productionmode, args))

        if self.productionmode == "ttH":
            if self.hffhypothesis != "Hff0+":
                raise ValueError("{} is not used to make templates {}!\n{}".format(self.hffhypothesis, args))
        else:
            if self.hffhypothesis is not None:
               raise ValueError("HffHypothesis {} provided for {}!\n{}".format(self.hffhypothesis, self.productionmode, args))
            dontcheck.append(HffHypothesis)

        super(IntTemplate, self).check(*args, dontcheck=dontcheck)

    def templatename(self):
        if self.productionmode in ("ggH", "ttH"):
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

        #cross talk - production discriminants for the wrong category don't make sense
        if self.category in ("VBFtagged", "VHHadrtagged") and self.productionmode in ("ggH", "ttH"):
            #ggH has no production information, and only using SM ttH, so mirror symmetric
            #note for ttH this is an approximation, since we could have H(0-)->2l2q tt->bbllnunu
            return {"type":"mirror", "antisymmetric":False, "axis":1}

        #Mirror antisymmetric for VH in VBF category and VBF in VH category
        #the interference templates are 0 to within error bars anyway,
        #except for some effects in ZH g13g41 D_CP_VBF which are antisymmetric

        #cross talk to the untagged category is exactly correct, since the decay is the same

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
            g1exp = giexp = 1
            hffhypothesis = "Hff0+" if self.productionmode == "ttH" else None
            multiplyby = getattr(ReweightingSample(self.productionmode, self.analysis.purehypotheses[1], hffhypothesis), self.analysis.couplingname)
        if self.interferencetype in ("g11gi3", "g12gi2", "g13gi1"):
            g1exp, giexp = (int(i) for i in str(self.interferencetype).replace("g1", "").replace("gi", ""))
            multiplyby = 1

        invertedmatrix = self.templatesfile.invertedmatrix
        vectoroftemplates = self.templatesfile.templates()  #ok technically it's a list not a vector
        rowofinvertedmatrix = giexp  #first row is labeled 0
        templatesum = []
        for j, template in enumerate(vectoroftemplates):
            templatesum.append({
                                "name": template.templatename(final=False),
                                "factor": invertedmatrix[rowofinvertedmatrix,j] * multiplyby,
                               })
        intjsn = [
                  {
                   "name": self.templatename(),
                   "templatesum":templatesum,
                   "postprocessing":[],
                  },
                 ]

        if self.domirror:
            intjsn[0]["postprocessing"].append(self.mirrorjsn)
        return intjsn

    @property
    def hascustomsmoothing(self):
        return False

    def docustomsmoothing(self):
        if not self.hascustomsmoothing: return
        assert False

class DataTree(MultiEnum):
    enums = [Channel, Production, Category, Analysis]
    enumname = "datatree"
    @property
    def originaltreefile(self):
        return Sample("data", self.production).withdiscriminantsfile()
    @property
    def treefile(self):
        return os.path.join(config.repositorydir, "step7_templates", "data_{}_{}_{}_{}.root".format(self.production, self.channel, self.category, self.analysis))
    def passescut(self, t):
        return abs(t.Z1Flav * t.Z2Flav) == self.channel.ZZFlav and config.m4lmin < t.ZZMass < config.m4lmax and config.unblindscans and getattr(t, self.analysis.categoryname) in self.category

datatrees = []
for channel in channels:
    for production in productions:
        for category in categories:
            for analysis in analyses:
                datatrees.append(DataTree(channel, production, category, analysis))



class SubtractProduction(MyEnum):
    enumname = "subtractproduction"
    enumitems = (
                 EnumItem("subtract160720"),
                )
    subtracttree = None
    def passescut(self, t):
        if self.subtracttree is None:
            self.subtracttree = tfiles[Sample("data", str(self).replace("subtract", "")).withdiscriminantsfile()].candTree
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
