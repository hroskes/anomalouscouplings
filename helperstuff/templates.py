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
from samples import ReweightingSample, ReweightingSamplePlus, Sample, SampleBasis
from utilities import cache, is_almost_integer, JsonDict, jsonloads, tfiles

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
            self.category = Category("Untagged")

        if self.analysis.isdecayonly:
            if self.category != "Untagged":
                raise ValueError("decay only analysis is only done for untagged!\n{}".format(args))
            if self.templategroup in ("vbf", "zh", "wh", "tth"):
                raise ValueError("decay only analysis is only done with decay information!\n{}".format(args))
        if config.LHE:
            if self.channel != "2e2mu":
                raise ValueError("LHE analysis is only done for 2e2mu for now!\n{}".format(args))

        super(TemplatesFile, self).check(*args, dontcheck=dontcheck)

        if not self.shapesystematic.appliesto(self.templategroup):
            raise ValueError("ShapeSystematic {} does not apply to {}\n{}".format(self.shapesystematic, self.templategroup, args))

    def jsonfile(self, iteration=None):
        folder = os.path.join(config.repositorydir, "step5_json")
        if iteration is not None:
            folder = os.path.join(folder, "bkp_iter{}".format(iteration))

        nameparts = ["templates", self.templategroup, self.analysis, self.channel, self.category, self.shapesystematic, self.production]

        nameparts = [str(x) for x in nameparts if x]
        result = os.path.join(folder, "_".join(x for x in nameparts if x) + ".json")

        return result

    def templatesfile(self, iteration=None, firststep=False):
        if self.copyfromothertemplatesfile is not None:
            return self.copyfromothertemplatesfile.templatesfile()

        folder = os.path.join(config.repositorydir, "step7_templates")
        if iteration is not None:
            folder = os.path.join(folder, "bkp_iter{}".format(iteration))
            if not os.path.exists(folder):
                raise IOError("No folder {}".format(folder))

        nameparts = ["templates", self.templategroup, self.analysis, self.channel, self.category, self.shapesystematic, self.production]
        if firststep: nameparts.append("firststep")

        nameparts = [str(x) for x in nameparts if x]
        result = os.path.join(folder, "_".join(x for x in nameparts if x) + ".root")

        return result

    def controlplotsdir(self, *args, **kwargs):
        relpath = os.path.relpath(self.templatesfile(*args, **kwargs), os.path.join(config.repositorydir, "step7_templates"))
        assert ".." not in relpath
        return os.path.join(config.plotsbasedir, "templateprojections", "controlplots", relpath.replace(".root", "").replace("bkp_", ""))

    def signalsamples(self):
        if self.templategroup == "ggh" and self.shapesystematic == "MINLO_SM":
            return [ReweightingSamplePlus("ggH", "0+", "MINLO")]

        elif self.templategroup in ("ggh", "tth"):
            if self.templategroup == "ggh":
                args = "ggH",
            elif self.templategroup == "tth":
                args = "ttH", "Hff0+"
            else:
                assert False, self.templategroup

            if self.analysis in ("fa3", "fa3_STXS"):
                reweightingsamples = [ReweightingSample("0+", *args), ReweightingSample("0-", *args), ReweightingSample("fa30.5", *args)]
            if self.analysis == "fa2":
                reweightingsamples = [ReweightingSample("0+", *args), ReweightingSample("a2", *args), ReweightingSample("fa2-0.5", *args)]
            if self.analysis == "fL1":
                reweightingsamples = [ReweightingSample("0+", *args), ReweightingSample("L1", *args), ReweightingSample("fL10.5", *args)]
            if self.analysis == "fL1Zg":
                reweightingsamples = [ReweightingSample("0+", *args), ReweightingSample("L1Zg", *args), ReweightingSample("fL1Zg-0.5", *args)]
            if self.analysis.isfL1fL1Zg:
                reweightingsamples = [ReweightingSample("0+", *args), ReweightingSample("L1", *args), ReweightingSample("L1Zg", *args), ReweightingSample("fL10.5", *args), ReweightingSample("fL1Zg0.5", *args), ReweightingSample("fL10.5fL1Zg0.5", *args)]
            if self.analysis == "fa2fa3fL1fL1Zg":
                hypotheses = ["0+", "0-", "a2", "L1", "L1Zg",
                              "fa30.5", "fa2-0.5", "fL10.5", "fL1Zg-0.5",
                              "fa30.5fa20.5", "fa30.5fL10.5", "fa30.5fL1Zg0.5",
                                              "fa20.5fL10.5", "fa20.5fL1Zg0.5",
                                                              "fL10.5fL1Zg0.5",
                             ]
                reweightingsamples = [
                  ReweightingSample(h, *args) for h in hypotheses
                ]

        elif self.templategroup in ("vbf", "zh", "wh"):
            p = str(self.templategroup).upper()
            if self.analysis in ("fa3", "fa3_STXS"):
                reweightingsamples = [ReweightingSample(p, "0+"), ReweightingSample(p, "0-"), ReweightingSample(p, "fa3prod0.5"), ReweightingSample(p, "fa3dec0.5"), ReweightingSample(p, "fa3proddec-0.5")]
            if self.analysis == "fa2":
                reweightingsamples = [ReweightingSample(p, "0+"), ReweightingSample(p, "a2"), ReweightingSample(p, "fa2prod0.5"), ReweightingSample(p, "fa2dec-0.5"), ReweightingSample(p, "fa2proddec-0.5")]
            if self.analysis == "fL1":
                reweightingsamples = [ReweightingSample(p, "0+"), ReweightingSample(p, "L1"), ReweightingSample(p, "fL1prod0.5"), ReweightingSample(p, "fL1dec0.5"), ReweightingSample(p, "fL1proddec-0.5")]
            if self.analysis == "fL1Zg":
                reweightingsamples = [ReweightingSample(p, "0+"), ReweightingSample(p, "L1Zg"), ReweightingSample(p, "fL1Zgprod0.5"), ReweightingSample(p, "fL1Zgdec0.5"), ReweightingSample(p, "fL1Zgproddec-0.5")]
            if self.analysis == "fa2fa3fL1fL1Zg":
                hypotheses = ["0+", "0-", "a2", "L1", "L1Zg",
                              "fa30.5",        "fa2-0.5",       "fL10.5",        "fL1Zg-0.5",
                              "fa3prod0.5",    "fa2prod0.5",    "fL1prod0.5",    "fL1Zgprod0.5",
                              "fa3proddec0.5", "fa2proddec0.5", "fL1proddec0.5", "fL1Zgproddec0.5",

                              "fa30.5fa20.5",                  "fa30.5fL10.5",                  "fa30.5fL1Zg0.5",
                                                               "fa20.5fL10.5",                  "fa20.5fL1Zg0.5",
                                                                                                "fL10.5fL1Zg0.5",
                              "fa3prod0.5fa2prod0.5",          "fa3prod0.5fL1prod0.5",          "fa3prod0.5fL1Zgprod0.5",
                                                               "fa2prod0.5fL1prod0.5",          "fa2prod0.5fL1Zgprod0.5",
                                                                                                "fL1prod0.5fL1Zgprod0.5",
                              "fa3proddec0.5fa2proddec-0.5",   "fa3proddec0.5fL1proddec-0.5",   "fa3proddec0.5fL1Zgproddec-0.5",
                                                               "fa2proddec0.5fL1proddec-0.5",   "fa2proddec0.5fL1Zgproddec-0.5",
                                                                                                "fL1proddec0.5fL1Zgproddec-0.5",

                              "fa30.33fa20.33",                "fa30.33fL10.33",                "fa30.33fL1Zg0.33",
                                                               "fa20.33fL10.33",                "fa20.33fL1Zg0.33",
                                                                                                "fL10.33fL1Zg0.33",
                              "fa3prod0.33fa2prod0.33",        "fa3prod0.33fL1prod0.33",        "fa3prod0.33fL1Zgprod0.33",
                                                               "fa2prod0.33fL1prod0.33",        "fa2prod0.33fL1Zgprod0.33",
                                                                                                "fL1prod0.33fL1Zgprod0.33",
                              "fa3proddec0.33fa2proddec-0.33", "fa3proddec0.33fL1proddec-0.33", "fa3proddec0.33fL1Zgproddec-0.33",
                                                               "fa2proddec0.33fL1proddec-0.33", "fa2proddec0.33fL1Zgproddec-0.33",
                                                                                                "fL1proddec0.33fL1Zgproddec-0.33",

                              "fa30.33fa20.33fL10.33",                         "fa30.33fa20.33fL1Zg0.33",
                              "fa30.33fL10.33fL1Zg0.33",                       "fa20.33fL10.33fL1Zg0.33",
                              "fa3prod0.33fa2prod0.33fL1prod0.33",             "fa3prod0.33fa2prod0.33fL1Zgprod0.33",
                              "fa3prod0.33fL1prod0.33fL1Zgprod0.33",           "fa2prod0.33fL1prod0.33fL1Zgprod0.33",
                              "fa3proddec0.33fa2proddec0.33fL1proddec-0.33",   "fa3proddec0.33fa2proddec0.33fL1Zgproddec-0.33",
                              "fa3proddec0.33fL1proddec0.33fL1Zgproddec-0.33", "fa2proddec0.33fL1proddec0.33fL1Zgproddec-0.33",

                              "fa3proddec0.25fa2proddec0.25fL1proddec0.25",   "fa3proddec0.25fa2proddec0.25fL1Zgproddec0.25",
                              "fa3proddec0.25fL1proddec0.25fL1Zgproddec0.25", "fa2proddec0.25fL1proddec0.25fL1Zgproddec0.25",

                              "fa3proddec0.25fa2proddec0.25fL1proddec0.25fL1Zgproddec0.25",
                             ]
                reweightingsamples = [
                  ReweightingSample(h, p) for h in hypotheses
                ]

        return reweightingsamples

    def templates(self):
        if self.templategroup in ["ggh", "vbf", "zh", "wh", "tth"]:
            return [Template(self, sample) for sample in self.signalsamples()]
        elif self.templategroup == "bkg":
            if config.LHE: return [Template(self, "qqZZ")]
            result = ["qqZZ", "ggZZ"]
            if config.usedata:
                result.append("ZX")
            result.append("VBF bkg")
            return [Template(self, productionmode) for productionmode in result]
        elif self.templategroup == "DATA":
            return [Template(self, "data")]
        assert False

    def inttemplates(self):
        if self.templategroup == "ggh" and self.shapesystematic != "MINLO_SM" or self.templategroup == "tth":
            h = str(self.templategroup).replace("h", "H")
            if self.analysis.dimensions == 4:
                return [
                    "g{}1g{}1".format(i1, i2)
                        for i1, i2 in itertools.combinations("1ijkl", 2):
                ]
            elif self.analysis.dimensions == 2:
                return [IntTemplate(self, h, _) for _ in ("g11gi1", "g11gj1", "gi1gj1")]
            elif self.analysis.dimensions == 1:
                return [IntTemplate(self, h, "g11gi1")]

        elif self.templategroup in ("vbf", "zh", "wh"):
            h = str(self.templategroup).upper()
            if self.analysis.dimensions == 4:
                return [
                    "g{}3g{}1".format(*indices)
                        for indices in itertools.combinations("1ijkl", 2)
                ] + [
                    "g{}2g{}2".format(*indices)
                        for indices in itertools.combinations("1ijkl", 2)
                ] + [
                    "g{}1g{}3".format(*indices)
                        for indices in itertools.combinations("1ijkl", 2)
                ] + [
                    "g{}2g{}1g{}1".format(*indices)
                        for indices in itertools.combinations("1ijkl", 3)
                ] + [
                    "g{}1g{}2g{}1".format(*indices)
                        for indices in itertools.combinations("1ijkl", 3)
                ] + [
                    "g{}1g{}1g{}2".format(*indices)
                        for indices in itertools.combinations("1ijkl", 3)
                ] + [
                    "g{}1g{}1g{}1g{}1".format(*indices)
                        for indices in itertools.combinations("1ijkl", 4)
                ]
            elif self.analysis.dimensions == 2:
                assert False
            elif self.analysis.dimensions == 1:
                return [IntTemplate(self, h, "g1{}gi{}".format(i, 4-i)) for i in (1, 2, 3)]

        elif self.templategroup in ("bkg", "DATA") or self.shapesystematic == "MINLO_SM":
            return []
        assert False, self

    @property
    def bkgdiscriminant(self):
        return self.shapesystematic.D_bkg()

    @property
    def purediscriminant(self):
        from discriminants import discriminant

        if self.shapesystematic in ("JECUp", "JECDn"):
            JECappend = "_{}".format(self.shapesystematic)
        else:
            JECappend = ""

        if self.category == "Untagged":
            if self.analysis == "fa3":
                return discriminant("D_0minus_decay")
            if self.analysis == "fa2":
                return discriminant("D_0hplus_decay")
            if self.analysis == "fL1":
                return discriminant("D_L1_decay")
            if self.analysis == "fL1Zg":
                return discriminant("D_L1Zg_decay")
            if self.analysis == "fL1fL1Zg_DL1_DL1L1Zgint":
                return discriminant("D_L1_decay")
            if self.analysis == "fL1fL1Zg_DL1_DL1Zgint":
                return discriminant("D_L1_decay")
            if self.analysis == "fL1fL1Zg_DeR_DeLeR":
                return discriminant("D_eR_decay")
            if self.analysis == "fL1fL1Zg_DeR_DeLint":
                return discriminant("D_eR_decay")
            if self.analysis == "fL1fL1Zg_m1_m2":
                return discriminant("Z1Mass")
            if self.analysis == "fL1fL1Zg_m1_phi":
                return discriminant("Z1Mass")
            if self.analysis == "fL1fL1Zg_m2_phi":
                return discriminant("Z2Mass")
            if self.analysis == "fa3_STXS":
                return discriminant("D_STXS_ggH_stage1"+JECappend)
            if self.analysis == "fa3fa2fL1fL1Zg":
                return disciminant("D_4couplings_decay")

        if self.category == "VBFtagged":
            if self.analysis == "fa3":
                return discriminant("D_0minus_VBFdecay"+JECappend)
            if self.analysis == "fa2":
                return discriminant("D_0hplus_VBFdecay"+JECappend)
            if self.analysis == "fL1":
                return discriminant("D_L1_VBFdecay"+JECappend)
            if self.analysis == "fL1Zg":
                return discriminant("D_L1Zg_VBFdecay"+JECappend)
            if self.analysis == "fa3_STXS":
                return discriminant("D_STXS_VBF_stage1"+JECappend)
            if self.analysis == "fa3fa2fL1fL1Zg":
                return disciminant("D_4couplings_VBFdecay"+JECappend)

        if self.category == "VHHadrtagged":
            if self.analysis == "fa3":
                return discriminant("D_0minus_HadVHdecay"+JECappend)
            if self.analysis == "fa2":
                return discriminant("D_0hplus_HadVHdecay"+JECappend)
            if self.analysis == "fL1":
                return discriminant("D_L1_HadVHdecay"+JECappend)
            if self.analysis == "fL1Zg":
                return discriminant("D_L1Zg_HadVHdecay"+JECappend)
            if self.analysis == "fa3_STXS":
                return discriminant("D_STXS_VBF_stage1"+JECappend)
            if self.analysis == "fa3fa2fL1fL1Zg":
                return disciminant("D_4couplings_HadVHdecay"+JECappend)

        assert False

    @property
    def mixdiscriminant(self):
        from discriminants import discriminant

        if self.category == "Untagged":
            if self.analysis == "fa3":
                return discriminant("D_CP_decay")
            if self.analysis == "fa2":
                return discriminant("D_int_decay")
            if self.analysis == "fL1":
                return discriminant("D_0hplus_decay")
            if self.analysis == "fL1Zg":
                return discriminant("D_0hplus_decay")
            if self.analysis == "fL1fL1Zg_DL1_DL1L1Zgint":
                return discriminant("D_L1L1Zgint_decay")
            if self.analysis == "fL1fL1Zg_DL1_DL1Zgint":
                return discriminant("D_L1Zgint_decay")
            if self.analysis == "fL1fL1Zg_DeR_DeLeR":
                return discriminant("D_eLeR_decay")
            if self.analysis == "fL1fL1Zg_DeR_DeLint":
                return discriminant("D_eLint_decay")
            if self.analysis == "fL1fL1Zg_m1_m2":
                return discriminant("Z2Mass")
            if self.analysis == "fL1fL1Zg_m1_phi":
                return discriminant("Phi")
            if self.analysis == "fL1fL1Zg_m2_phi":
                return discriminant("Phi")
            if self.analysis == "fa3_STXS":
                return discriminant("phistarZ2")
            if self.analysis == "fa3fa2fL1fL1Zg":
                return discriminant("D_CP_decay_2bins")

        if self.shapesystematic in ("JECUp", "JECDn"):
            JECappend = "_{}".format(self.shapesystematic)
        else:
            JECappend = ""

        if self.category == "VBFtagged":
            if self.analysis == "fa3":
                return discriminant("D_CP_VBF"+JECappend)
            if self.analysis == "fa2":
                return discriminant("D_int_VBF"+JECappend)
            if self.analysis == "fL1":
                return discriminant("D_0hplus_VBFdecay"+JECappend)
            if self.analysis == "fL1Zg":
                return discriminant("D_0hplus_VBFdecay"+JECappend)
            if self.analysis == "fa3_STXS":
                return discriminant("phistarZ2")
            if self.analysis == "fa3fa2fL1fL1Zg":
                return discriminant("D_CP_VBFdecay_2bins")

        if self.category == "VHHadrtagged":
            if self.analysis == "fa3":
                return discriminant("D_CP_HadVH"+JECappend)
            if self.analysis == "fa2":
                return discriminant("D_int_HadVH"+JECappend)
            if self.analysis == "fL1":
                return discriminant("D_0hplus_HadVHdecay"+JECappend)
            if self.analysis == "fL1Zg":
                return discriminant("D_0hplus_HadVHdecay"+JECappend)
            if self.analysis == "fa3_STXS":
                return discriminant("phistarZ2")
            if self.analysis == "fa3fa2fL1fL1Zg":
                return discriminant("D_CP_HadVHdecay_2bins")

        assert False

    @property
    def discriminants(self):
        if self.shapesystematic == "MINLO_SM": return (self.purediscriminant, self.mixdiscriminant)
        return (self.purediscriminant, self.mixdiscriminant, self.bkgdiscriminant)

    @property
    @cache
    def invertedmatrix(self):
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
            if self.analysis.dimensions == 1:
                #make sure the two pure templates can be used as is
                #these assertions should be equivalent to asserting that SMtemplate.g1 == anomaloustemplate.ganomalous == 1
                # and that the two pure samples are in the right places on the list
                assert invertedmatrix[0,0] == 1 and all(invertedmatrix[0,i] == 0 for i in range(1,5))
                assert invertedmatrix[4,1] == 1 and invertedmatrix[4,0] == 0 and all(invertedmatrix[4,i] == 0 for i in range(2,5))
            elif self.analysis.dimensions == 4:
                assert invertedmatrix[0, 0] == 1 and all(invertedmatrix[0, i] == 0 for i in             range(1, 70))
                assert invertedmatrix[4, 1] == 1 and all(invertedmatrix[4, i] == 0 for i in range(0, 1)+range(2, 70))
                assert invertedmatrix[14,2] == 1 and all(invertedmatrix[14,i] == 0 for i in range(0, 2)+range(3, 70))
                assert invertedmatrix[34,3] == 1 and all(invertedmatrix[34,i] == 0 for i in range(0, 3)+range(4, 70))
                assert invertedmatrix[69,4] == 1 and all(invertedmatrix[69,i] == 0 for i in range(0, 4)+range(5, 70))
            else:
                assert False

        if self.templategroup in ("ggh", "tth"):
            if self.analysis.dimensions == 4:
                assert invertedmatrix[0, 0] == 1 and all(invertedmatrix[0, i] == 0 for i in             range(1, 70))
                assert                               all(invertedmatrix[2, i] == 0 for i in range(0, 1)+range(2, 70))
                assert                               all(invertedmatrix[5, i] == 0 for i in range(0, 2)+range(3, 70))
                assert                               all(invertedmatrix[9, i] == 0 for i in range(0, 3)+range(4, 70))
                assert                               all(invertedmatrix[14,i] == 0 for i in range(0, 4)+range(5, 70))
            elif self.analysis.dimensions == 2:
                assert invertedmatrix[0,0] == 1 and all(invertedmatrix[0,i] == 0 for i in range(1,6))
                assert all(invertedmatrix[1,i] == 0 for i in (2, 4, 5))
                assert all(invertedmatrix[2,i] == 0 for i in (0, 2, 3, 4, 5))
                assert all(invertedmatrix[3,i] == 0 for i in (1, 3, 5))
                assert all(invertedmatrix[4,i] == 0 for i in (0, 3, 4))
                assert all(invertedmatrix[5,i] == 0 for i in (0, 1, 3, 4, 5))
                assert invertedmatrix[0,0] == 1 and all(invertedmatrix[0,i] == 0 for i in range(1,6))
            elif self.analysis.dimensions == 1:
                assert invertedmatrix[0,0] == 1 and all(invertedmatrix[0,i] == 0 for i in range(1,3))
                assert invertedmatrix[2,0] == 0 and invertedmatrix[2,2] == 0 #can't assert invertedmatrix[2,1] == 1 because different convention :(
            else:
                assert False

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
        return result

    def docustomsmoothing(self):
        if not self.hascustomsmoothing: return
        newf = ROOT.TFile(self.templatesfile(), "RECREATE")
        oldf = ROOT.TFile(self.templatesfile(firststep=True))
        newf.cd()

        controlplotsdir = newf.mkdir("controlPlots")
        controlplotsdir.cd()
        for key in oldf.controlPlots.GetListOfKeys():
            key.ReadObj().Write()

        for template in self.templates()+self.inttemplates():
            template.docustomsmoothing(newf, controlplotsdir)

    @property
    def copyfromothertemplatesfile(self):
        if self.production == "170203" and self.templategroup == "DATA":
          kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in self.enums}
          kwargs["production"] = "170222"
          return TemplatesFile(*kwargs.values())
        if self.production == "170222" and self.templategroup != "DATA":
          kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in self.enums}
          kwargs["production"] = "170203"
          return TemplatesFile(*kwargs.values())
        if self.production == "170712":
          kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in self.enums}
          kwargs["production"] = "170222"
          return TemplatesFile(*kwargs.values())
        return None

def listfromiterator(function):
    return list(function())

@listfromiterator
def templatesfiles():
    for channel in channels:
        if channel != "2e2mu" and config.LHE: continue
        for production in productions:
            for analysis in analyses:
                for category in categories:
                    if category != "Untagged" and analysis.isdecayonly: continue
                    for shapesystematic in treeshapesystematics:
                        if category != "Untagged" and shapesystematic in ("ScaleUp", "ScaleDown", "ResUp", "ResDown"): continue
                        if config.LHE and shapesystematic in ("ScaleUp", "ScaleDown", "ResUp", "ResDown"): continue
                        yield TemplatesFile(channel, shapesystematic, "ggh", analysis, production, category)
                        if shapesystematic in ("ScaleUp", "ScaleDown", "ResUp", "ResDown"): continue
                        yield TemplatesFile(channel, shapesystematic, "bkg", analysis, production, category)
                        if analysis.isdecayonly: continue
                        yield TemplatesFile(channel, shapesystematic, "vbf", analysis, production, category)
                        yield TemplatesFile(channel, shapesystematic, "zh", analysis, production, category)
                        yield TemplatesFile(channel, shapesystematic, "wh", analysis, production, category)
                        yield TemplatesFile(channel, shapesystematic, "tth", analysis, production, category)
                    if config.showblinddistributions:
                        yield TemplatesFile(channel, "DATA", analysis, production, category)
                    if category != "Untagged" and config.applyMINLOsystematics:
                        yield TemplatesFile(channel, "ggh", analysis, production, category, "MINLO_SM")


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
        result = sum(([d.bins, float(d.min), float(d.max)] for d in self.discriminants), [])
        return result

    @abc.abstractmethod
    def getjson(self):
        pass

    @abc.abstractproperty
    def customsmoothingkwargs(self):
        pass

    @abc.abstractmethod
    def hascustomsmoothing(self):
        pass

    def docustomsmoothing(self, newf, controlplotsdir):
        f = ROOT.TFile(self.templatesfile.templatesfile(firststep=True))
        h = getattr(f, self.templatename(final=False))
        rawprojections = [f.Get("controlPlots/control_{}_projAxis{}_afterNormalization".format(self.templatename(final=False), i)).GetListOfPrimitives()[0]
                              for i in range(3)]
        try:
            customsmoothing.customsmoothing(h, rawprojections, newf, controlplotsdir, **self.customsmoothingkwargs)
        except:
            print "Error while smoothing {}:".format(self.templatename(final=False))
            raise
        if self.domirror and not isinstance(self, IntTemplate):
            assert "axes" not in self.customsmoothingkwargs or 1 not in self.customsmoothingkwargs["axes"]
            hmirror = getattr(f, self.templatename(final=True))
            try:
                customsmoothing.customsmoothing(hmirror, rawprojections, newf, controlplotsdir, **self.customsmoothingkwargs)
            except:
                print "Error while custom smoothing {}:".format(self.templatename(final=True))
                raise

class Template(TemplateBase, MultiEnum):
    __metaclass__ = MultiEnumABCMeta
    enums = [TemplatesFile, ReweightingSamplePlus]
    enumname = "template"

    def __init__(self, *args, **kwargs):
        super(Template, self).__init__(*args, **kwargs)
        self.__smoothingparameters = None

    def initsmoothingparameters(self):
        if self.__smoothingparameters is None:
            self.__smoothingparameters = SmoothingParameters(self)

    def check(self, *args):
        if self.productionmode is None:
            raise ValueError("No option provided for productionmode\n{}".format(args))

        elif self.productionmode in ("ggH", "VBF", "ZH", "WH", "ttH"):
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.reweightingsampleplus not in self.templatesfile.signalsamples():
                print self.templatesfile.signalsamples()
                raise ValueError("{} is not used to make templates for {} {}!\n{}".format(self.reweightingsampleplus, self.templategroup, self.analysis, args))
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
            if self.hypothesis == "0+":
                name = "template0PlusAdapSmooth"
            elif self.hypothesis == "0-":
                name = "template0MinusAdapSmooth"
            elif self.hypothesis == "a2":
                name = "template0HPlusAdapSmooth"
            elif self.hypothesis == "L1":
                name = "template0L1AdapSmooth"
            elif self.hypothesis == "L1Zg":
                name = "template0L1ZgAdapSmooth"
            elif self.hypothesis == "fa20.5":
                name = "templateMixa1a2AdapSmooth"
            elif self.hypothesis == "fa30.5":
                name = "templateMixa1a3AdapSmooth"
            elif self.hypothesis == "fL10.5":
                name = "templateMixa1L1AdapSmooth"
            elif self.hypothesis == "fL1Zg0.5":
                name = "templateMixa1L1ZgAdapSmooth"
            elif self.hypothesis == "fa2-0.5":
                name = "templateMixa1a2PiAdapSmooth"
            elif self.hypothesis == "fa3-0.5":
                name = "templateMixa1a3PiAdapSmooth"
            elif self.hypothesis == "fL1-0.5":
                name = "templateMixa1L1PiAdapSmooth"
            elif self.hypothesis == "fL1Zg-0.5":
                name = "templateMixa1L1ZgPiAdapSmooth"
            elif self.hypothesis == "fL10.5fL1Zg0.5":
                name = "templateMixL1L1ZgAdapSmooth"
        elif self.productionmode in ("VBF", "ZH", "WH"):
            if self.hypothesis == "0+":
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
        if self.productionmode == "data":
            return None
        return self.reweightingsampleplus.weightname()

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
            if config.LHE:
                result = {Sample(self.production, self.productionmode, self.hypothesis)}
            elif self.shapesystematic == "MINLO_SM":
                result = {Sample(self.production, self.productionmode, self.hypothesis, "MINLO")}
            else:
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
        if self.productionmode == "ZX":
            result = {Sample(self.production, self.productionmode)}
        if self.productionmode == "qqZZ":
            result = {Sample(self.production, self.productionmode), Sample(self.production, self.productionmode, "ext")}
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
        if self.shapesystematic == "MINLO_SM":
            result = len(self.reweightfrom())
        elif self.productionmode in ("VBF", "ggH", "ZH", "WH"):
            result = ReweightingSample(self.productionmode, self.hypothesis).xsec / ReweightingSample(self.productionmode, "SM").xsec
        elif self.productionmode == "ttH":
            result = ReweightingSample(self.productionmode, self.hypothesis, self.hffhypothesis).xsec / ReweightingSample(self.productionmode, "SM", "Hff0+").xsec
        elif self.productionmode in ("VBF bkg", "ggZZ", "ZX", "data"):
            result = len(self.reweightfrom())
        elif self.productionmode == "qqZZ":
            result = 1
        result /= sum(
                      Sample.effectiveentries(
                                              reweightfrom=reweightfrom,
                                              reweightto=ReweightingSample(self.productionmode, self.hypothesis, self.hffhypothesis)
                                             )
                       for reweightfrom in self.reweightfrom()
                     )
        if isnan(result): result = 0
        if self.productionmode in ("VBF bkg", "ggZZ", "ZX", "data"):
            assert result == 1
        return result

    @property
    def domirror(self):
        if self.analysis not in ("fa3", "fa3_STXS"): return False
        if self.analysis == "fa3_STXS": return False #for now... could do it later
        if self.productionmode == "data": return False

        if self.hypothesis in ("fa30.5", "fa3prod0.5", "fa3proddec-0.5"): return False
        if self.hypothesis in ("0+", "0-"): return True
        if self.productionmode in ("qqZZ", "ggZZ", "VBF bkg", "ZX"): return True
        assert False

    @property
    def smoothingparameters(self):
      self.initsmoothingparameters()
      return self.__smoothingparameters.value

    @smoothingparameters.setter
    def smoothingparameters(self, value):
      self.initsmoothingparameters()
      self.__smoothingparameters.value = value

    @staticmethod
    def writedict():
      SmoothingParameters.writedict()

    @property
    def smoothentriesperbin(self):
      return self.smoothingparameters[0][0]

    @smoothentriesperbin.setter
    def smoothentriesperbin(self, value):
      self.smoothingparameters[0][0] = value

    @property
    def reweightaxes(self):
      return self.smoothingparameters[0][1]

    @reweightaxes.setter
    def reweightaxes(self, value):
      self.smoothingparameters[0][1] = value

    @property
    def reweightrebin(self):
      result = self.smoothingparameters[0][2]
      #validation
      if result is not None:
        if len(result) != len(self.discriminants):
          raise ValueError("len(reweightrebin) for {!r} != {}!\n{}".format(self, len(self.discriminants), reweightrebin))
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

    @reweightrebin.setter
    def reweightrebin(self, value):
        self.smoothingparameters[0][2] = value

    @property
    def customsmoothingkwargs(self):
      result = self.smoothingparameters[1].copy()
      if "othertemplateargs" in result:
        result["othertemplate"] = Template(*result["othertemplateargs"])
        del result["othertemplateargs"]
      return result

    @customsmoothingkwargs.setter
    def customsmoothingkwargs(self, value):
      result = self.smoothingparameters
      result[1] = value
      self.smoothingparameters = result

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
        result = ["ZZMass>{}".format(config.m4lmin), "ZZMass<{}".format(config.m4lmax)]
        if not config.LHE:
            result.append("Z1Flav*Z2Flav == {}".format(self.ZZFlav))
        if not self.analysis.isdecayonly:
            result.append("(" + " || ".join("{} == {}".format(self.categoryname, c) for c in self.category.idnumbers) + ")")
        return " && ".join(result)

    @property
    def ZZFlav(self):
        assert not config.LHE
        result = self.channel.ZZFlav
        if self.productionmode == "ZX": result *= -1
        return result

    @property
    def sample(self):
        return ReweightingSample(self.hypothesis, self.productionmode)

class SmoothingParameters(MultiEnum, JsonDict):
      __metaclass__ = MultiEnumABCMeta
      enums = Template.enums
      enumname = "smoothingparameters"
      dictfile = os.path.join(config.repositorydir, "step5_json", "smoothingparameters", "smoothingparameters.json")
      def check(self, *args):
          Template(*args)
      @property
      def keys(self):
        return (
                str(self.productionmode),
                str(self.category),
                str(self.channel),
                str(self.analysis),
                str(self.hypothesis),
                str(self.shapesystematic),
                str(self.production.productionforsmoothingparameters),
               )
      @property
      def default(self):
        return [None, None, None]

      @property
      def rawvalue(self):
        return super(SmoothingParameters, self).getvalue()

      def getvalue(self):
        result = self.rawvalue

        if self.shapesystematic != "" and result == [None, None, None]:
          kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in type(self).needenums}
          kwargs["shapesystematic"] = ""
          return type(self)(*kwargs.values()).getvalue()  #can't use actual kwargs

        if len(result) == 3:
          result = [result, None]
        if result[1] is None:
          result[1] = {}
        return result

class IntTemplate(TemplateBase, MultiEnum):
    __metaclass__ = MultiEnumABCMeta
    enumname = "inttemplate"
    class InterferenceType(MyEnum):
        enumname = "interferencetype"
        enumitems = (
                     EnumItem("g11gi1"),
                     EnumItem("g11gj1"),
                     EnumItem("gi1gj1"),
                     EnumItem("g11gi3"),
                     EnumItem("g12gi2"),
                     EnumItem("g13gi1"),
                    )
    enums = [TemplatesFile, ProductionMode, InterferenceType, HffHypothesis]

    def check(self, *args):

        dontcheck = []

        if self.interferencetype in ("g11gj1", "gi1gj1") and self.analysis.dimensions != 2:
            raise ValueError("Invalid interferencetype {} for 1D analysis\n{}".format(self.interferencetype, args))

        if self.productionmode in ("ggH", "ttH"):
            if self.interferencetype not in ("g11gi1", "g11gj1", "gi1gj1"):
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

    def templatename(self, final=True):
        if self.productionmode in ("ggH", "ttH"):
            if self.interferencetype == "g11gi1":
                if self.analysis in ("fa3", "fa3_STXS"):
                    result = "templatea1a3IntAdapSmooth"
                elif self.analysis == "fa2":
                    result = "templatea1a2IntAdapSmooth"
                elif self.analysis == "fL1" or self.analysis.isfL1fL1Zg:
                    result = "templatea1L1IntAdapSmooth"
                elif self.analysis == "fL1Zg":
                    result = "templatea1L1ZgIntAdapSmooth"
            elif self.interferencetype == "g11gj1":
                if self.analysis.isfL1fL1Zg:
                    result = "templatea1L1ZgIntAdapSmooth"
            elif self.interferencetype == "gi1gj1":
                if self.analysis.isfL1fL1Zg:
                    result = "templateL1L1ZgIntAdapSmooth"
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
        if self.analysis not in ("fa3", "fa3_STXS"): return None

        #cross talk - production discriminants for the wrong category don't make sense
        if self.category in ("VBFtagged", "VHHadrtagged") and self.productionmode in ("ggH", "ttH"):
            assert self.interferencetype == "g11gi1"
            #ggH has no production information, and only using SM ttH, so mirror antisymmetric
            #over the (pretend) D_CP_decay axis, which sets the whole thing to 0
            #note for ttH this is an approximation, since we could have H(0-)->2l2q tt->bbllnunu
            return {"type":"rescale", "factor":0}

        if self.analysis == "fa3_STXS" and self.interferencetype in ("g11gi1", "g11gi3", "g13gi1"):
            #same (antimirror over D_CP_whatever, which doesn't exist in STXS)
            return {"type":"rescale", "factor":0}

        if self.analysis == "fa3_STXS": return None #for now... could do it later

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

    @property
    @cache
    def templatesandfactors(self):
        if self.interferencetype == "g11gi1":
            g1exp = giexp = 1
            rowofinvertedmatrix = giexp #first row is labeled 0
            hffhypothesis = "Hff0+" if self.productionmode == "ttH" else None
            multiplyby = getattr(ReweightingSample(self.productionmode, self.analysis.purehypotheses[1], hffhypothesis), self.analysis.purehypotheses[1].couplingname)

        if self.interferencetype == "g11gj1":
            assert self.analysis.dimensions == 2
            rowofinvertedmatrix = 3 #first row is labeled 0
            hffhypothesis = "Hff0+" if self.productionmode == "ttH" else None
            multiplyby = getattr(ReweightingSample(self.productionmode, self.analysis.purehypotheses[2], hffhypothesis), self.analysis.purehypotheses[2].couplingname)

        if self.interferencetype == "gi1gj1":
            assert self.analysis.dimensions == 2
            rowofinvertedmatrix = 4 #first row is labeled 0
            hffhypothesis = "Hff0+" if self.productionmode == "ttH" else None
            multiplyby = (
                          getattr(ReweightingSample(self.productionmode, self.analysis.purehypotheses[1], hffhypothesis), self.analysis.purehypotheses[1].couplingname)
                         *
                          getattr(ReweightingSample(self.productionmode, self.analysis.purehypotheses[2], hffhypothesis), self.analysis.purehypotheses[2].couplingname)
                         )

        if self.interferencetype in ("g11gi3", "g12gi2", "g13gi1"):
            g1exp, giexp = (int(i) for i in str(self.interferencetype).replace("g1", "").replace("gi", ""))
            rowofinvertedmatrix = giexp #first row is labeled 0
            multiplyby = 1

        invertedmatrix = self.templatesfile.invertedmatrix
        vectoroftemplates = self.templatesfile.templates()  #ok technically it's a list not a vector
        templatesandfactors = []

        for j, template in enumerate(vectoroftemplates):
            templatesandfactors.append((template, invertedmatrix[rowofinvertedmatrix,j] * multiplyby))

        return templatesandfactors

    def getjson(self):
        templatesum = [{
                        "name": template.templatename(final=False),
                        "factor": factor,
                       } for template, factor in self.templatesandfactors]
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
        return any(template.hascustomsmoothing for template, factor in self.templatesandfactors if factor)

    @property
    def customsmoothingkwargs(self):
        return {
                "name": "redointerference",
                "newf": self.customsmoothing_newf,
                "templatesandfactors": self.templatesandfactors,
                "mirrorjsn": self.mirrorjsn,
               }

    def docustomsmoothing(self, newf, controlplotsdir):
        self.customsmoothing_newf = newf
        try:
            return super(IntTemplate, self).docustomsmoothing(newf, controlplotsdir)
        finally:
            del self.customsmoothing_newf

    def gettemplate(self, *args, **kwargs):
        result = super(IntTemplate, self).gettemplate(*args, **kwargs)
        if self.analysis == "fa3" and self.productionmode in ("ggH", "ttH") and self.interferencetype == "g11gi1" and self.category in ("VBFtagged", "VHHadrtagged"):
            assert self.domirror
            result.Scale(0)
        return result

class DataTree(MultiEnum):
    enums = [Channel, Production, Category, Analysis]
    enumname = "datatree"
    def check(self, *args):
        super(DataTree, self).check(*args)
        if self.analysis.isdecayonly:
            if self.category != "Untagged":
                raise ValueError("decayonly analysis is only done for untagged!\n{}".format(args))
        if config.LHE:
            if self.channel != "2e2mu":
                raise ValueError("LHE analysis is only done for 2e2mu for now!\n{}".format(args))

    @property
    def originaltreefile(self):
        return Sample("data", self.production).withdiscriminantsfile()
    @property
    def treefile(self):
        return os.path.join(config.repositorydir, "step7_templates", "data_{}_{}_{}_{}.root".format(self.production, self.channel, self.category, self.analysis))
    def passescut(self, t):
        return abs(t.Z1Flav * t.Z2Flav) == self.channel.ZZFlav and config.m4lmin < t.ZZMass < config.m4lmax and config.unblinddistributions and getattr(t, "category_"+self.analysis.categoryname) in self.category

@listfromiterator
def datatrees():
    for channel in channels:
        if channel != "2e2mu" and config.LHE: continue
        for production in productions:
            for category in categories:
                for analysis in analyses:
                    if category != "Untagged" and analysis.isdecayonly: continue
                    yield DataTree(channel, production, category, analysis)



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
