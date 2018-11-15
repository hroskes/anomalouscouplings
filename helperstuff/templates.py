#!/usr/bin/env python

import abc
from collections import Counter
from itertools import combinations, izip, product
import json
from math import isnan
import os
import re
import uncertainties

import numpy

import ROOT

import config
import customsmoothing
from enums import Analysis, analyses, Channel, channels, Category, categories, EnumItem, flavors, HffHypothesis, Hypothesis, MultiEnum, MultiEnumABCMeta, MyEnum, prodonlyhypotheses, Production, ProductionMode, productions, ShapeSystematic, shapesystematics, TemplateGroup, treeshapesystematics
from samples import ReweightingSample, ReweightingSamplePlus, ReweightingSampleWithPdf, Sample, SampleBasis
from utilities import cache, is_almost_integer, JsonDict, jsonloads, TFile, tfiles

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
        if self.production.LHE:
            if self.channel != "2e2mu":
                raise ValueError("LHE analysis is only done for 2e2mu for now!\n{}".format(args))

        super(TemplatesFile, self).check(*args, dontcheck=dontcheck)

        if self.shapesystematic not in self.allshapesystematics:
            raise ValueError("ShapeSystematic {} does not apply to {!r}\n{}".format(self.shapesystematic, self.nominal, args))

    def jsonfile(self, iteration=None):
        assert "Ulascan" not in str(self.production)
        folder = os.path.join(config.repositorydir, "step5_json", str(self.production))
        if iteration is not None:
            folder = os.path.join(folder, "bkp_iter{}".format(iteration))

        nameparts = ["templates", self.templategroup, self.analysis, self.channel, self.category, self.shapesystematic, self.production]

        nameparts = [str(x) for x in nameparts if x]
        result = os.path.join(folder, "_".join(x for x in nameparts if x) + ".json")

        return result

    def templatesfile(self, iteration=None, firststep=False, splitinterference=False):
        if self.copyfromothertemplatesfile is not None:
            return self.copyfromothertemplatesfile.templatesfile()

        folder = os.path.join(config.repositorydir, "step7_templates")
        if iteration is not None:
            folder = os.path.join(folder, "bkp_iter{}".format(iteration))
            if not os.path.exists(folder):
                raise IOError("No folder {}".format(folder))

        nameparts = ["templates", self.templategroup, self.analysis, self.channel, self.category, self.shapesystematic, self.production]
        if firststep: nameparts.append("firststep")
        if splitinterference: nameparts.append("splitint")

        nameparts = [str(x) for x in nameparts if x]
        result = os.path.join(folder, "_".join(x for x in nameparts if x) + ".root")

        return result

    def controlplotsdir(self, *args, **kwargs):
        relpath = os.path.relpath(self.templatesfile(*args, **kwargs), os.path.join(config.repositorydir, "step7_templates"))
        assert ".." not in relpath
        return os.path.join(config.plotsbasedir, "templateprojections", "controlplots", relpath.replace(".root", "").replace("bkp_", ""))

    @cache
    def signalsamples(self):
        if self.templategroup == "ggh" and self.shapesystematic == "MINLO_SM":
            return [ReweightingSamplePlus("ggH", "0+", "MINLO")]

        elif self.templategroup in ("ggh", "tth", "bbh"):
            if self.templategroup == "ggh":
                args = "ggH",
            elif self.templategroup == "tth":
                args = "ttH", "Hff0+"
            elif self.templategroup == "bbh":
                args = "bbH",
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
            if self.analysis == "fa3fa2fL1fL1Zg":
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
            if self.analysis == "fa3fa2fL1fL1Zg":
                hypotheses = ["0+", "0-", "a2", "L1", "L1Zg",
                              "fa30.5",         "fa2-0.5",        "fL10.5",         "fL1Zg-0.5",
                              "fa3prod0.5",     "fa2prod0.5",     "fL1prod0.5",     "fL1Zgprod0.5",
                              "fa3proddec-0.5", "fa2proddec-0.5", "fL1proddec-0.5", "fL1Zgproddec-0.5",

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
        if self.templategroup in ["ggh", "vbf", "zh", "wh", "tth", "bbh"]:
            return [Template(self, sample) for sample in self.signalsamples()]
        elif self.templategroup == "bkg":
            if self.production.LHE: return [Template(self, "qqZZ")]
            result = ["qqZZ", "ggZZ", "VBF bkg"]
            if config.usedata:
                result.append("ZX")
            if self.shapesystematic in ("ZXUp", "ZXDown"):
                result = [_ for _ in result if _ == "ZX"]
            return [Template(self, productionmode) for productionmode in result]
        elif self.templategroup == "DATA":
            return [Template(self, "data")]
        assert False

    def inttemplates(self):
        if self.templategroup == "ggh" and self.shapesystematic != "MINLO_SM" or self.templategroup in ("tth", "bbh"):
            h = str(self.templategroup).replace("h", "H")
            if self.analysis.dimensions == 4:
                return [IntTemplate(self, h, "g{}1g{}1".format(i1, i2))
                        for i1, i2 in combinations("1ijkl", 2)
                ]
            elif self.analysis.dimensions == 2:
                return [IntTemplate(self, h, _) for _ in ("g11gi1", "g11gj1", "gi1gj1")]
            elif self.analysis.dimensions == 1:
                return [IntTemplate(self, h, "g11gi1")]

        elif self.templategroup in ("vbf", "zh", "wh"):
            h = str(self.templategroup).upper()
            if self.analysis.dimensions == 4:
                return [IntTemplate(self, h, _) for _ in [
                    "g{}3g{}1".format(*indices)
                        for indices in combinations("1ijkl", 2)
                ] + [
                    "g{}2g{}2".format(*indices)
                        for indices in combinations("1ijkl", 2)
                ] + [
                    "g{}1g{}3".format(*indices)
                        for indices in combinations("1ijkl", 2)
                ] + [
                    "g{}2g{}1g{}1".format(*indices)
                        for indices in combinations("1ijkl", 3)
                ] + [
                    "g{}1g{}2g{}1".format(*indices)
                        for indices in combinations("1ijkl", 3)
                ] + [
                    "g{}1g{}1g{}2".format(*indices)
                        for indices in combinations("1ijkl", 3)
                ] + [
                    "g{}1g{}1g{}1g{}1".format(*indices)
                        for indices in combinations("1ijkl", 4)
                ]]
            elif self.analysis.dimensions == 2:
                assert False
            elif self.analysis.dimensions == 1:
                return [IntTemplate(self, h, "g1{}gi{}".format(i, 4-i)) for i in (1, 2, 3)]

        elif self.templategroup in ("bkg", "DATA") or self.shapesystematic == "MINLO_SM":
            return []
        assert False, self

    @property
    def bkgdiscriminant(self):
        from discriminants import discriminant

        name = "D_bkg"
        if self.production >= "180416" or self.production == "180224_newdiscriminants":
            if self.category == "Untagged": pass
            elif self.category == "VBFtagged": name += "_VBFdecay"
            elif self.category == "VHHadrtagged": name += "_HadVHdecay"
            else: assert False
        elif self.production.year == 2016: pass
        else: assert False

        name += self.shapesystematic.Dbkgappendname

        if self.analysis == "fa3fa2fL1fL1Zg":
            name += "_10bins"
        elif self.production >= "180416" or self.production in ("180224_newdiscriminants", "180224_10bins"):
            if self.category == "Untagged":
                name += "_20bins"
            else:
                name += "_10bins"

        if self.shapesystematic in ("JECUp", "JECDn") and ("VBF" in name or "HadVH" in name):
            name += "_{}".format(self.shapesystematic)

        return discriminant(name)

    @property
    def purediscriminant(self):
        from discriminants import discriminant

        if self.shapesystematic in ("JECUp", "JECDn"):
            JECappend = "_{}".format(self.shapesystematic)
        else:
            JECappend = ""

        if self.production >= "180416" or self.production in ("180224_10bins", "180224_newdiscriminants"):
            if self.category == "Untagged":
                binsappend = "_20bins"
            else:
                binsappend = "_10bins"
        elif self.production.year == 2016:
            binsappend = ""

        if self.category == "Untagged":
            if self.analysis == "fa3":
                return discriminant("D_0minus_decay"+binsappend)
            if self.analysis == "fa2":
                return discriminant("D_0hplus_decay"+binsappend)
            if self.analysis == "fL1":
                return discriminant("D_L1_decay"+binsappend)
            if self.analysis == "fL1Zg":
                return discriminant("D_L1Zg_decay"+binsappend)
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
                return discriminant("D_4couplings_decay")

        if self.category == "VBFtagged":
            if self.analysis == "fa3":
                return discriminant("D_0minus_VBFdecay"+binsappend+JECappend)
            if self.analysis == "fa2":
                return discriminant("D_0hplus_VBFdecay"+binsappend+JECappend)
            if self.analysis == "fL1":
                return discriminant("D_L1_VBFdecay"+binsappend+JECappend)
            if self.analysis == "fL1Zg":
                return discriminant("D_L1Zg_VBFdecay"+binsappend+JECappend)
            if self.analysis == "fa3_STXS":
                return discriminant("D_STXS_VBF_stage1"+JECappend)
            if self.analysis == "fa3fa2fL1fL1Zg":
                return discriminant("D_4couplings_VBFdecay"+JECappend)

        if self.category == "VHHadrtagged":
            if self.analysis == "fa3":
                return discriminant("D_0minus_HadVHdecay"+binsappend+JECappend)
            if self.analysis == "fa2":
                return discriminant("D_0hplus_HadVHdecay"+binsappend+JECappend)
            if self.analysis == "fL1":
                return discriminant("D_L1_HadVHdecay"+binsappend+JECappend)
            if self.analysis == "fL1Zg":
                return discriminant("D_L1Zg_HadVHdecay"+binsappend+JECappend)
            if self.analysis == "fa3_STXS":
                return discriminant("D_STXS_VBF_stage1"+JECappend)
            if self.analysis == "fa3fa2fL1fL1Zg":
                return discriminant("D_4couplings_HadVHdecay"+JECappend)

        assert False

    @property
    def mixdiscriminant(self):
        from discriminants import discriminant

        if self.production >= "180416" or self.production == "180224_newdiscriminants":
            if self.category == "Untagged":
                binsappend = "_new_20bins"
            else:
                binsappend = "_new_10bins"
            if self.analysis in ("fL1", "fL1Zg"): binsappend = binsappend.replace("_new", "")
        elif self.production == "180224_10bins":
            binsappend = "_10bins"
        elif self.production.year == 2016:
            binsappend = ""

        if self.category == "Untagged":
            if self.analysis == "fa3":
                return discriminant("D_CP_decay"+binsappend)
            if self.analysis == "fa2":
                return discriminant("D_int_decay"+binsappend)
            if self.analysis == "fL1":
                return discriminant("D_0hplus_decay"+binsappend)
            if self.analysis == "fL1Zg":
                return discriminant("D_0hplus_decay"+binsappend)
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
                return discriminant("D_CP_VBF"+binsappend+JECappend)
            if self.analysis == "fa2":
                return discriminant("D_int_VBF"+binsappend+JECappend)
            if self.analysis == "fL1":
                return discriminant("D_0hplus_VBFdecay"+binsappend+JECappend)
            if self.analysis == "fL1Zg":
                return discriminant("D_0hplus_VBFdecay"+binsappend+JECappend)
            if self.analysis == "fa3_STXS":
                return discriminant("phistarZ2")
            if self.analysis == "fa3fa2fL1fL1Zg":
                return discriminant("D_CP_VBF_2bins"+JECappend)

        if self.category == "VHHadrtagged":
            if self.analysis == "fa3":
                return discriminant("D_CP_HadVH"+binsappend+JECappend)
            if self.analysis == "fa2":
                return discriminant("D_int_HadVH"+binsappend+JECappend)
            if self.analysis == "fL1":
                return discriminant("D_0hplus_HadVHdecay"+binsappend+JECappend)
            if self.analysis == "fL1Zg":
                return discriminant("D_0hplus_HadVHdecay"+binsappend+JECappend)
            if self.analysis == "fa3_STXS":
                return discriminant("phistarZ2")
            if self.analysis == "fa3fa2fL1fL1Zg":
                return discriminant("D_CP_HadVH_2bins"+JECappend)

        assert False

    @property
    def discriminants(self):
        return (self.purediscriminant, self.mixdiscriminant, self.bkgdiscriminant)

    @property
    @cache
    def invertedmatrix(self):
        productionmode = str(self.templategroup).upper().replace("GGH", "ggH").replace("TTH", "ttH").replace("BBH", "bbH")
        basis = SampleBasis([template.hypothesis for template in self.templates()], productionmode, self.analysis)
        invertedmatrix = basis.invertedmatrix
        try:
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
                    def rows():
                        for matrixreturns in self.inttemplates():
                            j = matrixreturns.rowofinvertedmatrix
                            hypothesispowers = Counter({self.analysis.purehypotheses["1ijkl".index(k)]: v for k, v in matrixreturns.interferencetype.couplingpowers.iteritems()})
                            yield j, matrixreturns, hypothesispowers
                        for index, h in izip((0, 4, 14, 34, 69), self.analysis.purehypotheses):
                            yield index, h, Counter({h: 4})

                    for j, matrixreturns, hypothesispowers in rows():
                        divideby = 1e4**sum(hypothesispowers[Hypothesis(_)] for _ in ("L1", "L1Zg"))
                        for i, matrixmultiplies in enumerate(self.signalsamples()):
                            multiplyby = 1
                            if matrixmultiplies.hypothesis in ("L1", "L1Zg"): multiplyby = 1e16
                            threshold = 1e-10 * multiplyby / divideby
                            if abs(invertedmatrix[j,i]) < threshold: invertedmatrix[j,i] = 0
                            if abs(invertedmatrix[j,i] - 1) < 1e-10: invertedmatrix[j,i] = 1

                    assert invertedmatrix[0, 0] == 1 and all(invertedmatrix[0, i] == 0 for i in             range(1, 70))
                    assert invertedmatrix[4, 1] == 1 and all(invertedmatrix[4, i] == 0 for i in range(0, 1)+range(2, 70))
                    assert invertedmatrix[14,2] == 1 and all(invertedmatrix[14,i] == 0 for i in range(0, 2)+range(3, 70))
                    assert invertedmatrix[34,3] == 1 and all(invertedmatrix[34,i] == 0 for i in range(0, 3)+range(4, 70))
                    assert invertedmatrix[69,4] == 1 and all(invertedmatrix[69,i] == 0 for i in range(0, 4)+range(5, 70))
                else:
                    assert False

            if self.templategroup in ("ggh", "tth", "bbh"):
                if self.analysis.dimensions == 4:
                    assert invertedmatrix[0, 0] == 1 and all(invertedmatrix[0, i] == 0 for i in             range(1, 15))
                    assert                               all(invertedmatrix[2, i] == 0 for i in range(0, 1)+range(2, 15))
                    assert                               all(invertedmatrix[5, i] == 0 for i in range(0, 2)+range(3, 15))
                    assert                               all(invertedmatrix[9, i] == 0 for i in range(0, 3)+range(4, 15))
                    assert                               all(invertedmatrix[14,i] == 0 for i in range(0, 4)+range(5, 15))
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

        except:
            numpy.set_printoptions(edgeitems=35)
            print invertedmatrix
            raise

        return invertedmatrix

    def getjson(self):
        return {
                "inputDirectory": os.path.join("step3_withdiscriminants", str(self.production)),
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
        if self.production == "180224":
          kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in self.enums}
          kwargs["production"] = "170222"
          return TemplatesFile(*kwargs.values())

        if self.production == "180721":
          kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in self.enums}
          kwargs["production"] = "180530"
          return TemplatesFile(*kwargs.values())
        if self.production == "180722" and self.templategroup != "ggh":
          kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in self.enums}
          kwargs["production"] = "180531"
          return TemplatesFile(*kwargs.values())

        if self.production == "180530_Ulascan":
          kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in self.enums}
          kwargs["production"] = "180721_Ulascan"
          return TemplatesFile(*kwargs.values())
        if self.production == "180531_Ulascan":
          kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in self.enums}
          kwargs["production"] = "180722_Ulascan"
          return TemplatesFile(*kwargs.values())

        return None

    @property
    def actualtemplatesfile(self):
        if self.copyfromothertemplatesfile: return self.copyfromothertemplatesfile.actualtemplatesfile
        return self

    @property
    def nominal(self):
        return type(self)(*(getattr(self, enum.enumname) for enum in self.enums if enum != ShapeSystematic))

    @property
    def treeshapesystematics(self):
        for _ in self.allshapesystematics:
            if _ not in treeshapesystematics: continue
            if _ in ("ResUp", "ResDown", "ScaleUp", "ScaleDown") and self.templategroup != "ggh" and config.getm4lsystsfromggH: continue
            if config.getm4lsystsfromggHUntagged and self.category != "Untagged" and shapesystematic in ("ScaleUp", "ScaleDown", "ResUp", "ResDown"): continue
            yield _

    @property
    def allshapesystematics(self):
        for _ in shapesystematics:
            if not _.appliesto(self.templategroup): continue

            if _ in ("JECUp", "JECDn"):
                if self.templategroup not in ("ggh", "zh", "wh"): continue
                if self.templategroup in ("zh", "wh") and self.category == "VBFtagged": continue
                if self.category == "Untagged": continue

            yield _
            

def listfromiterator(function):
    return list(function())

@listfromiterator
def templatesfiles():
    for channel in channels:
        for production in productions:
            if channel != "2e2mu" and production.LHE: continue
            for analysis in analyses:
                for category in categories:
                    if category != "Untagged" and analysis.isdecayonly: continue
                    for templategroup in TemplateGroup.items():
                        nominal = TemplatesFile(channel, templategroup, analysis, production, category)
                        for shapesystematic in nominal.treeshapesystematics:
                            if config.getm4lsystsfromggHUntagged and category != "Untagged" and shapesystematic in ("ScaleUp", "ScaleDown", "ResUp", "ResDown"): continue
                            if (production.LHE or production.GEN) and shapesystematic != "": continue
                            if analysis.isdecayonly and templategroup not in ("bkg", "ggh"): continue
                            if category == "Untagged" and shapesystematic in ("JECUp", "JECDn", "MINLO_SM"): continue

                            yield TemplatesFile(channel, shapesystematic, templategroup, analysis, production, category)


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
            elif enumsdict[ProductionMode] == "bbH":
                enumsdict[TemplateGroup] = "bbh"
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

    @cache
    def gettemplate(self, *args, **kwargs):
        with TFile(self.templatefile(**kwargs)) as f:
            try:
                t = getattr(f, self.templatename())
            except AttributeError:
                raise IOError("No template {} in {}".format(self.templatename(), self.templatefile(**kwargs)))
            t.SetDirectory(0)
            return t

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

        elif self.productionmode in ("ggH", "VBF", "ZH", "WH", "ttH", "bbH"):
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
        if self.productionmode in ("ggH", "ttH", "bbH"):
            match1 = re.match("^f(a2|a3|L1|L1Zg)(?:dec)?(-?)0.5$", str(self.hypothesis))
            match2 = re.match("^f(a2|a3|L1|L1Zg)(?:dec)?0.5f(a2|a3|L1|L1Zg)(?:dec)?(-?)0.5$", str(self.hypothesis))
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
            elif match1:
                name = "templateMixa1{}{}AdapSmooth".format(match1.group(1), "Pi" if match1.group(2) else "")
            elif match2:
                name = "templateMix{}{}{}AdapSmooth".format(match2.group(1), match2.group(2), "Pi" if match2.group(3) else "")
            else:
                print self.hypothesis
        elif self.productionmode in ("VBF", "ZH", "WH"):
            match1 = re.match("^f(?P<ai>a2|a3|L1|L1Zg)(?P<proddec>(?:prod|dec|proddec)?)(?P<sign>-?)0[.]5$", str(self.hypothesis))
            match2 = re.match("^f(?P<ai>a2|a3|L1|L1Zg)(?P<proddec>(?:prod|dec|proddec)?)0[.]5f(?P<aj>a2|a3|L1|L1Zg)(?P=proddec)(?P<sign>-?)0[.]5$", str(self.hypothesis))
            match3 = re.match("^f(?P<ai>a2|a3|L1|L1Zg)(?P<proddec>(?:prod|dec|proddec)?)0[.]33f(?P<aj>a2|a3|L1|L1Zg)(?P=proddec)(?P<sign>-?)0[.]33$", str(self.hypothesis))
            match4 = re.match("^f(?P<ai>a2|a3|L1|L1Zg)(?P<proddec>(?:prod|dec|proddec)?)0[.]33f(?P<aj>a2|a3|L1|L1Zg)(?P=proddec)0[.]33f(?P<ak>a2|a3|L1|L1Zg)(?P=proddec)(?P<sign>-?)0[.]33$", str(self.hypothesis))
            match5 = re.match("^f(?P<ai>a2|a3|L1|L1Zg)(?P<proddec>(?:prod|dec|proddec)?)0[.]25f(?P<aj>a2|a3|L1|L1Zg)(?P=proddec)0[.]25f(?P<ak>a2|a3|L1|L1Zg)(?P=proddec)(?P<sign>-?)0[.]25$", str(self.hypothesis))
            match6 = re.match("^f(?P<ai>a2|a3|L1|L1Zg)(?P<proddec>(?:prod|dec|proddec)?)0[.]25f(?P<aj>a2|a3|L1|L1Zg)(?P=proddec)0[.]25f(?P<ak>a2|a3|L1|L1Zg)(?P=proddec)0[.]25f(?P<al>a2|a3|L1|L1Zg)(?P=proddec)(?P<sign>-?)0[.]25$", str(self.hypothesis))
            pddict = {"prod": "Prod", "dec": "Decay", "": "Decay", "proddec": "ProdDec"}
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
            elif match1:
                name = "templateMixa1{}{}{}AdapSmooth".format(
                    match1.group("ai"),
                    pddict[match1.group("proddec")],
                    "Pi" if match1.group("sign") else ""
                )
            elif match2:
                name = "templateMix{}{}{}{}AdapSmooth".format(
                    match2.group("ai"),
                    match2.group("aj"),
                    pddict[match2.group("proddec")],
                    "Pi" if match2.group("sign") else ""
                )
            elif match3:
                name = "templateMixa1{}{}{}{}AdapSmooth".format(
                    match3.group("ai"),
                    match3.group("aj"),
                    pddict[match3.group("proddec")],
                    "Pi" if match3.group("sign") else ""
                )
            elif match4:
                name = "templateMix{}{}{}{}{}AdapSmooth".format(
                    match4.group("ai"),
                    match4.group("aj"),
                    match4.group("ak"),
                    pddict[match4.group("proddec")],
                    "Pi" if match4.group("sign") else ""
                )
            elif match5:
                name = "templateMixa1{}{}{}{}{}AdapSmooth".format(
                    match5.group("ai"),
                    match5.group("aj"),
                    match5.group("ak"),
                    pddict[match5.group("proddec")],
                    "Pi" if match5.group("sign") else ""
                )
            elif match6:
                name = "templateMix{}{}{}{}{}{}AdapSmooth".format(
                    match6.group("ai"),
                    match6.group("aj"),
                    match6.group("ak"),
                    match6.group("al"),
                    pddict[match6.group("proddec")],
                    "Pi" if match6.group("sign") else ""
                )
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
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "ttH", "bbH"):
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
        result += self.shapesystematic.categoryappendname

        from treewrapper import TreeWrapper
        for categorization in TreeWrapper.categorizations:
            if categorization.category_function_name == result: return result
        assert False, "{} does not exist in TreeWrapper".format(result)

    def reweightfrom(self):
        if self.productionmode == "ggH":
            if self.production.LHE:
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
        if self.productionmode in ("ttH", "bbH"):
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
            if config.showblinddistributions:
                result = {Sample(self.production, self.productionmode)}
            else:
                result = {Sample("ggZZ", self.production, "4tau")}
        result = {sample for sample in result if tfiles[sample.withdiscriminantsfile()].candTree.GetEntries() != 0}
        assert result, self
        return result

    @property
    def scalefactor(self):
        if self.shapesystematic == "MINLO_SM":
            result = len(self.reweightfrom())
        elif self.productionmode in ("VBF", "ggH", "ZH", "WH", "bbH"):
            result = uncertainties.nominal_value(ReweightingSampleWithPdf(self.productionmode, self.hypothesis, self.production).xsec) / uncertainties.nominal_value(ReweightingSampleWithPdf(self.productionmode, "SM", self.production).xsec)
        elif self.productionmode == "ttH":
            result = uncertainties.nominal_value(ReweightingSampleWithPdf(self.productionmode, self.hypothesis, self.hffhypothesis, self.production).xsec) / uncertainties.nominal_value(ReweightingSampleWithPdf(self.productionmode, "SM", "Hff0+", self.production).xsec)
        elif self.productionmode == "data" and not config.showblinddistributions:
            return None  #this is to avoid calling effectiveentries
        elif self.productionmode in ("VBF bkg", "ggZZ", "ZX", "data"):
            result = len(self.reweightfrom())
        elif self.productionmode == "qqZZ":
            result = 1
        try:
            result /= sum(
                          Sample.effectiveentries(
                                                  reweightfrom=reweightfrom,
                                                  reweightto=ReweightingSample(self.productionmode, self.hypothesis, self.hffhypothesis)
                                                 )
                           for reweightfrom in self.reweightfrom()
                         )
        except ZeroDivisionError:
            result = 0
        if isnan(result): result = 0
        if self.productionmode in ("VBF bkg", "ggZZ", "ZX", "data"):
            assert result == 1
        return result

    @property
    def domirror(self):
        if self.analysis not in ("fa3", "fa3_STXS", "fa3fa2fL1fL1Zg"): return False
        if self.analysis == "fa3_STXS": return False #for now... could do it later
        if self.productionmode == "data": return False

        if self.hypothesis in ("fa30.5", "fa3prod0.5", "fa3proddec-0.5"): return False
        if self.hypothesis in ("0+", "0-", "a2", "L1", "L1Zg"): return True

        if self.hypothesis in ("fa20.5", "fa2-0.5", "fL10.5", "fL1Zg0.5", "fL1Zg-0.5"): return True
        if self.hypothesis in ("fa2prod0.5", "fa2prod-0.5", "fL1prod0.5", "fL1Zgprod0.5", "fL1Zgprod-0.5"): return True
        if self.hypothesis in ("fa2proddec0.5", "fa2proddec-0.5", "fL1proddec-0.5", "fL1Zgproddec0.5", "fL1Zgproddec-0.5"): return True

        if self.hypothesis in ("fa20.5fL10.5", "fa20.5fL1Zg0.5", "fL10.5fL1Zg0.5"): return True
        if self.hypothesis in ("fa2prod0.5fL1prod0.5", "fa2prod0.5fL1Zgprod0.5", "fL1prod0.5fL1Zgprod0.5"): return True
        if self.hypothesis in ("fa2proddec0.5fL1proddec-0.5", "fa2proddec0.5fL1Zgproddec-0.5", "fL1proddec0.5fL1Zgproddec-0.5"): return True

        if self.hypothesis in ("fa20.33fL10.33", "fa20.33fL1Zg0.33", "fL10.33fL1Zg0.33"): return True
        if self.hypothesis in ("fa2prod0.33fL1prod0.33", "fa2prod0.33fL1Zgprod0.33", "fL1prod0.33fL1Zgprod0.33"): return True
        if self.hypothesis in ("fa2proddec0.33fL1proddec-0.33", "fa2proddec0.33fL1Zgproddec-0.33", "fL1proddec0.33fL1Zgproddec-0.33"): return True

        if self.hypothesis in ("fa20.33fL10.33fL1Zg0.33",): return True
        if self.hypothesis in ("fa2prod0.33fL1prod0.33fL1Zgprod0.33",): return True
        if self.hypothesis in ("fa2proddec0.33fL1proddec0.33fL1Zgproddec-0.33",): return True

        if self.hypothesis in ("fa2proddec0.25fL1proddec0.25fL1Zgproddec0.25",): return True

        if self.hypothesis in ("fa30.5fa20.5", "fa30.5fL10.5", "fa30.5fL1Zg0.5"): return False
        if self.hypothesis in ("fa3prod0.5fa2prod0.5", "fa3prod0.5fL1prod0.5", "fa3prod0.5fL1Zgprod0.5"): return False
        if self.hypothesis in ("fa3proddec0.5fa2proddec-0.5", "fa3proddec0.5fL1proddec-0.5", "fa3proddec0.5fL1Zgproddec-0.5"): return False

        if self.hypothesis in ("fa30.33fa20.33", "fa30.33fL10.33", "fa30.33fL1Zg0.33"): return False
        if self.hypothesis in ("fa3prod0.33fa2prod0.33", "fa3prod0.33fL1prod0.33", "fa3prod0.33fL1Zgprod0.33"): return False
        if self.hypothesis in ("fa3proddec0.33fa2proddec-0.33", "fa3proddec0.33fL1proddec-0.33", "fa3proddec0.33fL1Zgproddec-0.33"): return False

        if self.hypothesis in ("fa30.33fa20.33fL10.33", "fa30.33fa20.33fL1Zg0.33", "fa30.33fL10.33fL1Zg0.33"): return False
        if self.hypothesis in ("fa3prod0.33fa2prod0.33fL1prod0.33", "fa3prod0.33fa2prod0.33fL1Zgprod0.33", "fa3prod0.33fL1prod0.33fL1Zgprod0.33"): return False
        if self.hypothesis in ("fa3proddec0.33fa2proddec0.33fL1proddec-0.33", "fa3proddec0.33fa2proddec0.33fL1Zgproddec-0.33", "fa3proddec0.33fL1proddec0.33fL1Zgproddec-0.33"): return False

        if self.hypothesis in ("fa3proddec0.25fa2proddec0.25fL1proddec0.25", "fa3proddec0.25fa2proddec0.25fL1Zgproddec0.25", "fa3proddec0.25fL1proddec0.25fL1Zgproddec0.25"): return False

        if self.hypothesis in ("fa3proddec0.25fa2proddec0.25fL1proddec0.25fL1Zgproddec0.25",): return False

        if self.productionmode in ("qqZZ", "ggZZ", "VBF bkg", "ZX"): return True

        assert False, self

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
        if self.copyfromothertemplate: return []
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
        if not self.production.LHE:
            result.append("Z1Flav*Z2Flav == {}".format(self.ZZFlav))
        if not self.analysis.isdecayonly:
            result.append("(" + " || ".join("{} == {}".format(self.categoryname, c) for c in self.category.idnumbers) + ")")
        if self.productionmode == "data" and not config.showblinddistributions:
            result.insert(0, "0")
        return " && ".join(result)

    @property
    def ZZFlav(self):
        assert not self.production.LHE
        result = self.channel.ZZFlav
        if self.productionmode == "ZX": result *= -1
        return result

    @property
    def sample(self):
        return ReweightingSample(self.hypothesis, self.productionmode)

    @property
    def copyfromothertemplate(self):
        if self.productionmode in ("ggZZ", "VBF bkg", "ZX") and self.production in ("180721", "180722"):
            return type(self)(self.productionmode, self.channel, self.shapesystematic, self.templategroup, self.analysis, str(self.production)+"_Ulascan", self.category)
        return None

    def gettemplate(self):
        if self.copyfromothertemplate: return self.copyfromothertemplate.gettemplate()
        return super(Template, self).gettemplate()

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
        enumitems = tuple(EnumItem(_) for _ in [
                    "g{}1g{}1".format(*indices)
                        for indices in combinations("1ijkl", 2)
                ] + [
                    "g{}3g{}1".format(*indices)
                        for indices in combinations("1ijkl", 2)
                ] + [
                    "g{}2g{}2".format(*indices)
                        for indices in combinations("1ijkl", 2)
                ] + [
                    "g{}1g{}3".format(*indices)
                        for indices in combinations("1ijkl", 2)
                ] + [
                    "g{}2g{}1g{}1".format(*indices)
                        for indices in combinations("1ijkl", 3)
                ] + [
                    "g{}1g{}2g{}1".format(*indices)
                        for indices in combinations("1ijkl", 3)
                ] + [
                    "g{}1g{}1g{}2".format(*indices)
                        for indices in combinations("1ijkl", 3)
                ] + [
                    "g{}1g{}1g{}1g{}1".format(*indices)
                        for indices in combinations("1ijkl", 4)
                ])

        @property
        def couplingpowers(self):
            return Counter({coupling: int(power) for coupling, power in re.findall("g(.)([1-9][0-9]*)", str(self))})
        @property
        def mindimensions(self):
            return max({"1": 0, "i": 1, "j": 2, "k": 3, "l": 4}[g] for g in self.couplingpowers.elements())
        @property
        def totalpower(self):
            return sum(self.couplingpowers.values())

    enums = [TemplatesFile, ProductionMode, InterferenceType, HffHypothesis]

    def check(self, *args):

        dontcheck = []

        if self.productionmode == "ttH":
            if self.hffhypothesis != "Hff0+":
                raise ValueError("{} is not used to make templates {}!\n{}".format(self.hffhypothesis, args))
        else:
            if self.hffhypothesis is not None:
               raise ValueError("HffHypothesis {} provided for {}!\n{}".format(self.hffhypothesis, self.productionmode, args))
            dontcheck.append(HffHypothesis)

        super(IntTemplate, self).check(*args, dontcheck=dontcheck)

        del dontcheck

        if self.interferencetype.mindimensions > self.analysis.dimensions:
            raise ValueError("Invalid interferencetype {} for {}D analysis\n{}".format(self.interferencetype, self.analysis.dimensions, args))

        if self.productionmode in ("ggH", "ttH", "bbH"):
            if self.interferencetype.totalpower != 2:
                raise ValueError("Invalid interferencetype {} for productionmode {}!\n{}".format(self.interferencetype, self.productionmode, args))
        elif self.productionmode in ("VBF", "ZH", "WH"):
            if self.interferencetype.totalpower != 4:
                raise ValueError("Invalid interferencetype {} for productionmode {}!\n{}".format(self.interferencetype, self.productionmode, args))
        else:
            raise ValueError("Invalid productionmode {}!\n{}".format(self.productionmode, args))

    def templatename(self, final=True):
        if self.analysis.dimensions == 4:
            result = "template"+str(self.interferencetype)+"IntAdapSmooth"
            result = result.replace("gi", self.analysis.couplingnames[0])
            result = result.replace("gj", self.analysis.couplingnames[1])
            result = result.replace("gk", self.analysis.couplingnames[2])
            result = result.replace("gl", self.analysis.couplingnames[3])
        elif self.productionmode in ("ggH", "ttH", "bbH"):
            if self.interferencetype == "g11gi1":
                if self.analysis in ("fa3", "fa3_STXS"):
                    result = "templatea1a3IntAdapSmooth"
                elif self.analysis == "fa2":
                    result = "templatea1a2IntAdapSmooth"
                elif self.analysis == "fL1" or self.analysis.isfL1fL1Zg:
                    result = "templatea1L1IntAdapSmooth"
                elif self.analysis == "fL1Zg":
                    result = "templatea1L1ZgIntAdapSmooth"
                else:
                    assert False
                if "170203" in self.templatesfile.templatesfile():
                    result = "templateIntAdapSmooth"
            elif self.interferencetype == "g11gj1":
                if self.analysis.isfL1fL1Zg:
                    result = "templatea1L1ZgIntAdapSmooth"
            elif self.interferencetype == "gi1gj1":
                if self.analysis.isfL1fL1Zg:
                    result = "templateL1L1ZgIntAdapSmooth"
        elif self.productionmode in ("VBF", "ZH", "WH"):
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
        if self.analysis not in ("fa3", "fa3_STXS", "fa3fa2fL1fL1Zg"): return None

        #cross talk - production discriminants for the wrong category don't make sense
        if self.category in ("VBFtagged", "VHHadrtagged") and self.productionmode in ("ggH", "ttH", "bbH"):
            if (self.interferencetype == "g11gi1"
                or self.analysis == "fa3fa2fL1fL1Zg" and self.interferencetype.couplingpowers["i"] == 1):
                #ggH has no production information, and only using SM ttH, so mirror antisymmetric
                #over the (pretend) D_CP_decay axis, which sets the whole thing to 0
                #note for ttH this is an approximation, since we could have H(0-)->2l2q tt->bbllnunu
                return {"type":"rescale", "factor":0}
            if self.analysis == "fa3fa2fL1fL1Zg" and self.interferencetype.couplingpowers["i"] == 0:
                return {"type":"mirror", "antisymmetric":False, "axis":1}
            assert False

        if self.analysis == "fa3_STXS" and self.interferencetype in ("g11gi1", "g11gi3", "g13gi1"):
            #same (antimirror over D_CP_whatever, which doesn't exist in STXS)
            return {"type":"rescale", "factor":0}

        if self.analysis == "fa3_STXS": return None #for now... could do it later

        #Mirror antisymmetric for VH in VBF category and VBF in VH category
        #the interference templates are 0 to within error bars anyway,
        #except for some effects in ZH g13g41 D_CP_VBF which are antisymmetric

        #cross talk to the untagged category is exactly correct, since the decay is the same

        if (self.interferencetype in ("g11gi1", "g11gi3", "g13gi1")
            or self.analysis == "fa3fa2fL1fL1Zg" and self.interferencetype.couplingpowers["i"] in (1, 3)):
            return {"type":"mirror", "antisymmetric":True, "axis":1}
        elif (self.interferencetype == "g12gi2"
              or self.analysis == "fa3fa2fL1fL1Zg" and self.interferencetype.couplingpowers["i"] in (0, 2, 4)):
            return {"type":"mirror", "antisymmetric":False, "axis":1}
        assert False

    @property
    def domirror(self):
        return bool(self.mirrorjsn)

    @property
    def rowofinvertedmatrix(self):
        index = 0
        maxpower = self.interferencetype.totalpower
        for l in range(maxpower+1):
            for k in range(maxpower+1-l):
                for j in range(maxpower+1-k-l):
                    for i in range(maxpower+1-j-k-l):
                        dct = Counter({"i": i, "j": j, "k": k, "l": l})
                        dct["1"] = maxpower - sum(dct.values())
                        for key, val in dct.items():
                            if not val:
                                del dct[key]
                        if dct == self.interferencetype.couplingpowers:
                            return index
                        index += 1
        assert False, self.interferencetype.couplingpowers

    @property
    @cache
    def templatesandfactors(self):
        if self.analysis.dimensions == 4:
            rowofinvertedmatrix = self.rowofinvertedmatrix
            hffhypothesis = "Hff0+" if self.productionmode == "ttH" else None
            multiplyby = 1

        elif self.interferencetype == "g11gi1":
            if self.analysis.dimensions in (1, 2):
                g1exp = giexp = 1
                rowofinvertedmatrix = giexp #first row is labeled 0
                hffhypothesis = "Hff0+" if self.productionmode == "ttH" else None
                multiplyby = getattr(ReweightingSample(self.productionmode, self.analysis.purehypotheses[1], hffhypothesis), self.analysis.purehypotheses[1].couplingname)
            else:
                assert False, analysis.dimensions

        elif self.interferencetype == "g11gj1":
            if self.analysis.dimensions == 2:
                rowofinvertedmatrix = 3 #first row is labeled 0
                hffhypothesis = "Hff0+" if self.productionmode == "ttH" else None
                multiplyby = getattr(ReweightingSample(self.productionmode, self.analysis.purehypotheses[2], hffhypothesis), self.analysis.purehypotheses[2].couplingname)
            else:
                assert False, analysis.dimensions

        elif self.interferencetype == "gi1gj1":
            if self.analysis.dimensions == 2:
                rowofinvertedmatrix = 4 #first row is labeled 0
                hffhypothesis = "Hff0+" if self.productionmode == "ttH" else None
                multiplyby = (
                              getattr(ReweightingSample(self.productionmode, self.analysis.purehypotheses[1], hffhypothesis), self.analysis.purehypotheses[1].couplingname)
                             *
                              getattr(ReweightingSample(self.productionmode, self.analysis.purehypotheses[2], hffhypothesis), self.analysis.purehypotheses[2].couplingname)
                             )
            else:
                assert False, analysis.dimensions

        elif self.interferencetype in ("g11gi3", "g12gi2", "g13gi1"):
            if self.analysis.dimensions == 1:
                g1exp, giexp = (int(i) for i in str(self.interferencetype).replace("g1", "").replace("gi", ""))
                rowofinvertedmatrix = giexp #first row is labeled 0
                multiplyby = 1
            else:
                assert False, analysis.dimensions

        invertedmatrix = self.templatesfile.invertedmatrix
        vectoroftemplates = self.templatesfile.templates()  #ok technically it's a list not a vector
        templatesandfactors = []

        for j, template in enumerate(vectoroftemplates):
            factor = invertedmatrix[rowofinvertedmatrix,j] * multiplyby
            if factor:
                templatesandfactors.append((template, factor))

        return templatesandfactors

    def getjson(self):
        templatesum = [{
                        "name": template.templatename(final=False),
                        "factor": factor,
                       } for template, factor in self.templatesandfactors if factor]
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

    @cache
    def gettemplate(self, *args, **kwargs):
        result = super(IntTemplate, self).gettemplate(*args, **kwargs)
        if self.analysis == "fa3" and self.productionmode in ("ggH", "ttH", "bbH") and self.interferencetype == "g11gi1" and self.category in ("VBFtagged", "VHHadrtagged"):
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
        if self.production.LHE:
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
        for production in productions:
            if channel != "2e2mu" and production.LHE: continue
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

if __name__ == "__main__":
    numpy.set_printoptions(edgeitems=35)
    tf = TemplatesFile("vbf", "2e2mu", "Untagged", "fa3fa2fL1fL1Zg")
    productionmode = str(tf.templategroup).upper().replace("GGH", "ggH").replace("TTH", "ttH").replace("BBH", "bbH")
    basis = SampleBasis([template.hypothesis for template in tf.templates()], productionmode, tf.analysis)
    print basis.matrix
    print numpy.linalg.det(basis.matrix)
    print basis.invertedmatrix
