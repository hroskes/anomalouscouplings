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
from enums import AlternateWeight, Analysis, analyses, Channel, channels, Category, categories, EnumItem, flavors, hffhypotheses, HffHypothesis, Hypothesis, MultiEnum, MultiEnumABCMeta, MyEnum, Production, ProductionMode, productions, ShapeSystematic, shapesystematics, TemplateGroup, templategroups, treeshapesystematics
from samples import ReweightingSample, ReweightingSamplePlus, ReweightingSampleWithPdf, Sample, SampleBasis, SumOfSamples
from utilities import cache, deprecate, is_almost_integer, JsonDict, jsonloads, TFile, withdiscriminantsfileisvalid

class TemplatesFile(MultiEnum):
    enumname = "templatesfile"
    enums = [Channel, ShapeSystematic, TemplateGroup, Analysis, Production, Category]

    def applysynonyms(self, enumsdict):
        if enumsdict[Production] is None and len(config.productionsforcombine) == 1:
            enumsdict[Production] = config.productionsforcombine[0]
        super(TemplatesFile, self).applysynonyms(enumsdict)

    def check(self, *args):
        dontcheck = [ShapeSystematic]

        if self.category is None:
            self.category = Category("Untagged")

        if self.analysis.isdecayonly:
            if self.category != "Untagged":
                raise ValueError("decay only analysis is only done for untagged!\n{}".format(args))
            if self.templategroup in ("vbf", "zh", "wh", "vh", "tth"):
                raise ValueError("decay only analysis is only done with decay information!\n{}".format(args))
        if not self.analysis.useboosted:
            if self.category == "Boosted":
                raise ValueError("{} doesn't use boosted\n{}".format(self.analysis, args))
        if not self.analysis.usemorecategories:
            if self.category in ("VBF1jtagged", "VHLepttagged"):
                raise ValueError("{} doesn't use {}\n{}".format(self.analysis, self.category, args))

        if self.production.LHE:
            if self.channel != "2e2mu":
                raise ValueError("LHE analysis is only done for 2e2mu for now!\n{}".format(args))

        super(TemplatesFile, self).check(*args, dontcheck=dontcheck)

        if None is not self.shapesystematic not in self.allshapesystematics:
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

        folder = os.path.join(config.repositorydir, "step7_templates", str(self.production))
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
                commonargs = "ggH",
                if self.category in ("VBFtagged", "VHHadrtagged"):
                    otherargslist = [["Hff0+", "JHUGen"], ["Hff0-", "JHUGen"]]
                else:
                    otherargslist = [["Hff0+"]]
            elif self.templategroup == "tth":
                commonargs = "ttH",
                otherargslist = [["Hff0+"], ["Hff0-"]]
            elif self.templategroup == "bbh":
                commonargs = "bbH",
                otherargslist = [[]]
            else:
                assert False, self.templategroup

            reweightingsamples = []

            for otherargs in otherargslist:
                args = commonargs + tuple(otherargs)
                if self.analysis.fais == ("fa3",):
                    reweightingsamples += [ReweightingSample("0+", *args), ReweightingSample("0-", *args), ReweightingSample("fa30.5", *args)]
                if self.analysis == "fa2":
                    reweightingsamples += [ReweightingSample("0+", *args), ReweightingSample("a2", *args), ReweightingSample("fa2-0.5", *args)]
                if self.analysis == "fL1":
                    reweightingsamples += [ReweightingSample("0+", *args), ReweightingSample("L1", *args), ReweightingSample("fL10.5", *args)]
                if self.analysis == "fL1Zg":
                    reweightingsamples += [ReweightingSample("0+", *args), ReweightingSample("L1Zg", *args), ReweightingSample("fL1Zg-0.5", *args)]
                if self.analysis.isfL1fL1Zg:
                    reweightingsamples += [ReweightingSample("0+", *args), ReweightingSample("L1", *args), ReweightingSample("L1Zg", *args), ReweightingSample("fL10.5", *args), ReweightingSample("fL1Zg0.5", *args), ReweightingSample("fL10.5fL1Zg0.5", *args)]
                if self.analysis.isfa3fa2fL1fL1Zg:
                    hypotheses = ["0+", "0-", "a2", "L1", "L1Zg",
                                  "fa30.5", "fa2-0.5", "fL10.5", "fL1Zg-0.5",
                                  "fa30.5fa20.5", "fa30.5fL10.5", "fa30.5fL1Zg0.5",
                                                  "fa20.5fL10.5", "fa20.5fL1Zg0.5",
                                                                  "fL10.5fL1Zg0.5",
                                 ]
                    reweightingsamples += [
                      ReweightingSamplePlus(h, *args) for h in hypotheses
                    ]

        elif self.templategroup in ("vbf", "zh", "wh", "vh"):
            p = str(self.templategroup).upper()
            if self.analysis.fais == ("fa3",):
                reweightingsamples = [ReweightingSample(p, "0+"), ReweightingSample(p, "0-"), ReweightingSample(p, "fa3prod0.5"), ReweightingSample(p, "fa3dec0.5"), ReweightingSample(p, "fa3proddec-0.5")]
            if self.analysis == "fa2":
                reweightingsamples = [ReweightingSample(p, "0+"), ReweightingSample(p, "a2"), ReweightingSample(p, "fa2prod0.5"), ReweightingSample(p, "fa2dec-0.5"), ReweightingSample(p, "fa2proddec-0.5")]
            if self.analysis == "fL1":
                reweightingsamples = [ReweightingSample(p, "0+"), ReweightingSample(p, "L1"), ReweightingSample(p, "fL1prod0.5"), ReweightingSample(p, "fL1dec0.5"), ReweightingSample(p, "fL1proddec-0.5")]
            if self.analysis == "fL1Zg":
                reweightingsamples = [ReweightingSample(p, "0+"), ReweightingSample(p, "L1Zg"), ReweightingSample(p, "fL1Zgprod0.5"), ReweightingSample(p, "fL1Zgdec0.5"), ReweightingSample(p, "fL1Zgproddec-0.5")]
            if self.analysis.isfa3fa2fL1fL1Zg:
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
                  ReweightingSamplePlus(h, p) for h in hypotheses
                ]

        return reweightingsamples

    def templates(self):
        if self.templategroup in ["ggh", "vbf", "zh", "wh", "vh", "tth", "bbh"]:
            return [Template(self, sample) for sample in self.signalsamples()]
        elif self.templategroup == "bkg":
            if self.production.LHE or self.production.GEN: return [Template(self, "qqZZ")]
            result = ["qqZZ", "ggZZ", "VBF bkg"]
            result.remove("VBF bkg")
            if config.usedata:
                result.append("ZX")

            if self.shapesystematic == "":
                pass
            elif self.shapesystematic in ("ZXUp", "ZXDown"):
                result = [_ for _ in result if _ == "ZX"]
            elif self.shapesystematic in ("JECUp", "JECDown"):
                result = [_ for _ in result if _ != "ZX"]
            else:
                assert False, self.shapesystematic

            return [Template(self, productionmode) for productionmode in result]
        elif self.templategroup == "DATA":
            return [Template(self, "data")]
        assert False

    def inttemplates(self):
        if self.templategroup == "ggh" and self.shapesystematic != "MINLO_SM" or self.templategroup in ("tth", "bbh"):
            hffhypotheses = [_ for _ in (None, "Hff0+", "Hff0-", "fCP0.5") if any(t.hffhypothesis == _ for t in self.templates())]
            h = str(self.templategroup).replace("h", "H")
            if self.analysis.dimensions == 4:
                return [IntTemplate(self, h, "g{}1g{}1".format(i1, i2), hh)
                        for i1, i2 in combinations("1ijkl", 2)
                        for hh in hffhypotheses
                ]
            elif self.analysis.dimensions == 2:
                return [IntTemplate(self, h, _, hh) for _ in ("g11gi1", "g11gj1", "gi1gj1") for hh in hffhypotheses]
            elif self.analysis.dimensions == 1:
                return [IntTemplate(self, h, "g11gi1", hh) for hh in self.hffhypotheses]

        elif self.templategroup in ("vbf", "zh", "wh", "vh"):
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
                ] if not ("gl3" in _ and self.templategroup == "wh")]
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

        if self.analysis in ("fa3_only6bins", "fa3fa2fL1fL1Zg_only6bins", "fa3_onlyDCP"): return discriminant("phistarZ2")

        name = "D_bkg"
        if self.category == "Untagged": pass
        elif self.category == "VBFtagged": name += "_VBFdecay"
        elif self.category == "VHHadrtagged": name += "_HadVHdecay"
        elif self.category == "Boosted": pass
        elif self.category == "VBF1jtagged": pass
        elif self.category == "VHLepttagged": pass
        else: assert False

        name += self.shapesystematic.Dbkgappendname

        if self.analysis == "fa3_multiparameter" or self.analysis.isfa3fa2fL1fL1Zg or self.analysis.isSTXS:
            name += "_3bins"
        elif self.category == "Untagged":
            name += "_20bins"
        else:
            name += "_10bins"

        if self.shapesystematic in ("JECUp", "JECDn") and ("VBF" in name or "HadVH" in name):
            name += "_{}".format(self.shapesystematic)

        return discriminant(name)

    @property
    def purediscriminant(self):
        from discriminants import discriminant

        if self.analysis == "fa3_onlyDbkg": return discriminant("phistarZ2")

        if self.shapesystematic in ("JECUp", "JECDn"):
            JECappend = "_{}".format(self.shapesystematic)
        else:
            JECappend = ""

        if self.category == "Untagged":
            binsappend = "_20bins"
        else:
            binsappend = "_10bins"

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
            if self.analysis.isSTXS:
                return discriminant("D_STXS_stage1p1_untagged"+JECappend)
            if self.analysis in ("fa3fa2fL1fL1Zg", "fa3fa2fL1fL1Zg_decay", "fa3fa2fL1fL1Zg_boosted", "fa3fa2fL1fL1Zg_morecategories", "fa3_multiparameter"):
                return discriminant("D_4couplings_decay")
            if self.analysis in ("fa3_only6bins", "fa3fa2fL1fL1Zg_only6bins"):
                return discriminant("D_0minus_decay_3bins")
            if self.analysis == "fa3_onlyDCP":
                return discriminant("phistarZ1")

        if self.category == "VBFtagged":
            if self.analysis == "fa3":
                return discriminant("D_0minus_VBFdecay"+binsappend+JECappend)
            if self.analysis == "fa2":
                return discriminant("D_0hplus_VBFdecay"+binsappend+JECappend)
            if self.analysis == "fL1":
                return discriminant("D_L1_VBFdecay"+binsappend+JECappend)
            if self.analysis == "fL1Zg":
                return discriminant("D_L1Zg_VBFdecay"+binsappend+JECappend)
            if self.analysis.isSTXS:
                return discriminant("D_STXS_stage1p1_VBF"+JECappend)
            if self.analysis in ("fa3fa2fL1fL1Zg", "fa3fa2fL1fL1Zg_boosted", "fa3fa2fL1fL1Zg_morecategories", "fa3_multiparameter"):
                return discriminant("D_4couplings_VBFdecay"+JECappend)
            if self.analysis in ("fa3_only6bins", "fa3fa2fL1fL1Zg_only6bins"):
                return discriminant("D_0minus_VBFdecay_3bins")
            if self.analysis == "fa3_onlyDCP":
                return discriminant("phistarZ1")

        if self.category == "VHHadrtagged":
            if self.analysis == "fa3":
                return discriminant("D_0minus_HadVHdecay"+binsappend+JECappend)
            if self.analysis == "fa2":
                return discriminant("D_0hplus_HadVHdecay"+binsappend+JECappend)
            if self.analysis == "fL1":
                return discriminant("D_L1_HadVHdecay"+binsappend+JECappend)
            if self.analysis == "fL1Zg":
                return discriminant("D_L1Zg_HadVHdecay"+binsappend+JECappend)
            if self.analysis.isSTXS:
                return discriminant("D_STXS_stage1p1_HadVH"+JECappend)
            if self.analysis in ("fa3fa2fL1fL1Zg", "fa3fa2fL1fL1Zg_boosted", "fa3fa2fL1fL1Zg_morecategories", "fa3_multiparameter"):
                return discriminant("D_4couplings_HadVHdecay"+JECappend)
            if self.analysis in ("fa3_only6bins", "fa3fa2fL1fL1Zg_only6bins"):
                return discriminant("D_0minus_HadVHdecay_3bins")
            if self.analysis == "fa3_onlyDCP":
                return discriminant("phistarZ1")

        if self.category == "Boosted":
            if self.analysis in ("fa3fa2fL1fL1Zg_boosted", "fa3fa2fL1fL1Zg_morecategories"):
                return discriminant("ZZPt_boosted")

        if self.category == "VBF1jtagged":
            if self.analysis == "fa3fa2fL1fL1Zg_morecategories":
                return discriminant("ZZPt_VBF1jtagged")

        if self.category == "VHLepttagged":
            if self.analysis == "fa3fa2fL1fL1Zg_morecategories":
                return discriminant("ZZPt_VHLepttagged")

        assert False

    @property
    def mixdiscriminant(self):
        from discriminants import discriminant

        if self.category == "Untagged":
            binsappend = "_new_20bins"
        else:
            binsappend = "_new_10bins"
        if self.analysis in ("fL1", "fL1Zg"): binsappend = binsappend.replace("_new", "")

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
            if self.analysis.isSTXS:
                return discriminant("phistarZ2")
            if self.analysis in ("fa3fa2fL1fL1Zg", "fa3fa2fL1fL1Zg_decay", "fa3fa2fL1fL1Zg_boosted", "fa3fa2fL1fL1Zg_morecategories", "fa3_multiparameter"):
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
            if self.analysis.isSTXS:
                return discriminant("phistarZ2")
            if self.analysis in ("fa3fa2fL1fL1Zg", "fa3fa2fL1fL1Zg_boosted", "fa3fa2fL1fL1Zg_morecategories", "fa3_multiparameter"):
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
            if self.analysis.isSTXS:
                return discriminant("phistarZ2")
            if self.analysis in ("fa3fa2fL1fL1Zg", "fa3fa2fL1fL1Zg_boosted", "fa3fa2fL1fL1Zg_morecategories", "fa3_multiparameter"):
                return discriminant("D_CP_HadVH_2bins"+JECappend)

        if self.category == "Boosted":
            if self.analysis in ("fa3fa2fL1fL1Zg_boosted", "fa3fa2fL1fL1Zg_morecategories"):
                return discriminant("phistarZ2")

        if self.category == "VBF1jtagged":
            if self.analysis == "fa3fa2fL1fL1Zg_morecategories":
                return discriminant("phistarZ2")

        if self.category == "VHLepttagged":
            if self.analysis == "fa3fa2fL1fL1Zg_morecategories":
                return discriminant("phistarZ2")

        assert False, self

    @property
    def discriminants(self):
        return (self.purediscriminant, self.mixdiscriminant, self.bkgdiscriminant)

    @property
    @cache
    def invertedmatrix(self):
        productionmode = str(self.templategroup).upper().replace("GGH", "ggH").replace("TTH", "ttH").replace("BBH", "bbH")
        basis = SampleBasis([template.hypothesis for template in self.templates() if template.hffhypothesis in (None, "Hff0+")], productionmode, self.analysis)
        invertedmatrix = basis.invertedmatrix
        try:
            if self.templategroup in ("vbf", "zh", "wh", "vh"):
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
                        divideby = 1
                        for i, matrixmultiplies in enumerate(self.signalsamples()):
                            multiplyby = 1
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
                    for j in xrange(invertedmatrix.shape[0]):
                        divideby = 1
                        for i in xrange(invertedmatrix.shape[0]):
                            multiplyby = 1
                            threshold = 1e-10 * multiplyby / divideby
                            if abs(invertedmatrix[j,i]) < threshold: invertedmatrix[j,i] = 0
                            if abs(invertedmatrix[j,i] - 1) < 1e-10: invertedmatrix[j,i] = 1

                    assert invertedmatrix[0, 0] == 1 and all(invertedmatrix[0, i] == 0 for i in             range(1, 15))
                    assert invertedmatrix[2, 1] == 1 and all(invertedmatrix[2, i] == 0 for i in range(0, 1)+range(2, 15))
                    assert invertedmatrix[5, 2] == 1 and all(invertedmatrix[5, i] == 0 for i in range(0, 2)+range(3, 15))
                    assert invertedmatrix[9, 3] == 1 and all(invertedmatrix[9, i] == 0 for i in range(0, 3)+range(4, 15))
                    assert invertedmatrix[14, 4] == 1 and all(invertedmatrix[14,i] == 0 for i in range(0, 4)+range(5, 15))
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
        c = Counter(t.templatename() for t in self.templates()+self.inttemplates())
        for k, v in c.iteritems():
            if v>1:
                raise ValueError("Multiple templates in {} with name {}".format(self, k))
        result = {
          "outputFile": self.templatesfile(firststep=False),
          "templates": sum((_.getjson() for _ in self.templates()+self.inttemplates()), []),
          "constraints": self.constraints,
        }
        if any(sum(len(filelist) for filelist in template["files"]) > 15 for template in result["templates"]): result["maxthreads"] = 1
        return result

    @property
    def copyfromothertemplatesfile(self):
        kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in type(self).needenums}
        if self.templategroup == "bkg" and self.shapesystematic in ("JEC0PMUp", "JEC0MUp", "JEC0PHUp", "JEC0L1Up", "JEC0L1ZgUp"):
            kwargs["shapesystematic"] = "JECUp"
        elif self.templategroup == "bkg" and self.shapesystematic in ("JEC0PMDn", "JEC0MDn", "JEC0PHDn", "JEC0L1Dn", "JEC0L1ZgDn"):
            kwargs["shapesystematic"] = "JECDn"
#        elif self.production == "190821_2016":
#            kwargs["production"] = "190703_2016"
#        elif self.production == "190821_2017":
#            kwargs["production"] = "190703_2017"
#        elif self.production == "190821_2018":
#            kwargs["production"] = "190703_2018"
        else:
            return None
        return type(self)(*kwargs.itervalues())

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
            yield _

    @property
    def allshapesystematics(self):
        for _ in shapesystematics:
            if not _.appliesto(self.templategroup): continue

            if _ in ("JECUp", "JECDn", "JEC0PMUp", "JEC0PMDn", "JEC0MUp", "JEC0MDn", "JEC0PHUp", "JEC0PHDn", "JEC0L1Up", "JEC0L1Dn", "JEC0L1ZgUp", "JEC0L1ZgDn"):
                if self.category not in ("VBFtagged", "VHHadrtagged"): continue

            if _.isTHUggH:
                STXSuncertainties = "Mu", "Res", "Mig01", "Mig12", "VBF2j", "VBF3j", "PT60", "PT120", "qmtop"
                if self.category == "VBFtagged":
                    indices = 1, 3, 4, 5, 6, 7, 8
                elif self.category == "VHHadrtagged":
                    indices = 3, 6, 7, 8
                elif self.category == "VBF1jtagged":
                    indices = 2, 3, 6, 7
                elif self.category == "VHLepttagged":
                    indices = 1, 2, 3, 6, 7, 8
                elif self.category == "Boosted":
                    indices = 3, 7, 8
                elif self.category == "Untagged":
                    indices = ()
                else:
                    assert False, self
                if not any(_ == "THU_ggH_" + STXSuncertainties[index] + direction for index in indices for direction in ("Up", "Dn", "0PMUp", "0PMDn")):
                    continue
            yield _

    @property
    def constraints(self):
        if self.templategroup in ("background", "DATA"): return []
        if self.analysis == "fa3_STXS": return []
        if (
          self.templategroup in ("ggh", "tth", "bbh")
          and self.category in ("VBFtagged", "VHHadrtagged")
          and self.analysis.fais == Analysis("fa3").fais
        ): return []
        if "" != self.shapesystematic is not None: return []

        productionmode = {_.productionmode for _ in self.signalsamples()}
        assert len(productionmode) == 1
        productionmode = productionmode.pop()
        templates2 = []

        if self.analysis.dimensions == 1:
            if self.templategroup in ("ggh", "tth", "bbh"):
                constrainttype = "oneparameterHVV"
                templates = [
                    Template(self, productionmode, self.analysis.purehypotheses[0]),
                    IntTemplate(self, productionmode, "g11gi1"),
                    Template(self, productionmode, self.analysis.purehypotheses[1]),
                ]
            if self.templategroup in ("vbf", "zh", "wh", "vh"):
                constrainttype = "oneparameterVVHVV"
                templates = [
                    Template(self, productionmode, self.analysis.purehypotheses[0]),
                    IntTemplate(self, productionmode, "g13gi1"),
                    IntTemplate(self, productionmode, "g12gi2"),
                    IntTemplate(self, productionmode, "g11gi3"),
                    Template(self, productionmode, self.analysis.purehypotheses[1]),
                ]
        elif self.analysis.dimensions == 4:
            if self.templategroup in ("ggh", "tth", "bbh"):
                assert self.analysis.isfa3fa2fL1fL1Zg
                if self.category == "Untagged" and not self.analysis.isSTXS:
                    constrainttype = "fourparameterHVV"
                    templates = [
                        Template(self, productionmode, self.analysis.purehypotheses[0]),
                        IntTemplate(self, productionmode, "g11gi1"),
                        IntTemplate(self, productionmode, "g11gj1"),
                        IntTemplate(self, productionmode, "g11gk1"),
                        IntTemplate(self, productionmode, "g11gl1"),
                        Template(self, productionmode, self.analysis.purehypotheses[1]),
                        IntTemplate(self, productionmode, "gi1gj1"),
                        IntTemplate(self, productionmode, "gi1gk1"),
                        IntTemplate(self, productionmode, "gi1gl1"),
                        Template(self, productionmode, self.analysis.purehypotheses[2]),
                        IntTemplate(self, productionmode, "gj1gk1"),
                        IntTemplate(self, productionmode, "gj1gl1"),
                        Template(self, productionmode, self.analysis.purehypotheses[3]),
                        IntTemplate(self, productionmode, "gk1gl1"),
                        Template(self, productionmode, self.analysis.purehypotheses[4]),
                    ]
                elif self.category in ("VBFtagged", "VHHadrtagged", "Boosted", "VBF1jtagged", "VHLepttagged") or self.analysis.isSTXS:
                    #leave out fa3, because those interferences are 0
                    constrainttype = "threeparameterHVV"
                    if self.templategroup == "tth" or self.category in ("VBFtagged", "VHHadrtagged") and self.templategroup == "ggh":
                        generator = None
                        if self.templategroup == "ggh": generator = "JHUGen"
                        templates = [
                            Template(self, productionmode, self.analysis.purehypotheses[0], "Hff0+", generator),
                            IntTemplate(self, productionmode, "g11gj1", "Hff0+"),
                            IntTemplate(self, productionmode, "g11gk1", "Hff0+"),
                            IntTemplate(self, productionmode, "g11gl1", "Hff0+"),
                            Template(self, productionmode, self.analysis.purehypotheses[2], "Hff0+", generator),
                            IntTemplate(self, productionmode, "gj1gk1", "Hff0+"),
                            IntTemplate(self, productionmode, "gj1gl1", "Hff0+"),
                            Template(self, productionmode, self.analysis.purehypotheses[3], "Hff0+", generator),
                            IntTemplate(self, productionmode, "gk1gl1", "Hff0+"),
                            Template(self, productionmode, self.analysis.purehypotheses[4], "Hff0+", generator),
                        ]
                        templates2 = [
                            Template(self, productionmode, self.analysis.purehypotheses[0], "Hff0-", generator),
                            IntTemplate(self, productionmode, "g11gj1", "Hff0-"),
                            IntTemplate(self, productionmode, "g11gk1", "Hff0-"),
                            IntTemplate(self, productionmode, "g11gl1", "Hff0-"),
                            Template(self, productionmode, self.analysis.purehypotheses[2], "Hff0-", generator),
                            IntTemplate(self, productionmode, "gj1gk1", "Hff0-"),
                            IntTemplate(self, productionmode, "gj1gl1", "Hff0-"),
                            Template(self, productionmode, self.analysis.purehypotheses[3], "Hff0-", generator),
                            IntTemplate(self, productionmode, "gk1gl1", "Hff0-"),
                            Template(self, productionmode, self.analysis.purehypotheses[4], "Hff0-", generator),
                        ]
                    else:
                        templates = [
                            Template(self, productionmode, self.analysis.purehypotheses[0]),
                            IntTemplate(self, productionmode, "g11gj1"),
                            IntTemplate(self, productionmode, "g11gk1"),
                            IntTemplate(self, productionmode, "g11gl1"),
                            Template(self, productionmode, self.analysis.purehypotheses[2]),
                            IntTemplate(self, productionmode, "gj1gk1"),
                            IntTemplate(self, productionmode, "gj1gl1"),
                            Template(self, productionmode, self.analysis.purehypotheses[3]),
                            IntTemplate(self, productionmode, "gk1gl1"),
                            Template(self, productionmode, self.analysis.purehypotheses[4]),
                        ]
            if self.templategroup in ("vbf", "zh", "wh", "vh"):
                constrainttype = "fourparameterVVHVV"
                templates = [
                    Template(self, productionmode, self.analysis.purehypotheses[0]),
                    IntTemplate(self, productionmode, "g13gi1"),
                    IntTemplate(self, productionmode, "g13gj1"),
                    IntTemplate(self, productionmode, "g13gk1"),
                    IntTemplate(self, productionmode, "g13gl1"),
                    IntTemplate(self, productionmode, "g12gi2"),
                    IntTemplate(self, productionmode, "g12gi1gj1"),
                    IntTemplate(self, productionmode, "g12gi1gk1"),
                    IntTemplate(self, productionmode, "g12gi1gl1"),
                    IntTemplate(self, productionmode, "g12gj2"),
                    IntTemplate(self, productionmode, "g12gj1gk1"),
                    IntTemplate(self, productionmode, "g12gj1gl1"),
                    IntTemplate(self, productionmode, "g12gk2"),
                    IntTemplate(self, productionmode, "g12gk1gl1"),
                    IntTemplate(self, productionmode, "g12gl2"),
                    IntTemplate(self, productionmode, "g11gi3"),
                    IntTemplate(self, productionmode, "g11gi2gj1"),
                    IntTemplate(self, productionmode, "g11gi2gk1"),
                    IntTemplate(self, productionmode, "g11gi2gl1"),
                    IntTemplate(self, productionmode, "g11gi1gj2"),
                    IntTemplate(self, productionmode, "g11gi1gj1gk1"),
                    IntTemplate(self, productionmode, "g11gi1gj1gl1"),
                    IntTemplate(self, productionmode, "g11gi1gk2"),
                    IntTemplate(self, productionmode, "g11gi1gk1gl1"),
                    IntTemplate(self, productionmode, "g11gi1gl2"),
                    IntTemplate(self, productionmode, "g11gj3"),
                    IntTemplate(self, productionmode, "g11gj2gk1"),
                    IntTemplate(self, productionmode, "g11gj2gl1"),
                    IntTemplate(self, productionmode, "g11gj1gk2"),
                    IntTemplate(self, productionmode, "g11gj1gk1gl1"),
                    IntTemplate(self, productionmode, "g11gj1gl2"),
                    IntTemplate(self, productionmode, "g11gk3"),
                    IntTemplate(self, productionmode, "g11gk2gl1"),
                    IntTemplate(self, productionmode, "g11gk1gl2"),
                    IntTemplate(self, productionmode, "g11gl3"),
                    Template(self, productionmode, self.analysis.purehypotheses[1]),
                    IntTemplate(self, productionmode, "gi3gj1"),
                    IntTemplate(self, productionmode, "gi3gk1"),
                    IntTemplate(self, productionmode, "gi3gl1"),
                    IntTemplate(self, productionmode, "gi2gj2"),
                    IntTemplate(self, productionmode, "gi2gj1gk1"),
                    IntTemplate(self, productionmode, "gi2gj1gl1"),
                    IntTemplate(self, productionmode, "gi2gk2"),
                    IntTemplate(self, productionmode, "gi2gk1gl1"),
                    IntTemplate(self, productionmode, "gi2gl2"),
                    IntTemplate(self, productionmode, "gi1gj3"),
                    IntTemplate(self, productionmode, "gi1gj2gk1"),
                    IntTemplate(self, productionmode, "gi1gj2gl1"),
                    IntTemplate(self, productionmode, "gi1gj1gk2"),
                    IntTemplate(self, productionmode, "gi1gj1gk1gl1"),
                    IntTemplate(self, productionmode, "gi1gj1gl2"),
                    IntTemplate(self, productionmode, "gi1gk3"),
                    IntTemplate(self, productionmode, "gi1gk2gl1"),
                    IntTemplate(self, productionmode, "gi1gk1gl2"),
                    IntTemplate(self, productionmode, "gi1gl3"),
                    Template(self, productionmode, self.analysis.purehypotheses[2]),
                    IntTemplate(self, productionmode, "gj3gk1"),
                    IntTemplate(self, productionmode, "gj3gl1"),
                    IntTemplate(self, productionmode, "gj2gk2"),
                    IntTemplate(self, productionmode, "gj2gk1gl1"),
                    IntTemplate(self, productionmode, "gj2gl2"),
                    IntTemplate(self, productionmode, "gj1gk3"),
                    IntTemplate(self, productionmode, "gj1gk2gl1"),
                    IntTemplate(self, productionmode, "gj1gk1gl2"),
                    IntTemplate(self, productionmode, "gj1gl3"),
                    Template(self, productionmode, self.analysis.purehypotheses[3]),
                    IntTemplate(self, productionmode, "gk3gl1"),
                    IntTemplate(self, productionmode, "gk2gl2"),
                    IntTemplate(self, productionmode, "gk1gl3"),
                    Template(self, productionmode, self.analysis.purehypotheses[4]),
                ]
                if self.analysis.isSTXS or self.category in ("Boosted", "VBF1jtagged", "VHLepttagged"):
                    for i, _ in reversed(list(enumerate(templates[:]))):
                        if isinstance(_, IntTemplate) and _.interferencetype.couplingpowers["i"] in (1, 3):
                            del templates[i]
                    constrainttype = "fourparameterVVHVV_nog4int"
            if self.templategroup == ("wh"):
                constrainttype = "fourparameterWWHVV"
                templates = [
                    Template(self, productionmode, self.analysis.purehypotheses[0]),
                    IntTemplate(self, productionmode, "g13gi1"),
                    IntTemplate(self, productionmode, "g13gj1"),
                    IntTemplate(self, productionmode, "g13gk1"),
                    IntTemplate(self, productionmode, "g13gl1"),
                    IntTemplate(self, productionmode, "g12gi2"),
                    IntTemplate(self, productionmode, "g12gi1gj1"),
                    IntTemplate(self, productionmode, "g12gi1gk1"),
                    IntTemplate(self, productionmode, "g12gi1gl1"),
                    IntTemplate(self, productionmode, "g12gj2"),
                    IntTemplate(self, productionmode, "g12gj1gk1"),
                    IntTemplate(self, productionmode, "g12gj1gl1"),
                    IntTemplate(self, productionmode, "g12gk2"),
                    IntTemplate(self, productionmode, "g12gk1gl1"),
                    IntTemplate(self, productionmode, "g12gl2"),
                    IntTemplate(self, productionmode, "g11gi3"),
                    IntTemplate(self, productionmode, "g11gi2gj1"),
                    IntTemplate(self, productionmode, "g11gi2gk1"),
                    IntTemplate(self, productionmode, "g11gi2gl1"),
                    IntTemplate(self, productionmode, "g11gi1gj2"),
                    IntTemplate(self, productionmode, "g11gi1gj1gk1"),
                    IntTemplate(self, productionmode, "g11gi1gj1gl1"),
                    IntTemplate(self, productionmode, "g11gi1gk2"),
                    IntTemplate(self, productionmode, "g11gi1gk1gl1"),
                    IntTemplate(self, productionmode, "g11gi1gl2"),
                    IntTemplate(self, productionmode, "g11gj3"),
                    IntTemplate(self, productionmode, "g11gj2gk1"),
                    IntTemplate(self, productionmode, "g11gj2gl1"),
                    IntTemplate(self, productionmode, "g11gj1gk2"),
                    IntTemplate(self, productionmode, "g11gj1gk1gl1"),
                    IntTemplate(self, productionmode, "g11gj1gl2"),
                    IntTemplate(self, productionmode, "g11gk3"),
                    IntTemplate(self, productionmode, "g11gk2gl1"),
                    IntTemplate(self, productionmode, "g11gk1gl2"),
                    Template(self, productionmode, self.analysis.purehypotheses[1]),
                    IntTemplate(self, productionmode, "gi3gj1"),
                    IntTemplate(self, productionmode, "gi3gk1"),
                    IntTemplate(self, productionmode, "gi3gl1"),
                    IntTemplate(self, productionmode, "gi2gj2"),
                    IntTemplate(self, productionmode, "gi2gj1gk1"),
                    IntTemplate(self, productionmode, "gi2gj1gl1"),
                    IntTemplate(self, productionmode, "gi2gk2"),
                    IntTemplate(self, productionmode, "gi2gk1gl1"),
                    IntTemplate(self, productionmode, "gi2gl2"),
                    IntTemplate(self, productionmode, "gi1gj3"),
                    IntTemplate(self, productionmode, "gi1gj2gk1"),
                    IntTemplate(self, productionmode, "gi1gj2gl1"),
                    IntTemplate(self, productionmode, "gi1gj1gk2"),
                    IntTemplate(self, productionmode, "gi1gj1gk1gl1"),
                    IntTemplate(self, productionmode, "gi1gj1gl2"),
                    IntTemplate(self, productionmode, "gi1gk3"),
                    IntTemplate(self, productionmode, "gi1gk2gl1"),
                    IntTemplate(self, productionmode, "gi1gk1gl2"),
                    Template(self, productionmode, self.analysis.purehypotheses[2]),
                    IntTemplate(self, productionmode, "gj3gk1"),
                    IntTemplate(self, productionmode, "gj3gl1"),
                    IntTemplate(self, productionmode, "gj2gk2"),
                    IntTemplate(self, productionmode, "gj2gk1gl1"),
                    IntTemplate(self, productionmode, "gj2gl2"),
                    IntTemplate(self, productionmode, "gj1gk3"),
                    IntTemplate(self, productionmode, "gj1gk2gl1"),
                    IntTemplate(self, productionmode, "gj1gk1gl2"),
                    Template(self, productionmode, self.analysis.purehypotheses[3]),
                    IntTemplate(self, productionmode, "gk3gl1"),
                    IntTemplate(self, productionmode, "gk2gl2"),
                ]
                if self.analysis.isSTXS or self.category in ("Boosted", "VBF1jtagged", "VHLepttagged"):
                    for i, _ in reversed(list(enumerate(templates[:]))):
                        if isinstance(_, IntTemplate) and _.interferencetype.couplingpowers["i"] in (1, 3):
                            del templates[i]
                    constrainttype = "fourparameterWWHVV_nog4int"

        for _ in templates+templates2:
            if _ not in {Template: self.templates(), IntTemplate: self.inttemplates()}[type(_)]:
                raise ValueError("{} is not a template in {}. Templates:\n  {}".format(_, self, "\n  ".join(str(t) for t in self.templates() + self.inttemplates())))
        return [
            {
                "type": constrainttype,
                "templates": [_.templatename() for _ in templateslist],
            } for templateslist in (templates, templates2) if templateslist
        ]

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
                    if category == "Boosted" and not analysis.useboosted: continue
                    if category in ("VBF1jtagged", "VHLepttagged") and not analysis.usemorecategories: continue
                    for templategroup in templategroups:
                        if analysis.isdecayonly and templategroup not in ("bkg", "ggh", "DATA"): continue
                        nominal = TemplatesFile(channel, templategroup, analysis, production, category)
                        for shapesystematic in nominal.treeshapesystematics:
                            if (production.LHE or production.GEN) and shapesystematic != "": continue
                            if category not in ("VBFtagged", "VHHadrtagged") and shapesystematic in ("JECUp", "JECDn", "MINLO_SM"): continue

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
            elif enumsdict[ProductionMode] == "VH":
                enumsdict[TemplateGroup] = "vh"
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
            if enumsdict[ProductionMode] in ("ggH", "ttH"):
                enumsdict[HffHypothesis] = "Hff0+"

        if enumsdict[ShapeSystematic] is None:
            enumsdict[ShapeSystematic] = ""

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
    def bkgdiscriminant(self):
        return self.templatesfile.bkgdiscriminant
    @property
    def purediscriminant(self):
        return self.templatesfile.purediscriminant
    @property
    def mixdiscriminant(self):
        return self.templatesfile.mixdiscriminant
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

        elif self.productionmode in ("ggH", "VBF", "ZH", "WH", "VH", "ttH", "bbH"):
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
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

        if self.productionmode in ("ggH", "ttH"):
            if self.hffhypothesis == "Hff0-" and self.category not in ("VBFtagged", "VHHadrtagged") and self.productionmode != "ttH":
                raise ValueError("{} {} is not used to make templates for {}!\n{}".format(self.productionmode, self.hffhypothesis, self.category, args))
            if self.hffhypothesis not in ("Hff0+", "Hff0-"):
                raise ValueError("{} {} is not used to make templates!\n{}".format(self.productionmode, self.hffhypothesis, args))
        else:
            if self.hffhypothesis is not None:
               raise ValueError("HffHypothesis {} provided for {}!\n{}".format(self.hffhypothesis, self.productionmode, args))

    def templatename(self, final=True):
        if self.productionmode in ("ggH", "ttH", "bbH"):
            match1 = re.match("^f(a2|a3|L1|L1Zg)(?:dec)?(-?)0.5$", str(self.hypothesis))
            match2 = re.match("^f(a2|a3|L1|L1Zg)(?:dec)?0.5f(a2|a3|L1|L1Zg)(?:dec)?(-?)0.5$", str(self.hypothesis))
            if self.hypothesis == "0+":
                name = "template0Plus"
            elif self.hypothesis == "0-":
                name = "template0Minus"
            elif self.hypothesis == "a2":
                name = "template0HPlus"
            elif self.hypothesis == "L1":
                name = "template0L1"
            elif self.hypothesis == "L1Zg":
                name = "template0L1Zg"
            elif match1:
                name = "templateMixa1{}{}".format(match1.group(1), "Pi" if match1.group(2) else "")
            elif match2:
                name = "templateMix{}{}{}".format(match2.group(1), match2.group(2), "Pi" if match2.group(3) else "")
            else:
                print self.hypothesis
            if self.hffhypothesis == "Hff0-":
                name = name.replace("template", "templateHff0MinusHVV")
        elif self.productionmode in ("VBF", "ZH", "WH", "VH"):
            match1 = re.match("^f(?P<ai>a2|a3|L1|L1Zg)(?P<proddec>(?:prod|dec|proddec)?)(?P<sign>-?)0[.]5$", str(self.hypothesis))
            match2 = re.match("^f(?P<ai>a2|a3|L1|L1Zg)(?P<proddec>(?:prod|dec|proddec)?)0[.]5f(?P<aj>a2|a3|L1|L1Zg)(?P=proddec)(?P<sign>-?)0[.]5$", str(self.hypothesis))
            match3 = re.match("^f(?P<ai>a2|a3|L1|L1Zg)(?P<proddec>(?:prod|dec|proddec)?)0[.]33f(?P<aj>a2|a3|L1|L1Zg)(?P=proddec)(?P<sign>-?)0[.]33$", str(self.hypothesis))
            match4 = re.match("^f(?P<ai>a2|a3|L1|L1Zg)(?P<proddec>(?:prod|dec|proddec)?)0[.]33f(?P<aj>a2|a3|L1|L1Zg)(?P=proddec)0[.]33f(?P<ak>a2|a3|L1|L1Zg)(?P=proddec)(?P<sign>-?)0[.]33$", str(self.hypothesis))
            match5 = re.match("^f(?P<ai>a2|a3|L1|L1Zg)(?P<proddec>(?:prod|dec|proddec)?)0[.]25f(?P<aj>a2|a3|L1|L1Zg)(?P=proddec)0[.]25f(?P<ak>a2|a3|L1|L1Zg)(?P=proddec)(?P<sign>-?)0[.]25$", str(self.hypothesis))
            match6 = re.match("^f(?P<ai>a2|a3|L1|L1Zg)(?P<proddec>(?:prod|dec|proddec)?)0[.]25f(?P<aj>a2|a3|L1|L1Zg)(?P=proddec)0[.]25f(?P<ak>a2|a3|L1|L1Zg)(?P=proddec)0[.]25f(?P<al>a2|a3|L1|L1Zg)(?P=proddec)(?P<sign>-?)0[.]25$", str(self.hypothesis))
            pddict = {"prod": "Prod", "dec": "Decay", "": "Decay", "proddec": "ProdDec"}
            if self.hypothesis == "0+":
                name = "template0Plus"
            elif self.hypothesis == "0-":
                name = "template0Minus"
            elif self.hypothesis == "a2":
                name = "template0HPlus"
            elif self.hypothesis == "L1":
                name = "template0L1"
            elif self.hypothesis == "L1Zg":
                name = "template0L1Zg"
            elif match1:
                name = "templateMixa1{}{}{}".format(
                    match1.group("ai"),
                    pddict[match1.group("proddec")],
                    "Pi" if match1.group("sign") else ""
                )
            elif match2:
                name = "templateMix{}{}{}{}".format(
                    match2.group("ai"),
                    match2.group("aj"),
                    pddict[match2.group("proddec")],
                    "Pi" if match2.group("sign") else ""
                )
            elif match3:
                name = "templateMixa1{}{}{}{}".format(
                    match3.group("ai"),
                    match3.group("aj"),
                    pddict[match3.group("proddec")],
                    "Pi" if match3.group("sign") else ""
                )
            elif match4:
                name = "templateMix{}{}{}{}{}".format(
                    match4.group("ai"),
                    match4.group("aj"),
                    match4.group("ak"),
                    pddict[match4.group("proddec")],
                    "Pi" if match4.group("sign") else ""
                )
            elif match5:
                name = "templateMixa1{}{}{}{}{}".format(
                    match5.group("ai"),
                    match5.group("aj"),
                    match5.group("ak"),
                    pddict[match5.group("proddec")],
                    "Pi" if match5.group("sign") else ""
                )
            elif match6:
                name = "templateMix{}{}{}{}{}{}".format(
                    match6.group("ai"),
                    match6.group("aj"),
                    match6.group("ak"),
                    match6.group("al"),
                    pddict[match6.group("proddec")],
                    "Pi" if match6.group("sign") else ""
                )
        elif self.productionmode == "ZX" and not config.usedata:
            name = "templateqqZZ"
        elif self.productionmode in ("ggZZ", "qqZZ", "VBF bkg", "ZX"):
            name = "template{}".format(self.productionmode)
        elif self.productionmode == "data":
            name = "datatemplate"

        name

        if self.domirror and final:
            name += "Mirror"

        return name

    def title(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "VH", "ttH", "bbH"):
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
        if self.productionmode == "VH":
            samples = [ReweightingSample("ZH", self.hypothesis)]
            if self.hypothesis != "L1Zg":
                samples.append(ReweightingSample("WH", self.hypothesis))
        elif self.productionmode == "ggZZ":
            samples = [self.reweightingsampleplus]*6
        else:
            samples = [self.reweightingsampleplus]

        multiplyweight = ""
        if self.multiplyweight is not None:
            multiplyweight = " * (" + self.multiplyweight + ")"

        return ["MC_weight_nominal * (" + sample.MC_weight + ")" + multiplyweight for sample in samples]

    @property
    def multiplyweight(self):
        if self.shapesystematic in ("", "JECUp", "JECDn", "ScaleUp", "ScaleDn", "ResUp", "ResDn"):
            return None
        if self.shapesystematic.isTHUggH:
            return AlternateWeight(str(self.shapesystematic)).weightname
        assert False, self

    @property
    def categoryname(self):
        result = "category_"
        result += self.analysis.categoryname
        if self.shapesystematic is not None: result += self.shapesystematic.categoryappendname

        from treewrapper import TreeWrapper
        for categorization in TreeWrapper.categorizations:
            if categorization.category_function_name == result: return result
        assert False, "{} does not exist in TreeWrapper".format(result)

    @cache
    def reweightfrom(self):
        if self.productionmode == "ggH":
            if self.production.LHE:
                result = [{Sample(self.production, self.productionmode, self.hypothesis)}]
            elif self.shapesystematic == "MINLO_SM":
                result = [{Sample(self.production, self.productionmode, self.hypothesis, "MINLO")}]
            else:
                if self.category in ("VBFtagged", "VHHadrtagged"):
                    generators = ["JHUGen"]
                else:
                    generators = [None]
                result=[{
                        Sample(self.production, self.productionmode, hypothesis, ext, generator, hffhypothesis)
                            for hypothesis in self.productionmode.generatedhypotheses(self.production)
                            for hffhypothesis in hffhypotheses
                            for ext in (None, "ext1")
                            for generator in generators
                            if (hypothesis == "0+" or generator == None)
                            and (hffhypothesis == "Hff0+" or generator == "JHUGen")
                       }]
        if self.productionmode in ["VBF", "ZH"]:
            result=[{
                    Sample(self.production, self.productionmode, hypothesis, ext)
                    for hypothesis in self.productionmode.generatedhypotheses(self.production)
                    for ext in (None, "ext1")
                   }]
            if self.productionmode == "VBF" and self.reweightingsample.hasZZ:
              result = [{s for s in st if s.hasZZ} for st in result]
            if self.productionmode == "ZH" and (self.reweightingsample.hasZZ or self.reweightingsample.hasZg):
              result = [{s for s in st if s.hasZZ or s.hasZg} for st in result]
        if self.productionmode == "WH":
            result=[{
                    Sample(self.production, self.productionmode, hypothesis, ext)
                    for hypothesis in self.productionmode.generatedhypotheses(self.production)
                    for ext in (None, "ext1")
                   }]
        if self.productionmode == "VH":
            result = [
                {
                    Sample(self.production, "ZH", hypothesis, ext)
                    for hypothesis in ProductionMode("ZH").generatedhypotheses(self.production)
                    for ext in (None, "ext1")
                }, {
                    Sample(self.production, "WH", hypothesis, ext)
                    for hypothesis in ProductionMode("WH").generatedhypotheses(self.production)
                    for ext in (None, "ext1")
                }
            ]
            if self.hypothesis == "L1Zg":
                del result[1]
        if self.productionmode in ("ttH", "bbH"):
            result=[{
                    Sample(self.production, self.productionmode, "0+", self.hffhypothesis),
                    Sample(self.production, self.productionmode, "0+", self.hffhypothesis, "ext1"),
                   }]
        if self.productionmode == "ZX":
            result = [{Sample(self.production, self.productionmode)}]
        if self.productionmode == "qqZZ":
            result = [{Sample(self.production, self.productionmode, ext) for ext in (None, "ext", "ext1", "ext2")}]
        if self.productionmode == "ggZZ":
            result = [{Sample(self.production, self.productionmode, flavor)} for flavor in flavors]
        if self.productionmode == "VBF bkg":
            result = [{Sample(self.production, self.productionmode, flavor)} for flavor in flavors if not flavor.hastaus]
        if self.productionmode == "data":
            if config.showblinddistributions:
                result = [{Sample(self.production, self.productionmode)}]
            elif self.production.GEN or self.production.LHE:
                result = [{Sample("WplusH", "POWHEG", "ext", self.production, "0+")}]
            else:
                result = [{Sample("ggZZ", self.production, "4tau")}]

        for st in result:
            for sample in set(st):
                if not withdiscriminantsfileisvalid(sample.withdiscriminantsfile()):
                    st.remove(sample)
        assert result and all(result), self
        return result

    @property
    def scalefactor(self): return None

    @property
    def domirror(self):
        if "fa3" not in self.analysis.fais: return False
        if self.analysis.isSTXS or self.category in ("Boosted", "VBF1jtagged", "VHLepttagged"): return False
        if self.productionmode == "data": return False

        assert "D_CP" in self.mixdiscriminant.name, (self, self.mixdiscriminant.name)

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

    @smoothingparameters.deleter
    def smoothingparameters(self):
      self.initsmoothingparameters()
      self.__smoothingparameters.delvalue()

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
    def postprocessingjson(self):
      result = []
      if self.scalefactor is not None:
        result.append({"type": "rescale", "factor": self.scalefactor})
      result.append({"type": "floor", "floorvalue": 1e-10})
      if self.domirror:
        result.append({"type": "mirror", "antisymmetric": False})
      return result

    def getjson(self):
        import datetime; print "   ", self, datetime.datetime.now()
        if self.copyfromothertemplate: return []
        if self.hypothesis is not None and not self.hypothesis.ispure: return []
        if self.shapesystematic == "":
            pass
        elif self.shapesystematic in ("JECUp", "JECDn"):
            pass
        elif self.shapesystematic in ("ScaleUp", "ScaleDn", "ResUp", "ResDn") or self.shapesystematic.isTHUggH:
            if self.hypothesis != "0+" or None is not self.hffhypothesis != "Hff0+": return []
        else:
            assert False, self
        jsn = [
               {
                "name": self.templatename(final=True),
                "files": [sorted([sample.withdiscriminantsfile() for sample in st]) for st in self.reweightfrom()],
                "tree": "candTree",
                "variables": [d.formula for d in self.discriminants],
                "weight": self.weightname(),
                "selection": self.selection,
                "binning": {
                  "bins": self.binning,
                },
                "postprocessing": self.postprocessingjson,
               },
              ]

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
            idnumbers = self.category.idnumbers
            if not self.analysis.useboosted and self.category == "Untagged":
                idnumbers |= Category("Boosted").idnumbers
            if not self.analysis.usemorecategories and self.category == "Untagged":
                idnumbers |= Category("VBF1jtagged").idnumbers
                idnumbers |= Category("VHLepttagged").idnumbers
            result.append("(" + " || ".join("{} == {}".format(self.categoryname, c) for c in idnumbers) + ")")
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

        if self.productionmode in ("ggH", "ttH"):
            if self.hffhypothesis == "Hff0-" and self.category not in ("VBFtagged", "VHHadrtagged") and self.productionmode != "ttH":
                raise ValueError("{} {} is not used to make templates for {}!\n{}".format(self.productionmode, self.hffhypothesis, self.category, args))
            if self.hffhypothesis not in ("Hff0+", "Hff0-"):
                raise ValueError("{} {} is not used to make templates!\n{}".format(self.productionmode, self.hffhypothesis, args))
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
        elif self.productionmode in ("VBF", "ZH", "WH", "VH"):
            if self.interferencetype.totalpower != 4:
                raise ValueError("Invalid interferencetype {} for productionmode {}!\n{}".format(self.interferencetype, self.productionmode, args))
        else:
            raise ValueError("Invalid productionmode {}!\n{}".format(self.productionmode, args))

    def templatename(self, final=True):
        if self.analysis.dimensions == 4:
            result = "template"+str(self.interferencetype)+"Int"
            result = result.replace("gi", self.analysis.couplingnames[0])
            result = result.replace("gj", self.analysis.couplingnames[1])
            result = result.replace("gk", self.analysis.couplingnames[2])
            result = result.replace("gl", self.analysis.couplingnames[3])
        elif self.productionmode in ("ggH", "ttH", "bbH"):
            if self.interferencetype == "g11gi1":
                if self.analysis.fais == ("fa3",):
                    result = "templatea1a3Int"
                elif self.analysis == "fa2":
                    result = "templatea1a2Int"
                elif self.analysis == "fL1" or self.analysis.isfL1fL1Zg:
                    result = "templatea1L1Int"
                elif self.analysis == "fL1Zg":
                    result = "templatea1L1ZgInt"
                else:
                    assert False
            elif self.interferencetype == "g11gj1":
                if self.analysis.isfL1fL1Zg:
                    result = "templatea1L1ZgInt"
            elif self.interferencetype == "gi1gj1":
                if self.analysis.isfL1fL1Zg:
                    result = "templateL1L1ZgInt"
        elif self.productionmode in ("VBF", "ZH", "WH", "VH"):
            if self.interferencetype == "g11gi3":
                result = "templateg11{}3".format(self.analysis.couplingname)
            if self.interferencetype == "g12gi2":
                result = "templateg12{}2".format(self.analysis.couplingname)
            if self.interferencetype == "g13gi1":
                result = "templateg13{}1".format(self.analysis.couplingname)

        if self.hffhypothesis == "Hff0-":
            result = result.replace("template", "templateHff0MinusHVV")

        if self.domirror:
            result += "Mirror"
        return result

    @property
    def mirrorjsn(self):
        if "fa3" not in self.analysis.fais: return None

        #cross talk - production discriminants for the wrong category don't make sense
        if self.category in ("VBFtagged", "VHHadrtagged") and self.productionmode in ("ggH", "ttH", "bbH"):
            if (self.interferencetype == "g11gi1"
                or self.analysis.isfa3fa2fL1fL1Zg and self.interferencetype.couplingpowers["i"] == 1):
                #ggH has no production information, and only using SM ttH, so mirror antisymmetric
                #over the (pretend) D_CP_decay axis, which sets the whole thing to 0
                #note for ttH this is an approximation, since we could have H(0-)->2l2q tt->bbllnunu
                return {"type":"rescale", "factor":0}
            if self.analysis.isfa3fa2fL1fL1Zg and self.interferencetype.couplingpowers["i"] == 0:
                if self.analysis.isSTXS or self.category in ("Boosted", "VBF1jtagged", "VHLepttagged"): return None
                return {"type":"mirror", "antisymmetric":False, "axis":1}
            assert False

        if (self.analysis.isSTXS or self.category in ("Boosted", "VBF1jtagged", "VHLepttagged")) and self.interferencetype.couplingpowers["i"] in (1, 3):
            #same (antimirror over D_CP_whatever, which doesn't exist in STXS)
            assert "fa3" == self.analysis.fais[0]
            return {"type":"rescale", "factor":0}

        if self.analysis.isSTXS or self.category in ("Boosted", "VBF1jtagged", "VHLepttagged"): return None

        #Mirror antisymmetric for VH in VBF category and VBF in VH category
        #the interference templates are 0 to within error bars anyway,
        #except for some effects in ZH g13g41 D_CP_VBF which are antisymmetric

        #cross talk to the untagged category is exactly correct, since the decay is the same

        if (self.interferencetype in ("g11gi1", "g11gi3", "g13gi1")
            or self.analysis.isfa3fa2fL1fL1Zg and self.interferencetype.couplingpowers["i"] in (1, 3)):
            return {"type":"mirror", "antisymmetric":True, "axis":1}
        elif (self.interferencetype == "g12gi2"
              or self.analysis.isfa3fa2fL1fL1Zg and self.interferencetype.couplingpowers["i"] in (0, 2, 4)):
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
                if self.analysis.dimensions == 1:
                    #using new PDF class that uses a1=ai=1 for interference
                    multiplyby = 1
                else:
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
        vectoroftemplates = [t for t in self.templatesfile.templates() if t.hffhypothesis == self.hffhypothesis]  #ok technically it's a list not a vector
        templatesandfactors = []

        for j, template in enumerate(vectoroftemplates):
            factor = invertedmatrix[rowofinvertedmatrix,j] * multiplyby
            if factor:
                templatesandfactors.append((template, factor))

        return templatesandfactors

    def getjson(self):
        import datetime; print "   ", self, datetime.datetime.now()
        if self.shapesystematic != "": return []
        intjsn = [
          {
            "name": self.templatename(final=True),
            "files": [sorted([sample.withdiscriminantsfile() for sample in st]) for st in self.reweightfrom()],
            "tree": "candTree",
            "variables": [d.formula for d in self.discriminants],
            "weight": self.weightname(),
            "selection": self.selection,
            "binning": {
              "bins": self.binning,
            },
            "postprocessing": [],
          },
        ]

        if self.domirror:
          mirrorjsn = self.mirrorjsn
          if "axis" in mirrorjsn: del mirrorjsn["axis"]
          intjsn[0]["postprocessing"].append(mirrorjsn)
        return intjsn

    @cache
    def gettemplate(self, *args, **kwargs):
        result = super(IntTemplate, self).gettemplate(*args, **kwargs)
        if self.analysis == "fa3" and self.productionmode in ("ggH", "ttH", "bbH") and self.interferencetype == "g11gi1" and self.category in ("VBFtagged", "VHHadrtagged"):
            assert self.domirror
            result.Scale(0)
        return result

    @property
    def sumsofsamples(self):
        if self.productionmode == "VH":
            result = []
            kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in type(self).needenums}
            kwargs["productionmode"] = "ZH"
            kwargs["templategroup"] = "zh"
            ZHtemplate = type(self)(*kwargs.itervalues())
            result.append(sum((t.reweightingsampleplus*factor for t, factor in ZHtemplate.templatesandfactors), SumOfSamples()))
            
            if not self.analysis.isfa3fa2fL1fL1Zg: assert False
            if self.interferencetype.couplingpowers["l"] < 3:
                kwargs["productionmode"] = "WH"
                kwargs["templategroup"] = "wh"
                WHtemplate = type(self)(*kwargs.itervalues())
                result.append(sum((t.reweightingsampleplus*factor for t, factor in WHtemplate.templatesandfactors), SumOfSamples()))
            return result
        return [sum((t.reweightingsampleplus*factor for t, factor in self.templatesandfactors), SumOfSamples())]

    def weightname(self):
        return ["MC_weight_nominal * (" + sumofsamples.MC_weight + ")" for sumofsamples in self.sumsofsamples]

    @cache
    def reweightfrom(self):
        if self.productionmode == "VH":
            if not self.analysis.isfa3fa2fL1fL1Zg: assert False
            kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in self.enums}
            del kwargs["interferencetype"]; kwargs["hypothesis"] = "0+"
            result = Template(*kwargs.itervalues()).reweightfrom()[:]
            if self.interferencetype.couplingpowers["l"] >= 3:
                del result[1]
        else:
            others = [t.reweightfrom() for t, factor in self.templatesandfactors]
            assert all(len(other) == 1 for other in others)
            result = [set.intersection(*(other[0] for other in others))]
        assert result and all(result), result
        return result

    @property
    def selection(self):
        if self.productionmode == "VH":
            kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in self.enums}
            del kwargs["interferencetype"]; kwargs["hypothesis"] = "0+"
            return Template(*kwargs.itervalues()).selection
        result = {t.selection for t, factor in self.templatesandfactors}
        assert len(result) == 1, result
        return result.pop()

class DataTree(MultiEnum):
    enums = [Channel, Production, Category, Analysis]
    enumname = "datatree"
    def check(self, *args):
        super(DataTree, self).check(*args)
        if self.analysis.isdecayonly:
            if self.category != "Untagged":
                raise ValueError("decayonly analysis is only done for untagged!\n{}".format(args))
        if not self.analysis.useboosted:
            if self.category == "Boosted":
                raise ValueError("{} doesn't use boosted\n{}".format(self.analysis, args))
        if not self.analysis.usemorecategories:
            if self.category in ("VBF1jtagged", "VHLepttagged"):
                raise ValueError("{} doesn't use {}\n{}".format(self.analysis, self.category, args))
        if self.production.LHE:
            if self.channel != "2e2mu":
                raise ValueError("LHE analysis is only done for 2e2mu for now!\n{}".format(args))

    @property
    def originaltreefile(self):
        return Sample("data", self.production).withdiscriminantsfile()
    @property
    def treefile(self):
        return os.path.join(config.repositorydir, "step7_templates", str(self.production), "data_{}_{}_{}_{}.root".format(self.production, self.channel, self.category, self.analysis))
    def passescut(self, t):
        return (
          abs(t.Z1Flav * t.Z2Flav) == self.channel.ZZFlav
          and config.m4lmin < t.ZZMass < config.m4lmax
          and config.unblinddistributions
          and (
            getattr(t, "category_"+self.analysis.categoryname) in self.category
            or (
              getattr(t, "category_"+self.analysis.categoryname) in Category("Boosted")
              and not self.analysis.useboosted
              and self.category == "Untagged"
            ) or (
              (
                getattr(t, "category_"+self.analysis.categoryname) in Category("VBF1jtagged")
                or getattr(t, "category_"+self.analysis.categoryname) in Category("VHLepttagged")
              )
              and not self.analysis.usemorecategories
              and self.category == "Untagged"
            )
          )
        )

@listfromiterator
def datatrees():
    for channel in channels:
        for production in productions:
            if channel != "2e2mu" and production.LHE: continue
            for category in categories:
                for analysis in analyses:
                    if category != "Untagged" and analysis.isdecayonly: continue
                    if category == "Boosted" and not analysis.useboosted: continue
                    if category in ("VBF1jtagged", "VHLepttagged") and not analysis.usemorecategories: continue
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
