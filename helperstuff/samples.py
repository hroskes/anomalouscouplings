import config
from enums import BlindStatus, Flavor, hypotheses, Hypothesis, MultiEnum, ProductionMode, Production
import os
import ROOT

class ReweightingSample(MultiEnum):
    enumname = "reweightingsample"
    enums = [ProductionMode, Hypothesis, Flavor]

    def check(self, *args):
        if self.productionmode is None:
            raise ValueError("No option provided for productionmode\n{}".format(args))
        elif self.productionmode == "ggH":
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for ggH productionmode\n{}".format(args))
            if self.flavor is not None:
                raise ValueError("Flavor provided for ggH productionmode\n{}".format(args))
        elif self.productionmode in ("VBF", "ZH", "WplusH", "WminusH", "ttH"):
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hypothesis != "0+":
                raise ValueError("Bad hypothesis {} provided for {} productionmode\n{}".format(self.hypothesis, self.productionmode, args))
            if self.flavor is not None:
                raise ValueError("Flavor provided for {} productionmode\n{}".format(self.productionmode, args))
        elif self.productionmode == "ggZZ":
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for ggZZ productionmode\n{}".format(args))
            if self.flavor is None:
                raise ValueError("No flavor provided for ggZZ productionmode\n{}".format(args))
        elif self.productionmode == "qqZZ" or self.productionmode == "ZX" or self.productionmode == "data":
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.flavor is not None:
                raise ValueError("Flavor provided for {} productionmode\n{}".format(self.productionmode, args))
        else:
            raise ValueError("Bad productionmode {}\n{}".format(self.productionmode, args))

    def ValueError(self, functionname):
        return ValueError("invalid sample {} in function {}".format(self, functionname))

    def reweightingsamples(self):
        if self.productionmode == "ggH":
            return [ReweightingSample("ggH", hypothesis) for hypothesis in hypotheses]
        elif self.productionmode in ("ggZZ", "qqZZ", "ZX", "VBF", "ZH", "WplusH", "WminusH", "ttH"):
            return [self]
        elif self.productionmode == "data":
            return []
        raise self.ValueError("reweightingsamples")

    def isbkg(self):
        if self.productionmode in ("ggH", "data", "VBF", "ZH", "WplusH", "WminusH", "ttH"):
            return False
        elif self.productionmode in ("ggZZ", "qqZZ", "ZX"):
            return True
        raise self.ValueError("isbkg")

    def isZX(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "data", "VBF", "ZH", "WplusH", "WminusH", "ttH"):
            return False
        elif self.productionmode == "ZX":
            return True
        raise self.ValueError("isZX")

    def isdata(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "ZX", "VBF", "ZH", "WplusH", "WminusH", "ttH"):
            return False
        elif self.productionmode == "data":
            return True
        raise self.ValueError("isdata")
        
    def weightname(self):
        if self.productionmode == "ggH":
            if self.hypothesis == "0+":
                return "MC_weight_ggH_g1"
            elif self.hypothesis == "0-":
                return "MC_weight_ggH_g4"
            elif self.hypothesis == "fa30.5":
                return "MC_weight_ggH_g1g4"
            elif self.hypothesis == "a2":
                return "MC_weight_ggH_g2"
            elif self.hypothesis == "fa20.5":
                return "MC_weight_ggH_g1g2"
            elif self.hypothesis == "L1":
                return "MC_weight_ggH_g1prime2"
            elif self.hypothesis == "fL10.5":
                return "MC_weight_ggH_g1g1prime2"
        elif self.productionmode == "VBF":
            return "MC_weight_VBF"
        elif self.productionmode == "ZH":
            return "MC_weight_ZH"
        elif self.productionmode in ("WplusH", "WminusH"):
            return "MC_weight_WH"
        elif self.productionmode == "ttH":
            return "MC_weight_ttH"
        elif self.productionmode == "ggZZ":
            return "MC_weight_ggZZ"
        elif self.productionmode == "qqZZ":
            return "MC_weight_qqZZ"
        elif self.productionmode == "ZX":
            return "MC_weight_ZX"
        raise self.ValueError("weightname")

    def TDirectoryname(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "data", "VBF", "ZH", "WplusH", "WminusH", "ttH"):
            return "ZZTree"
        if self.productionmode == "ZX":
            return "CRZLLTree"
        raise self.ValueError("TDirectoryname")

    def color(self):
        if self.productionmode == "ggH":
            if self.hypothesis == "0+":
                return 1
            elif self.hypothesis == "0-":
                return 2
            elif self.hypothesis == "fa30.5":
                return ROOT.kRed+3
            elif self.hypothesis == "a2":
                return ROOT.kCyan
            elif self.hypothesis == "fa20.5":
                return 4
            elif self.hypothesis == "L1":
                return 3
            elif self.hypothesis == "fL10.5":
                return ROOT.kGreen+3
        raise self.ValueError("color")

    def onlyweights(self):
        """True if this sample is not needed for making templates,
           and only the weight and ZZMass should be recorded in the tree"""
        if self.productionmode in ("VBF", "ZH", "WplusH", "WminusH", "ttH"):
            return True
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "ZX", "data"):
            return False
        raise self.ValueError("onlyweights")

    def weightingredients(self):
        """only needs to be defined if self.onlyweights() is True.
           gives a list of variable names that are used to calculate
           the weight"""
        if self.productionmode in ("VBF", "ZH", "WplusH", "WminusH", "ttH"):
            return "overallEventWeight", "xsec"
        raise self.ValueError("weightingredients")

class Sample(ReweightingSample):
    enums = [ReweightingSample, Production, BlindStatus]

    def check(self, *args):
        if self.production is None:
            raise ValueError("No option provided for production\n{}".format(args))
        if self.blindstatus is None and self.productionmode == "data":
            raise ValueError("No blindstatus provided for data sample!\n{}".format(args))
        if self.blindstatus is not None and self.productionmode != "data":
            raise ValueError("blindstatus provided for MC sample!\n{}".format(args))
        super(Sample, self).check(*args)

    def CJLSTmaindir(self):
        if self.productionmode == "ggH":
            return self.production.CJLSTdir_anomalous()
        if self.productionmode in ("data", "ZX"):
            return self.production.CJLSTdir_data()
        return self.production.CJLSTdir()

    def CJLSTdirname(self):
        if self.productionmode == "ggH" and self.production == "160225":
            if self.hypothesis == "0+": return "0PM_v2"
            if self.hypothesis == "a2": return "0PH_v2"
            if self.hypothesis == "0-": return "0M_v1"
            if self.hypothesis == "L1": return "0L1_v2"
            if self.hypothesis == "fa20.5": return "0PHf05ph0_v2"
            if self.hypothesis == "fa30.5": return "0Mf05ph0_v2"
            if self.hypothesis == "fL10.5": return "0L1f05ph0_v2"
        if self.productionmode == "ggH" and self.production == "160624":
            if self.hypothesis == "0+": return "0PM"
            if self.hypothesis == "a2": return "0PH"
            if self.hypothesis == "0-": return "0M"
            if self.hypothesis == "L1": return "0L1"
            if self.hypothesis == "fa20.5": return "0PHf05ph0"
            if self.hypothesis == "fa30.5": return "0Mf05ph0"
            if self.hypothesis == "fL10.5": return "0L1f05ph0"
        if self.productionmode == "ggH" and self.production in ("160714", "160720", "160725"):
            if self.hypothesis == "0+": return "ggH0PM_M125"
            if self.hypothesis == "a2": return "ggH0PH_M125"
            if self.hypothesis == "0-": return "ggH0M_M125"
            if self.hypothesis == "L1": return "ggH0L1_M125"
            if self.hypothesis == "fa20.5": return "ggH0PHf05ph0_M125"
            if self.hypothesis == "fa30.5": return "ggH0Mf05ph0_M125"
            if self.hypothesis == "fL10.5": return "ggH0L1f05ph0_M125"
        if self.productionmode in ("VBF", "ZH", "WplusH", "WminusH", "ttH"):
            name = str(self.productionmode)
            if self.productionmode == "VBF":
                name += "H"
            name += "125"
            return name
        if self.productionmode == "ggZZ":
            if self.production == "160714" and self.flavor in ("4e", "4mu", "4tau", "2mu2tau") or self.production in ("160720", "160725"):
                return "ggTo{}_Contin_MCFM701".format(self.flavor)
            elif self.production in ("160225", "160624", "160714"):
                return "ggZZ{}".format(self.flavor)
        if self.productionmode == "qqZZ":
            return "ZZTo4l"
        if self.productionmode == "ZX" or self.productionmode == "data":
            return "AllData"
        raise self.ValueError("CJLSTdirname")

    def CJLSTfile(self):
        return os.path.join(self.CJLSTmaindir(), self.CJLSTdirname(), "ZZ4lAnalysis.root")

    def withdiscriminantsfile(self):
        result = os.path.join(config.repositorydir, "step3_withdiscriminants", "{}.root".format(self).replace(" ", ""))
        return result
        raise self.ValueError("withdiscriminantsfile")

    @property
    def useMELAv2(self):
        return self.production.useMELAv2

    @property
    def unblind(self):
        if self.productionmode == "data":
            return self.blindstatus == "unblind"
        raise self.ValueError("unblind")
