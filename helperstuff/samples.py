import config
from enums import Flavor, hypotheses, Hypothesis, MultiEnum, ProductionMode, Release
import os
import ROOT

class ReweightingSample(MultiEnum):
    enumname = "reweightingsample"
    enums = [ProductionMode, Hypothesis, Flavor]

    def __init__(self, *args, **kwargs):
        super(ReweightingSample, self).__init__(*args, **kwargs)
        if self.productionmode == "ZX":
            import ZX

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
    enums = [ReweightingSample, Release]

    def check(self, *args):
        if self.release is None:
            raise ValueError("No option provided for release\n{}".format(args))
        super(Sample, self).check(*args)

    def CJLSTmaindir(self):
        if self.productionmode == "ggH" and self.release == "76X":
            return config.CJLSTmaindir_76Xanomalous
        return self.release.CJLSTdir()

    def CJLSTdirname(self):
        if self.productionmode == "ggH" and self.release == "76X":
            if self.hypothesis == "0+": return "0PM_v2"
            if self.hypothesis == "a2": return "0PH_v2"
            if self.hypothesis == "0-": return "0M_v1"
            if self.hypothesis == "L1": return "0L1_v2"
            if self.hypothesis == "fa20.5": return "0PHf05ph0_v2"
            if self.hypothesis == "fa30.5": return "0Mf05ph0_v2"
            if self.hypothesis == "fL10.5": return "0L1f05ph0_v2"
        if self.productionmode == "ggH" and self.release == "80X":
            if self.hypothesis == "0+": return "0PM"
            if self.hypothesis == "a2": return "0PH"
            if self.hypothesis == "0-": return "0M"
            if self.hypothesis == "L1": return "0L1"
            if self.hypothesis == "fa20.5": return "0PHf05ph0"
            if self.hypothesis == "fa30.5": return "0Mf05ph0"
            if self.hypothesis == "fL10.5": return "0L1f05ph0"
        if self.productionmode in ("VBF", "ZH", "WplusH", "WminusH", "ttH"):
            name = str(self.productionmode)
            if self.productionmode == "VBF":
                name += "H"
            return name
        if self.productionmode == "ggZZ":
            return "ggZZ{}".format(self.flavor)
        if self.productionmode == "qqZZ":
            return "ZZTo4l"
        if self.productionmode == "ZX" or self.productionmode == "data":
            return "AllData"
        raise self.ValueError("CJLSTdirname")

    def CJLSTfile(self):
        return os.path.join(self.CJLSTmaindir(), self.CJLSTdirname(), "ZZ4lAnalysis.root")
        if self.productionmode == "ggH":
            dirnames = {Hypothesis(k): v for k, v in config.CJLSTdirnames_ggH.iteritems()}
            if self.release == "76X":
                return os.path.join(config.CJLSTmaindir_76Xanomalous, dirnames[self.hypothesis], "ZZ4lAnalysis.root")
            else:
                return os.path.join(self.release.CJLSTdir(), dirnames[self.hypothesis], "ZZ4lAnalysis.root")
        elif self.productionmode in ("VBF", "ZH", "WplusH", "WminusH", "ttH"):
            name = str(self.productionmode)
            if self.productionmode == "VBF":
                name += "H"
            return os.path.join(self.release.CJLSTdir(), name+"125", "ZZ4lAnalysis.root")
        elif self.productionmode == "ggZZ":
            return os.path.join(self.release.CJLSTdir(), "ggZZ{}".format(self.flavor), "ZZ4lAnalysis.root")
        elif self.productionmode == "qqZZ":
            return os.path.join(self.release.CJLSTdir(), "ZZTo4l", "ZZ4lAnalysis.root")
        elif self.productionmode == "ZX" or self.productionmode == "data":
            return os.path.join(self.release.CJLSTdir(), "AllData", "ZZ4lAnalysis.root")
        raise self.ValueError("CJLSTfile")

    def withdiscriminantsfile(self):
        return os.path.join(config.repositorydir, "step3_withdiscriminants", "{}.root".format(self).replace(" ", ""))
        raise self.ValueError("withdiscriminantsfile")
