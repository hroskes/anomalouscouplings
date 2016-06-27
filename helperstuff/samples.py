import config
from enums import Flavor, hypotheses, Hypothesis, MultiEnum, ProductionMode
import os
import ROOT

class Sample(MultiEnum):
    enums = [ProductionMode, Hypothesis, Flavor]

    def __init__(self, *args, **kwargs):
        super(Sample, self).__init__(*args, **kwargs)
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

    def CJLSTfile(self):
        if self.productionmode == "ggH":
            dirnames = {Hypothesis(k): v for k, v in config.CJLSTdirnames_ggH.iteritems()}
            return os.path.join(config.CJLSTmaindir, dirnames[self.hypothesis], "ZZ4lAnalysis.root")
####################################################################################################
#Temporary while 80X samples are not there yet
#remove this section once the samples are there
        elif self.productionmode == "ggZZ" and self.flavor in ["2e2mu", "2e2tau"]:
            return "root://lxcms03//data3/Higgs/160225/ggZZ{}/ZZ4lAnalysis.root".format(self.flavor)
        elif self.productionmode in ("WplusH", "WminusH"):
            return "root://lxcms03//data3/Higgs/160225/{}125/ZZ4lAnalysis.root".format(self.productionmode)
#####################################################################################################
        elif self.productionmode in ("VBF", "ZH", "WplusH", "WminusH", "ttH"):
            name = str(self.productionmode)
            if self.productionmode == "VBF":
                name += "H"
            return "root://lxcms03//data3/Higgs/160624/{}125/ZZ4lAnalysis.root".format(name)
        elif self.productionmode == "ggZZ":
            return "root://lxcms03//data3/Higgs/160624/ggZZ{}/ZZ4lAnalysis.root".format(self.flavor)
        elif self.productionmode == "qqZZ":
            return "root://lxcms03//data3/Higgs/160624/ZZTo4l/ZZ4lAnalysis.root"
        elif self.productionmode == "ZX" or self.productionmode == "data":
            return "root://lxcms03//data3/Higgs/160624/AllData/ZZ4lAnalysis.root"
        raise self.ValueError("CJLSTfile")

    def withdiscriminantsfile(self):
        return os.path.join(config.repositorydir, "step3_withdiscriminants", "{}.root".format(self).replace(" ", ""))
        raise self.ValueError("withdiscriminantsfile")

    def reweightingsamples(self):
        if self.productionmode == "ggH":
            return [Sample("ggH", hypothesis) for hypothesis in hypotheses]
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
