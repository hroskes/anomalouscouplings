import config
from enums import Flavor, hypotheses, Hypothesis, MultiEnum, ProductionMode
import os
import ROOT

class Sample(MultiEnum):
    enums = [ProductionMode, Hypothesis, Flavor]

    def check(self, *args):
        if self.productionmode is None:
            raise ValueError("No option provided for productionmode\n{}".format(args))
        elif self.productionmode == "ggH":
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for ggH productionmode\n{}".format(args))
            if self.flavor is not None:
                raise ValueError("Flavor provided for ggH productionmode\n{}".format(args))
        elif self.productionmode == "ggZZ":
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for ggZZ productionmode\n{}".format(args))
            if self.flavor is None:
                raise ValueError("No flavor provided for ggZZ productionmode\n{}".format(args))
        elif self.productionmode == "qqZZ":
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for qqZZ productionmode\n{}".format(args))
            if self.flavor is not None:
                raise ValueError("Flavor provided for qqZZ productionmode\n{}".format(args))
        else:
            raise ValueError("Bad productionmode {}\n{}".format(self.productionmode, args))

    def ValueError(self, functionname):
        return ValueError("invalid sample {} in function {}".format(self, functionname))

    def CJLSTfile(self):
        if self.productionmode == "ggH":
            dirnames = {Hypothesis(k): v for k, v in config.CJLSTdirnames_ggH.iteritems()}
            return os.path.join(config.CJLSTmaindir, dirnames[self.hypothesis], "ZZ4lAnalysis.root")
        elif self.productionmode == "ggZZ":
            return "root://lxcms03//data3/Higgs/160225/ggZZ{}/ZZ4lAnalysis.root".format(self.flavor)
        elif self.productionmode == "qqZZ":
            return "root://lxcms03//data3/Higgs/160225/ZZTo4l/ZZ4lAnalysis.root"
        raise self.ValueError("CJLSTfile")

    def withdiscriminantsfile(self):
        return os.path.join(config.repositorydir, "step3_withdiscriminants", "{}.root".format(self).replace(" ", ""))
        raise self.ValueError("withdiscriminantsfile")

    def reweightingsamples(self):
        if self.productionmode == "ggH":
            return [Sample("ggH", hypothesis) for hypothesis in hypotheses]
        elif self.productionmode == "ggZZ":
            return [self]
        elif self.productionmode == "qqZZ":
            return [self]
        raise self.ValueError("reweightingsamples")

    def isbkg(self):
        if self.productionmode == "ggH":
            return False
        elif self.productionmode == "ggZZ":
            return True
        elif self.productionmode == "qqZZ":
            return True
        raise self.ValueError("isbkg")

    def isdata(self):
        if self.productionmode == "ggH":
            return False
        elif self.productionmode == "ggZZ":
            return False
        elif self.productionmode == "qqZZ":
            return False
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
        elif self.productionmode == "ggZZ":
            return "MC_weight_ggZZ"
        elif self.productionmode == "qqZZ":
            return "MC_weight_qqZZ"
        raise self.ValueError("weightname")

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
