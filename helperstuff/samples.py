import config
from enums import hypotheses, Hypothesis, MultiEnum, ProductionMode
import os
import ROOT

class Sample(MultiEnum):
    enums = [ProductionMode, Hypothesis]

    def ValueError(self, functionname):
        return ValueError("invalid sample {} in function {}".format(self, functionname))

    def CJLSTfile(self):
        if self.productionmode == "ggH":
            dirnames = {Hypothesis(k): v for k, v in config.CJLSTdirnames_ggH.iteritems()}
            return os.path.join(config.CJLSTmaindir, dirnames[self.hypothesis], "ZZ4lAnalysis.root")
        raise self.ValueError("CJLSTfile")

    def withdiscriminantsfile(self):
        return os.path.join(config.repositorydir, "step3_withdiscriminants", "{}.root".format(self).replace(" ", ""))

    def reweightingsamples(self):
        if self.productionmode == "ggH":
            return [Sample("ggH", hypothesis) for hypothesis in hypotheses]
        raise self.ValueError("reweightingsamples")

    def isbkg(self):
        if self.productionmode == "ggH":
            return False
        raise self.ValueError("isbkg")

    def isdata(self):
        if self.productionmode == "ggH":
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
                return "MC_weight_ggH_g1g2"
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
