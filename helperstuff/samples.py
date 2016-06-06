import config
from enums import hypotheses, Hypothesis, MultiEnum, ProductionMode
import os

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
        
