from enums import hypotheses, Hypothesis, MultiEnum, ProductionMode
import os

class Sample(MultiEnum):
    enums = [ProductionMode, Hypothesis]

    def ValueError(self, functionname):
        return ValueError("invalid sample {} in function {}".format(self, functionname))

    def CJLSTfile(self):
        if self.productionmode == "ggH":
            maindir = "/afs/cern.ch/work/h/hroskes/reweighting_CJLST/CMSSW_7_6_3_patch2/src/ZZAnalysis/AnalysisStep/test/prod/AnomalousCouplingsReweighting/PT13TeV"
            dirnames = {
                        "0+": "0PM_v2",
                        "a2": "0PH_v2",
                        "0-": "0M_v1",
                        "L1": "0L1_v2",
                        "fa20.5": "0PHf05ph0_v2",
                        "fa30.5": "0Mf05ph0_v2",
                        "fL10.5": "0L1f05ph0_v2",
                       }
            dirnames = {Hypothesis(k): v for k, v in dirnames.iteritems()}
            return os.path.join(maindir, dirnames[self.hypothesis], "ZZ4lAnalysis.root")
        raise self.ValueError("CJLSTfile")

    def withdiscriminantsfile(self):
        maindir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings/withdiscriminants"
        return os.path.join(maindir, "{}.root".format(self).replace(" ", ""))

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
        
