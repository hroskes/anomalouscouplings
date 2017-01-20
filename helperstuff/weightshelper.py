#!/usr/bin/env python
from enums import Analysis, MultiEnum, ProductionMode

class WeightsHelper(MultiEnum):
    enums = (ProductionMode, Analysis)

    def __init__(self, *args, **kwargs):
        from samples import SampleBase
        if len(args) == 1 and not kwargs and isinstance(args[0], SampleBase):
            sample = args[0]
            if   sample.g2 == sample.g1prime2 == sample.ghzgs1prime2 == 0: args = (sample.productionmode, "fa3")
            elif sample.g4 == sample.g1prime2 == sample.ghzgs1prime2 == 0: args = (sample.productionmode, "fa2")
            elif sample.g2 == sample.g4       == sample.ghzgs1prime2 == 0: args = (sample.productionmode, "fL1")
            elif sample.g2 == sample.g4       == sample.ghz1prime2   == 0: args = (sample.productionmode, "fL1Zgs")
            else: assert False
        return super(WeightsHelper, self).__init__(*args, **kwargs)


    def check(self, *args):
        if self.productionmode in ("WplusH", "WminusH"): self.productionmode = ProductionMode("WH")
        if self.productionmode not in ("ggH", "VBF", "ZH", "WH", "ttH", "HJJ"):
            raise ValueError("Has to be a signal productionmode, not {}!\n{}".format(self.productionmode, args))

    @property
    def weightdecaystring(self):
        if self.productionmode == "ggH":
            return "GG"
        if self.productionmode in ("VBF", "ZH", "WH", "ttH", "HJJ"):
            return "Dec"
    @property
    def weightprodstring(self):
        if self.productionmode == "ggH":
            return None
        if self.productionmode == "ttH": #production reweighting for ttH doesn't work unfortunately
            return None
        if self.productionmode in ("VBF", "ZH", "WH", "ttH", "HJJ"):
            return str(self.productionmode)

    @property
    def decaySMcouplingname(self):
        return "ghz1"
    @property
    def decayBSMcouplingname(self):
        if self.analysis == "fL1Zg": return self.analysis.couplingname
        return self.analysis.couplingname.replace("g", "ghz")
    @property
    def decayBSMcouplingvalue(self):
        if self.analysis == "fL1": return "1E4"
        else: return "1"

    @property
    def prodSMcouplingname(self):
        if self.productionmode in ("ggH", "HJJ"): return "ghg2"
        if self.productionmode == "VBF": return "ghv1"
        if self.productionmode == "WH": return "ghw1"
        if self.productionmode == "ZH": return "ghz1"
        if self.productionmode == "ttH": return "kappa"
    @property
    def prodBSMcouplingname(self):
        if self.productionmode == "ttH": return "kappa_tilde"
        if self.productionmode == "HJJ": return "ghg4"
        if self.productionmode == "ggH": return None
        if self.analysis == "fL1Zg": return self.analysis.couplingname
        if self.productionmode == "VBF": return self.analysis.couplingname.replace("g", "ghv")
        if self.productionmode == "WH": return self.analysis.couplingname.replace("g", "ghw")
        if self.productionmode == "ZH": return self.analysis.couplingname.replace("g", "ghz")
    @property
    def prodBSMcouplingvalue(self):
        if self.analysis in ["fL1", "fL1Zg"] and self.productionmode in ("VBF", "WH", "ZH"): return "1E4"
        else: return "1"

    @property
    def formatdict(self):
        return {
                "decay": self.weightdecaystring,
                "prod": self.weightprodstring,
                "decaySM": self.decaySMcouplingname,
                "prodSM": self.prodSMcouplingname,
                "decayBSM": self.decayBSMcouplingname,
                "prodBSM": self.prodBSMcouplingname,
                "decayBSMvalue": self.decayBSMcouplingvalue,
                "prodBSMvalue": self.prodBSMcouplingvalue,
               }

    @property
    def decayweightSM(self):
        result = "p_Gen_{decay}_SIG_"
        if self.productionmode=="ggH":
            result += "{prodSM}_1_"
        result += "{decaySM}_1_JHUGen"
        return result.format(**self.formatdict)
    @property
    def decayweightBSM(self):
        result = "p_Gen_{decay}_SIG_"
        if self.productionmode=="ggH":
            result += "{prodSM}_1_"
        result += "{decayBSM}_{decayBSMvalue}_JHUGen"
        return result.format(**self.formatdict)
    @property
    def decayweightmix(self):
        result = "p_Gen_{decay}_SIG_"
        if self.productionmode=="ggH":
            result += "{prodSM}_1_"
        result += "{decaySM}_1_{decayBSM}_{decayBSMvalue}_JHUGen".format(**self.formatdict)
        return result.format(**self.formatdict)
    @property
    def decaySMcoupling(self):
        return "g1"
    @property
    def decayBSMcoupling(self):
        return self.analysis.couplingname

    @property
    def prodweightSM(self):
        return "p_Gen_{prod}_SIG_{prodSM}_1_JHUGen".format(**self.formatdict)
    @property
    def prodweightBSM(self):
        return "p_Gen_{prod}_SIG_{prodBSM}_{prodBSMvalue}_JHUGen".format(**self.formatdict)
    @property
    def prodweightmix(self):
        return "p_Gen_{prod}_SIG_{prodSM}_1_{prodBSM}_{prodBSMvalue}_JHUGen".format(**self.formatdict)
    @property
    def prodSMcoupling(self):
        if self.productionmode in ("VBF", "ZH", "WH"): return "g1"
        if self.productionmode == "HJJ":               return "ghg2"
        if self.productionmode == "ttH":               return "kappa"
    @property
    def prodBSMcoupling(self):
        if self.productionmode in ("VBF", "ZH", "WH"): return self.analysis.couplingname
        if self.productionmode == "HJJ":               return "ghg4"
        if self.productionmode == "ttH":               return "kappa_tilde"

if __name__ == "__main__":
    from samples import ArbitraryCouplingsSample, ReweightingSample
    print ArbitraryCouplingsSample("ttH", g1=1, g2=0, g4=0, g1prime2=12345, ghzgs1prime2=0, kappa=1, kappa_tilde=4).MC_weight
