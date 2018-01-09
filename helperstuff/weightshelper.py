#!/usr/bin/env python
from enums import Analysis, MultiEnum, ProductionMode

class WeightsHelper(MultiEnum):
    enums = (ProductionMode,)

    def __init__(self, *args, **kwargs):
        from samples import SampleBase
        if len(args) == 1 and not kwargs and isinstance(args[0], SampleBase):
            sample = args[0]
            args = sample.productionmode
        return super(WeightsHelper, self).__init__(*args, **kwargs)

    def check(self, *args):
        if self.productionmode in ("WplusH", "WminusH"): self.productionmode = ProductionMode("WH")
        if self.productionmode not in ("ggH", "VBF", "ZH", "WH", "ttH", "HJJ"):
            raise ValueError("Has to be a signal productionmode, not {}!\n{}".format(self.productionmode, args))

    def weightstring(self, prodordec):
        if prodordec == "dec":
            if self.productionmode == "ggH":
                return "GG"
            if self.productionmode in ("VBF", "ZH", "WH", "ttH", "HJJ"):
                return "Dec"
        if prodordec == "prod":
            if self.productionmode == "ggH":
                return None
            if self.productionmode == "ttH": #production reweighting for ttH doesn't work unfortunately
                return None
            if self.productionmode in ("VBF", "ZH", "WH", "ttH", "HJJ"):
                return str(self.productionmode)

    def useproddec(self, prodordec):
      return bool(self.weightstring(prodordec))
      assert False

    def allcouplings(self, prodordec):
      if not self.useproddec(prodordec): return None
      if prodordec == "dec" or prodordec == "prod" and self.productionmode == "ZH":
        return "ghz1", "ghz1_prime2", "ghz2", "ghz4", "ghza1_prime2"
      if prodordec == "prod" and self.productionmode == "VBF":
        return "ghv1", "ghv1_prime2", "ghv2", "ghv4", "ghza1_prime2"
      if prodordec == "prod" and self.productionmode == "WH":
        return "ghw1", "ghw1_prime2", "ghw2", "ghw4"
      if prodordec == "prod" and self.productionmode == "ttH":
        return "kappa", "kappa_tilde"
      if prodordec == "prod" and self.productionmode == "HJJ":
        return "ghg2", "ghg4"
      assert False

    @staticmethod
    def couplingvalue(coupling):
      if "_prime2" in coupling: return "1E4"
      return "1"

    def couplingsandweights(self, prodordec, mix):
      if not mix:
        for coupling in self.allcouplings(prodordec):
          couplingvalue = self.couplingvalue(coupling)
          dct = {
            "weightstring": self.weightstring(prodordec),
            "coupling": coupling,
            "couplingvalue": couplingvalue,
          }
        yield coupling, couplingvalue, self.weight(prodordec).format(**dct)
      else:
        for (coupling1, coupling1value, weight1), (coupling2, coupling2value, weight2) in combinations(self.couplingsandweights(prodordec, False), 2):
          dct = {
            "weightstring": self.weightstring(prodordec),
            "coupling1": coupling1,
            "coupling1value": coupling1value,
            "coupling2": coupling2,
            "coupling2value": coupling2value,
          }
          yield (coupling1, coupling1value, weight1), (coupling2, coupling2value, weight2), self.weightmix(prodordec).format(**dct)

    @property
    def weight(self):
        result = "p_Gen_{weightstring}_SIG_"
        if prodordec == "dec" and self.productionmode=="ggH":
            result += "ghg2_1_"
        result += "{coupling}_{couplingvalue}_JHUGen"
        return result
    @property
    def weightmix(self):
        result = "p_Gen_{weightstring}_SIG_"
        if prodordec == "dec" and self.productionmode=="ggH":
            result += "ghg2_1_"
        result += "{coupling1}_{coupling1value}_{coupling2}_{coupling2value}_JHUGen".format(**self.formatdict)
        return result

if __name__ == "__main__":
    from samples import ArbitraryCouplingsSample, ReweightingSample
    print ArbitraryCouplingsSample("ttH", g1=1, g2=0, g4=0, g1prime2=12345, ghzgs1prime2=0, kappa=1, kappa_tilde=4).MC_weight
    print ArbitraryCouplingsSample("ggH", g1=0, g2=0, g4=0, g1prime2=12345, ghzgs1prime2=23456, ghg2=1, ghg4=0).MC_weight

