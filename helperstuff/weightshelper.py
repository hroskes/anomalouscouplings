#!/usr/bin/env python
from itertools import combinations

from helperstuff import config
from helperstuff.enums import Analysis, MultiEnum, ProductionMode

class WeightsHelper(MultiEnum):
    enums = (ProductionMode,)

    def __init__(self, *args, **kwargs):
        from samples import SampleBase
        if len(args) == 1 and not kwargs and isinstance(args[0], SampleBase):
            sample = args[0]
            args = sample.productionmode,
            if hasattr(sample, "alternategenerator") and sample.alternategenerator in ("JHUGen", "MCatNLO") and "useHJJ" not in kwargs:
                kwargs["useHJJ"] = True

        self.__useHJJ = kwargs.pop("useHJJ", False)

        super(WeightsHelper, self).__init__(*args, **kwargs)

        if self.__useHJJ and self.productionmode != "ggH":
            raise ValueError("useHJJ only makes sense for ggH")

    def check(self, *args):
        if self.productionmode in ("WplusH", "WminusH"): self.productionmode = ProductionMode("WH")
        if self.productionmode not in ("ggH", "VBF", "ZH", "WH", "ttH", "bbH", "tqH"):
            raise ValueError("Has to be a signal productionmode, not {}!\n{}".format(self.productionmode, args))

    def weightstring(self, prodordec):
        if prodordec == "dec":
            if self.productionmode == "ggH" and not self.__useHJJ:
                return "GG"
            if self.productionmode in ("VBF", "ZH", "WH", "ttH", "bbH", "ggH"):
                return "Dec"
        if prodordec == "prod":
            if self.productionmode == "ggH" and self.__useHJJ:
                return "HJJ"
            if self.productionmode in ("ggH", "bbH"):
                return None
            if self.productionmode == "ttH":
                return "ttH"
            if self.productionmode in ("VBF", "ZH", "WH", "ttH"):
                return str(self.productionmode)

        assert False, (self, prodordec)

    def useproddec(self, prodordec):
      return bool(self.weightstring(prodordec))
      assert False

    def allcouplings(self, prodordec):
      if not self.useproddec(prodordec): return None
      if prodordec == "dec" or prodordec == "prod" and self.productionmode == "ZH":
        return "ghz1", "ghz1prime2", "ghz2", "ghz4", "ghza1prime2", "ghza2", "ghza4", "gha2", "gha4"
      if prodordec == "prod" and self.productionmode == "VBF":
        if config.separateZZWWVBFweights:
          return "ghv1", "ghz1prime2", "ghw1prime2", "ghz2", "ghw2", "ghz4", "ghw4", "ghza1prime2", "ghza2", "ghza4", "gha2", "gha4"
        else:
          return "ghv1", "ghv1prime2", "ghv2", "ghv4", "ghza1prime2", "ghza2", "ghza4", "gha2", "gha4"
      if prodordec == "prod" and self.productionmode == "WH":
        return "ghw1", "ghw1prime2", "ghw2", "ghw4"
      if prodordec == "prod" and self.productionmode == "ttH":
        return "kappa", "kappa_tilde"
      if prodordec == "prod" and self.productionmode == "ggH":
        return "ghg2", "ghg4"
      assert False

    @staticmethod
    def couplingname(coupling):
      if coupling == "ghza1prime2": return "ghzgs1prime2"
      if coupling == "ghza2": return "ghzgs2"
      if coupling == "ghza4": return "ghzgs4"
      if coupling == "gha2": return "ghgsgs2"
      if coupling == "gha4": return "ghgsgs4"
      return coupling.replace("ghv", "g")

    @staticmethod
    def couplingvalue(coupling):
      if "prime2" in coupling: return "1E4"
      return "1"

    def couplingsandweights(self, prodordec, mix, __forrecursivecall=False):
      if not mix:
        for coupling in self.allcouplings(prodordec):
          couplingvalue = self.couplingvalue(coupling)
          dct = {
            "weightstring": self.weightstring(prodordec),
            "coupling": coupling,
            "couplingvalue": couplingvalue,
          }
          if not __forrecursivecall: coupling = self.couplingname(coupling)
          yield coupling, couplingvalue, self.weight(prodordec).format(**dct)
      else:
        assert not __forrecursivecall
        for (coupling1, coupling1value, weight1), (coupling2, coupling2value, weight2) in combinations(self.couplingsandweights(prodordec, False, True), 2):
          dct = {
            "weightstring": self.weightstring(prodordec),
            "coupling1": coupling1,
            "coupling1value": coupling1value,
            "coupling2": coupling2,
            "coupling2value": coupling2value,
          }
          if coupling1.replace("z", "w") == coupling2 and coupling1value == coupling2value:
            dct["coupling"] = dct["coupling1"].replace("z", "v")
            dct["couplingvalue"] = dct["coupling1value"]
            wt = self.weight(prodordec).format(**dct)
          else:
            wt = self.weightmix(prodordec).format(**dct)

          coupling1 = self.couplingname(coupling1)
          coupling2 = self.couplingname(coupling2)
          yield (coupling1, coupling1value, weight1), (coupling2, coupling2value, weight2), wt

    def weight(self, prodordec):
        result = "p_Gen_{weightstring}_SIG_"
        if prodordec == "dec" and self.productionmode=="ggH" and not self.__useHJJ:
            result += "ghg2_1_"
        result += "{coupling}_{couplingvalue}_JHUGen"
        return result
    def weightmix(self, prodordec):
        result = "p_Gen_{weightstring}_SIG_"
        if prodordec == "dec" and self.productionmode=="ggH" and not self.__useHJJ:
            result += "ghg2_1_"
        result += "{coupling1}_{coupling1value}_{coupling2}_{coupling2value}_JHUGen"
        return result

if __name__ == "__main__":
    from helperstuff.samples import ArbitraryCouplingsSample, ReweightingSample
    from pprint import pprint
    pprint(ArbitraryCouplingsSample("VBF", g1=1, g2=0, g4=0, g1prime2=0, ghzgs1prime2=0, ghzgs2=1, ghzgs4=0, ghgsgs2=0, ghgsgs4=1).MC_weight_terms)
