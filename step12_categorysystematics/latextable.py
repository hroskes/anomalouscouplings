#!/usr/bin/env python

from itertools import product

from helperstuff import config

from helperstuff.enums import categories, Category, flavors, JECSystematic, ProductionMode
from helperstuff.samples import ReweightingSample, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import keydefaultdict, MultiplyCounter, tfiles

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

def gettree(productionmode):
    t = ROOT.TChain("candTree")
    for sample in productionmode.allsamples(production):
        t.Add(sample.withdiscriminantsfile())
    return t

trees = KeyDefaultDict(gettree)

class Row(MultiEnum):
    enums = (Productionmode, Hypothesis, Analysis, Channel)
    def check(self, *args):
        dontcheck = []
        if self.productionmode.isbkg:
            dontcheck.append(Hypothesis)
        super(Row, self).check(*args, dontcheck=dontcheck)

    @property
    def tree(self): return trees[self.productionmode]
    @property
    def reweightingsample(self): return ReweightingSample(self.productionmode, self.hypothesis)

    def scalefactor(self):
        result = 1
        if self.productionmode in ("VBF", "ZH", "WH"):
            result *= (
                        (ReweightingSample(self.productionmode, self.hypothesis).xsec 
                            / sum(ReweightingSample(_, self.hypothesis).xsec for _ in ("VBF", "ZH", "WH")))
                       /
                        (ReweightingSample(self.productionmode, "SM"           ).xsec 
                            / sum(ReweightingSample(_, "SM"           ).xsec for _ in ("VBF", "ZH", "WH")))
                      )
        if self.productionmode in ("VBF", "ZH", "WH", "ggH", "ttH"):
            result /= sum(
                      Sample.effectiveentries(
                                              reweightfrom=reweightfrom
                                              reweightto=self.reweightingsample
                                             )
                       for reweightfrom in self.productionmode.allsamples(production)
                     )
        return result

    def getcategorydistribution(self):
        t = self.tree
        weightparts = [
                       "ZZMass>{}"
        t.Draw("1>>h", "" && ZZMass<{} && Z1Flav*Z2Flav == {}".format(config.m4lmin, config.m4lmax, self.ZZFl

    @property
    def ZZFlav(self):
        result = self.channel.ZZFlav
        if self.productionmode == "ZX": result *= -1
        return result
