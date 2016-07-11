from collections import Counter, Iterator
import config
import constants
import ROOT
from samples import ReweightingSample

class TreeWrapper(Iterator):

    def __init__(self, tree, treesample, Counters, Counters_reweighted, minevent=0, maxevent=None, isdummy=False):
        """
        tree - a TTree object
        treesample - which sample the TTree was created from
        Counters_reweighted - from the CJLST file
        """
        self.tree = tree
        self.treesample = treesample
        self.productionmode = str(treesample.productionmode)
        self.hypothesis = str(treesample.hypothesis)
        self.isbkg = treesample.isbkg()
        self.isdata = treesample.isdata()
        self.isZX = treesample.isZX()
        self.isdummy = isdummy
        self.baseweight = None
        if not self.isdata:
            self.baseweight = self.getbaseweightfunction()
        self.weightfunctions = [self.getweightfunction(sample) for sample in treesample.reweightingsamples()]
        self.nevents = self.nevents2L2l = None
        if Counters is not None:
            self.nevents = Counters.GetBinContent(1)
        if Counters_reweighted is not None:
            self.nevents2L2l = [
                                Counters_reweighted.GetBinContent(4, i) #2e2mu
                              + Counters_reweighted.GetBinContent(8, i) #2e2tau+2mu2tau
                                  for i, sample in enumerate(treesample.reweightingsamples(), start=1)
                               ]
        self.minevent = minevent
        if self.isdata and not config.usedata or self.isdummy:
            self.length = 0
        elif maxevent is None or maxevent >= tree.GetEntries():
            self.length = tree.GetEntries() - minevent
        else:
            self.length = maxevent - minevent + 1

        self.initlists()
        if treesample.onlyweights(): self.onlyweights()
        tree.GetEntry(0)
        self.xsec = tree.xsec * 1000 #pb to fb
        self.checkfunctions()

    def __iter__(self):
        self.__i = 0                               #at the beginning of next self.__i and self.__treeentry are
        self.__treeentry = self.minevent-1         # incremented, so they start at 1 and self.minevent
        return super(TreeWrapper, self).__iter__() # for the first entry

    def next(self):
        while True:
            self.__i += 1
            self.__treeentry += 1
            i, t = self.__i, self.tree
            t.GetEntry(self.__treeentry)
            if i > self.length:
                raise StopIteration
            if i % 10000 == 0 or i == self.length:
                print i, "/", self.length

            if self.isdata:
                self.MC_weight = 1
            elif self.isZX:
                self.MC_weight = 1
            else:
                self.MC_weight = t.overallEventWeight
            if self.productionmode == "ggH":
                self.reweightingweights = t.reweightingweights
            isSelected = bool(self.MC_weight)

            self.flavor = self.flavordict[abs(t.Z1Flav*t.Z2Flav)]

            if isSelected:
                break

        #I prefer this to defining __getattr__ because it's faster
        self.p0plus_VAJHU = t.p0plus_VAJHU
        self.p0minus_VAJHU = t.p0minus_VAJHU
        self.pg1g4_VAJHU = t.pg1g4_VAJHU
        self.p0plus_m4l = t.p0plus_m4l
        self.bkg_VAMCFM = t.bkg_VAMCFM
        self.bkg_m4l = t.bkg_m4l
        for a in "p0plus", "bkg":
            for b in "Scale", "Res":
                for c in "Up", "Down":
                    attr = "{}_m4l_{}{}".format(a, b, c)
                    setattr(self, attr, getattr(t, attr))

        self.p0hplus_VAJHU = t.p0hplus_VAJHU
        self.pg1g2_VAJHU = t.pg1g2_VAJHU
        self.p0_g1prime2_VAJHU = t.p0_g1prime2_VAJHU
        self.pg1g1prime2_VAJHU = t.pg1g1prime2_VAJHU

        self.vbf_p0plus_VAJHU = t.pvbf_VAJHU_new
        self.hjj_p0plus_VAJHU = t.phjj_VAJHU_new

        #express in terms of |M|^2, this will make life easier
        self.M2g1_decay   = self.p0plus_VAJHU
        self.M2g4_decay   = self.p0minus_VAJHU / constants.CJLSTg4decay_pure[self.flavor]**2
        self.M2g1g4_decay = self.pg1g4_VAJHU / constants.CJLSTg4decay_mix
        self.M2g2_decay   = self.p0hplus_VAJHU / constants.CJLSTg2decay_pure[self.flavor]**2
        self.M2g1g2_decay = self.pg1g2_VAJHU / constants.CJLSTg2decay_mix
        self.M2g1prime2_decay   = self.p0_g1prime2_VAJHU / constants.CJLSTg1prime2decay_pure**2
        self.M2g1g1prime2_decay = self.pg1g1prime2_VAJHU / constants.CJLSTg1prime2decay_mix

        self.notdijet = self.vbf_p0plus_VAJHU == 0 or self.hjj_p0plus_VAJHU == 0
        return self

##########################
#background discriminants#
##########################

    def D_bkg_0plus(self):
        return self.p0plus_VAJHU*self.p0plus_m4l / (self.p0plus_VAJHU*self.p0plus_m4l  + self.bkg_VAMCFM*self.bkg_m4l)
    def D_bkg_0plus_ResUp(self):
        return self.p0plus_VAJHU*self.p0plus_m4l_ResUp / (self.p0plus_VAJHU*self.p0plus_m4l_ResUp  + self.bkg_VAMCFM*self.bkg_m4l_ResUp)
    def D_bkg_0plus_ResDown(self):
        return self.p0plus_VAJHU*self.p0plus_m4l_ResDown / (self.p0plus_VAJHU*self.p0plus_m4l_ResDown  + self.bkg_VAMCFM*self.bkg_m4l_ResDown)
    def D_bkg_0plus_ScaleUp(self):
        return self.p0plus_VAJHU*self.p0plus_m4l_ScaleUp / (self.p0plus_VAJHU*self.p0plus_m4l_ScaleUp  + self.bkg_VAMCFM*self.bkg_m4l_ScaleUp)
    def D_bkg_0plus_ScaleDown(self):
        return self.p0plus_VAJHU*self.p0plus_m4l_ScaleDown / (self.p0plus_VAJHU*self.p0plus_m4l_ScaleDown  + self.bkg_VAMCFM*self.bkg_m4l_ScaleDown)
    def D_bkg_0minus(self):
        return self.p0minus_VAJHU*self.p0plus_m4l / (self.p0minus_VAJHU*self.p0plus_m4l  + self.bkg_VAMCFM*self.bkg_m4l)

    def D_jet_0plus(self):
        if self.notdijet: return -999
        return self.vbf_p0plus_VAJHU / (self.vbf_p0plus_VAJHU + self.hjj_p0plus_VAJHU)

###################################
#anomalous couplings discriminants#
###################################

    def D_0minus_decay(self):
        return self.M2g1_decay / (self.M2g1_decay + self.M2g4_decay*constants.g4decay**2)
    def D_CP_decay(self):
        return self.M2g1g4_decay*constants.g4decay / (self.M2g1_decay + self.M2g4_decay*constants.g4decay**2)
    def D_g2_decay(self):
        return self.M2g1_decay / (self.M2g1_decay + self.M2g2_decay*constants.g2decay**2)
    def D_g1g2_decay(self):
        return self.M2g1g2_decay*constants.g2decay / (self.M2g1_decay + self.M2g2_decay*constants.g2decay**2)
    def D_g1prime2_decay(self):
        return self.M2g1_decay / (self.M2g1_decay + self.M2g1prime2_decay*constants.g1prime2decay_reco**2)
    def D_g1g1prime2_decay(self):
        return self.M2g1g1prime2_decay*constants.g1prime2decay_reco / (self.M2g1_decay + self.M2g1prime2_decay*constants.g1prime2decay_reco**2)

#####################
#Reweighting weights#
#####################

    def MC_weight_ggH(self, index):
        return self.MC_weight * self.reweightingweights[index] * constants.SMXS2L2l / self.nevents2L2l[index]
    def MC_weight_ggH_g1(self):
        return self.MC_weight_ggH(0)
    def MC_weight_ggH_g2(self):
        return self.MC_weight_ggH(1)
    def MC_weight_ggH_g4(self):
        return self.MC_weight_ggH(2)
    def MC_weight_ggH_g1prime2(self):
        return self.MC_weight_ggH(3)
    def MC_weight_ggH_g1g2(self):
        return self.MC_weight_ggH(4)
    def MC_weight_ggH_g1g4(self):
        return self.MC_weight_ggH(5)
    def MC_weight_ggH_g1g1prime2(self):
        return self.MC_weight_ggH(6)

    def MC_weight_ggZZ(self):
        return self.MC_weight * self.xsec * self.tree.KFactorggZZ / self.nevents
    def MC_weight_qqZZ(self):
        return self.MC_weight * self.xsec * self.tree.KFactorEWKqqZZ * self.tree.KFactorQCDqqZZ_M / self.nevents
    def MC_weight_ZX(self):
        self.LepPt, self.LepEta, self.LepLepId = self.tree.LepPt, self.tree.LepEta, self.tree.LepLepId
        return ROOT.fakeRate13TeV(self.LepPt[2],self.LepEta[2],self.LepLepId[2]) * ROOT.fakeRate13TeV(self.LepPt[3],self.LepEta[3],self.LepLepId[3])

    def MC_weight_plain(self):
        return self.MC_weight * self.xsec / self.nevents
    MC_weight_VBF = MC_weight_ZH = MC_weight_WH = MC_weight_ttH = MC_weight_plain

    def getweightfunction(self, sample):
        return getattr(self, sample.weightname())
        raise RuntimeError("{} does not work!".format(sample))

    def getbaseweightfunction(self):
        return self.getweightfunction(self.treesample)

    def initlists(self):
        self.toaddtotree = [
            "D_bkg_0plus",
            "D_bkg_0plus_ResUp",
            "D_bkg_0plus_ResDown",
            "D_bkg_0plus_ScaleUp",
            "D_bkg_0plus_ScaleDown",
            "D_bkg_0minus",
            "D_jet_0plus",
            "D_0minus_decay",
            "D_CP_decay",
            "D_g2_decay",
            "D_g1g2_decay",
            "D_g1prime2_decay",
            "D_g1g1prime2_decay",
        ]

        self.exceptions = [
            "baseweight",
            "checkfunctions",
            "exceptions",
            "flavordict",
            "getbaseweightfunction",
            "getweightfunction",
            "hypothesis",
            "initlists",
            "isbkg",
            "isdata",
            "isdummy",
            "isZX",
            "length",
            "MC_weight_ggH",
            "MC_weight_plain",
            "minevent",
            "nevents",
            "nevents2L2l",
            "next",
            "onlyweights",
            "productionmode",
            "toaddtotree",
            "tree",
            "treesample",
            "weightfunctions",
            "xsec",
        ]

        allsamples = [    #all samples that have weight functions defined in this class
            ReweightingSample("ggH", "0+"),
            ReweightingSample("ggH", "a2"),
            ReweightingSample("ggH", "0-"),
            ReweightingSample("ggH", "L1"),
            ReweightingSample("ggH", "fa20.5"),
            ReweightingSample("ggH", "fa30.5"),
            ReweightingSample("ggH", "fL10.5"),
            ReweightingSample("VBF", "0+"),
            ReweightingSample("ZH", "0+"),
            ReweightingSample("WplusH", "0+"),
            ReweightingSample("WminusH", "0+"),
            ReweightingSample("ttH", "0+"),
            ReweightingSample("ggZZ", "2e2mu"),  #flavor doesn't matter
            ReweightingSample("qqZZ"),
            ReweightingSample("ZX"),
        ]
        reweightingweightnames = [sample.weightname() for sample in self.treesample.reweightingsamples()]
        allreweightingweightnames = [sample.weightname() for sample in allsamples]
        for name in reweightingweightnames:
            if name not in allreweightingweightnames:
                raise ValueError("{} not in allreweightingweightnames!".format(name))
        for sample in allsamples:
            if sample.weightname() in self.toaddtotree or sample.weightname() in self.exceptions: continue
            if sample.weightname() in reweightingweightnames:
                self.toaddtotree.append(sample.weightname())
            else:
                self.exceptions.append(sample.weightname())

        self.flavordict = {
                           13*13*13*13: 0,
                           11*11*11*11: 1,
                           11*11*13*13: 2,
                          }

    def onlyweights(self):
        """Call this to only add the weights and ZZMass to the new tree"""
        #only want the weight, and ZZMass for the range
        reweightingweightnames = [sample.weightname() for sample in self.treesample.reweightingsamples()]
        for name in self.toaddtotree[:]:
            if name not in reweightingweightnames:
                self.exceptions.append(name)
                self.toaddtotree.remove(name)
        self.tree.SetBranchStatus("*", 0)
        self.tree.SetBranchStatus("ZZMass", 1)
        self.tree.SetBranchStatus("Z*Flav", 1)
        for variable in self.treesample.weightingredients():
            self.tree.SetBranchStatus(variable, 1)

    def checkfunctions(self):

        #some cross checking in case of stupid mistakes
        #if a function is added in the class but not added to toaddtotree
        #all member variables, unless they have __, should be added to either toaddtotree or exceptions
        notanywhere, inboth, nonexistent, multipletimes = [], [], [], []
        for key in set(list(type(self).__dict__) + list(self.__dict__) + self.toaddtotree + self.exceptions):
            if key.startswith("__"): continue
            if key.startswith("_abc"): continue
            if key not in self.exceptions and key not in self.toaddtotree and (key in self.__dict__ or key in type(self).__dict__):
                notanywhere.append(key)
            if key in self.toaddtotree and key in self.exceptions:
                inboth.append(key)
            if key not in type(self).__dict__ and key not in self.__dict__:
                nonexistent.append(key)
        for key, occurences in Counter(self.toaddtotree + self.exceptions).iteritems():
            if occurences >= 2 and key not in inboth or occurences >= 3: multipletimes.append(key)
        error = ""
        if notanywhere: error += "the following items are not in toaddtotree or exceptions! " + ", ".join(notanywhere) + "\n"
        if inboth: error += "the following items are in both toaddtotree and exceptions! " + ", ".join(inboth) + "\n"
        if nonexistent: error += "the following items are in toaddtotree or exceptions, but don't exist! " + ", ".join(nonexistent) + "\n"
        if multipletimes: error += "the following items appear multiple times in toaddtotree or exceptions! " + ", ".join(multipletimes) + "\n"
        if error:
            raise SyntaxError(error)


if __name__ == '__main__':
    class DummyTree(object):
        def GetEntries(self):
            return 0
        def GetEntry(self, entry): pass
        xsec = 0
    class DummySample(object):
        productionmode = "graviton fusion"
        hypothesis = "spin 3"
        def isbkg(self): return False
        def isdata(self): return False
        def isZX(self): return False
        def onlyweights(self): return False
        def reweightingsamples(self): return []
        def weightname(self): return "__init__"
    TreeWrapper(DummyTree(), DummySample(), None, None)
    print "You are good to go!"
