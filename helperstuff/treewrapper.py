#!/usr/bin/env python
from array import array
import CJLSTscripts
from collections import Counter, Iterator
import config
import constants
import resource
from samples import ReweightingSample
import sys
import ZX

resource.setrlimit(resource.RLIMIT_STACK, (2**29,-1))
sys.setrecursionlimit(10000)

#to pass to the category code when there are no jets
dummyfloatstar = array('f', [0])

class TreeWrapper(Iterator):

    def __init__(self, tree, treesample, Counters, Counters_reweighted, couplings, minevent=0, maxevent=None, isdummy=False):
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
        self.isPOWHEG = treesample.alternategenerator == "POWHEG"
        self.isdummy = isdummy
        if self.isdata:
            self.unblind = treesample.unblind
        else:
            self.unblind = True

        self.weightfunctions = [self.getweightfunction(sample) for sample in treesample.reweightingsamples()]

        self.nevents = self.nevents2L2l = None
        if Counters is not None:
            self.nevents = Counters.GetBinContent(40)
        if Counters_reweighted is not None:
            self.nevents2L2l = [
                                Counters_reweighted.GetBinContent(4, i) #2e2mu
                              + Counters_reweighted.GetBinContent(8, i) #2e2tau+2mu2tau
                                  for i, sample in enumerate(treesample.reweightingsamples(), start=1)
                               ]

        self.minevent = minevent
        if (
               self.isdata and (self.unblind and not config.unblinddistributions or not config.usedata)
            or self.isZX and not config.usedata
            or self.isdummy
           ):
            self.__length = 0
        elif maxevent is None or maxevent >= self.tree.GetEntries():
            self.__length = self.tree.GetEntries() - minevent
        else:
            self.__length = maxevent - minevent + 1

        if self.isZX:
            ZX.setup(treesample.production)

        self.initlists()
        if treesample.onlyweights(): self.onlyweights()

        tree.GetEntry(0)
        if self.isdata or self.isZX:
            self.xsec = 0
        else:
            self.xsec = tree.xsec * 1000 #pb to fb

        self.cconstantforDbkg = self.cconstantforD2jet = None
        self.checkfunctions(couplings)

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
            if i > len(self):
                raise StopIteration
            if i % 10000 == 0 or i == len(self):
                print i, "/", len(self)
                #raise StopIteration

            if self.isdata:
                self.MC_weight = 1
            elif self.isZX:
                CRflag = t.CRflag
                if CRflag and ZX.test_bit(CRflag, ZX.CRZLLss):
                    self.MC_weight = 1
                else:
                    self.MC_weight = 0
            else:
                self.MC_weight = t.overallEventWeight
            if self.productionmode in ("ggH", "VBF", "ZH", "WH") and not self.isPOWHEG:
                self.reweightingweights = t.reweightingweights
            isSelected = bool(self.MC_weight)

            self.flavor = abs(t.Z1Flav*t.Z2Flav)

            if isSelected:
                break

        #I prefer this to defining __getattr__ because it's faster
        self.ZZMass = t.ZZMass

        #self.cconstantforDbkgkin = CJLSTscripts.getDbkgkinConstant(self.flavor, self.ZZMass)
        self.cconstantforDbkg = CJLSTscripts.getDbkgConstant(self.flavor, self.ZZMass)
        self.cconstantforD2jet = CJLSTscripts.getDVBF2jetsConstant(self.ZZMass)

        self.p_m4l_BKG = t.p_m4l_BKG
        self.p_m4l_SIG = t.p_m4l_SIG

        for a in "SIG", "BKG":
            for b in "Scale", "Res":
                for c in "Up", "Down":
                    attr = "p_m4l_{} {}{}".format(a, b, c)
                    setattr(self, attr, getattr(t, attr))
                    attr = "p_m4l_{}_{}{}".format(a, b, c)
                    setattr(self, attr, getattr(t, attr))

        #express in terms of |M|^2, this will make life easier
        self.M2qqZZ = t.p_QQB_BKG_MCFM

        self.M2g1_decay         = t.p_GG_SIG_ghg2_1_ghz1_1_JHUGen
        self.M2g4_decay         = t.p_GG_SIG_ghg2_1_ghz4_1_JHUGen
        self.M2g1g4_decay       = t.p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen
        self.M2g2_decay         = t.p_GG_SIG_ghg2_1_ghz2_1_JHUGen
        self.M2g1g2_decay       = t.p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen
        self.M2g1prime2_decay   = t.p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen / 1e4
        self.M2g1g1prime2_decay = t.p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen / 1e4

        self.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal = \
        self.M2g1_VBF         = t.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal
        self.M2g4_VBF         = t.p_JJVBF_SIG_ghv4_1_JHUGen_JECNominal
        self.M2g1g4_VBF       = t.p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECNominal
        self.M2g2_VBF         = t.p_JJVBF_SIG_ghv2_1_JHUGen_JECNominal
        self.M2g1g2_VBF       = t.p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECNominal
        self.M2g1prime2_VBF   = t.p_JJVBF_SIG_ghv1prime2_1_JHUGen_JECNominal
        self.M2g1g1prime2_VBF = t.p_JJVBF_SIG_ghv1_1_ghv1prime2_1_JHUGen_JECNominal

        self.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal = \
        self.M2g2_HJJ   = t.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal
        self.M2g4_HJJ   = t.p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal
        self.M2g2g4_HJJ = t.p_JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_JECNominal

        self.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal = \
        self.M2g1_HadZH         = t.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal
        self.M2g4_HadZH         = t.p_HadZH_SIG_ghz4_1_JHUGen_JECNominal
        self.M2g1g4_HadZH       = t.p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECNominal
        self.M2g2_HadZH         = t.p_HadZH_SIG_ghz2_1_JHUGen_JECNominal
        self.M2g1g2_HadZH       = t.p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECNominal
        self.M2g1prime2_HadZH   = t.p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECNominal / 1e4
        self.M2g1g1prime2_HadZH = t.p_HadZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen_JECNominal / 1e4

        self.p_HadWH_SIG_ghz1_1_JHUGen_JECNominal = \
        self.M2g1_HadWH         = t.p_HadWH_SIG_ghz1_1_JHUGen_JECNominal
        self.M2g4_HadWH         = t.p_HadWH_SIG_ghz4_1_JHUGen_JECNominal
        self.M2g1g4_HadWH       = t.p_HadWH_SIG_ghz1_1_ghz4_1_JHUGen_JECNominal
        self.M2g2_HadWH         = t.p_HadWH_SIG_ghz2_1_JHUGen_JECNominal
        self.M2g1g2_HadWH       = t.p_HadWH_SIG_ghz1_1_ghz2_1_JHUGen_JECNominal
        self.M2g1prime2_HadWH   = t.p_HadWH_SIG_ghz1prime2_1E4_JHUGen_JECNominal / 1e4
        self.M2g1g1prime2_HadWH = t.p_HadWH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen_JECNominal / 1e4

        #category variables
        self.nExtraLep = t.nExtraLep
        self.nExtraZ = t.nExtraZ
        self.nCleanedJetsPt30 = t.nCleanedJetsPt30
        self.nCleanedJetsPt30BTagged = t.nCleanedJetsPt30BTagged
        self.jetQGLikelihood = t.JetQGLikelihood.data()
        self.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal = t.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal
        self.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal = t.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal
        self.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal = t.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal
        self.jetPhi = t.JetPhi.data()

        if self.nCleanedJetsPt30 == 0:
            self.jetQGLikelihood = self.jetPhi = dummyfloatstar

        if self.isdata and not self.unblind and not self.passesblindcut():
            return next(self)

        self.notdijet = self.M2g1_VBF < 0 or self.M2g2_HJJ < 0
        return self

    def __len__(self):
        return self.__length

##########################
#background discriminants#
##########################

    def D_bkg_0plus(self):
        try:
            return self.p0plus_VAJHU*self.p_m4l_SIG / (self.p0plus_VAJHU*self.p_m4l_SIG  + self.M2qqZZ*self.p_m4l_BKG*self.cconstantforDbkg)
        except ZeroDivisionError:
            return 0
    def D_bkg_0plus_ResUp(self):
        try:
            return self.p0plus_VAJHU*self.p_m4l_SIG_ResUp / (self.p0plus_VAJHU*self.p_m4l_SIG_ResUp  + self.M2qqZZ*self.p_m4l_BKG_ResUp*self.cconstantforDbkg)
        except ZeroDivisionError:
            return 0
    def D_bkg_0plus_ResDown(self):
        try:
            return self.p0plus_VAJHU*self.p_m4l_SIG_ResDown / (self.p0plus_VAJHU*self.p_m4l_SIG_ResDown  + self.M2qqZZ*self.p_m4l_BKG_ResDown*self.cconstantforDbkg)
        except ZeroDivisionError:
            return 0
    def D_bkg_0plus_ScaleUp(self):
        try:
            return self.p0plus_VAJHU*self.p_m4l_SIG_ScaleUp / (self.p0plus_VAJHU*self.p_m4l_SIG_ScaleUp  + self.M2qqZZ*self.p_m4l_BKG_ScaleUp*self.cconstantforDbkg)
        except ZeroDivisionError:
            return 0
    def D_bkg_0plus_ScaleDown(self):
        try:
            return self.p0plus_VAJHU*self.p_m4l_SIG_ScaleDown / (self.p0plus_VAJHU*self.p_m4l_SIG_ScaleDown  + self.M2qqZZ*self.p_m4l_BKG_ScaleDown*self.cconstantforDbkg)
        except ZeroDivisionError:
            return 0

    def D_2jet_0plus(self):
        if self.notdijet: return -999
        return self.M2g1_VBF / (self.M2g1_VBF + self.M2g2_HJJ*self.cconstantforD2jet)
    def D_2jet_0minus(self):
        if self.notdijet: return -999
        return self.M2g4_VBF / (self.M2g4_VBF + self.M2g2_HJJ*self.cconstantforD2jet)

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

#######################################
#VBF anomalous couplings discriminants#
#######################################

    def D_0minus_VBF(self):
        if self.notdijet: return -999
        return self.M2g1_VBF / (self.M2g1_VBF + self.M2g4_VBF*constants.g4VBF**2)
    def D_CP_VBF(self):
        if self.notdijet: return -999
        return self.M2g1g4_VBF*constants.g4VBF / (self.M2g1_VBF + self.M2g4_VBF*constants.g4VBF**2)
    def D_g2_VBF(self):
        if self.notdijet: return -999
        return self.M2g1_VBF / (self.M2g1_VBF + self.M2g2_VBF*constants.g2VBF**2)
    def D_g1g2_VBF(self):
        if self.notdijet: return -999
        return self.M2g1g2_VBF*constants.g2VBF / (self.M2g1_VBF + self.M2g2_VBF*constants.g2VBF**2)
    def D_g1prime2_VBF(self):
        if self.notdijet: return -999
        return self.M2g1_VBF / (self.M2g1_VBF + self.M2g1prime2_VBF*constants.g1prime2VBF_reco**2)
    def D_g1g1prime2_VBF(self):
        if self.notdijet: return -999
        return self.M2g1g1prime2_VBF*constants.g1prime2VBF_reco / (self.M2g1_VBF + self.M2g1prime2_VBF*constants.g1prime2VBF_reco**2)

###############################################
#ZH hadronic anomalous couplings discriminants#
###############################################

    def D_0minus_HadZH(self):
        if self.notdijet: return -999
        return self.M2g1_HadZH / (self.M2g1_HadZH + self.M2g4_HadZH*constants.g4ZH**2)
    def D_CP_HadZH(self):
        if self.notdijet: return -999
        return self.M2g1g4_HadZH*constants.g4ZH / (self.M2g1_HadZH + self.M2g4_HadZH*constants.g4ZH**2)
    def D_g2_HadZH(self):
        if self.notdijet: return -999
        return self.M2g1_HadZH / (self.M2g1_HadZH + self.M2g2_HadZH*constants.g2ZH**2)
    def D_g1g2_HadZH(self):
        if self.notdijet: return -999
        return self.M2g1g2_HadZH*constants.g2ZH / (self.M2g1_HadZH + self.M2g2_HadZH*constants.g2ZH**2)
    def D_g1prime2_HadZH(self):
        if self.notdijet: return -999
        return self.M2g1_HadZH / (self.M2g1_HadZH + self.M2g1prime2_HadZH*constants.g1prime2ZH_reco**2)
    def D_g1g1prime2_HadZH(self):
        if self.notdijet: return -999
        return self.M2g1g1prime2_HadZH*constants.g1prime2ZH_reco / (self.M2g1_HadZH + self.M2g1prime2_HadZH*constants.g1prime2ZH_reco**2)

###############################################
#WH hadronic anomalous couplings discriminants#
###############################################

    def D_0minus_HadWH(self):
        if self.notdijet: return -999
        return self.M2g1_HadWH / (self.M2g1_HadWH + self.M2g4_HadWH*constants.g4WH**2)
    def D_CP_HadWH(self):
        if self.notdijet: return -999
        return self.M2g1g4_HadWH*constants.g4WH / (self.M2g1_HadWH + self.M2g4_HadWH*constants.g4WH**2)
    def D_g2_HadWH(self):
        if self.notdijet: return -999
        return self.M2g1_HadWH / (self.M2g1_HadWH + self.M2g2_HadWH*constants.g2WH**2)
    def D_g1g2_HadWH(self):
        if self.notdijet: return -999
        return self.M2g1g2_HadWH*constants.g2WH / (self.M2g1_HadWH + self.M2g2_HadWH*constants.g2WH**2)
    def D_g1prime2_HadWH(self):
        if self.notdijet: return -999
        return self.M2g1_HadWH / (self.M2g1_HadWH + self.M2g1prime2_HadWH*constants.g1prime2WH_reco**2)
    def D_g1g1prime2_HadWH(self):
        if self.notdijet: return -999
        return self.M2g1g1prime2_HadWH*constants.g1prime2WH_reco / (self.M2g1_HadWH + self.M2g1prime2_HadWH*constants.g1prime2WH_reco**2)

###############################################
#VH hadronic anomalous couplings discriminants#
###############################################

    def D_0minus_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1_HadWH*constants.JHUXSWHa1 + self.M2g1_HadZH*constants.JHUXSZHa1)
               /
                 (
                   (self.M2g1_HadWH + self.M2g4_HadWH*constants.g4WH**2)*constants.JHUXSWHa1
                 +
                   (self.M2g1_HadZH + self.M2g4_HadZH*constants.g4ZH**2)*constants.JHUXSZHa1
                 )
               )
    def D_CP_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1g4_HadWH*constants.g4WH*constants.JHUXSWHa1 + self.M2g1g4_HadZH*constants.g4ZH*constants.JHUXSZHa1)
               /
                 (
                   (self.M2g1_HadWH + self.M2g4_HadWH*constants.g4WH**2)*constants.JHUXSWHa1
                 +
                   (self.M2g1_HadZH + self.M2g4_HadZH*constants.g4ZH**2)*constants.JHUXSZHa1
                 )
               )
    def D_g2_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1_HadWH*constants.JHUXSWHa1 + self.M2g1_HadZH*constants.JHUXSZHa1)
               /
                 (
                   (self.M2g1_HadWH + self.M2g2_HadWH*constants.g2WH**2)*constants.JHUXSWHa1
                 +
                   (self.M2g1_HadZH + self.M2g2_HadZH*constants.g2ZH**2)*constants.JHUXSZHa1
                 )
               )
    def D_g1g2_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1g2_HadWH*constants.g2WH*constants.JHUXSWHa1 + self.M2g1g2_HadZH*constants.g2ZH*constants.JHUXSZHa1)
               /
                 (
                   (self.M2g1_HadWH + self.M2g2_HadWH*constants.g2WH**2)*constants.JHUXSWHa1
                 +
                   (self.M2g1_HadZH + self.M2g2_HadZH*constants.g2ZH**2)*constants.JHUXSZHa1
                 )
               )
    def D_g1prime2_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1_HadWH*constants.JHUXSWHa1 + self.M2g1_HadZH*constants.JHUXSZHa1)
               /
                 (
                   (self.M2g1_HadWH + self.M2g1prime2_HadWH*constants.g1prime2WH_reco**2)*constants.JHUXSWHa1
                 +
                   (self.M2g1_HadZH + self.M2g1prime2_HadZH*constants.g1prime2ZH_reco**2)*constants.JHUXSZHa1
                 )
               )
    def D_g1g1prime2_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1g1prime2_HadWH*constants.g1prime2WH_reco*constants.JHUXSWHa1 + self.M2g1g1prime2_HadZH*constants.g1prime2ZH_reco*constants.JHUXSZHa1)
               /
                 (
                   (self.M2g1_HadWH + self.M2g1prime2_HadWH*constants.g1prime2WH_reco**2)*constants.JHUXSWHa1
                 +
                   (self.M2g1_HadZH + self.M2g1prime2_HadZH*constants.g1prime2ZH_reco**2)*constants.JHUXSZHa1
                 )
               )

############################################
#VBFdecay anomalous couplings discriminants#
############################################

    def D_0minus_VBFdecay(self):
        if self.notdijet: return -999
        return self.M2g1_VBF*self.M2g1_decay / (self.M2g1_VBF*self.M2g1_decay + self.M2g4_VBF*self.M2g4_decay*(constants.g4VBF*constants.g4decay)**2)
    def D_g2_VBFdecay(self):
        if self.notdijet: return -999
        return self.M2g1_VBF*self.M2g1_decay / (self.M2g1_VBF*self.M2g1_decay + self.M2g2_VBF*self.M2g2_decay * (constants.g2VBF*constants.g2decay)**2)
    def D_g1prime2_VBFdecay(self):
        if self.notdijet: return -999
        return self.M2g1_VBF*self.M2g1_decay / (self.M2g1_VBF*self.M2g1_decay + self.M2g1prime2_VBF*self.M2g1prime2_decay * (constants.g1prime2VBF_reco*constants.g1prime2decay_reco)**2)

####################################################
#ZHdecay hadronic anomalous couplings discriminants#
####################################################

    def D_0minus_HadZHdecay(self):
        if self.notdijet: return -999
        return self.M2g1_HadZH*self.M2g1_decay / (self.M2g1_HadZH*self.M2g1_decay + self.M2g4_HadZH*self.M2g4_decay*(constants.g4ZH*constants.g4decay)**2)
    def D_g2_HadZHdecay(self):
        if self.notdijet: return -999
        return self.M2g1_HadZH*self.M2g1_decay / (self.M2g1_HadZH*self.M2g1_decay + self.M2g2_HadZH*self.M2g2_decay * (constants.g2ZH*constants.g2decay)**2)
    def D_g1prime2_HadZHdecay(self):
        if self.notdijet: return -999
        return self.M2g1_HadZH*self.M2g1_decay / (self.M2g1_HadZH*self.M2g1_decay + self.M2g1prime2_HadZH*self.M2g1prime2_decay * (constants.g1prime2ZH_reco*constants.g1prime2decay_reco)**2)

####################################################
#WHdecay hadronic anomalous couplings discriminants#
####################################################

    def D_0minus_HadWHdecay(self):
        if self.notdijet: return -999
        return self.M2g1_HadWH*self.M2g1_decay / (self.M2g1_HadWH*self.M2g1_decay + self.M2g4_HadWH*self.M2g4_decay*(constants.g4WH*constants.g4decay)**2)
    def D_g2_HadWHdecay(self):
        if self.notdijet: return -999
        return self.M2g1_HadWH*self.M2g1_decay / (self.M2g1_HadWH*self.M2g1_decay + self.M2g2_HadWH*self.M2g2_decay * (constants.g2WH*constants.g2decay)**2)
    def D_g1prime2_HadWHdecay(self):
        if self.notdijet: return -999
        return self.M2g1_HadWH*self.M2g1_decay / (self.M2g1_HadWH*self.M2g1_decay + self.M2g1prime2_HadWH*self.M2g1prime2_decay * (constants.g1prime2WH_reco*constants.g1prime2decay_reco)**2)

####################################################
#VHdecay hadronic anomalous couplings discriminants#
####################################################

    def D_0minus_HadVHdecay(self):
        if self.notdijet: return -999
        return (
                 ((self.M2g1_HadWH*constants.JHUXSWHa1+self.M2g1_HadZH*constants.JHUXSZHa1)
                    *self.M2g1_decay)
               /
                 (
                   (self.M2g1_HadWH*constants.JHUXSWHa1+self.M2g1_HadZH*constants.JHUXSZHa1)
                        *self.M2g1_decay
                 + (self.M2g4_HadWH*constants.g4WH**2*constants.JHUXSWHa1+self.M2g4_HadZH*constants.g4ZH**2*constants.JHUXSZHa1)
                        *self.M2g4_decay*constants.g4decay**2
                 )
               )
    def D_g2_HadVHdecay(self):
        if self.notdijet: return -999
        return (
                 ((self.M2g1_HadWH*constants.JHUXSWHa1+self.M2g1_HadZH*constants.JHUXSZHa1)
                    *self.M2g1_decay)
               /
                 (
                   (self.M2g1_HadWH*constants.JHUXSWHa1+self.M2g1_HadZH*constants.JHUXSZHa1)
                        *self.M2g1_decay
                 + (self.M2g2_HadWH*constants.g2WH**2*constants.JHUXSWHa1+self.M2g2_HadZH*constants.g2ZH**2*constants.JHUXSZHa1)
                        *self.M2g2_decay*constants.g2decay**2
                 )
               )
    def D_g1prime2_HadVHdecay(self):
        if self.notdijet: return -999
        return (
                 ((self.M2g1_HadWH*constants.JHUXSWHa1+self.M2g1_HadZH*constants.JHUXSZHa1)
                    *self.M2g1_decay)
               /
                 (
                   (self.M2g1_HadWH*constants.JHUXSWHa1+self.M2g1_HadZH*constants.JHUXSZHa1)
                        *self.M2g1_decay
                 + (self.M2g1prime2_HadWH*constants.g1prime2WH_reco**2*constants.JHUXSWHa1+self.M2g1prime2_HadZH*constants.g1prime2ZH_reco**2*constants.JHUXSZHa1)
                        *self.M2g1prime2_decay*constants.g1prime2decay_reco**2
                 )
               )

######################################
#VBFdecay g1^x gi^(4-x) discriminants#
######################################

    ####
    #g4#
    ####

    def D_g14_g40_VBFdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g1_VBF*self.M2g1_decay
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g4VBF**2*self.M2g4_VBF * constants.g4decay**2*self.M2g4_decay)
               )
    def D_g13_g41_VBFdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_VBF * self.M2g1g4_decay*constants.g4decay + self.M2g1g4_VBF*constants.g4VBF * self.M2g1_decay)
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g4VBF**2*self.M2g4_VBF * constants.g4decay**2*self.M2g4_decay)
               )
    def D_g12_g42_VBFdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_VBF * self.M2g4_decay*constants.g4decay**2 + self.M2g4_VBF*constants.g4VBF**2 * self.M2g1_decay
                   + self.M2g1g4_VBF*constants.g4VBF * self.M2g1g4_decay*constants.g4decay)
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g4VBF**2*self.M2g4_VBF * constants.g4decay**2*self.M2g4_decay)
               )
    def D_g11_g43_VBFdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g4_VBF*constants.g4VBF**2 * self.M2g1g4_decay*constants.g4decay + self.M2g1g4_VBF*constants.g4VBF * self.M2g4_decay*constants.g4decay**2)
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g4VBF**2*self.M2g4_VBF * constants.g4decay**2*self.M2g4_decay)
               )
    def D_g10_g44_VBFdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g4_VBF*constants.g4VBF**2 * self.M2g4_decay*constants.g4decay**2
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g4VBF**2*self.M2g4_VBF * constants.g4decay**2*self.M2g4_decay)
               )

    ####
    #g2#
    ####

    def D_g14_g20_VBFdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g1_VBF*self.M2g1_decay
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g2VBF**2*self.M2g2_VBF * constants.g2decay**2*self.M2g2_decay)
               )
    def D_g13_g21_VBFdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_VBF * self.M2g1g2_decay*constants.g2decay + self.M2g1g2_VBF*constants.g2VBF * self.M2g1_decay)
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g2VBF**2*self.M2g2_VBF * constants.g2decay**2*self.M2g2_decay)
               )
    def D_g12_g22_VBFdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_VBF * self.M2g2_decay*constants.g2decay**2 + self.M2g2_VBF*constants.g2VBF**2 * self.M2g1_decay
                   + self.M2g1g2_VBF*constants.g2VBF * self.M2g1g2_decay*constants.g2decay)
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g2VBF**2*self.M2g2_VBF * constants.g2decay**2*self.M2g2_decay)
               )
    def D_g11_g23_VBFdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g2_VBF*constants.g2VBF**2 * self.M2g1g2_decay*constants.g2decay + self.M2g1g2_VBF*constants.g2VBF * self.M2g2_decay*constants.g2decay**2)
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g2VBF**2*self.M2g2_VBF * constants.g2decay**2*self.M2g2_decay)
               )
    def D_g10_g24_VBFdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g2_VBF*constants.g2VBF**2 * self.M2g2_decay*constants.g2decay**2
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g2VBF**2*self.M2g2_VBF * constants.g2decay**2*self.M2g2_decay)
               )

    ##########
    #g1prime2#
    ##########

    def D_g14_g1prime20_VBFdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g1_VBF*self.M2g1_decay
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g1prime2VBF_reco**2*self.M2g1prime2_VBF * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )
    def D_g13_g1prime21_VBFdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_VBF * self.M2g1g1prime2_decay*constants.g1prime2decay_reco + self.M2g1g1prime2_VBF*constants.g1prime2VBF_reco * self.M2g1_decay)
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g1prime2VBF_reco**2*self.M2g1prime2_VBF * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )
    def D_g12_g1prime22_VBFdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_VBF * self.M2g1prime2_decay*constants.g1prime2decay_reco**2 + self.M2g1prime2_VBF*constants.g1prime2VBF_reco**2 * self.M2g1_decay
                   + self.M2g1g1prime2_VBF*constants.g1prime2VBF_reco * self.M2g1g1prime2_decay*constants.g1prime2decay_reco)
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g1prime2VBF_reco**2*self.M2g1prime2_VBF * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )
    def D_g11_g1prime23_VBFdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1prime2_VBF*constants.g1prime2VBF_reco**2 * self.M2g1g1prime2_decay*constants.g1prime2decay_reco + self.M2g1g1prime2_VBF*constants.g1prime2VBF_reco * self.M2g1prime2_decay*constants.g1prime2decay_reco**2)
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g1prime2VBF_reco**2*self.M2g1prime2_VBF * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )
    def D_g10_g1prime24_VBFdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g1prime2_VBF*constants.g1prime2VBF_reco**2 * self.M2g1prime2_decay*constants.g1prime2decay_reco**2
                 /
                (self.M2g1_VBF * self.M2g1_decay + constants.g1prime2VBF_reco**2*self.M2g1prime2_VBF * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )

############################################
#VBFdecay g1^x gi^(4-x) discriminants prime#
############################################

    ####
    #g4#
    ####

    def D_g14_g40_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g1_VBF*self.M2g1_decay
                 /
                ((self.M2g1_VBF + constants.g4VBF**2*self.M2g4_VBF) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )
    def D_g13_g41_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_VBF * self.M2g1g4_decay*constants.g4decay + self.M2g1g4_VBF*constants.g4VBF * self.M2g1_decay)
                 /
                ((self.M2g1_VBF + constants.g4VBF**2*self.M2g4_VBF) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )
    def D_g12_g42_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_VBF * self.M2g4_decay*constants.g4decay**2 + self.M2g4_VBF*constants.g4VBF**2 * self.M2g1_decay
                   + self.M2g1g4_VBF*constants.g4VBF * self.M2g1g4_decay*constants.g4decay)
                 /
                ((self.M2g1_VBF + constants.g4VBF**2*self.M2g4_VBF) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )
    def D_g11_g43_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g4_VBF*constants.g4VBF**2 * self.M2g1g4_decay*constants.g4decay + self.M2g1g4_VBF*constants.g4VBF * self.M2g4_decay*constants.g4decay**2)
                 /
                ((self.M2g1_VBF + constants.g4VBF**2*self.M2g4_VBF) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )
    def D_g10_g44_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g4_VBF*constants.g4VBF**2 * self.M2g4_decay*constants.g4decay**2
                 /
                ((self.M2g1_VBF + constants.g4VBF**2*self.M2g4_VBF) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )

    ####
    #g2#
    ####

    def D_g14_g20_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g1_VBF*self.M2g1_decay
                 /
                ((self.M2g1_VBF + constants.g2VBF**2*self.M2g2_VBF) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )
    def D_g13_g21_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_VBF * self.M2g1g2_decay*constants.g2decay + self.M2g1g2_VBF*constants.g2VBF * self.M2g1_decay)
                 /
                ((self.M2g1_VBF + constants.g2VBF**2*self.M2g2_VBF) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )
    def D_g12_g22_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_VBF * self.M2g2_decay*constants.g2decay**2 + self.M2g2_VBF*constants.g2VBF**2 * self.M2g1_decay
                   + self.M2g1g2_VBF*constants.g2VBF * self.M2g1g2_decay*constants.g2decay)
                 /
                ((self.M2g1_VBF + constants.g2VBF**2*self.M2g2_VBF) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )
    def D_g11_g23_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g2_VBF*constants.g2VBF**2 * self.M2g1g2_decay*constants.g2decay + self.M2g1g2_VBF*constants.g2VBF * self.M2g2_decay*constants.g2decay**2)
                 /
                ((self.M2g1_VBF + constants.g2VBF**2*self.M2g2_VBF) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )
    def D_g10_g24_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g2_VBF*constants.g2VBF**2 * self.M2g2_decay*constants.g2decay**2
                 /
                ((self.M2g1_VBF + constants.g2VBF**2*self.M2g2_VBF) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )

    ##########
    #g1prime2#
    ##########

    def D_g14_g1prime20_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g1_VBF*self.M2g1_decay
                 /
                ((self.M2g1_VBF + constants.g1prime2VBF_reco**2*self.M2g1prime2_VBF) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )
    def D_g13_g1prime21_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_VBF * self.M2g1g1prime2_decay*constants.g1prime2decay_reco + self.M2g1g1prime2_VBF*constants.g1prime2VBF_reco * self.M2g1_decay)
                 /
                ((self.M2g1_VBF + constants.g1prime2VBF_reco**2*self.M2g1prime2_VBF) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )
    def D_g12_g1prime22_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_VBF * self.M2g1prime2_decay*constants.g1prime2decay_reco**2 + self.M2g1prime2_VBF*constants.g1prime2VBF_reco**2 * self.M2g1_decay
                   + self.M2g1g1prime2_VBF*constants.g1prime2VBF_reco * self.M2g1g1prime2_decay*constants.g1prime2decay_reco)
                 /
                ((self.M2g1_VBF + constants.g1prime2VBF_reco**2*self.M2g1prime2_VBF) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )
    def D_g11_g1prime23_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1prime2_VBF*constants.g1prime2VBF_reco**2 * self.M2g1g1prime2_decay*constants.g1prime2decay_reco + self.M2g1g1prime2_VBF*constants.g1prime2VBF_reco * self.M2g1prime2_decay*constants.g1prime2decay_reco**2)
                 /
                ((self.M2g1_VBF + constants.g1prime2VBF_reco**2*self.M2g1prime2_VBF) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )
    def D_g10_g1prime24_VBFdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g1prime2_VBF*constants.g1prime2VBF_reco**2 * self.M2g1prime2_decay*constants.g1prime2decay_reco**2
                 /
                ((self.M2g1_VBF + constants.g1prime2VBF_reco**2*self.M2g1prime2_VBF) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )

##############################################
#ZHdecay hadronic g1^x gi^(4-x) discriminants#
##############################################

    ####
    #g4#
    ####

    def D_g14_g40_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g1_HadZH*self.M2g1_decay
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g4ZH**2*self.M2g4_HadZH * constants.g4decay**2*self.M2g4_decay)
               )
    def D_g13_g41_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadZH * self.M2g1g4_decay*constants.g4decay + self.M2g1g4_HadZH*constants.g4ZH * self.M2g1_decay)
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g4ZH**2*self.M2g4_HadZH * constants.g4decay**2*self.M2g4_decay)
               )
    def D_g12_g42_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadZH * self.M2g4_decay*constants.g4decay**2 + self.M2g4_HadZH*constants.g4ZH**2 * self.M2g1_decay
                   + self.M2g1g4_HadZH*constants.g4ZH * self.M2g1g4_decay*constants.g4decay)
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g4ZH**2*self.M2g4_HadZH * constants.g4decay**2*self.M2g4_decay)
               )
    def D_g11_g43_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g4_HadZH*constants.g4ZH**2 * self.M2g1g4_decay*constants.g4decay + self.M2g1g4_HadZH*constants.g4ZH * self.M2g4_decay*constants.g4decay**2)
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g4ZH**2*self.M2g4_HadZH * constants.g4decay**2*self.M2g4_decay)
               )
    def D_g10_g44_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g4_HadZH*constants.g4ZH**2 * self.M2g4_decay*constants.g4decay**2
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g4ZH**2*self.M2g4_HadZH * constants.g4decay**2*self.M2g4_decay)
               )

    ####
    #g2#
    ####

    def D_g14_g20_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g1_HadZH*self.M2g1_decay
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g2ZH**2*self.M2g2_HadZH * constants.g2decay**2*self.M2g2_decay)
               )
    def D_g13_g21_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadZH * self.M2g1g2_decay*constants.g2decay + self.M2g1g2_HadZH*constants.g2ZH * self.M2g1_decay)
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g2ZH**2*self.M2g2_HadZH * constants.g2decay**2*self.M2g2_decay)
               )
    def D_g12_g22_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadZH * self.M2g2_decay*constants.g2decay**2 + self.M2g2_HadZH*constants.g2ZH**2 * self.M2g1_decay
                   + self.M2g1g2_HadZH*constants.g2ZH * self.M2g1g2_decay*constants.g2decay)
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g2ZH**2*self.M2g2_HadZH * constants.g2decay**2*self.M2g2_decay)
               )
    def D_g11_g23_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g2_HadZH*constants.g2ZH**2 * self.M2g1g2_decay*constants.g2decay + self.M2g1g2_HadZH*constants.g2ZH * self.M2g2_decay*constants.g2decay**2)
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g2ZH**2*self.M2g2_HadZH * constants.g2decay**2*self.M2g2_decay)
               )
    def D_g10_g24_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g2_HadZH*constants.g2ZH**2 * self.M2g2_decay*constants.g2decay**2
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g2ZH**2*self.M2g2_HadZH * constants.g2decay**2*self.M2g2_decay)
               )

    ##########
    #g1prime2#
    ##########

    def D_g14_g1prime20_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g1_HadZH*self.M2g1_decay
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g1prime2ZH_reco**2*self.M2g1prime2_HadZH * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )
    def D_g13_g1prime21_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadZH * self.M2g1g1prime2_decay*constants.g1prime2decay_reco + self.M2g1g1prime2_HadZH*constants.g1prime2ZH_reco * self.M2g1_decay)
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g1prime2ZH_reco**2*self.M2g1prime2_HadZH * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )
    def D_g12_g1prime22_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadZH * self.M2g1prime2_decay*constants.g1prime2decay_reco**2 + self.M2g1prime2_HadZH*constants.g1prime2ZH_reco**2 * self.M2g1_decay
                   + self.M2g1g1prime2_HadZH*constants.g1prime2ZH_reco * self.M2g1g1prime2_decay*constants.g1prime2decay_reco)
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g1prime2ZH_reco**2*self.M2g1prime2_HadZH * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )
    def D_g11_g1prime23_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1prime2_HadZH*constants.g1prime2ZH_reco**2 * self.M2g1g1prime2_decay*constants.g1prime2decay_reco + self.M2g1g1prime2_HadZH*constants.g1prime2ZH_reco * self.M2g1prime2_decay*constants.g1prime2decay_reco**2)
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g1prime2ZH_reco**2*self.M2g1prime2_HadZH * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )
    def D_g10_g1prime24_HadZHdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g1prime2_HadZH*constants.g1prime2ZH_reco**2 * self.M2g1prime2_decay*constants.g1prime2decay_reco**2
                 /
                (self.M2g1_HadZH * self.M2g1_decay + constants.g1prime2ZH_reco**2*self.M2g1prime2_HadZH * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )

####################################################
#ZHdecay hadronic g1^x gi^(4-x) discriminants prime#
####################################################

    ####
    #g4#
    ####

    def D_g14_g40_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g1_HadZH*self.M2g1_decay
                 /
                ((self.M2g1_HadZH + constants.g4ZH**2*self.M2g4_HadZH) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )
    def D_g13_g41_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadZH * self.M2g1g4_decay*constants.g4decay + self.M2g1g4_HadZH*constants.g4ZH * self.M2g1_decay)
                 /
                ((self.M2g1_HadZH + constants.g4ZH**2*self.M2g4_HadZH) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )
    def D_g12_g42_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadZH * self.M2g4_decay*constants.g4decay**2 + self.M2g4_HadZH*constants.g4ZH**2 * self.M2g1_decay
                   + self.M2g1g4_HadZH*constants.g4ZH * self.M2g1g4_decay*constants.g4decay)
                 /
                ((self.M2g1_HadZH + constants.g4ZH**2*self.M2g4_HadZH) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )
    def D_g11_g43_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g4_HadZH*constants.g4ZH**2 * self.M2g1g4_decay*constants.g4decay + self.M2g1g4_HadZH*constants.g4ZH * self.M2g4_decay*constants.g4decay**2)
                 /
                ((self.M2g1_HadZH + constants.g4ZH**2*self.M2g4_HadZH) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )
    def D_g10_g44_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g4_HadZH*constants.g4ZH**2 * self.M2g4_decay*constants.g4decay**2
                 /
                ((self.M2g1_HadZH + constants.g4ZH**2*self.M2g4_HadZH) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )

    ####
    #g2#
    ####

    def D_g14_g20_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g1_HadZH*self.M2g1_decay
                 /
                ((self.M2g1_HadZH + constants.g2ZH**2*self.M2g2_HadZH) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )
    def D_g13_g21_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadZH * self.M2g1g2_decay*constants.g2decay + self.M2g1g2_HadZH*constants.g2ZH * self.M2g1_decay)
                 /
                ((self.M2g1_HadZH + constants.g2ZH**2*self.M2g2_HadZH) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )
    def D_g12_g22_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadZH * self.M2g2_decay*constants.g2decay**2 + self.M2g2_HadZH*constants.g2ZH**2 * self.M2g1_decay
                   + self.M2g1g2_HadZH*constants.g2ZH * self.M2g1g2_decay*constants.g2decay)
                 /
                ((self.M2g1_HadZH + constants.g2ZH**2*self.M2g2_HadZH) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )
    def D_g11_g23_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g2_HadZH*constants.g2ZH**2 * self.M2g1g2_decay*constants.g2decay + self.M2g1g2_HadZH*constants.g2ZH * self.M2g2_decay*constants.g2decay**2)
                 /
                ((self.M2g1_HadZH + constants.g2ZH**2*self.M2g2_HadZH) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )
    def D_g10_g24_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g2_HadZH*constants.g2ZH**2 * self.M2g2_decay*constants.g2decay**2
                 /
                ((self.M2g1_HadZH + constants.g2ZH**2*self.M2g2_HadZH) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )

    ##########
    #g1prime2#
    ##########

    def D_g14_g1prime20_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g1_HadZH*self.M2g1_decay
                 /
                ((self.M2g1_HadZH + constants.g1prime2ZH_reco**2*self.M2g1prime2_HadZH) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )
    def D_g13_g1prime21_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadZH * self.M2g1g1prime2_decay*constants.g1prime2decay_reco + self.M2g1g1prime2_HadZH*constants.g1prime2ZH_reco * self.M2g1_decay)
                 /
                ((self.M2g1_HadZH + constants.g1prime2ZH_reco**2*self.M2g1prime2_HadZH) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )
    def D_g12_g1prime22_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadZH * self.M2g1prime2_decay*constants.g1prime2decay_reco**2 + self.M2g1prime2_HadZH*constants.g1prime2ZH_reco**2 * self.M2g1_decay
                   + self.M2g1g1prime2_HadZH*constants.g1prime2ZH_reco * self.M2g1g1prime2_decay*constants.g1prime2decay_reco)
                 /
                ((self.M2g1_HadZH + constants.g1prime2ZH_reco**2*self.M2g1prime2_HadZH) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )
    def D_g11_g1prime23_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1prime2_HadZH*constants.g1prime2ZH_reco**2 * self.M2g1g1prime2_decay*constants.g1prime2decay_reco + self.M2g1g1prime2_HadZH*constants.g1prime2ZH_reco * self.M2g1prime2_decay*constants.g1prime2decay_reco**2)
                 /
                ((self.M2g1_HadZH + constants.g1prime2ZH_reco**2*self.M2g1prime2_HadZH) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )
    def D_g10_g1prime24_HadZHdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g1prime2_HadZH*constants.g1prime2ZH_reco**2 * self.M2g1prime2_decay*constants.g1prime2decay_reco**2
                 /
                ((self.M2g1_HadZH + constants.g1prime2ZH_reco**2*self.M2g1prime2_HadZH) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )

##############################################
#WHdecay hadronic g1^x gi^(4-x) discriminants#
##############################################

    ####
    #g4#
    ####

    def D_g14_g40_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g1_HadWH*self.M2g1_decay
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g4WH**2*self.M2g4_HadWH * constants.g4decay**2*self.M2g4_decay)
               )
    def D_g13_g41_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadWH * self.M2g1g4_decay*constants.g4decay + self.M2g1g4_HadWH*constants.g4WH * self.M2g1_decay)
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g4WH**2*self.M2g4_HadWH * constants.g4decay**2*self.M2g4_decay)
               )
    def D_g12_g42_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadWH * self.M2g4_decay*constants.g4decay**2 + self.M2g4_HadWH*constants.g4WH**2 * self.M2g1_decay
                   + self.M2g1g4_HadWH*constants.g4WH * self.M2g1g4_decay*constants.g4decay)
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g4WH**2*self.M2g4_HadWH * constants.g4decay**2*self.M2g4_decay)
               )
    def D_g11_g43_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g4_HadWH*constants.g4WH**2 * self.M2g1g4_decay*constants.g4decay + self.M2g1g4_HadWH*constants.g4WH * self.M2g4_decay*constants.g4decay**2)
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g4WH**2*self.M2g4_HadWH * constants.g4decay**2*self.M2g4_decay)
               )
    def D_g10_g44_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g4_HadWH*constants.g4WH**2 * self.M2g4_decay*constants.g4decay**2
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g4WH**2*self.M2g4_HadWH * constants.g4decay**2*self.M2g4_decay)
               )

    ####
    #g2#
    ####

    def D_g14_g20_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g1_HadWH*self.M2g1_decay
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g2WH**2*self.M2g2_HadWH * constants.g2decay**2*self.M2g2_decay)
               )
    def D_g13_g21_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadWH * self.M2g1g2_decay*constants.g2decay + self.M2g1g2_HadWH*constants.g2WH * self.M2g1_decay)
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g2WH**2*self.M2g2_HadWH * constants.g2decay**2*self.M2g2_decay)
               )
    def D_g12_g22_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadWH * self.M2g2_decay*constants.g2decay**2 + self.M2g2_HadWH*constants.g2WH**2 * self.M2g1_decay
                   + self.M2g1g2_HadWH*constants.g2WH * self.M2g1g2_decay*constants.g2decay)
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g2WH**2*self.M2g2_HadWH * constants.g2decay**2*self.M2g2_decay)
               )
    def D_g11_g23_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g2_HadWH*constants.g2WH**2 * self.M2g1g2_decay*constants.g2decay + self.M2g1g2_HadWH*constants.g2WH * self.M2g2_decay*constants.g2decay**2)
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g2WH**2*self.M2g2_HadWH * constants.g2decay**2*self.M2g2_decay)
               )
    def D_g10_g24_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g2_HadWH*constants.g2WH**2 * self.M2g2_decay*constants.g2decay**2
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g2WH**2*self.M2g2_HadWH * constants.g2decay**2*self.M2g2_decay)
               )

    ##########
    #g1prime2#
    ##########

    def D_g14_g1prime20_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g1_HadWH*self.M2g1_decay
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g1prime2WH_reco**2*self.M2g1prime2_HadWH * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )
    def D_g13_g1prime21_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadWH * self.M2g1g1prime2_decay*constants.g1prime2decay_reco + self.M2g1g1prime2_HadWH*constants.g1prime2WH_reco * self.M2g1_decay)
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g1prime2WH_reco**2*self.M2g1prime2_HadWH * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )
    def D_g12_g1prime22_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadWH * self.M2g1prime2_decay*constants.g1prime2decay_reco**2 + self.M2g1prime2_HadWH*constants.g1prime2WH_reco**2 * self.M2g1_decay
                   + self.M2g1g1prime2_HadWH*constants.g1prime2WH_reco * self.M2g1g1prime2_decay*constants.g1prime2decay_reco)
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g1prime2WH_reco**2*self.M2g1prime2_HadWH * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )
    def D_g11_g1prime23_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                (self.M2g1prime2_HadWH*constants.g1prime2WH_reco**2 * self.M2g1g1prime2_decay*constants.g1prime2decay_reco + self.M2g1g1prime2_HadWH*constants.g1prime2WH_reco * self.M2g1prime2_decay*constants.g1prime2decay_reco**2)
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g1prime2WH_reco**2*self.M2g1prime2_HadWH * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )
    def D_g10_g1prime24_HadWHdecay(self):
        if self.notdijet: return -999
        return (
                self.M2g1prime2_HadWH*constants.g1prime2WH_reco**2 * self.M2g1prime2_decay*constants.g1prime2decay_reco**2
                 /
                (self.M2g1_HadWH * self.M2g1_decay + constants.g1prime2WH_reco**2*self.M2g1prime2_HadWH * constants.g1prime2decay_reco**2*self.M2g1prime2_decay)
               )

####################################################
#WHdecay hadronic g1^x gi^(4-x) discriminants prime#
####################################################

    ####
    #g4#
    ####

    def D_g14_g40_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g1_HadWH*self.M2g1_decay
                 /
                ((self.M2g1_HadWH + constants.g4WH**2*self.M2g4_HadWH) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )
    def D_g13_g41_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadWH * self.M2g1g4_decay*constants.g4decay + self.M2g1g4_HadWH*constants.g4WH * self.M2g1_decay)
                 /
                ((self.M2g1_HadWH + constants.g4WH**2*self.M2g4_HadWH) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )
    def D_g12_g42_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadWH * self.M2g4_decay*constants.g4decay**2 + self.M2g4_HadWH*constants.g4WH**2 * self.M2g1_decay
                   + self.M2g1g4_HadWH*constants.g4WH * self.M2g1g4_decay*constants.g4decay)
                 /
                ((self.M2g1_HadWH + constants.g4WH**2*self.M2g4_HadWH) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )
    def D_g11_g43_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g4_HadWH*constants.g4WH**2 * self.M2g1g4_decay*constants.g4decay + self.M2g1g4_HadWH*constants.g4WH * self.M2g4_decay*constants.g4decay**2)
                 /
                ((self.M2g1_HadWH + constants.g4WH**2*self.M2g4_HadWH) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )
    def D_g10_g44_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g4_HadWH*constants.g4WH**2 * self.M2g4_decay*constants.g4decay**2
                 /
                ((self.M2g1_HadWH + constants.g4WH**2*self.M2g4_HadWH) * (self.M2g1_decay + constants.g4decay**2*self.M2g4_decay))
               )

    ####
    #g2#
    ####

    def D_g14_g20_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g1_HadWH*self.M2g1_decay
                 /
                ((self.M2g1_HadWH + constants.g2WH**2*self.M2g2_HadWH) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )
    def D_g13_g21_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadWH * self.M2g1g2_decay*constants.g2decay + self.M2g1g2_HadWH*constants.g2WH * self.M2g1_decay)
                 /
                ((self.M2g1_HadWH + constants.g2WH**2*self.M2g2_HadWH) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )
    def D_g12_g22_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadWH * self.M2g2_decay*constants.g2decay**2 + self.M2g2_HadWH*constants.g2WH**2 * self.M2g1_decay
                   + self.M2g1g2_HadWH*constants.g2WH * self.M2g1g2_decay*constants.g2decay)
                 /
                ((self.M2g1_HadWH + constants.g2WH**2*self.M2g2_HadWH) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )
    def D_g11_g23_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g2_HadWH*constants.g2WH**2 * self.M2g1g2_decay*constants.g2decay + self.M2g1g2_HadWH*constants.g2WH * self.M2g2_decay*constants.g2decay**2)
                 /
                ((self.M2g1_HadWH + constants.g2WH**2*self.M2g2_HadWH) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )
    def D_g10_g24_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g2_HadWH*constants.g2WH**2 * self.M2g2_decay*constants.g2decay**2
                 /
                ((self.M2g1_HadWH + constants.g2WH**2*self.M2g2_HadWH) * (self.M2g1_decay + constants.g2decay**2*self.M2g2_decay))
               )

    ##########
    #g1prime2#
    ##########

    def D_g14_g1prime20_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g1_HadWH*self.M2g1_decay
                 /
                ((self.M2g1_HadWH + constants.g1prime2WH_reco**2*self.M2g1prime2_HadWH) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )
    def D_g13_g1prime21_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadWH * self.M2g1g1prime2_decay*constants.g1prime2decay_reco + self.M2g1g1prime2_HadWH*constants.g1prime2WH_reco * self.M2g1_decay)
                 /
                ((self.M2g1_HadWH + constants.g1prime2WH_reco**2*self.M2g1prime2_HadWH) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )
    def D_g12_g1prime22_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1_HadWH * self.M2g1prime2_decay*constants.g1prime2decay_reco**2 + self.M2g1prime2_HadWH*constants.g1prime2WH_reco**2 * self.M2g1_decay
                   + self.M2g1g1prime2_HadWH*constants.g1prime2WH_reco * self.M2g1g1prime2_decay*constants.g1prime2decay_reco)
                 /
                ((self.M2g1_HadWH + constants.g1prime2WH_reco**2*self.M2g1prime2_HadWH) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )
    def D_g11_g1prime23_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                (self.M2g1prime2_HadWH*constants.g1prime2WH_reco**2 * self.M2g1g1prime2_decay*constants.g1prime2decay_reco + self.M2g1g1prime2_HadWH*constants.g1prime2WH_reco * self.M2g1prime2_decay*constants.g1prime2decay_reco**2)
                 /
                ((self.M2g1_HadWH + constants.g1prime2WH_reco**2*self.M2g1prime2_HadWH) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )
    def D_g10_g1prime24_HadWHdecay_prime(self):
        if self.notdijet: return -999
        return (
                self.M2g1prime2_HadWH*constants.g1prime2WH_reco**2 * self.M2g1prime2_decay*constants.g1prime2decay_reco**2
                 /
                ((self.M2g1_HadWH + constants.g1prime2WH_reco**2*self.M2g1prime2_HadWH) * (self.M2g1_decay + constants.g1prime2decay_reco**2*self.M2g1prime2_decay))
               )

##########
#Category#
##########

    def category(self):
        return CJLSTscripts.categoryIchep16(
              nExtraLep,
              nExtraZ,
              nCleanedJetsPt30,
              nCleanedJetsPt30BTagged,
              jetQGLikelihood,
              p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
              p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
              p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
              p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
              pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
              p_HadWH_SIG_ghv1_1_JHUGen_JECNominal,
              p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
              jetPhi,
              ZZMass,
              useQGTagging,
             )

#####################
#Reweighting weights#
#####################

    def MC_weight_plain_xsec(self):
        return self.MC_weight * self.xsec / self.nevents
    def MC_weight_plain_noxsec(self):
        return self.MC_weight / self.nevents
    MC_weight_WplusH_g1 = MC_weight_WminusH_g1 = \
    MC_weight_ttH_kappatilde = MC_weight_ttH_kappakappatilde = \
    MC_weight_HJJ_g2    = MC_weight_HJJ_g4         = MC_weight_HJJ_g2g4 = \
        MC_weight_plain_noxsec

    def MC_weight_ggH(self, index):
        return self.MC_weight * self.reweightingweights[index] * constants.SMXSggH2L2l / self.nevents2L2l[index]
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

    def MC_weight_VBF(self, index):
        return self.MC_weight * self.reweightingweights[index] * constants.SMXSVBF2L2l / self.nevents2L2l[index]
    def MC_weight_VBF_linearcombination(self, indicesandfactors):
        return self.MC_weight * constants.SMXSVBF2L2l * sum(factor*self.reweightingweights[index] for index, factor in indicesandfactors) / sum(factor * self.nevents2L2l[index] for index, factor in indicesandfactors)
    def MC_weight_VBF_g1(self):
        if self.isPOWHEG: return self.MC_weight_plain_xsec()
        return self.MC_weight_VBF(0)
    def MC_weight_VBF_g2(self):
        return self.MC_weight_VBF(1)
    def MC_weight_VBF_g4(self):
        return self.MC_weight_VBF(2)
    def MC_weight_VBF_g1prime2(self):
        return self.MC_weight_VBF(3)
    def MC_weight_VBF_g1g2_dec(self):
        return self.MC_weight_VBF(4)
    def MC_weight_VBF_g1g4_dec(self):
        return self.MC_weight_VBF(5)
    def MC_weight_VBF_g1g1prime2_dec(self):
        return self.MC_weight_VBF(6)
    def MC_weight_VBF_g1g2_prod(self):
        return self.MC_weight_VBF(7)
    def MC_weight_VBF_g1g4_prod(self):
        return self.MC_weight_VBF(8)
    def MC_weight_VBF_g1g1prime2_prod(self):
        return self.MC_weight_VBF(9)
    def MC_weight_VBF_g1g2_proddec_pi(self):
        return self.MC_weight_VBF(10)
    def MC_weight_VBF_g1g4_proddec_pi(self):
        return self.MC_weight_VBF(11)
    def MC_weight_VBF_g1g1prime2_proddec(self):
        return self.MC_weight_VBF(12)

    MC_weight_VBF_g1g2_dec_pi = ReweightingSample("VBF", "fa2dec-0.5").get_MC_weight_function()
    MC_weight_VBF_g1g2_prod_pi = ReweightingSample("VBF", "fa2prod-0.5").get_MC_weight_function()
    MC_weight_VBF_g1g2_proddec = ReweightingSample("VBF", "fa2proddec0.5").get_MC_weight_function()

    def MC_weight_ZH(self, index):
        return self.MC_weight * self.reweightingweights[index] * constants.SMXSZH2L2l / self.nevents2L2l[index]
    def MC_weight_ZH_g1(self):
        if self.isPOWHEG: return self.MC_weight_plain_xsec()
        return self.MC_weight_ZH(0)
    def MC_weight_ZH_g2(self):
        return self.MC_weight_ZH(1)
    def MC_weight_ZH_g4(self):
        return self.MC_weight_ZH(2)
    def MC_weight_ZH_g1prime2(self):
        return self.MC_weight_ZH(3)
    def MC_weight_ZH_g1g2_dec(self):
        return self.MC_weight_ZH(4)
    def MC_weight_ZH_g1g4_dec(self):
        return self.MC_weight_ZH(5)
    def MC_weight_ZH_g1g1prime2_dec(self):
        return self.MC_weight_ZH(6)
    def MC_weight_ZH_g1g2_prod(self):
        return self.MC_weight_ZH(7)
    def MC_weight_ZH_g1g4_prod(self):
        return self.MC_weight_ZH(8)
    def MC_weight_ZH_g1g1prime2_prod(self):
        return self.MC_weight_ZH(9)
    def MC_weight_ZH_g1g2_proddec_pi(self):
        return self.MC_weight_ZH(10)
    def MC_weight_ZH_g1g4_proddec_pi(self):
        return self.MC_weight_ZH(11)
    def MC_weight_ZH_g1g1prime2_proddec(self):
        return self.MC_weight_ZH(12)

    def MC_weight_WH(self, index):
        return self.MC_weight * self.reweightingweights[index] * constants.SMXSWH2L2l / self.nevents2L2l[index]
    def MC_weight_WH_g1(self):
        return self.MC_weight_WH(0)
    def MC_weight_WH_g2(self):
        return self.MC_weight_WH(1)
    def MC_weight_WH_g4(self):
        return self.MC_weight_WH(2)
    def MC_weight_WH_g1prime2(self):
        return self.MC_weight_WH(3)
    def MC_weight_WH_g1g2_dec(self):
        return self.MC_weight_WH(4)
    def MC_weight_WH_g1g4_dec(self):
        return self.MC_weight_WH(5)
    def MC_weight_WH_g1g1prime2_dec(self):
        return self.MC_weight_WH(6)
    def MC_weight_WH_g1g2_prod(self):
        return self.MC_weight_WH(7)
    def MC_weight_WH_g1g4_prod(self):
        return self.MC_weight_WH(8)
    def MC_weight_WH_g1g1prime2_prod(self):
        return self.MC_weight_WH(9)
    def MC_weight_WH_g1g2_proddec_pi(self):
        return self.MC_weight_WH(10)
    def MC_weight_WH_g1g4_proddec_pi(self):
        return self.MC_weight_WH(11)
    def MC_weight_WH_g1g1prime2_proddec(self):
        return self.MC_weight_WH(12)

    def MC_weight_ttH_kappa(self):
        if self.isPOWHEG: return self.MC_weight_plain_xsec()
        return self.MC_weight_plain_noxsec()

    def MC_weight_ggZZ(self):
        KFactor = self.tree.KFactor_QCD_ggZZ_Nominal
        return self.MC_weight * self.xsec * KFactor / self.nevents
    def MC_weight_qqZZ(self):
        KFactor = self.tree.KFactor_EW_qqZZ * self.tree.KFactor_QCD_qqZZ_M
        return self.MC_weight * self.xsec * KFactor / self.nevents
    MC_weight_VBFbkg = MC_weight_plain_xsec
    def MC_weight_ZX(self):
        self.LepPt, self.LepEta, self.LepLepId = self.tree.LepPt, self.tree.LepEta, self.tree.LepLepId
        return ZX.fakeRate13TeV(self.LepPt[2],self.LepEta[2],self.LepLepId[2]) * ZX.fakeRate13TeV(self.LepPt[3],self.LepEta[3],self.LepLepId[3])

    def getweightfunction(self, sample):
        return getattr(self, sample.weightname())
        raise RuntimeError("{} does not work!".format(sample))

    def getmainweightfunction(self):
        return getattr(self, "MC_weight_{}".format(self.productionmode))

    def initlists(self):
        self.toaddtotree = [
            "D_bkg_0plus",
            "D_bkg_0plus_ResUp",
            "D_bkg_0plus_ResDown",
            "D_bkg_0plus_ScaleUp",
            "D_bkg_0plus_ScaleDown",
            "D_2jet_0plus",
            "D_2jet_0minus",
            "D_0minus_decay",
            "D_CP_decay",
            "D_g2_decay",
            "D_g1g2_decay",
            "D_g1prime2_decay",
            "D_g1g1prime2_decay",
        ]

        self.exceptions = [
            "cconstantforDbkg",
            "cconstantforD2jet",
            "checkfunctions",
            "exceptions",
            "getweightfunction",
            "getmainweightfunction",
            "hypothesis",
            "initlists",
            "isbkg",
            "isdata",
            "isdummy",
            "isPOWHEG",
            "isZX",
            "MC_weight_ggH",
            "MC_weight_VBF",
            "MC_weight_VBF_linearcombination",
            "MC_weight_WH",
            "MC_weight_ZH",
            "MC_weight_plain_xsec",
            "MC_weight_plain_noxsec",
            "minevent",
            "nevents",
            "nevents2L2l",
            "next",
            "onlyweights",
            "passesblindcut",
            "productionmode",
            "toaddtotree",
            "toaddtotree_int",
            "tree",
            "treesample",
            "unblind",
            "weightfunctions",
            "xsec",
        ]

        proddiscriminants = [
            "D_0minus_{prod}{suffix}",
            "D_CP_{prod}{suffix}",
            "D_g2_{prod}{suffix}",
            "D_g1g2_{prod}{suffix}",
            "D_g1prime2_{prod}{suffix}",
            "D_g1g1prime2_{prod}{suffix}",
            "D_0minus_{prod}decay{suffix}",
            "D_g2_{prod}decay{suffix}",
            "D_g1prime2_{prod}decay{suffix}",
        ]
        prodcomponentdiscriminants = [
            "D_g1{}_{}{}_{{prod}}decay{{suffix}}{}".format(i, gj, 4-i, prime)
                for prime in ("", "_prime")
                for gj in ("g4", "g2", "g1prime2")
                for i in range(5)
        ]
        for prod, suffix in (("VBF", ""), ("ZH", "_hadronic"), ("WH", "_hadronic"), ("VH", "_hadronic")):
            self.toaddtotree += [_.format(prod=prod, suffix=suffix) for _ in proddiscriminants]
            if prod != "VH":
                self.exceptions += [_.format(prod=prod, suffix=suffix) for _ in prodcomponentdiscriminants]

        self.toaddtotree_int = [
            "category",
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
            ReweightingSample("VBF", "a2"),
            ReweightingSample("VBF", "0-"),
            ReweightingSample("VBF", "L1"),
            ReweightingSample("VBF", "fa2dec0.5"),
            ReweightingSample("VBF", "fa3dec0.5"),
            ReweightingSample("VBF", "fL1dec0.5"),
            ReweightingSample("VBF", "fa2prod0.5"),
            ReweightingSample("VBF", "fa3prod0.5"),
            ReweightingSample("VBF", "fL1prod0.5"),
            ReweightingSample("VBF", "fa2proddec-0.5"),
            ReweightingSample("VBF", "fa3proddec-0.5"),
            ReweightingSample("VBF", "fL1proddec0.5"),
            ReweightingSample("VBF", "fa2dec-0.5"),
            ReweightingSample("VBF", "fa2prod-0.5"),
            ReweightingSample("VBF", "fa2proddec0.5"),
            ReweightingSample("ZH", "0+"),
            ReweightingSample("ZH", "a2"),
            ReweightingSample("ZH", "0-"),
            ReweightingSample("ZH", "L1"),
            ReweightingSample("ZH", "fa2dec0.5"),
            ReweightingSample("ZH", "fa3dec0.5"),
            ReweightingSample("ZH", "fL1dec0.5"),
            ReweightingSample("ZH", "fa2prod0.5"),
            ReweightingSample("ZH", "fa3prod0.5"),
            ReweightingSample("ZH", "fL1prod0.5"),
            ReweightingSample("ZH", "fa2proddec-0.5"),
            ReweightingSample("ZH", "fa3proddec-0.5"),
            ReweightingSample("ZH", "fL1proddec0.5"),
            ReweightingSample("WH", "0+"),
            ReweightingSample("WH", "a2"),
            ReweightingSample("WH", "0-"),
            ReweightingSample("WH", "L1"),
            ReweightingSample("WH", "fa2dec0.5"),
            ReweightingSample("WH", "fa3dec0.5"),
            ReweightingSample("WH", "fL1dec0.5"),
            ReweightingSample("WH", "fa2prod0.5"),
            ReweightingSample("WH", "fa3prod0.5"),
            ReweightingSample("WH", "fL1prod0.5"),
            ReweightingSample("WH", "fa2proddec-0.5"),
            ReweightingSample("WH", "fa3proddec-0.5"),
            ReweightingSample("WH", "fL1proddec0.5"),
            ReweightingSample("ttH", "0+"),
            ReweightingSample("ttH", "0-"),
            ReweightingSample("ttH", "fCP0.5"),
            ReweightingSample("HJJ", "0+"),
            ReweightingSample("HJJ", "0-"),
            ReweightingSample("HJJ", "fCP0.5"),
            ReweightingSample("WplusH", "0+"),
            ReweightingSample("WminusH", "0+"),
            ReweightingSample("ggZZ", "2e2mu"),  #flavor doesn't matter
            ReweightingSample("qqZZ"),
            ReweightingSample("VBF bkg", "2e2mu"),  #flavor doesn't matter
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

    def onlyweights(self):
        """Call this to only add the weights and ZZMass to the new tree"""
        #only want the weight, and ZZMass for the range
        reweightingweightnames = [sample.weightname() for sample in self.treesample.reweightingsamples()]

        categoryingredients = [
            "nExtraLep",
            "nExtraZ",
            "nCleanedJetsPt30",
            "nCleanedJetsPt30BTagged",
            "jetQGLikelihood",
            "p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal",
            "p_JQCD_SIG_ghg2_1_JHUGen_JECNominal",
            "p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal",
            "p_JVBF_SIG_ghv1_1_JHUGen_JECNominal",
            "pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal",
            "p_HadWH_SIG_ghv1_1_JHUGen_JECNominal",
            "p_HadZH_SIG_ghz1_1_JHUGen_JECNominal",
            "jetPhi",
            "ZZMass",
            "useQGTagging",
        ]

        for lst in (self.toaddtotree, self.toaddtotree_int):
            for name in lst[:]:
                if name not in reweightingweightnames + ["category"]:
                    self.exceptions.append(name)
                    lst.remove(name)
        self.tree.GetEntry(0)  #to initialize the vectors to something other than (vector<float>*)0 which gives a segfault in python
        self.tree.SetBranchStatus("*", 0)
        self.tree.SetBranchStatus("ZZMass", 1)
        self.tree.SetBranchStatus("Z*Flav", 1)
        for variable in categoryingredients:
            self.tree.SetBranchStatus(variable, 1)
        for variable in self.treesample.weightingredients():
            self.tree.SetBranchStatus(variable, 1)

    def checkfunctions(self, couplings=None):

        #some cross checking in case of stupid mistakes
        #if a function is added in the class but not added to toaddtotree
        #all member variables, unless they have __, should be added to either toaddtotree or exceptions
        notanywhere, inboth, nonexistent, multipletimes = [], [], [], []
        for key in set(list(type(self).__dict__) + list(self.__dict__) + self.toaddtotree+self.toaddtotree_int + self.exceptions):
            if key.startswith("__"): continue
            if key.startswith("_abc"): continue
            if any(key.startswith("_{}__".format(cls.__name__)) for cls in type(self).__mro__):
                continue
            if key not in self.exceptions and key not in self.toaddtotree+self.toaddtotree_int and (key in self.__dict__ or key in type(self).__dict__):
                notanywhere.append(key)
            if key in self.toaddtotree+self.toaddtotree_int and key in self.exceptions:
                inboth.append(key)
            if key not in type(self).__dict__ and key not in self.__dict__:
                nonexistent.append(key)
        for key, occurences in Counter(self.toaddtotree+self.toaddtotree_int + self.exceptions).iteritems():
            if occurences >= 2 and key not in inboth or occurences >= 3: multipletimes.append(key)
        error = ""
        if notanywhere: error += "the following items are not in toaddtotree or exceptions! " + ", ".join(notanywhere) + "\n"
        if inboth: error += "the following items are in both toaddtotree and exceptions! " + ", ".join(inboth) + "\n"
        if nonexistent: error += "the following items are in toaddtotree or exceptions, but don't exist! " + ", ".join(nonexistent) + "\n"
        if multipletimes: error += "the following items appear multiple times in toaddtotree or exceptions! " + ", ".join(multipletimes) + "\n"
        if error:
            raise SyntaxError(error)

        #cross checking that the order of weights is defined right
        if (
                self
            and self.treesample.productionmode in ("ggH", "VBF", "ZH", "WH")
            and self.treesample.alternategenerator is None
           ):
            if couplings is None:
                raise SyntaxError("No couplings tree for {}".format(self.treesample))

            for entry in self:
                if self.getmainweightfunction()(0) != 0:
                    break
            else:
                raise SyntaxError("All weights are 0?!")

            if len(self.treesample.reweightingsamples()) != len(self.weightfunctions):
                raise SyntaxError("Something is very very wrong.  {} {}".format(len(self.treesample.reweightingsamples), len(self.weightfunctions)))

            couplings.GetEntry(0)
            for i, (sample, function) in enumerate(zip(self.treesample.directreweightingsamples(), self.weightfunctions)):
                if sample.indexinreweightingweights != i:
                    raise SyntaxError("{!r}.indexinreweightingweights == {} when it should be {}!".format(sample, sample.indexinreweightingweights, i))
                if self.getweightfunction(sample)() != self.getmainweightfunction()(i):
                    raise SyntaxError("{}() == {}, but {}({}) == {}!\nCheck the order of reweightingsamples or the weight functions!".format(self.getweightfunction(sample).__name__, self.getweightfunction(sample)(), self.getmainweightfunction().__name__, i, self.getmainweightfunction()(i)))
                print couplings.ghz1Re[i], couplings.ghz2Re[i], couplings.ghz4Re[i], couplings.ghz1_prime2Re[i], sample.g1, sample.g2, sample.g4, sample.g1prime2
                if sample.productionmode != "ggH" and self.treesample.production <= "160729":
                    continue

                gs = g1, g2, g4, g1prime2 = sample.g1, sample.g2, sample.g4, sample.g1prime2
                if len([g for g in gs if g!=0]) == 1 and self.productionmode == "ggH":
                    #to set gi=1 for the pure samples to compare to the couplings tree
                    #should just check ratios, but that's annoying when most of them are 0
                    #can only do this for ggH, because VBF and VH have the linear combination weight branches
                    # and they have to be the correct couplings not just ratios
                    gs = g1, g2, g4, g1prime2 = [1 if g else 0 for g in gs]

                if not (couplings.spin[i] == 0 and couplings.ghz1Re[i] == g1 and couplings.ghz2Re[i] == g2 and couplings.ghz4Re[i] == g4 and couplings.ghz1_prime2Re[i] == g1prime2):
                    raise SyntaxError("Order of reweightingsamples or the weight functions seems wrong!  Check entry {} in couplings with respect to {}.".format(i, sample))


    passesblindcut = config.blindcut

if __name__ == '__main__':
    class DummyTree(object):
        def GetEntries(self):
            return 0
        def GetEntry(self, entry): pass
        xsec = 0
    class DummySample(object):
        productionmode = "graviton fusion"
        hypothesis = "spin 3"
        alternategenerator = "magic"
        def isbkg(self): return False
        def isdata(self): return False
        def isZX(self): return False
        def onlyweights(self): return False
        def reweightingsamples(self): return []
        def weightname(self): return "__init__"
    TreeWrapper(DummyTree(), DummySample(), None, None, None)
    print "You are good to go!"

