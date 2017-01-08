#!/usr/bin/env python
from array import array
import CJLSTscripts
from collections import Counter, Iterator
import config
import constants
import numpy
import resource
from samples import ReweightingSample
import sys
from utilities import callclassinitfunctions
import ZX

resource.setrlimit(resource.RLIMIT_STACK, (2**29,-1))
sys.setrecursionlimit(10000)

#to pass to the category code when there are no jets
dummyfloatstar = array('f', [0])

@callclassinitfunctions("initweightfunctions")
class TreeWrapper(Iterator):

    def __init__(self, tree, treesample, Counters, minevent=0, maxevent=None, isdummy=False):
        """
        tree - a TTree object
        treesample - which sample the TTree was created from
        Counters - from the CJLST file
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

        self.nevents = self.nevents2L2l = self.cutoffs = None
        if Counters is not None:
            self.nevents = Counters.GetBinContent(40)

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
        self.checkfunctions()

        self.preliminaryloop()

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
                raise StopIteration

            if self.isdata:
                self.overallEventWeight = 1
            elif self.isZX:
                CRflag = t.CRflag
                if CRflag and ZX.test_bit(CRflag, ZX.CRZLLss):
                    self.overallEventWeight = 1
                else:
                    self.overallEventWeight = 0
            else:
                self.overallEventWeight = t.overallEventWeight

            self.flavor = abs(t.Z1Flav*t.Z2Flav)

            if self.overallEventWeight:
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
                    attr = "p_m4l_{}_{}{}".format(a, b, c)
                    setattr(self, attr, getattr(t, attr))

        #express in terms of |M|^2, this will make life easier
        self.M2qqZZ = t.p_QQB_BKG_MCFM

        self.M2g1_decay         = t.p_GG_SIG_ghg2_1_ghz1_1_JHUGen
        self.M2g4_decay         = t.p_GG_SIG_ghg2_1_ghz4_1_JHUGen
        self.M2g1g4_decay       = t.p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen
        self.M2g2_decay         = t.p_GG_SIG_ghg2_1_ghz2_1_JHUGen
        self.M2g1g2_decay       = t.p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen
        self.M2g1prime2_decay   = t.p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen / 1e4**2
        self.M2g1g1prime2_decay = t.p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen / 1e4

        self.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal = \
        self.M2g1_VBF         = t.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal
        self.M2g4_VBF         = t.p_JJVBF_SIG_ghv4_1_JHUGen_JECNominal
        self.M2g1g4_VBF       = t.p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECNominal
        self.M2g2_VBF         = t.p_JJVBF_SIG_ghv2_1_JHUGen_JECNominal
        self.M2g1g2_VBF       = t.p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECNominal
        self.M2g1prime2_VBF   = t.p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal / 1e4**2
        self.M2g1g1prime2_VBF = t.p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_JECNominal / 1e4

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
        self.M2g1prime2_HadZH   = t.p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECNominal / 1e4**2
        self.M2g1g1prime2_HadZH = t.p_HadZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen_JECNominal / 1e4

        self.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal = \
        self.M2g1_HadWH         = t.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal
        self.M2g4_HadWH         = t.p_HadWH_SIG_ghw4_1_JHUGen_JECNominal
        self.M2g1g4_HadWH       = t.p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECNominal
        self.M2g2_HadWH         = t.p_HadWH_SIG_ghw2_1_JHUGen_JECNominal
        self.M2g1g2_HadWH       = t.p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECNominal
        self.M2g1prime2_HadWH   = t.p_HadWH_SIG_ghw1prime2_1E4_JHUGen_JECNominal / 1e4**2
        self.M2g1g1prime2_HadWH = t.p_HadWH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen_JECNominal / 1e4

        #Gen MEs
        for weightname in self.genMEs:
            setattr(self, weightname, getattr(t, weightname))

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

        self.notdijet = self.M2g1_VBF == 0
        return self

    def __len__(self):
        return self.__length

    def Show(self, *args, **kwargs):
        self.tree.Show(*args, **kwargs)

##########################
#background discriminants#
##########################

    def D_bkg_0plus(self):
        try:
            return self.M2g1_decay*self.p_m4l_SIG / (self.M2g1_decay*self.p_m4l_SIG  + self.M2qqZZ*self.p_m4l_BKG*self.cconstantforDbkg)
        except ZeroDivisionError:
            return 0
    def D_bkg_0plus_ResUp(self):
        try:
            return self.M2g1_decay*self.p_m4l_SIG_ResUp / (self.M2g1_decay*self.p_m4l_SIG_ResUp  + self.M2qqZZ*self.p_m4l_BKG_ResUp*self.cconstantforDbkg)
        except ZeroDivisionError:
            return 0
    def D_bkg_0plus_ResDown(self):
        try:
            return self.M2g1_decay*self.p_m4l_SIG_ResDown / (self.M2g1_decay*self.p_m4l_SIG_ResDown  + self.M2qqZZ*self.p_m4l_BKG_ResDown*self.cconstantforDbkg)
        except ZeroDivisionError:
            return 0
    def D_bkg_0plus_ScaleUp(self):
        try:
            return self.M2g1_decay*self.p_m4l_SIG_ScaleUp / (self.M2g1_decay*self.p_m4l_SIG_ScaleUp  + self.M2qqZZ*self.p_m4l_BKG_ScaleUp*self.cconstantforDbkg)
        except ZeroDivisionError:
            return 0
    def D_bkg_0plus_ScaleDown(self):
        try:
            return self.M2g1_decay*self.p_m4l_SIG_ScaleDown / (self.M2g1_decay*self.p_m4l_SIG_ScaleDown  + self.M2qqZZ*self.p_m4l_BKG_ScaleDown*self.cconstantforDbkg)
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
              self.nExtraLep,
              self.nExtraZ,
              self.nCleanedJetsPt30,
              self.nCleanedJetsPt30BTagged,
              self.jetQGLikelihood,
              self.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
              self.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
              self.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
              self.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
              self.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
              self.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
              self.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
              self.jetPhi,
              self.ZZMass,
              config.useQGTagging,
             )

#####################
#Reweighting weights#
#####################

    def getweightfunction(self, sample):
        return getattr(self, sample.weightname())

    @classmethod
    def initweightfunctions(cls):
        for sample in cls.allsamples:
            setattr(cls, sample.weightname(), sample.get_MC_weight_function())

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
            "allsamples",
            "cconstantforDbkg",
            "cconstantforD2jet",
            "checkfunctions",
            "cutoffs",
            "exceptions",
            "genMEs",
            "getweightfunction",
            "hypothesis",
            "initlists",
            "initweightfunctions",
            "isbkg",
            "isdata",
            "isdummy",
            "isPOWHEG",
            "isZX",
            "minevent",
            "nevents",
            "nevents2L2l",
            "next",
            "onlyweights",
            "passesblindcut",
            "productionmode",
            "preliminaryloop",
            "toaddtotree",
            "toaddtotree_int",
            "tree",
            "treesample",
            "Show",
            "unblind",
            "weightfunctions",
            "xsec",
        ]

        proddiscriminants = [
            "D_0minus_{prod}",
            "D_CP_{prod}",
            "D_g2_{prod}",
            "D_g1g2_{prod}",
            "D_g1prime2_{prod}",
            "D_g1g1prime2_{prod}",
            "D_0minus_{prod}decay",
            "D_g2_{prod}decay",
            "D_g1prime2_{prod}decay",
        ]
        prodcomponentdiscriminants = [
            "D_g1{}_{}{}_{{prod}}decay{}".format(i, gj, 4-i, prime)
                for prime in ("", "_prime")
                for gj in ("g4", "g2", "g1prime2")
                for i in range(5)
        ]
        for prod in ("VBF", "HadZH", "HadWH", "HadVH"):
            self.toaddtotree += [_.format(prod=prod) for _ in proddiscriminants]
            if prod != "HadVH":
                self.exceptions += [_.format(prod=prod) for _ in prodcomponentdiscriminants]

        self.toaddtotree_int = [
            "category",
        ]

        reweightingweightnames = [sample.weightname() for sample in self.treesample.reweightingsamples()]
        allreweightingweightnames = [sample.weightname() for sample in self.allsamples]
        for name in reweightingweightnames:
            if name not in allreweightingweightnames:
                raise ValueError("{} not in allreweightingweightnames!".format(name))
        for sample in self.allsamples:
            if sample.weightname() in self.toaddtotree or sample.weightname() in self.exceptions: continue
            if sample.weightname() in reweightingweightnames:
                self.toaddtotree.append(sample.weightname())
            else:
                self.exceptions.append(sample.weightname())

        self.genMEs = []
        if not self.isbkg:
            for sample in self.treesample.reweightingsamples():
                for factor in sample.MC_weight_terms:
                    for weightname, couplingsq in factor:
                        self.genMEs.append(weightname)

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
            "p_HadWH_SIG_ghw1_1_JHUGen_JECNominal",
            "p_HadZH_SIG_ghz1_1_JHUGen_JECNominal",
            "jetPhi",
            "ZZMass",
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

    def checkfunctions(self):

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

    def preliminaryloop(self):
        """
        Do the initial loops through the tree to find, for each hypothesis,
        the cutoff and then the sum of weights for 2L2l
        """
        if self.isbkg: return
        print "Doing initial loop through tree"
        self.tree.SetBranchStatus("*", 0)
        for weightname in self.genMEs:
            self.tree.SetBranchStatus(weightname, 1)
        self.tree.SetBranchStatus("GenZ*Flav", 1)

        functionsandarrays = {sample: (sample.get_MC_weight_function(reweightingonly=True), []) for sample in self.treesample.reweightingsamples()}
        is2L2l = []
        flavs2L2l = {11*11*13*13, 11*11*15*15, 13*13*15*15}
        values = functionsandarrays.values()
        #will fail if multiple have the same str() which would make no sense
        assert len(functionsandarrays) == len(self.treesample.reweightingsamples())

        length = self.tree.GetEntries()
        for i, entry in enumerate(self.tree, start=1):
            is2L2l.append(entry.GenZ1Flav * entry.GenZ2Flav in flavs2L2l)
            for function, array in values:
                array.append(function(entry))
            if i % 10000 == 0 or i == length:
                print i, "/", length, "   (preliminary run)"
                break

        self.cutoffs = {}
        self.nevents2L2l = {}

        for sample, (function, array) in functionsandarrays.iteritems():
            array = numpy.array(array)
            cutoff = self.cutoffs[str(sample)] = numpy.percentile(array, 99.99)
            self.nevents2L2l[str(sample)] = sum(
                                                weight if weight < cutoff else cutoff**2/weight
                                                     for weight, isthis2L2l in zip(array, is2L2l)
                                                     if isthis2L2l
                                               )

        self.tree.SetBranchStatus("*", 1)

        print "Cutoffs:"
        for sample, cutoff in self.cutoffs.iteritems():
             print "    {:15} {}".format(sample, cutoff)
        print "nevents2L2l:"
        for sample, nevents in self.nevents2L2l.iteritems():
             print "    {:15} {}".format(sample, nevents)

    passesblindcut = config.blindcut

    allsamples = [    #all samples that should have weight functions defined in this class
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

