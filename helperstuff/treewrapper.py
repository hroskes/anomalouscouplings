#!/usr/bin/env python
from abc import abstractmethod, abstractproperty
from array import array
from collections import Counter, Iterator
import inspect
from itertools import chain, izip, izip_longest
from math import sqrt
import resource
import sys

import numpy
import ROOT

import categorization
import CJLSTscripts
import config
import constants
import discriminants
import enums
import STXS
import xrd
import ZX
from gconstants import gconstant
from makesystematics import MakeBtagSystematics, MakeJECSystematics, MakeSystematics
from samples import ReweightingSample, ReweightingSamplePlus, Sample
from utilities import cache_instancemethod, callclassinitfunctions, deprecate, Fake_LSF_creating, getmembernames, MultiplyCounter, product, TFile, tlvfromptetaphim

resource.setrlimit(resource.RLIMIT_STACK, (2**29,-1))
sys.setrecursionlimit(100000)

#to pass to the category code when there are no jets
dummyfloatstar = array('f', [0])

class TreeWrapperBase(Iterator):
    def __init__(self, treesample, minevent=0, maxevent=None):
        self.treesample = treesample
        self.productionmode = str(treesample.productionmode)
        self.hypothesis = str(treesample.hypothesis)

        self.year = treesample.production.year

        if self.isdata:
            self.unblind = config.unblinddistributions
        else:
            self.unblind = True

        self.printevery = 10000
        if self.isZX:
            self.printevery = 1000

        self.cconstantforDbkg = self.cconstantforD2jet = self.cconstantforDHadWH = self.cconstantforDHadZH = None

        self.minevent = minevent
        self.maxevent = maxevent

        self.GEN = self.treesample.production.GEN

        self.initlists()
        self.checkfunctions()

    @abstractmethod
    def __len__(self): pass
    @abstractmethod
    def initlists(self): pass
    @abstractmethod
    def Show(self): pass

    @property
    @cache_instancemethod
    def isdummy(self):
        if self.isZX and not config.usedata: return True
        if self.isdata and not config.showblinddistributions: return True
        return False

    @property
    @cache_instancemethod
    def isdata(self): return self.treesample.isdata()
    @property
    @cache_instancemethod
    def isbkg(self): return not self.isdata and self.treesample.isbkg
    @property
    @cache_instancemethod
    def isZX(self): return self.treesample.isZX()
    @property
    @cache_instancemethod
    def isalternate(self): return self.treesample.alternategenerator in ("POWHEG", "MINLO", "NNLOPS") or self.treesample.pythiasystematic is not None


    def checkfunctions(self):
        #some cross checking in case of stupid mistakes
        #if a function is added in the class but not added to toaddtotree
        #all member variables, unless they have __, should be added to either toaddtotree or exceptions
        notanywhere, inboth, nonexistent, multipletimes = [], [], [], []
        toaddtotree = self.toaddtotree+self.toaddtotree_int+self.toaddtotree_float
        for key in set(getmembernames(self) + toaddtotree + self.exceptions):
            if key.startswith("__"): continue
            if key.startswith("_abc"): continue
            if any(key.startswith("_{}__".format(cls.__name__)) for cls in type(self).__mro__):
                continue
            if key not in self.exceptions and key not in toaddtotree and key in getmembernames(self):
                notanywhere.append(key)
            if key in toaddtotree and key in self.exceptions:
                inboth.append(key)
            if key not in getmembernames(self):
                nonexistent.append(key)
        for key, occurences in Counter(toaddtotree + self.exceptions).iteritems():
            if occurences >= 2 and key not in inboth or occurences >= 3: multipletimes.append(key)
        error = ""
        if notanywhere: error += "the following items are not in toaddtotree or exceptions! " + ", ".join(notanywhere) + "\n"
        if inboth: error += "the following items are in both toaddtotree and exceptions! " + ", ".join(inboth) + "\n"
        if nonexistent: error += "the following items are in toaddtotree or exceptions, but don't exist! " + ", ".join(nonexistent) + "\n"
        if multipletimes: error += "the following items appear multiple times in toaddtotree or exceptions! " + ", ".join(multipletimes) + "\n"
        if error:
            raise SyntaxError(error)

    def per_event_scale_factor(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "ttH", "bbH", "tqH", "WplusH", "WminusH"): return 1
        if self.productionmode == "ggZZ": return self.KFactor_QCD_ggZZ_Nominal
        if self.productionmode == "qqZZ":
            if self.GEN: return 1
            return self.KFactor_EW_qqZZ * self.KFactor_QCD_qqZZ_M
        assert False, self

    def MC_weight_nominal(self):
        if self.productionmode == "ZX" or self.productionmode == "data":
            assert self.overallEventWeight == 1, self.overallEventWeight
            if self.productionmode == "data": return 1
            LepPt, LepEta, LepLepId = self.tree.LepPt, self.tree.LepEta, self.tree.LepLepId
            return ZX.normalizeZX(self.year, self.tree.Z1Flav, self.tree.Z2Flav) * ZX.getfakerate(self.year, LepPt[2], LepEta[2], LepLepId[2]) * ZX.getfakerate(self.year, LepPt[3], LepEta[3], LepLepId[3])

        pb_to_fb = 1000
        return self.overallEventWeight * self.genxsec * self.genBR * pb_to_fb * self.per_event_scale_factor() / self.nevents

##########################
#background discriminants#
##########################

    def D_bkg(self):
        if self.p_m4l_SIG <= 0: return -999
        return self.M2g1_decay*self.p_m4l_SIG / (self.M2g1_decay*self.p_m4l_SIG  + self.M2qqZZ*self.p_m4l_BKG*self.cconstantforDbkg)
    def D_bkg_ResUp(self):
        if self.p_m4l_SIG_ResUp <= 0: return -999
        return self.M2g1_decay*self.p_m4l_SIG_ResUp / (self.M2g1_decay*self.p_m4l_SIG_ResUp  + self.M2qqZZ*self.p_m4l_BKG_ResUp*self.cconstantforDbkg)
    def D_bkg_ResDown(self):
        if self.p_m4l_SIG_ResDown <= 0: return -999
        return self.M2g1_decay*self.p_m4l_SIG_ResDown / (self.M2g1_decay*self.p_m4l_SIG_ResDown  + self.M2qqZZ*self.p_m4l_BKG_ResDown*self.cconstantforDbkg)
    def D_bkg_ScaleUp(self):
        if self.p_m4l_SIG_ScaleUp <= 0: return -999
        return self.M2g1_decay*self.p_m4l_SIG_ScaleUp / (self.M2g1_decay*self.p_m4l_SIG_ScaleUp  + self.M2qqZZ*self.p_m4l_BKG_ScaleUp*self.cconstantforDbkg)
    def D_bkg_ScaleDown(self):
        if self.p_m4l_SIG_ScaleDown <= 0: return -999
        return self.M2g1_decay*self.p_m4l_SIG_ScaleDown / (self.M2g1_decay*self.p_m4l_SIG_ScaleDown  + self.M2qqZZ*self.p_m4l_BKG_ScaleDown*self.cconstantforDbkg)

    def D_2jet_0plus(self):
        if self.notdijet: return -999
        return self.M2g1_VBF / (self.M2g1_VBF + self.M2g2_HJJ*self.cconstantforD2jet)
    def D_2jet_0minus(self):
        if self.notdijet: return -999
        return self.M2g4_VBF*self.g4VBF_m4l**2 / (self.M2g4_VBF*self.g4VBF_m4l**2 + self.M2g2_HJJ*self.cconstantforD2jet)
    def D_2jet_a2(self):
        if self.notdijet: return -999
        return self.M2g2_VBF*self.g2VBF_m4l**2 / (self.M2g2_VBF*self.g2VBF_m4l**2 + self.M2g2_HJJ*self.cconstantforD2jet)
    def D_2jet_L1(self):
        if self.notdijet: return -999
        return self.M2g1prime2_VBF*self.g1prime2VBF_m4l**2 / (self.M2g1prime2_VBF*self.g1prime2VBF_m4l**2 + self.M2g2_HJJ*self.cconstantforD2jet)
    def D_2jet_L1Zg(self):
        if self.notdijet: return -999
        return self.M2ghzgs1prime2_VBF*self.ghzgs1prime2VBF_m4l**2 / (self.M2ghzgs1prime2_VBF*self.ghzgs1prime2VBF_m4l**2 + self.M2g2_HJJ*self.cconstantforD2jet)

    def D_HadWH_0plus(self):
        if self.notdijet: return -999
        return (self.M2g1_HadWH / (self.M2g1_HadWH + self.M2g2_HJJ*self.cconstantforDHadWH / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal))
    def D_HadWH_0minus(self):
        if self.notdijet: return -999
        return (self.M2g4_HadWH*self.g4WH_m4l**2 / (self.M2g4_HadWH*self.g4WH_m4l**2 + self.M2g2_HJJ*self.cconstantforDHadWH / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal))
    def D_HadWH_a2(self):
        if self.notdijet: return -999
        return (self.M2g2_HadWH*self.g2WH_m4l**2 / (self.M2g2_HadWH*self.g2WH_m4l**2 + self.M2g2_HJJ*self.cconstantforDHadWH / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal))
    def D_HadWH_L1(self):
        if self.notdijet: return -999
        return (self.M2g1prime2_HadWH*self.g1prime2WH_m4l**2 / (self.M2g1prime2_HadWH*self.g1prime2WH_m4l**2 + self.M2g2_HJJ*self.cconstantforDHadWH / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal))
    def D_HadWH_L1Zg(self):
        if self.notdijet: return -999
        #return (self.M2ghzgs1prime2_HadWH*self.ghzgs1prime2WH_m4l**2 / (self.M2ghzgs1prime2_HadWH*self.ghzgs1prime2WH_m4l**2 + self.M2g2_HJJ*self.cconstantforDHadWH / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal))
        return 0

    def D_HadZH_0plus(self):
        if self.notdijet: return -999
        return self.M2g1_HadZH / (self.M2g1_HadZH + self.M2g2_HJJ*self.cconstantforDHadZH / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal)
    def D_HadZH_0minus(self):
        if self.notdijet: return -999
        return self.M2g4_HadZH*self.g4ZH_m4l**2 / (self.M2g4_HadZH*self.g4ZH_m4l**2 + self.M2g2_HJJ*self.cconstantforDHadZH / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal)
    def D_HadZH_a2(self):
        if self.notdijet: return -999
        return self.M2g2_HadZH*self.g2ZH_m4l**2 / (self.M2g2_HadZH*self.g2ZH_m4l**2 + self.M2g2_HJJ*self.cconstantforDHadZH / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal)
    def D_HadZH_L1(self):
        if self.notdijet: return -999
        return self.M2g1prime2_HadZH*self.g1prime2ZH_m4l**2 / (self.M2g1prime2_HadZH*self.g1prime2ZH_m4l**2 + self.M2g2_HJJ*self.cconstantforDHadZH / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal)
    def D_HadZH_L1Zg(self):
        if self.notdijet: return -999
        return self.M2ghzgs1prime2_HadZH*self.ghzgs1prime2ZH_m4l**2 / (self.M2ghzgs1prime2_HadZH*self.ghzgs1prime2ZH_m4l**2 + self.M2g2_HJJ*self.cconstantforDHadZH / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal)

    @MakeJECSystematics
    def D_bkg_kin_VBFdecay(self):
        if self.notdijet: return -999
        return CJLSTscripts.D_bkg_VBFdec(
          self.p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
          self.p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
          self.p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
          self.p_JJVBF_BKG_MCFM_JECNominal,
          self.p_HadZH_BKG_MCFM_JECNominal,
          self.p_HadWH_BKG_MCFM_JECNominal,
          self.p_JJQCD_BKG_MCFM_JECNominal,
          self.p_HadZH_mavjj_JECNominal,
          self.p_HadZH_mavjj_true_JECNominal,
          self.p_HadWH_mavjj_JECNominal,
          self.p_HadWH_mavjj_true_JECNominal,
          self.pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
          self.pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
          self.pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
          self.pConst_JJVBF_BKG_MCFM_JECNominal,
          self.pConst_HadZH_BKG_MCFM_JECNominal,
          self.pConst_HadWH_BKG_MCFM_JECNominal,
          self.pConst_JJQCD_BKG_MCFM_JECNominal,
	  self.flavor,
          self.ZZMass,
        )

    @MakeJECSystematics
    def D_bkg_VBFdecay(self):
        if self.notdijet or self.p_m4l_SIG <= 0: return -999

        #result = signal / (signal+bkg) = 1 / (1+bkg/signal)
        #1/result - 1 = bkg/signal

        result = self.D_bkg_kin_VBFdecay()

        result = 1/result - 1
        result *= self.p_m4l_BKG / self.p_m4l_SIG * self.cconstantforDbkg / self.cconstantforDbkgkin
        result = 1/(1+result)

        return result

    def D_bkg_VBFdecay_ScaleUp(self):
        if self.notdijet or self.p_m4l_SIG_ScaleUp <= 0: return -999

        #result = signal / (signal+bkg) = 1 / (1+bkg/signal)
        #1/result - 1 = bkg/signal

        result = self.D_bkg_kin_VBFdecay()

        result = 1/result - 1
        result *= self.p_m4l_BKG_ScaleUp / self.p_m4l_SIG_ScaleUp * self.cconstantforDbkg / self.cconstantforDbkgkin
        result = 1/(1+result)

        return result

    def D_bkg_VBFdecay_ScaleDown(self):
        if self.notdijet or self.p_m4l_SIG_ScaleDown <= 0: return -999

        #result = signal / (signal+bkg) = 1 / (1+bkg/signal)
        #1/result - 1 = bkg/signal

        result = self.D_bkg_kin_VBFdecay()

        result = 1/result - 1
        result *= self.p_m4l_BKG_ScaleDown / self.p_m4l_SIG_ScaleDown * self.cconstantforDbkg / self.cconstantforDbkgkin
        result = 1/(1+result)

        return result

    def D_bkg_VBFdecay_ResUp(self):
        if self.notdijet or self.p_m4l_SIG_ResUp <= 0: return -999

        #result = signal / (signal+bkg) = 1 / (1+bkg/signal)
        #1/result - 1 = bkg/signal

        result = self.D_bkg_kin_VBFdecay()

        result = 1/result - 1
        result *= self.p_m4l_BKG_ResUp / self.p_m4l_SIG_ResUp * self.cconstantforDbkg / self.cconstantforDbkgkin
        result = 1/(1+result)

        return result

    def D_bkg_VBFdecay_ResDown(self):
        if self.notdijet or self.p_m4l_SIG_ResDown <= 0: return -999

        #result = signal / (signal+bkg) = 1 / (1+bkg/signal)
        #1/result - 1 = bkg/signal

        result = self.D_bkg_kin_VBFdecay()

        result = 1/result - 1
        result *= self.p_m4l_BKG_ResDown / self.p_m4l_SIG_ResDown * self.cconstantforDbkg / self.cconstantforDbkgkin
        result = 1/(1+result)

        return result

    @MakeJECSystematics
    def D_bkg_kin_HadVHdecay(self):
        if self.notdijet: return -999
        return CJLSTscripts.D_bkg_VHdec(
          self.p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
          self.p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
          self.p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
          self.p_JJVBF_BKG_MCFM_JECNominal,
          self.p_HadZH_BKG_MCFM_JECNominal,
          self.p_HadWH_BKG_MCFM_JECNominal,
          self.p_JJQCD_BKG_MCFM_JECNominal,
          self.p_HadZH_mavjj_JECNominal,
          self.p_HadZH_mavjj_true_JECNominal,
          self.p_HadWH_mavjj_JECNominal,
          self.p_HadWH_mavjj_true_JECNominal,
          self.pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
          self.pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
          self.pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
          self.pConst_JJVBF_BKG_MCFM_JECNominal,
          self.pConst_HadZH_BKG_MCFM_JECNominal,
          self.pConst_HadWH_BKG_MCFM_JECNominal,
          self.pConst_JJQCD_BKG_MCFM_JECNominal,
	  self.flavor,
          self.ZZMass,
        )

    @MakeJECSystematics
    def D_bkg_HadVHdecay(self):
        if self.notdijet or self.p_m4l_SIG <= 0: return -999

        #result = signal / (signal+bkg) = 1 / (1+bkg/signal)
        #1/result - 1 = bkg/signal

        result = self.D_bkg_kin_HadVHdecay()

        result = 1/result - 1
        result *= self.p_m4l_BKG / self.p_m4l_SIG * self.cconstantforDbkg / self.cconstantforDbkgkin
        result = 1/(1+result)

        return result

    def D_bkg_HadVHdecay_ScaleUp(self):
        if self.notdijet or self.p_m4l_SIG_ScaleUp <= 0: return -999

        #result = signal / (signal+bkg) = 1 / (1+bkg/signal)
        #1/result - 1 = bkg/signal

        result = self.D_bkg_kin_HadVHdecay()

        result = 1/result - 1
        result *= self.p_m4l_BKG_ScaleUp / self.p_m4l_SIG_ScaleUp * self.cconstantforDbkg / self.cconstantforDbkgkin
        result = 1/(1+result)

        return result

    def D_bkg_HadVHdecay_ScaleDown(self):
        if self.notdijet or self.p_m4l_SIG_ScaleDown <= 0: return -999

        #result = signal / (signal+bkg) = 1 / (1+bkg/signal)
        #1/result - 1 = bkg/signal

        result = self.D_bkg_kin_HadVHdecay()

        result = 1/result - 1
        result *= self.p_m4l_BKG_ScaleDown / self.p_m4l_SIG_ScaleDown * self.cconstantforDbkg / self.cconstantforDbkgkin
        result = 1/(1+result)

        return result

    def D_bkg_HadVHdecay_ResUp(self):
        if self.notdijet or self.p_m4l_SIG_ResUp <= 0: return -999

        #result = signal / (signal+bkg) = 1 / (1+bkg/signal)
        #1/result - 1 = bkg/signal

        result = self.D_bkg_kin_HadVHdecay()

        result = 1/result - 1
        result *= self.p_m4l_BKG_ResUp / self.p_m4l_SIG_ResUp * self.cconstantforDbkg / self.cconstantforDbkgkin
        result = 1/(1+result)

        return result

    def D_bkg_HadVHdecay_ResDown(self):
        if self.notdijet or self.p_m4l_SIG_ResDown <= 0: return -999

        #result = signal / (signal+bkg) = 1 / (1+bkg/signal)
        #1/result - 1 = bkg/signal

        result = self.D_bkg_kin_HadVHdecay()

        result = 1/result - 1
        result *= self.p_m4l_BKG_ResDown / self.p_m4l_SIG_ResDown * self.cconstantforDbkg / self.cconstantforDbkgkin
        result = 1/(1+result)

        return result

###################################
#anomalous couplings discriminants#
###################################

    def D_0minus_decay(self):
        return self.M2g1_decay / (self.M2g1_decay + self.M2g4_decay*self.g4HZZ_m4l**2)
    def D_CP_decay(self):
        return self.M2g1g4_decay*self.g4HZZ_m4l / (self.M2g1_decay + self.M2g4_decay*self.g4HZZ_m4l**2)
    def D_CP_decay_new(self):
        return self.M2g1g4_decay*self.g4HZZ_m4l / (2 * sqrt(self.M2g1_decay * self.M2g4_decay*self.g4HZZ_m4l**2))
    def D_0hplus_decay(self):
        return self.M2g1_decay / (self.M2g1_decay + self.M2g2_decay*self.g2HZZ_m4l**2)
    def D_int_decay(self):
        return self.M2g1g2_decay*self.g2HZZ_m4l / (self.M2g1_decay + self.M2g2_decay*self.g2HZZ_m4l**2)
    def D_int_decay_new(self):
        return self.M2g1g2_decay*self.g2HZZ_m4l / (2 * sqrt(self.M2g1_decay * self.M2g2_decay*self.g2HZZ_m4l**2))
    def D_L1_decay(self):
        return self.M2g1_decay / (self.M2g1_decay + self.M2g1prime2_decay*self.g1prime2HZZ_m4l**2)
    def D_L1int_decay(self):
        return self.M2g1g1prime2_decay*self.g1prime2HZZ_m4l / (self.M2g1_decay + self.M2g1prime2_decay*self.g1prime2HZZ_m4l**2)
    def D_L1int_decay_new(self):
        return self.M2g1g1prime2_decay*self.g1prime2HZZ_m4l / (2 * sqrt(self.M2g1_decay * self.M2g1prime2_decay*self.g1prime2HZZ_m4l**2))
    def D_L1Zg_decay(self):
        return self.M2g1_decay / (self.M2g1_decay + self.M2ghzgs1prime2_decay*self.ghzgs1prime2HZZ_m4l**2)
    def D_L1Zgint_decay(self):
        return self.M2g1ghzgs1prime2_decay*self.ghzgs1prime2HZZ_m4l / (self.M2g1_decay + self.M2ghzgs1prime2_decay*self.ghzgs1prime2HZZ_m4l**2)
    def D_L1Zgint_decay_new(self):
        return self.M2g1ghzgs1prime2_decay*self.ghzgs1prime2HZZ_m4l / (2 * sqrt(self.M2g1_decay * self.M2ghzgs1prime2_decay*self.ghzgs1prime2HZZ_m4l**2))
    def D_L1L1Zg_decay(self):
        return self.M2g1prime2_decay*self.g1prime2HZZ_m4l**2 / (self.M2g1prime2_decay*self.g1prime2HZZ_m4l**2 + self.M2ghzgs1prime2_decay*self.ghzgs1prime2HZZ_m4l**2)
    def D_L1L1Zgint_decay(self):
        return self.M2g1prime2ghzgs1prime2_decay*self.g1prime2HZZ_m4l*self.ghzgs1prime2HZZ_m4l / (self.M2g1prime2_decay*self.g1prime2HZZ_m4l**2 + self.M2ghzgs1prime2_decay*self.ghzgs1prime2HZZ_m4l**2)
    def D_L1L1Zgint_decay_new(self):
        return self.M2g1prime2ghzgs1prime2_decay*self.g1prime2HZZ_m4l*self.ghzgs1prime2HZZ_m4l / (2 * sqrt(self.M2g1prime2_decay*self.g1prime2HZZ_m4l**2 * self.M2ghzgs1prime2_decay*self.ghzgs1prime2HZZ_m4l**2))

############################
#contact term discriminants#
############################

    def D_eL_decay(self):
        return self.M2g1_decay / (self.M2g1_decay + self.M2eL_decay*constants.eLHZZ**2)
    def D_eLint_decay(self):
        return self.M2g1eL_decay*constants.eLHZZ / (self.M2g1_decay + self.M2eL_decay*constants.eLHZZ**2)
    def D_eLint_decay_new(self):
        return self.M2g1eL_decay*constants.eLHZZ / (2 * sqrt(self.M2g1_decay * self.M2eL_decay*constants.eLHZZ**2))
    def D_eR_decay(self):
        return self.M2g1_decay / (self.M2g1_decay + self.M2eR_decay*constants.eRHZZ**2)
    def D_eRint_decay(self):
        return self.M2g1eR_decay*constants.eRHZZ / (self.M2g1_decay + self.M2eR_decay*constants.eRHZZ**2)
    def D_eRint_decay_new(self):
        return self.M2g1eR_decay*constants.eRHZZ / (2 * sqrt(self.M2g1_decay * self.M2eR_decay*constants.eRHZZ**2))
    def D_eLeR_decay(self):
        return self.M2eL_decay*constants.eLHZZ**2 / (self.M2eR_decay*constants.eRHZZ**2 + self.M2eL_decay*constants.eLHZZ**2)
    def D_eLeRint_decay(self):
        return self.M2eLeR_decay*constants.eLHZZ*constants.eRHZZ / (self.M2eR_decay*constants.eRHZZ**2 + self.M2eL_decay*constants.eLHZZ**2)
    def D_eLeRint_decay_new(self):
        return self.M2eLeR_decay*constants.eLHZZ*constants.eRHZZ / (2 * sqrt(self.M2eR_decay*constants.eRHZZ**2 * self.M2eL_decay*constants.eLHZZ**2))

#######################################
#VBF anomalous couplings discriminants#
#######################################

    @MakeJECSystematics
    def D_0minus_VBF(self):
        if self.notdijet: return -999
        return self.M2g1_VBF / (self.M2g1_VBF + self.M2g4_VBF*self.g4VBF_m4l**2)
    @MakeJECSystematics
    def D_CP_VBF(self):
        if self.notdijet: return -999
        return self.M2g1g4_VBF*self.g4VBF_m4l / (self.M2g1_VBF + self.M2g4_VBF*self.g4VBF_m4l**2)
    @MakeJECSystematics
    def D_CP_VBF_new(self):
        if self.notdijet: return -999
        return self.M2g1g4_VBF*self.g4VBF_m4l / (2 * sqrt(self.M2g1_VBF * self.M2g4_VBF*self.g4VBF_m4l**2))
    @MakeJECSystematics
    def D_0hplus_VBF(self):
        if self.notdijet: return -999
        return self.M2g1_VBF / (self.M2g1_VBF + self.M2g2_VBF*self.g2VBF_m4l**2)
    @MakeJECSystematics
    def D_int_VBF(self):
        if self.notdijet: return -999
        return self.M2g1g2_VBF*self.g2VBF_m4l / (self.M2g1_VBF + self.M2g2_VBF*self.g2VBF_m4l**2)
    @MakeJECSystematics
    def D_int_VBF_new(self):
        if self.notdijet: return -999
        return self.M2g1g2_VBF*self.g2VBF_m4l / (2 * sqrt(self.M2g1_VBF * self.M2g2_VBF*self.g2VBF_m4l**2))
    @MakeJECSystematics
    def D_L1_VBF(self):
        if self.notdijet: return -999
        return self.M2g1_VBF / (self.M2g1_VBF + self.M2g1prime2_VBF*self.g1prime2VBF_m4l**2)
    @MakeJECSystematics
    def D_L1int_VBF(self):
        if self.notdijet: return -999
        return self.M2g1g1prime2_VBF*self.g1prime2VBF_m4l / (self.M2g1_VBF + self.M2g1prime2_VBF*self.g1prime2VBF_m4l**2)
    @MakeJECSystematics
    def D_L1int_VBF_new(self):
        if self.notdijet: return -999
        return self.M2g1g1prime2_VBF*self.g1prime2VBF_m4l / (2 * sqrt(self.M2g1_VBF * self.M2g1prime2_VBF*self.g1prime2VBF_m4l**2))
    @MakeJECSystematics
    def D_L1Zg_VBF(self):
        if self.notdijet: return -999
        return self.M2g1_VBF / (self.M2g1_VBF + self.M2ghzgs1prime2_VBF*self.ghzgs1prime2VBF_m4l**2)
    @MakeJECSystematics
    def D_L1Zgint_VBF(self):
        if self.notdijet: return -999
        return self.M2g1ghzgs1prime2_VBF*self.ghzgs1prime2VBF_m4l / (self.M2g1_VBF + self.M2ghzgs1prime2_VBF*self.ghzgs1prime2VBF_m4l**2)
    @MakeJECSystematics
    def D_L1Zgint_VBF_new(self):
        if self.notdijet: return -999
        return self.M2g1ghzgs1prime2_VBF*self.ghzgs1prime2VBF_m4l / (2 * sqrt(self.M2g1_VBF * self.M2ghzgs1prime2_VBF*self.ghzgs1prime2VBF_m4l**2))

###############################################
#VH hadronic anomalous couplings discriminants#
###############################################

    @MakeJECSystematics
    def D_0minus_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1_HadWH + self.M2g1_HadZH)
               /
                 (
                   (self.M2g1_HadWH + self.M2g1_HadZH)
                 +
                   (self.M2g4_HadWH + self.M2g4_HadZH)*self.g4VH_m4l**2
                 )
               )
    @MakeJECSystematics
    def D_CP_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1g4_HadWH + self.M2g1g4_HadZH)*self.g4VH_m4l
               /
                 (
                   (self.M2g1_HadWH + self.M2g1_HadZH)
                 +
                   (self.M2g4_HadWH + self.M2g4_HadZH)*self.g4VH_m4l**2
                 )
               )
    @MakeJECSystematics
    def D_CP_HadVH_new(self):
        if self.notdijet: return -999
        return .5 * (
                     self.M2g1g4_HadWH / (2 * sqrt(self.M2g1_HadWH * self.M2g4_HadWH))
                    +
                     self.M2g1g4_HadZH / (2 * sqrt(self.M2g1_HadZH * self.M2g4_HadZH))
                    )
    @MakeJECSystematics
    def D_0hplus_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1_HadWH + self.M2g1_HadZH)
               /
                 (
                   (self.M2g1_HadWH + self.M2g1_HadZH)
                 +
                   (self.M2g2_HadWH + self.M2g2_HadZH)*self.g2VH_m4l**2
                 )
               )
    @MakeJECSystematics
    def D_int_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1g2_HadWH + self.M2g1g2_HadZH)*self.g2VH_m4l
               /
                 (
                   (self.M2g1_HadWH + self.M2g1_HadZH)
                 +
                   (self.M2g2_HadWH + self.M2g2_HadZH)*self.g2VH_m4l**2
                 )
               )
    @MakeJECSystematics
    def D_int_HadVH_new(self):
        if self.notdijet: return -999
        return .5 * (
                     self.M2g1g2_HadWH / (2 * sqrt(self.M2g1_HadWH * self.M2g2_HadWH))
                    +
                     self.M2g1g2_HadZH / (2 * sqrt(self.M2g1_HadZH * self.M2g2_HadZH))
                    )
    @MakeJECSystematics
    def D_L1_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1_HadWH + self.M2g1_HadZH)
               /
                 (
                   (self.M2g1_HadWH + self.M2g1_HadZH)
                 +
                   (self.M2g1prime2_HadWH + self.M2g1prime2_HadZH)*self.g1prime2VH_m4l**2
                 )
               )
    @MakeJECSystematics
    def D_L1int_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1g1prime2_HadWH + self.M2g1g1prime2_HadZH)*self.g1prime2VH_m4l
               /
                 (
                   (self.M2g1_HadWH + self.M2g1_HadZH)
                 +
                   (self.M2g1prime2_HadWH + self.M2g1prime2_HadZH)*self.g1prime2VH_m4l**2
                 )
               )
    @MakeJECSystematics
    def D_L1int_HadVH_new(self):
        if self.notdijet: return -999
        return .5 * (
                     self.M2g1g1prime2_HadWH / (2 * sqrt(self.M2g1_HadWH * self.M2g1prime2_HadWH))
                    +
                     self.M2g1g1prime2_HadZH / (2 * sqrt(self.M2g1_HadZH * self.M2g1prime2_HadZH))
                    )
    @MakeJECSystematics
    def D_L1Zg_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1_HadWH + self.M2g1_HadZH)
               /
                 (
                   (self.M2g1_HadWH + self.M2g1_HadZH)
                 +
                   (self.M2ghzgs1prime2_HadWH + self.M2ghzgs1prime2_HadZH)
                          * self.ghzgs1prime2VH_m4l**2
                 )
               )
    @MakeJECSystematics
    def D_L1Zgint_HadVH(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1ghzgs1prime2_HadWH + self.M2g1ghzgs1prime2_HadZH)*self.ghzgs1prime2VH_m4l
               /
                 (
                   (self.M2g1_HadWH + self.M2g1_HadZH)
                 +
                   (self.M2ghzgs1prime2_HadWH + self.M2ghzgs1prime2_HadZH)*self.ghzgs1prime2VH_m4l**2
                 )
               )
    @MakeJECSystematics
    def D_L1Zgint_HadVH_new(self):
        if self.notdijet: return -999
        return self.M2g1ghzgs1prime2_HadZH / (2 * sqrt(self.M2g1_HadZH * self.M2ghzgs1prime2_HadZH))
        """
        return .5 * (
                     self.M2g1ghzgs1prime2_HadWH / (2 * sqrt(self.M2g1_HadWH * self.M2ghzgs1prime2_HadWH))
                    +
                     self.M2g1ghzgs1prime2_HadZH / (2 * sqrt(self.M2g1_HadZH * self.M2ghzgs1prime2_HadZH))
                    )
        """

############################################
#VBFdecay anomalous couplings discriminants#
############################################

    @MakeJECSystematics
    def D_0minus_VBFdecay(self):
        if self.notdijet: return -999
        return self.M2g1_VBF*self.M2g1_decay / (self.M2g1_VBF*self.M2g1_decay + self.M2g4_VBF*self.M2g4_decay*(self.g4VBF_m4l*self.g4HZZ_m4l)**2)
    @MakeJECSystematics
    def D_0hplus_VBFdecay(self):
        if self.notdijet: return -999
        return self.M2g1_VBF*self.M2g1_decay / (self.M2g1_VBF*self.M2g1_decay + self.M2g2_VBF*self.M2g2_decay * (self.g2VBF_m4l*self.g2HZZ_m4l)**2)
    @MakeJECSystematics
    def D_L1_VBFdecay(self):
        if self.notdijet: return -999
        return self.M2g1_VBF*self.M2g1_decay / (self.M2g1_VBF*self.M2g1_decay + self.M2g1prime2_VBF*self.M2g1prime2_decay * (self.g1prime2VBF_m4l*self.g1prime2HZZ_m4l)**2)
    @MakeJECSystematics
    def D_L1Zg_VBFdecay(self):
        if self.notdijet: return -999
        return self.M2g1_VBF*self.M2g1_decay / (self.M2g1_VBF*self.M2g1_decay + self.M2ghzgs1prime2_VBF*self.M2ghzgs1prime2_decay * (self.ghzgs1prime2VBF_m4l*self.ghzgs1prime2HZZ_m4l)**2)

####################################################
#VHdecay hadronic anomalous couplings discriminants#
####################################################

    @MakeJECSystematics
    def D_0minus_HadVHdecay(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1_HadWH + self.M2g1_HadZH)*self.M2g1_decay
               /
                 (
                   (self.M2g1_HadWH + self.M2g1_HadZH)
                        *self.M2g1_decay
                 + (self.M2g4_HadWH + self.M2g4_HadZH)*self.g4VH_m4l**2
                        *self.M2g4_decay*self.g4HZZ_m4l**2
                 )
               )
    @MakeJECSystematics
    def D_0hplus_HadVHdecay(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1_HadWH + self.M2g1_HadZH)*self.M2g1_decay
               /
                 (
                   (self.M2g1_HadWH + self.M2g1_HadZH)
                        *self.M2g1_decay
                 + (self.M2g2_HadWH + self.M2g2_HadZH)*self.g2VH_m4l**2
                        *self.M2g2_decay*self.g2HZZ_m4l**2
                 )
               )
    @MakeJECSystematics
    def D_L1_HadVHdecay(self):
        if self.notdijet: return -999
        return (
                 (self.M2g1_HadWH + self.M2g1_HadZH)*self.M2g1_decay
               /
                 (
                   (self.M2g1_HadWH + self.M2g1_HadZH)
                        *self.M2g1_decay
                 + (self.M2g1prime2_HadWH + self.M2g1prime2_HadZH)*self.g1prime2VH_m4l**2
                        *self.M2g1prime2_decay*self.g1prime2HZZ_m4l**2
                 )
               )
    @MakeJECSystematics
    def D_L1Zg_HadVHdecay(self):
        if self.notdijet: return -999
        return (
                 ((self.M2g1_HadWH + self.M2g1_HadZH)
                    *self.M2g1_decay)
               /
                 (
                   (self.M2g1_HadWH + self.M2g1_HadZH)
                        *self.M2g1_decay
                 + (self.M2ghzgs1prime2_HadWH + self.M2ghzgs1prime2_HadZH)*self.ghzgs1prime2VH_m4l**2
                        *self.M2ghzgs1prime2_decay*self.ghzgs1prime2HZZ_m4l**2
                 )
               )

######
#STXS#
######

    @MakeBtagSystematics
    @MakeJECSystematics
    def D_STXS_stage1p1(self):
        category = CJLSTscripts.categoryMor18(
          self.nExtraLep,
          self.nExtraZ,
          self.nCleanedJetsPt30,
          self.nCleanedJetsPt30BTagged_bTagSF,
          self.jetQGLikelihood,
          self.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
          self.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
          self.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
          self.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
          self.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
          self.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
          self.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
          self.p_HadWH_mavjj_JECNominal,
          self.p_HadWH_mavjj_true_JECNominal,
          self.p_HadZH_mavjj_JECNominal,
          self.p_HadZH_mavjj_true_JECNominal,
          self.jetPhi,
          self.ZZMass,
          self.PFMET,
          config.useVHMETTagged,
          config.useQGTagging,
        )
        return STXS.stage1_reco_1p1(self.nCleanedJetsPt30, self.DiJetMass, self.ZZPt, category, self.HjjPt)

#########################
#4 coupling discriminant#
#########################
    def D_4couplings_general_raw(self, *variables_and_bins):
        """
        variables_and_bins is something like:
          ("D_0minus_decay", [.333, .667]), ("D_CP_decay", [0]), ...
          do not include either endpoint in the binning.
        """
        result = 0
        for variablename, binning in variables_and_bins:
          result *= len(binning)+1
          variable = getattr(self, variablename)()
          for bin in binning:
            if variable > bin:
              result += 1
        return result

    def D_4couplings_general(self, variables_and_bins, foldbins):
      result = self.D_4couplings_general_raw(*variables_and_bins)
      if result in foldbins: result = product(len(bins)+1 for variable, bins in variables_and_bins)
      for foldbin in foldbins:
          if result > foldbin:
              result -= 1
      return result

    #here we specify the discriminants and bin separations.
    #Note D_CP is NOT here.  It's used as the third dimension
    #so that mirroring is easier.
    binning_4couplings_decay = (
      ("D_0minus_decay", [.333, .667]),
      ("D_0hplus_decay", [.5, .7]),
      ("D_L1_decay", [.55, .8]),
      ("D_L1Zg_decay", [.4, .55]),
      ("D_int_decay", [.8]),
    )
    foldbins_4couplings_decay = tuple(sorted((
      #bins with very few entries, summed over all production modes
      #see misc/foldbins
      15, 17, 22, 33, 34, 35, 36, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 68, 69, 70, 71, 89, 101, 103, 105, 106, 107, 114, 116, 118, 120, 121, 122, 123, 124, 125, 140, 141, 143, 155, 157, 159, 160, 161,
    ), reverse=True))
    assert 162 - len(foldbins_4couplings_decay) == discriminants.discriminant("D_4couplings_decay").bins == discriminants.discriminant("D_4couplings_decay").max, (162 - len(foldbins_4couplings_decay), discriminants.discriminant("D_4couplings_decay").bins)

    binning_4couplings_VBFdecay = (
      ("D_0minus_VBFdecay", [.1, .9]),
      ("D_0hplus_VBFdecay", [.1, .9]),
      ("D_L1_VBFdecay", [.1, .9]),
      ("D_L1Zg_VBFdecay", [.1, .8]),
      ("D_int_VBF", [0]),
    )
    binning_4couplings_VBFdecay_JECUp = tuple((name+"_JECUp", bins) for name, bins in binning_4couplings_VBFdecay)
    binning_4couplings_VBFdecay_JECDn = tuple((name+"_JECDn", bins) for name, bins in binning_4couplings_VBFdecay)

    foldbins_4couplings_VBFdecay = tuple(sorted((
      #bins with very few entries, summed over all production modes
      #see misc/foldbins
      4, 5, 12, 13, 14, 15, 17, 22, 23, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 58, 59, 66, 67, 68, 69, 70, 71, 76, 77, 90, 91, 92, 93, 94, 95, 102, 103, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 138, 139, 144, 145, 146, 147, 148, 149, 156, 157,
    ), reverse=True))
    assert 162 - len(foldbins_4couplings_VBFdecay) == discriminants.discriminant("D_4couplings_VBFdecay").bins == discriminants.discriminant("D_4couplings_VBFdecay").max, (162 - len(foldbins_4couplings_VBFdecay), discriminants.discriminant("D_4couplings_VBFdecay").bins)

    binning_4couplings_HadVHdecay = (
      ("D_0minus_HadVHdecay", [.2, .8]),
      ("D_0hplus_HadVHdecay", [.333, .667]),
      ("D_L1_HadVHdecay", [.333, .667]),
      ("D_L1Zg_HadVHdecay", [.1, .9]),
      ("D_int_HadVH", [-.6]),
    )
    binning_4couplings_HadVHdecay_JECUp = tuple((name+"_JECUp", bins) for name, bins in binning_4couplings_HadVHdecay)
    binning_4couplings_HadVHdecay_JECDn = tuple((name+"_JECDn", bins) for name, bins in binning_4couplings_HadVHdecay)

    foldbins_4couplings_HadVHdecay = tuple(sorted((
      #bins with very few entries, summed over all production modes
      #see misc/foldbins
      0, 4, 6, 18, 22, 23, 24, 25, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 58, 72, 76, 77, 78, 90, 91, 92, 94, 95, 96, 97, 100, 108, 110, 112, 113, 114, 120, 121, 125, 126, 127, 128, 130, 131, 132, 133, 137, 144, 145, 146, 147, 148, 149, 150, 151, 154, 155,
    ), reverse=True))
    assert 162 - len(foldbins_4couplings_HadVHdecay) == discriminants.discriminant("D_4couplings_HadVHdecay").bins == discriminants.discriminant("D_4couplings_HadVHdecay").max, (162 - len(foldbins_4couplings_HadVHdecay), discriminants.discriminant("D_4couplings_HadVHdecay").bins)

    def D_4couplings_decay_raw(self):
      return self.D_4couplings_general_raw(*self.binning_4couplings_decay)
    @MakeJECSystematics
    def D_4couplings_VBFdecay_raw(self):
      if self.notdijet: return -999
      return self.D_4couplings_general_raw(*self.binning_4couplings_VBFdecay)
    @MakeJECSystematics
    def D_4couplings_HadVHdecay_raw(self):
      if self.notdijet: return -999
      return self.D_4couplings_general_raw(*self.binning_4couplings_HadVHdecay)

    def D_4couplings_decay(self):
      return self.D_4couplings_general(self.binning_4couplings_decay, self.foldbins_4couplings_decay)
    @MakeJECSystematics
    def D_4couplings_VBFdecay(self):
      if self.notdijet: return -999
      return self.D_4couplings_general(self.binning_4couplings_VBFdecay, self.foldbins_4couplings_VBFdecay)
    @MakeJECSystematics
    def D_4couplings_HadVHdecay(self):
      if self.notdijet: return -999
      return self.D_4couplings_general(self.binning_4couplings_HadVHdecay, self.foldbins_4couplings_HadVHdecay)

class FixWrongWH(object):
  def __init__(self, tree):
    self.tree = tree
  def __getattr__(self, attr):
    try:
      return getattr(self.tree, attr)
    except AttributeError:
      return getattr(self.tree, attr.replace("ghw", "ghv"))

@callclassinitfunctions("initcategoryfunctions", "initsystematics")
class TreeWrapper(TreeWrapperBase):

    definitelyexists = Sample("ggH", "0+", config.productionsforcombine[0])
    if not xrd.exists(definitelyexists.CJLSTfile()):
        raise ValueError("{} does not exist!".format(definitelyexists.CJLSTfile()))

    def __init__(self, treesample, minevent=0, maxevent=None, LSF=Fake_LSF_creating()):
        """
        treesample - which sample the TTree was created from
        """
        self.treesample = treesample

        filename = treesample.CJLSTfile()
        if xrd.exists(filename):
            filename = LSF.basename(treesample.CJLSTfile())
            self.f = TFile(filename, contextmanager=False)

        self.filename = filename

        if self.isdummy:
            print "{} does not exist or is bad, using {}".format(filename, self.definitelyexists.CJLSTfile())
            #give it a tree so that it can get the format, but not fill any entries
            filename = self.definitelyexists.CJLSTfile()
        self.f = TFile(filename, contextmanager=False)

        self.counters = None

        if not self.isdata and not self.isZX:
            self.counters = self.f.Get("{}/Counters".format(treesample.TDirectoryname()))
            if not self.counters:
                raise ValueError("No Counters in file "+filename)

        self.tree = self.f.Get("{}/candTree".format(treesample.TDirectoryname()))
        self.failedtree = self.f.Get("{}/candTree_failed".format(treesample.TDirectoryname()))
        if not self.failedtree: self.failedtree = None

        self.nevents = self.nevents2e2mu = self.cutoffs = None
        self.xsec = None

        super(TreeWrapper, self).__init__(treesample, minevent, maxevent)

        if self.counters is not None:
            self.nevents = self.counters.GetBinContent(40)
            #========================
            if "GEN" in str(self.treesample.production): assert self.treesample.production == "GEN_181119", "remove this section"
            if self.nevents == 0:
                self.nevents = self.counters.GetBinContent(41)
            assert self.nevents != 0
            #========================

        self.tree.GetEntry(0)
        if not self.isdata and not self.isZX:
            self.xsec = self.tree.xsec * 1000 #pb to fb
            self.genxsec = self.tree.genxsec
            self.genBR = self.tree.genBR

    @property
    @cache_instancemethod
    def isdummy(self):
        if not xrd.exists(self.filename):
            return True
        else:
            if not self.f.Get("{}/candTree".format(self.treesample.TDirectoryname())):
                return True
        return super(TreeWrapper, self).isdummy

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
            if i % self.printevery == 0 or i == len(self):
                print i, "/", len(self)
                #raise StopIteration

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
        self.ZZPt = t.ZZPt
        self.ZZEta = t.ZZEta
        Hvector = tlvfromptetaphim(self.ZZPt, self.ZZEta, t.ZZPhi, self.ZZMass)

        self.cconstantforDbkgkin = CJLSTscripts.getDbkgkinConstant(self.flavor, self.ZZMass)
        self.cconstantforDbkg = CJLSTscripts.getDbkgConstant(self.flavor, self.ZZMass)
        self.cconstantforD2jet = CJLSTscripts.getDVBF2jetsConstant_shiftWP(self.ZZMass, config.useQGTagging, 0.5)
        self.cconstantforDHadWH = CJLSTscripts.getDWHhConstant_shiftWP(self.ZZMass, config.useQGTagging, 0.5)
        self.cconstantforDHadZH = CJLSTscripts.getDZHhConstant_shiftWP(self.ZZMass, config.useQGTagging, 0.5)

        self.p_m4l_BKG = t.p_m4l_BKG
        self.p_m4l_SIG = t.p_m4l_SIG

        for a in "SIG", "BKG":
            for b in "Scale", "Res":
                for c in "Up", "Down":
                    attr = "p_m4l_{}_{}{}".format(a, b, c)
                    setattr(self, attr, getattr(t, attr))

        #express in terms of |M|^2, this will make life easier
        self.M2qqZZ = t.p_QQB_BKG_MCFM

        self.M2g1_decay                   = t.p_GG_SIG_ghg2_1_ghz1_1_JHUGen
        self.M2g4_decay                   = t.p_GG_SIG_ghg2_1_ghz4_1_JHUGen
        self.M2g1g4_decay                 = t.p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen
        self.M2g2_decay                   = t.p_GG_SIG_ghg2_1_ghz2_1_JHUGen
        self.M2g1g2_decay                 = t.p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen
        self.M2g1prime2_decay             = t.p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen / 1e4**2
        self.M2g1g1prime2_decay           = t.p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen / 1e4
        self.M2ghzgs1prime2_decay         = t.p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen / 1e4**2
        self.M2g1ghzgs1prime2_decay       = t.p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen / 1e4
        self.M2g1prime2ghzgs1prime2_decay = t.p_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen / 1e4**2

        self.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal               = self.M2g1_VBF                     = t.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal
        self.p_JJVBF_SIG_ghv4_1_JHUGen_JECNominal               = self.M2g4_VBF                     = t.p_JJVBF_SIG_ghv4_1_JHUGen_JECNominal
        self.p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECNominal        = self.M2g1g4_VBF                   = t.p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECNominal
        self.p_JJVBF_SIG_ghv2_1_JHUGen_JECNominal               = self.M2g2_VBF                     = t.p_JJVBF_SIG_ghv2_1_JHUGen_JECNominal
        self.p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECNominal        = self.M2g1g2_VBF                   = t.p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECNominal
        self.p_JJVBF_SIG_ghv1prime2_1_JHUGen_JECNominal         = self.M2g1prime2_VBF               = t.p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal / 1e4**2
        self.p_JJVBF_SIG_ghv1_1_ghv1prime2_1_JHUGen_JECNominal  = self.M2g1g1prime2_VBF             = t.p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_JECNominal / 1e4
        self.p_JJVBF_SIG_ghza1prime2_1_JHUGen_JECNominal        = self.M2ghzgs1prime2_VBF           = t.p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_JECNominal / 1e4**2
        self.p_JJVBF_SIG_ghv1_1_ghza1prime2_1_JHUGen_JECNominal = self.M2g1ghzgs1prime2_VBF         = t.p_JJVBF_SIG_ghv1_1_ghza1prime2_1E4_JHUGen_JECNominal / 1e4

        self.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal               = self.M2g2_HJJ                     = t.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal
        self.p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal                                                   = t.p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal
        self.p_JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_JECNominal                                            = t.p_JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_JECNominal

        self.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal               = self.M2g1_HadZH                   = t.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal
        self.p_HadZH_SIG_ghz4_1_JHUGen_JECNominal               = self.M2g4_HadZH                   = t.p_HadZH_SIG_ghz4_1_JHUGen_JECNominal
        self.p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECNominal        = self.M2g1g4_HadZH                 = t.p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECNominal
        self.p_HadZH_SIG_ghz2_1_JHUGen_JECNominal               = self.M2g2_HadZH                   = t.p_HadZH_SIG_ghz2_1_JHUGen_JECNominal
        self.p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECNominal        = self.M2g1g2_HadZH                 = t.p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECNominal
        self.p_HadZH_SIG_ghz1prime2_1_JHUGen_JECNominal         = self.M2g1prime2_HadZH             = t.p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECNominal / 1e4**2
        self.p_HadZH_SIG_ghz1_1_ghz1prime2_1_JHUGen_JECNominal  = self.M2g1g1prime2_HadZH           = t.p_HadZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen_JECNominal / 1e4
        self.p_HadZH_SIG_ghza1prime2_1_JHUGen_JECNominal        = self.M2ghzgs1prime2_HadZH         = t.p_HadZH_SIG_ghza1prime2_1E4_JHUGen_JECNominal / 1e4**2
        self.p_HadZH_SIG_ghz1_1_ghza1prime2_1_JHUGen_JECNominal = self.M2g1ghzgs1prime2_HadZH       = t.p_HadZH_SIG_ghz1_1_ghza1prime2_1E4_JHUGen_JECNominal / 1e4

        self.p_HadZH_mavjj_JECNominal                                                               = t.p_HadZH_mavjj_JECNominal
        self.p_HadZH_mavjj_true_JECNominal                                                          = t.p_HadZH_mavjj_true_JECNominal
        self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal                                              = t.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal

        self.M2g1_HadZH                   *= self.p_HadZH_mavjj_JECNominal / self.p_HadZH_mavjj_true_JECNominal / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal
        self.M2g4_HadZH                   *= self.p_HadZH_mavjj_JECNominal / self.p_HadZH_mavjj_true_JECNominal / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal
        self.M2g1g4_HadZH                 *= self.p_HadZH_mavjj_JECNominal / self.p_HadZH_mavjj_true_JECNominal / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal
        self.M2g2_HadZH                   *= self.p_HadZH_mavjj_JECNominal / self.p_HadZH_mavjj_true_JECNominal / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal
        self.M2g1g2_HadZH                 *= self.p_HadZH_mavjj_JECNominal / self.p_HadZH_mavjj_true_JECNominal / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal
        self.M2g1prime2_HadZH             *= self.p_HadZH_mavjj_JECNominal / self.p_HadZH_mavjj_true_JECNominal / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal
        self.M2g1g1prime2_HadZH           *= self.p_HadZH_mavjj_JECNominal / self.p_HadZH_mavjj_true_JECNominal / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal
        self.M2ghzgs1prime2_HadZH         *= self.p_HadZH_mavjj_JECNominal / self.p_HadZH_mavjj_true_JECNominal / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal
        self.M2g1ghzgs1prime2_HadZH       *= self.p_HadZH_mavjj_JECNominal / self.p_HadZH_mavjj_true_JECNominal / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal

        self.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal               = self.M2g1_HadWH                   = t.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal
        self.p_HadWH_SIG_ghw4_1_JHUGen_JECNominal               = self.M2g4_HadWH                   = t.p_HadWH_SIG_ghw4_1_JHUGen_JECNominal
        self.p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECNominal        = self.M2g1g4_HadWH                 = t.p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECNominal
        self.p_HadWH_SIG_ghw2_1_JHUGen_JECNominal               = self.M2g2_HadWH                   = t.p_HadWH_SIG_ghw2_1_JHUGen_JECNominal
        self.p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECNominal        = self.M2g1g2_HadWH                 = t.p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECNominal
        self.p_HadWH_SIG_ghw1prime2_1_JHUGen_JECNominal         = self.M2g1prime2_HadWH             = t.p_HadWH_SIG_ghw1prime2_1E4_JHUGen_JECNominal / 1e4**2
        self.p_HadWH_SIG_ghw1_1_ghw1prime2_1_JHUGen_JECNominal  = self.M2g1g1prime2_HadWH           = t.p_HadWH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen_JECNominal / 1e4
        self.p_HadWH_SIG_ghza1prime2_1_JHUGen_JECNominal        = self.M2ghzgs1prime2_HadWH         = 0
        self.p_HadWH_SIG_ghz1_1_ghza1prime2_1_JHUGen_JECNominal = self.M2g1ghzgs1prime2_HadWH       = 0

        self.p_HadWH_mavjj_JECNominal                                                               = t.p_HadWH_mavjj_JECNominal
        self.p_HadWH_mavjj_true_JECNominal                                                          = t.p_HadWH_mavjj_true_JECNominal
        self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal                                              = t.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal

        self.M2g1_HadWH                   *= self.p_HadWH_mavjj_JECNominal / self.p_HadWH_mavjj_true_JECNominal / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal
        self.M2g4_HadWH                   *= self.p_HadWH_mavjj_JECNominal / self.p_HadWH_mavjj_true_JECNominal / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal
        self.M2g1g4_HadWH                 *= self.p_HadWH_mavjj_JECNominal / self.p_HadWH_mavjj_true_JECNominal / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal
        self.M2g2_HadWH                   *= self.p_HadWH_mavjj_JECNominal / self.p_HadWH_mavjj_true_JECNominal / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal
        self.M2g1g2_HadWH                 *= self.p_HadWH_mavjj_JECNominal / self.p_HadWH_mavjj_true_JECNominal / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal
        self.M2g1prime2_HadWH             *= self.p_HadWH_mavjj_JECNominal / self.p_HadWH_mavjj_true_JECNominal / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal
        self.M2g1g1prime2_HadWH           *= self.p_HadWH_mavjj_JECNominal / self.p_HadWH_mavjj_true_JECNominal / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal
        self.M2ghzgs1prime2_HadWH         *= self.p_HadWH_mavjj_JECNominal / self.p_HadWH_mavjj_true_JECNominal / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal
        self.M2g1ghzgs1prime2_HadWH       *= self.p_HadWH_mavjj_JECNominal / self.p_HadWH_mavjj_true_JECNominal / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal

        self.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal                                                    = t.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal
        self.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal                                                    = t.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal
        self.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal                                                 = t.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal

        #for Dbkgs
        self.p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal = t.p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal
        self.p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal = t.p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal
        self.p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal = t.p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal
        self.p_JJVBF_BKG_MCFM_JECNominal = t.p_JJVBF_BKG_MCFM_JECNominal
        self.p_HadZH_BKG_MCFM_JECNominal = t.p_HadZH_BKG_MCFM_JECNominal
        self.p_HadWH_BKG_MCFM_JECNominal = t.p_HadWH_BKG_MCFM_JECNominal
        self.p_JJQCD_BKG_MCFM_JECNominal = t.p_JJQCD_BKG_MCFM_JECNominal
        self.p_HadZH_mavjj_JECNominal = t.p_HadZH_mavjj_JECNominal
        self.p_HadZH_mavjj_true_JECNominal = t.p_HadZH_mavjj_true_JECNominal
        self.p_HadWH_mavjj_JECNominal = t.p_HadWH_mavjj_JECNominal
        self.p_HadWH_mavjj_true_JECNominal = t.p_HadWH_mavjj_true_JECNominal
        self.pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal = t.pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal
        self.pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal = t.pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal
        self.pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal = t.pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal
        self.pConst_JJVBF_BKG_MCFM_JECNominal = t.pConst_JJVBF_BKG_MCFM_JECNominal
        self.pConst_HadZH_BKG_MCFM_JECNominal = t.pConst_HadZH_BKG_MCFM_JECNominal
        self.pConst_HadWH_BKG_MCFM_JECNominal = t.pConst_HadWH_BKG_MCFM_JECNominal
        self.pConst_JJQCD_BKG_MCFM_JECNominal = t.pConst_JJQCD_BKG_MCFM_JECNominal

        self.g2HZZ_m4l           = gconstant("HZZ2e2mu", "g2", self.ZZMass)
        self.g4HZZ_m4l           = gconstant("HZZ2e2mu", "g4", self.ZZMass)
        self.g1prime2HZZ_m4l     = gconstant("HZZ2e2mu", "L1", self.ZZMass)
        self.ghzgs1prime2HZZ_m4l = gconstant("HZZ2e2mu", "L1Zg", self.ZZMass)

        self.g2VBF_m4l           = gconstant("VBF",      "g2", self.ZZMass)
        self.g4VBF_m4l           = gconstant("VBF",      "g4", self.ZZMass)
        self.g1prime2VBF_m4l     = gconstant("VBF",      "L1", self.ZZMass)
        self.ghzgs1prime2VBF_m4l = gconstant("VBF",      "L1Zg", self.ZZMass)

        self.g2VH_m4l            = gconstant("VH",       "g2", self.ZZMass)
        self.g4VH_m4l            = gconstant("VH",       "g4", self.ZZMass)
        self.g1prime2VH_m4l      = gconstant("VH",       "L1", self.ZZMass)
        self.ghzgs1prime2VH_m4l  = gconstant("VH",       "L1Zg", self.ZZMass)

        self.g2ZH_m4l            = gconstant("ZH",       "g2", self.ZZMass)
        self.g4ZH_m4l            = gconstant("ZH",       "g4", self.ZZMass)
        self.g1prime2ZH_m4l      = gconstant("ZH",       "L1", self.ZZMass)
        self.ghzgs1prime2ZH_m4l  = gconstant("ZH",       "L1Zg", self.ZZMass)

        self.g2WH_m4l            = gconstant("WH",       "g2", self.ZZMass)
        self.g4WH_m4l            = gconstant("WH",       "g4", self.ZZMass)
        self.g1prime2WH_m4l      = gconstant("WH",       "L1", self.ZZMass)
        #self.ghzgs1prime2WH_m4l  = gconstant("WH",       "L1Zg", self.ZZMass)

        #k factors
        for kfactor in self.kfactors:
            setattr(self, kfactor, getattr(t, kfactor))

        #category variables
        self.nExtraLep = t.nExtraLep
        self.nExtraZ = t.nExtraZ
        self.PFMET = t.PFMET

        nCleanedJets = t.nCleanedJets

        self.nCleanedJetsPt30 = t.nCleanedJetsPt30
        if not self.GEN:
            self.nCleanedJetsPt30BTagged_bTagSF = t.nCleanedJetsPt30BTagged_bTagSF
            self.nCleanedJetsPt30BTagged_bTagSFUp = t.nCleanedJetsPt30BTagged_bTagSFUp
            self.nCleanedJetsPt30BTagged_bTagSFDn = t.nCleanedJetsPt30BTagged_bTagSFDn
            self.jetQGLikelihood = t.JetQGLikelihood.data()
        else:
            assert not hasattr(t, "nCleanedJetsPt30BTagged_bTagSF")
            self.nCleanedJetsPt30BTagged_bTagSF = 0
            self.jetQGLikelihood = dummyfloatstar
        self.jetPhi = t.JetPhi.data()
        jets = [tlvfromptetaphim(pt, eta, phi, m) for pt, eta, phi, m in zip(t.JetPt, t.JetEta, t.JetPhi, t.JetMass)[:2]]

        if nCleanedJets == 0:
            self.jetQGLikelihood = self.jetPhi = dummyfloatstar
            self.Jet1Pt = -999
        else:
            self.Jet1Pt = t.JetPt[0]

        self.DiJetMass = t.DiJetMass
        self.DiJetDEta = t.DiJetDEta
        self.HjjPt = (Hvector + jets[0] + jets[1]).Pt() if nCleanedJets >= 2 else -99

        #LHE info, needed for photon cut
        try:
            self.LHEDaughterId = t.LHEDaughterId
            self.LHEDaughterPt = t.LHEDaughterPt
            self.LHEDaughterEta = t.LHEDaughterEta
            self.LHEDaughterPhi = t.LHEDaughterPhi
            self.LHEDaughterMass = t.LHEDaughterMass
            self.LHEAssociatedParticlePt = t.LHEAssociatedParticlePt
            self.LHEAssociatedParticleEta = t.LHEAssociatedParticleEta
            self.LHEAssociatedParticlePhi = t.LHEAssociatedParticlePhi
            self.LHEAssociatedParticleMass = t.LHEAssociatedParticleMass
        except AttributeError:
            pass

        if self.isdata and not self.unblind and not self.passesblindcut():
            return next(self)

        self.notdijet = self.M2g1_VBF <= 0

        if self.GEN: return self

        #JECUp
        self.p_JJVBF_SIG_ghv1_1_JHUGen_JECUp                    = self.M2g1_VBF_JECUp               = t.p_JJVBF_SIG_ghv1_1_JHUGen_JECUp
        self.p_JJVBF_SIG_ghv4_1_JHUGen_JECUp                    = self.M2g4_VBF_JECUp               = t.p_JJVBF_SIG_ghv4_1_JHUGen_JECUp
        self.p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECUp             = self.M2g1g4_VBF_JECUp             = t.p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECUp
        self.p_JJVBF_SIG_ghv2_1_JHUGen_JECUp                    = self.M2g2_VBF_JECUp               = t.p_JJVBF_SIG_ghv2_1_JHUGen_JECUp
        self.p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECUp             = self.M2g1g2_VBF_JECUp             = t.p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECUp
        self.p_JJVBF_SIG_ghv1prime2_1_JHUGen_JECUp              = self.M2g1prime2_VBF_JECUp         = t.p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECUp / 1e4**2
        self.p_JJVBF_SIG_ghv1_1_ghv1prime2_1_JHUGen_JECUp       = self.M2g1g1prime2_VBF_JECUp       = t.p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_JECUp / 1e4
        self.p_JJVBF_SIG_ghza1prime2_1_JHUGen_JECUp             = self.M2ghzgs1prime2_VBF_JECUp     = t.p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_JECUp / 1e4**2
        self.p_JJVBF_SIG_ghv1_1_ghza1prime2_1_JHUGen_JECUp      = self.M2g1ghzgs1prime2_VBF_JECUp   = t.p_JJVBF_SIG_ghv1_1_ghza1prime2_1E4_JHUGen_JECUp / 1e4

        self.p_JJQCD_SIG_ghg2_1_JHUGen_JECUp                    = self.M2g2_HJJ_JECUp               = t.p_JJQCD_SIG_ghg2_1_JHUGen_JECUp
        self.p_JJQCD_SIG_ghg4_1_JHUGen_JECUp                                                        = t.p_JJQCD_SIG_ghg4_1_JHUGen_JECUp
        self.p_JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_JECUp                                                 = t.p_JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_JECUp

        self.p_HadZH_SIG_ghz1_1_JHUGen_JECUp                    = self.M2g1_HadZH_JECUp             = t.p_HadZH_SIG_ghz1_1_JHUGen_JECUp
        self.p_HadZH_SIG_ghz4_1_JHUGen_JECUp                    = self.M2g4_HadZH_JECUp             = t.p_HadZH_SIG_ghz4_1_JHUGen_JECUp
        self.p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECUp             = self.M2g1g4_HadZH_JECUp           = t.p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECUp
        self.p_HadZH_SIG_ghz2_1_JHUGen_JECUp                    = self.M2g2_HadZH_JECUp             = t.p_HadZH_SIG_ghz2_1_JHUGen_JECUp
        self.p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECUp             = self.M2g1g2_HadZH_JECUp           = t.p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECUp
        self.p_HadZH_SIG_ghz1prime2_1_JHUGen_JECUp              = self.M2g1prime2_HadZH_JECUp       = t.p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECUp / 1e4**2
        self.p_HadZH_SIG_ghz1_1_ghz1prime2_1_JHUGen_JECUp       = self.M2g1g1prime2_HadZH_JECUp     = t.p_HadZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen_JECUp / 1e4
        self.p_HadZH_SIG_ghza1prime2_1_JHUGen_JECUp             = self.M2ghzgs1prime2_HadZH_JECUp   = t.p_HadZH_SIG_ghza1prime2_1E4_JHUGen_JECUp / 1e4**2
        self.p_HadZH_SIG_ghz1_1_ghza1prime2_1_JHUGen_JECUp      = self.M2g1ghzgs1prime2_HadZH_JECUp = t.p_HadZH_SIG_ghz1_1_ghza1prime2_1E4_JHUGen_JECUp / 1e4

        self.p_HadZH_mavjj_JECUp                                                                    = t.p_HadZH_mavjj_JECUp
        self.p_HadZH_mavjj_true_JECUp                                                               = t.p_HadZH_mavjj_true_JECUp
        self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECUp                                                   = t.pConst_HadZH_SIG_ghz1_1_JHUGen_JECUp

        self.M2g1_HadZH_JECUp             *= self.p_HadZH_mavjj_JECUp / self.p_HadZH_mavjj_true_JECUp / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECUp
        self.M2g4_HadZH_JECUp             *= self.p_HadZH_mavjj_JECUp / self.p_HadZH_mavjj_true_JECUp / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECUp
        self.M2g1g4_HadZH_JECUp           *= self.p_HadZH_mavjj_JECUp / self.p_HadZH_mavjj_true_JECUp / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECUp
        self.M2g2_HadZH_JECUp             *= self.p_HadZH_mavjj_JECUp / self.p_HadZH_mavjj_true_JECUp / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECUp
        self.M2g1g2_HadZH_JECUp           *= self.p_HadZH_mavjj_JECUp / self.p_HadZH_mavjj_true_JECUp / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECUp
        self.M2g1prime2_HadZH_JECUp       *= self.p_HadZH_mavjj_JECUp / self.p_HadZH_mavjj_true_JECUp / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECUp
        self.M2g1g1prime2_HadZH_JECUp     *= self.p_HadZH_mavjj_JECUp / self.p_HadZH_mavjj_true_JECUp / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECUp
        self.M2ghzgs1prime2_HadZH_JECUp   *= self.p_HadZH_mavjj_JECUp / self.p_HadZH_mavjj_true_JECUp / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECUp
        self.M2g1ghzgs1prime2_HadZH_JECUp *= self.p_HadZH_mavjj_JECUp / self.p_HadZH_mavjj_true_JECUp / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECUp

        self.p_HadWH_SIG_ghw1_1_JHUGen_JECUp                    = self.M2g1_HadWH_JECUp             = t.p_HadWH_SIG_ghw1_1_JHUGen_JECUp
        self.p_HadWH_SIG_ghw4_1_JHUGen_JECUp                    = self.M2g4_HadWH_JECUp             = t.p_HadWH_SIG_ghw4_1_JHUGen_JECUp
        self.p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECUp             = self.M2g1g4_HadWH_JECUp           = t.p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECUp
        self.p_HadWH_SIG_ghw2_1_JHUGen_JECUp                    = self.M2g2_HadWH_JECUp             = t.p_HadWH_SIG_ghw2_1_JHUGen_JECUp
        self.p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECUp             = self.M2g1g2_HadWH_JECUp           = t.p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECUp
        self.p_HadWH_SIG_ghw1prime2_1_JHUGen_JECUp              = self.M2g1prime2_HadWH_JECUp       = t.p_HadWH_SIG_ghw1prime2_1E4_JHUGen_JECUp / 1e4**2
        self.p_HadWH_SIG_ghw1_1_ghw1prime2_1_JHUGen_JECUp       = self.M2g1g1prime2_HadWH_JECUp     = t.p_HadWH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen_JECUp / 1e4
        self.p_HadWH_SIG_ghza1prime2_1_JHUGen_JECUp             = self.M2ghzgs1prime2_HadWH_JECUp   = 0
        self.p_HadWH_SIG_ghz1_1_ghza1prime2_1_JHUGen_JECUp      = self.M2g1ghzgs1prime2_HadWH_JECUp = 0

        self.p_HadWH_mavjj_JECUp                                                               = t.p_HadWH_mavjj_JECUp
        self.p_HadWH_mavjj_true_JECUp                                                          = t.p_HadWH_mavjj_true_JECUp
        self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECUp                                                   = t.pConst_HadWH_SIG_ghw1_1_JHUGen_JECUp

        self.M2g1_HadWH_JECUp             *= self.p_HadWH_mavjj_JECUp / self.p_HadWH_mavjj_true_JECUp / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECUp
        self.M2g4_HadWH_JECUp             *= self.p_HadWH_mavjj_JECUp / self.p_HadWH_mavjj_true_JECUp / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECUp
        self.M2g1g4_HadWH_JECUp           *= self.p_HadWH_mavjj_JECUp / self.p_HadWH_mavjj_true_JECUp / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECUp
        self.M2g2_HadWH_JECUp             *= self.p_HadWH_mavjj_JECUp / self.p_HadWH_mavjj_true_JECUp / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECUp
        self.M2g1g2_HadWH_JECUp           *= self.p_HadWH_mavjj_JECUp / self.p_HadWH_mavjj_true_JECUp / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECUp
        self.M2g1prime2_HadWH_JECUp       *= self.p_HadWH_mavjj_JECUp / self.p_HadWH_mavjj_true_JECUp / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECUp
        self.M2g1g1prime2_HadWH_JECUp     *= self.p_HadWH_mavjj_JECUp / self.p_HadWH_mavjj_true_JECUp / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECUp
        self.M2ghzgs1prime2_HadWH_JECUp   *= self.p_HadWH_mavjj_JECUp / self.p_HadWH_mavjj_true_JECUp / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECUp
        self.M2g1ghzgs1prime2_HadWH_JECUp *= self.p_HadWH_mavjj_JECUp / self.p_HadWH_mavjj_true_JECUp / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECUp

        self.p_JQCD_SIG_ghg2_1_JHUGen_JECUp                                                         = t.p_JQCD_SIG_ghg2_1_JHUGen_JECUp
        self.p_JVBF_SIG_ghv1_1_JHUGen_JECUp                                                         = t.p_JVBF_SIG_ghv1_1_JHUGen_JECUp
        self.pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp                                                      = t.pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp

        self.p_JJVBF_S_SIG_ghv1_1_MCFM_JECUp = t.p_JJVBF_S_SIG_ghv1_1_MCFM_JECUp
        self.p_HadZH_S_SIG_ghz1_1_MCFM_JECUp = t.p_HadZH_S_SIG_ghz1_1_MCFM_JECUp
        self.p_HadWH_S_SIG_ghw1_1_MCFM_JECUp = t.p_HadWH_S_SIG_ghw1_1_MCFM_JECUp
        self.p_JJVBF_BKG_MCFM_JECUp = t.p_JJVBF_BKG_MCFM_JECUp
        self.p_HadZH_BKG_MCFM_JECUp = t.p_HadZH_BKG_MCFM_JECUp
        self.p_HadWH_BKG_MCFM_JECUp = t.p_HadWH_BKG_MCFM_JECUp
        self.p_JJQCD_BKG_MCFM_JECUp = t.p_JJQCD_BKG_MCFM_JECUp
        self.p_HadZH_mavjj_JECUp = t.p_HadZH_mavjj_JECUp
        self.p_HadZH_mavjj_true_JECUp = t.p_HadZH_mavjj_true_JECUp
        self.p_HadWH_mavjj_JECUp = t.p_HadWH_mavjj_JECUp
        self.p_HadWH_mavjj_true_JECUp = t.p_HadWH_mavjj_true_JECUp
        self.pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECUp = t.pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECUp
        self.pConst_HadZH_S_SIG_ghz1_1_MCFM_JECUp = t.pConst_HadZH_S_SIG_ghz1_1_MCFM_JECUp
        self.pConst_HadWH_S_SIG_ghw1_1_MCFM_JECUp = t.pConst_HadWH_S_SIG_ghw1_1_MCFM_JECUp
        self.pConst_JJVBF_BKG_MCFM_JECUp = t.pConst_JJVBF_BKG_MCFM_JECUp
        self.pConst_HadZH_BKG_MCFM_JECUp = t.pConst_HadZH_BKG_MCFM_JECUp
        self.pConst_HadWH_BKG_MCFM_JECUp = t.pConst_HadWH_BKG_MCFM_JECUp
        self.pConst_JJQCD_BKG_MCFM_JECUp = t.pConst_JJQCD_BKG_MCFM_JECUp

        self.PFMET_jesUp = t.PFMET_jesUp
        self.nCleanedJetsPt30_jecUp = t.nCleanedJetsPt30_jecUp
        ptUp = sorted(((pt*(1+sigma), i) for i, (pt, sigma) in enumerate(izip(t.JetPt, t.JetSigma))), reverse=True)
        jecUpIndices = [i for pt, i in ptUp if pt>30]

        self.Jet1Pt_jecUp = ptUp[0][0] if ptUp else None
        if self.nCleanedJetsPt30_jecUp < 2:
            self.DiJetMass_jecUp = self.DiJetDEta_jecUp = self.HjjPt_jecUp = -99
        else:
            jets_jecUp = [tlvfromptetaphim(pt, t.JetEta[i], t.JetPhi[i], t.JetMass[i]) for pt, i in ptUp[:2]]
            self.DiJetMass_jecUp = (jets_jecUp[0]+jets_jecUp[1]).M()
            self.DiJetDEta_jecUp = jets_jecUp[0].Eta() - jets_jecUp[1].Eta()
            self.HjjPt_jecUp = (Hvector + jets_jecUp[0] + jets_jecUp[1]).Pt() if nCleanedJets >= 2 else -99

        if self.nCleanedJetsPt30 == self.nCleanedJetsPt30_jecUp:
            self.nCleanedJetsPt30BTagged_bTagSF_jecUp = self.nCleanedJetsPt30BTagged_bTagSF
            self.jetQGLikelihood_jecUp = self.jetQGLikelihood
            self.jetPhi_jecUp = self.jetPhi
            #note that this automatically takes care of the dummyfloatstar case, since if there are no jets
            # nCleanedJetsPt30 will be the same as nCleanedJetsPt30_jecUp
        else:
            self.nCleanedJetsPt30BTagged_bTagSF_jecUp = int(sum(t.JetIsBtaggedWithSF[i] for i in jecUpIndices))
            if jecUpIndices == range(self.nCleanedJetsPt30_jecUp):
                #should happen almost always
                self.jetQGLikelihood_jecUp = self.jetQGLikelihood
                self.jetPhi_jecUp = self.jetPhi
            else:
                self.jetQGLikelihood_jecUp = array("f", [t.JetQGLikelihood[i] for i in jecUpIndices])
                self.jetPhi_jecUp = array("f", [t.JetPhi[i] for i in jecUpIndices])

        self.notdijet_JECUp = self.M2g1_VBF_JECUp <= 0

        #JECDn
        self.p_JJVBF_SIG_ghv1_1_JHUGen_JECDn                    = self.M2g1_VBF_JECDn               = t.p_JJVBF_SIG_ghv1_1_JHUGen_JECDn
        self.p_JJVBF_SIG_ghv4_1_JHUGen_JECDn                    = self.M2g4_VBF_JECDn               = t.p_JJVBF_SIG_ghv4_1_JHUGen_JECDn
        self.p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECDn             = self.M2g1g4_VBF_JECDn             = t.p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECDn
        self.p_JJVBF_SIG_ghv2_1_JHUGen_JECDn                    = self.M2g2_VBF_JECDn               = t.p_JJVBF_SIG_ghv2_1_JHUGen_JECDn
        self.p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECDn             = self.M2g1g2_VBF_JECDn             = t.p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECDn
        self.p_JJVBF_SIG_ghv1prime2_1_JHUGen_JECDn              = self.M2g1prime2_VBF_JECDn         = t.p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECDn / 1e4**2
        self.p_JJVBF_SIG_ghv1_1_ghv1prime2_1_JHUGen_JECDn       = self.M2g1g1prime2_VBF_JECDn       = t.p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_JECDn / 1e4
        self.p_JJVBF_SIG_ghza1prime2_1_JHUGen_JECDn             = self.M2ghzgs1prime2_VBF_JECDn     = t.p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_JECDn / 1e4**2
        self.p_JJVBF_SIG_ghv1_1_ghza1prime2_1_JHUGen_JECDn      = self.M2g1ghzgs1prime2_VBF_JECDn   = t.p_JJVBF_SIG_ghv1_1_ghza1prime2_1E4_JHUGen_JECDn / 1e4

        self.p_JJQCD_SIG_ghg2_1_JHUGen_JECDn                    = self.M2g2_HJJ_JECDn               = t.p_JJQCD_SIG_ghg2_1_JHUGen_JECDn
        self.p_JJQCD_SIG_ghg4_1_JHUGen_JECDn                                                        = t.p_JJQCD_SIG_ghg4_1_JHUGen_JECDn
        self.p_JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_JECDn                                                 = t.p_JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_JECDn

        self.p_HadZH_SIG_ghz1_1_JHUGen_JECDn                    = self.M2g1_HadZH_JECDn             = t.p_HadZH_SIG_ghz1_1_JHUGen_JECDn
        self.p_HadZH_SIG_ghz4_1_JHUGen_JECDn                    = self.M2g4_HadZH_JECDn             = t.p_HadZH_SIG_ghz4_1_JHUGen_JECDn
        self.p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECDn             = self.M2g1g4_HadZH_JECDn           = t.p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECDn
        self.p_HadZH_SIG_ghz2_1_JHUGen_JECDn                    = self.M2g2_HadZH_JECDn             = t.p_HadZH_SIG_ghz2_1_JHUGen_JECDn
        self.p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECDn             = self.M2g1g2_HadZH_JECDn           = t.p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECDn
        self.p_HadZH_SIG_ghz1prime2_1_JHUGen_JECDn              = self.M2g1prime2_HadZH_JECDn       = t.p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECDn / 1e4**2
        self.p_HadZH_SIG_ghz1_1_ghz1prime2_1_JHUGen_JECDn       = self.M2g1g1prime2_HadZH_JECDn     = t.p_HadZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen_JECDn / 1e4
        self.p_HadZH_SIG_ghza1prime2_1_JHUGen_JECDn             = self.M2ghzgs1prime2_HadZH_JECDn   = t.p_HadZH_SIG_ghza1prime2_1E4_JHUGen_JECDn / 1e4**2
        self.p_HadZH_SIG_ghz1_1_ghza1prime2_1_JHUGen_JECDn      = self.M2g1ghzgs1prime2_HadZH_JECDn = t.p_HadZH_SIG_ghz1_1_ghza1prime2_1E4_JHUGen_JECDn / 1e4

        self.p_HadZH_mavjj_JECDn                                                                    = t.p_HadZH_mavjj_JECDn
        self.p_HadZH_mavjj_true_JECDn                                                               = t.p_HadZH_mavjj_true_JECDn
        self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECDn                                                   = t.pConst_HadZH_SIG_ghz1_1_JHUGen_JECDn

        self.M2g1_HadZH_JECDn             *= self.p_HadZH_mavjj_JECDn / self.p_HadZH_mavjj_true_JECDn / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECDn
        self.M2g4_HadZH_JECDn             *= self.p_HadZH_mavjj_JECDn / self.p_HadZH_mavjj_true_JECDn / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECDn
        self.M2g1g4_HadZH_JECDn           *= self.p_HadZH_mavjj_JECDn / self.p_HadZH_mavjj_true_JECDn / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECDn
        self.M2g2_HadZH_JECDn             *= self.p_HadZH_mavjj_JECDn / self.p_HadZH_mavjj_true_JECDn / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECDn
        self.M2g1g2_HadZH_JECDn           *= self.p_HadZH_mavjj_JECDn / self.p_HadZH_mavjj_true_JECDn / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECDn
        self.M2g1prime2_HadZH_JECDn       *= self.p_HadZH_mavjj_JECDn / self.p_HadZH_mavjj_true_JECDn / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECDn
        self.M2g1g1prime2_HadZH_JECDn     *= self.p_HadZH_mavjj_JECDn / self.p_HadZH_mavjj_true_JECDn / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECDn
        self.M2ghzgs1prime2_HadZH_JECDn   *= self.p_HadZH_mavjj_JECDn / self.p_HadZH_mavjj_true_JECDn / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECDn
        self.M2g1ghzgs1prime2_HadZH_JECDn *= self.p_HadZH_mavjj_JECDn / self.p_HadZH_mavjj_true_JECDn / self.pConst_HadZH_SIG_ghz1_1_JHUGen_JECDn

        self.p_HadWH_SIG_ghw1_1_JHUGen_JECDn                    = self.M2g1_HadWH_JECDn             = t.p_HadWH_SIG_ghw1_1_JHUGen_JECDn
        self.p_HadWH_SIG_ghw4_1_JHUGen_JECDn                    = self.M2g4_HadWH_JECDn             = t.p_HadWH_SIG_ghw4_1_JHUGen_JECDn
        self.p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECDn             = self.M2g1g4_HadWH_JECDn           = t.p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECDn
        self.p_HadWH_SIG_ghw2_1_JHUGen_JECDn                    = self.M2g2_HadWH_JECDn             = t.p_HadWH_SIG_ghw2_1_JHUGen_JECDn
        self.p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECDn             = self.M2g1g2_HadWH_JECDn           = t.p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECDn
        self.p_HadWH_SIG_ghw1prime2_1_JHUGen_JECDn              = self.M2g1prime2_HadWH_JECDn       = t.p_HadWH_SIG_ghw1prime2_1E4_JHUGen_JECDn / 1e4**2
        self.p_HadWH_SIG_ghw1_1_ghw1prime2_1_JHUGen_JECDn       = self.M2g1g1prime2_HadWH_JECDn     = t.p_HadWH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen_JECDn / 1e4
        self.p_HadWH_SIG_ghza1prime2_1_JHUGen_JECDn             = self.M2ghzgs1prime2_HadWH_JECDn   = 0
        self.p_HadWH_SIG_ghz1_1_ghza1prime2_1_JHUGen_JECDn      = self.M2g1ghzgs1prime2_HadWH_JECDn = 0

        self.p_HadWH_mavjj_JECDn                                                                    = t.p_HadWH_mavjj_JECDn
        self.p_HadWH_mavjj_true_JECDn                                                               = t.p_HadWH_mavjj_true_JECDn
        self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECDn                                                   = t.pConst_HadWH_SIG_ghw1_1_JHUGen_JECDn

        self.M2g1_HadWH_JECDn             *= self.p_HadWH_mavjj_JECDn / self.p_HadWH_mavjj_true_JECDn / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECDn
        self.M2g4_HadWH_JECDn             *= self.p_HadWH_mavjj_JECDn / self.p_HadWH_mavjj_true_JECDn / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECDn
        self.M2g1g4_HadWH_JECDn           *= self.p_HadWH_mavjj_JECDn / self.p_HadWH_mavjj_true_JECDn / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECDn
        self.M2g2_HadWH_JECDn             *= self.p_HadWH_mavjj_JECDn / self.p_HadWH_mavjj_true_JECDn / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECDn
        self.M2g1g2_HadWH_JECDn           *= self.p_HadWH_mavjj_JECDn / self.p_HadWH_mavjj_true_JECDn / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECDn
        self.M2g1prime2_HadWH_JECDn       *= self.p_HadWH_mavjj_JECDn / self.p_HadWH_mavjj_true_JECDn / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECDn
        self.M2g1g1prime2_HadWH_JECDn     *= self.p_HadWH_mavjj_JECDn / self.p_HadWH_mavjj_true_JECDn / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECDn
        self.M2ghzgs1prime2_HadWH_JECDn   *= self.p_HadWH_mavjj_JECDn / self.p_HadWH_mavjj_true_JECDn / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECDn
        self.M2g1ghzgs1prime2_HadWH_JECDn *= self.p_HadWH_mavjj_JECDn / self.p_HadWH_mavjj_true_JECDn / self.pConst_HadWH_SIG_ghw1_1_JHUGen_JECDn

        self.p_JQCD_SIG_ghg2_1_JHUGen_JECDn                                                         = t.p_JQCD_SIG_ghg2_1_JHUGen_JECDn
        self.p_JVBF_SIG_ghv1_1_JHUGen_JECDn                                                         = t.p_JVBF_SIG_ghv1_1_JHUGen_JECDn
        self.pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn                                                      = t.pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn

        self.p_JJVBF_S_SIG_ghv1_1_MCFM_JECDn = t.p_JJVBF_S_SIG_ghv1_1_MCFM_JECDn
        self.p_HadZH_S_SIG_ghz1_1_MCFM_JECDn = t.p_HadZH_S_SIG_ghz1_1_MCFM_JECDn
        self.p_HadWH_S_SIG_ghw1_1_MCFM_JECDn = t.p_HadWH_S_SIG_ghw1_1_MCFM_JECDn
        self.p_JJVBF_BKG_MCFM_JECDn = t.p_JJVBF_BKG_MCFM_JECDn
        self.p_HadZH_BKG_MCFM_JECDn = t.p_HadZH_BKG_MCFM_JECDn
        self.p_HadWH_BKG_MCFM_JECDn = t.p_HadWH_BKG_MCFM_JECDn
        self.p_JJQCD_BKG_MCFM_JECDn = t.p_JJQCD_BKG_MCFM_JECDn
        self.p_HadZH_mavjj_JECDn = t.p_HadZH_mavjj_JECDn
        self.p_HadZH_mavjj_true_JECDn = t.p_HadZH_mavjj_true_JECDn
        self.p_HadWH_mavjj_JECDn = t.p_HadWH_mavjj_JECDn
        self.p_HadWH_mavjj_true_JECDn = t.p_HadWH_mavjj_true_JECDn
        self.pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECDn = t.pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECDn
        self.pConst_HadZH_S_SIG_ghz1_1_MCFM_JECDn = t.pConst_HadZH_S_SIG_ghz1_1_MCFM_JECDn
        self.pConst_HadWH_S_SIG_ghw1_1_MCFM_JECDn = t.pConst_HadWH_S_SIG_ghw1_1_MCFM_JECDn
        self.pConst_JJVBF_BKG_MCFM_JECDn = t.pConst_JJVBF_BKG_MCFM_JECDn
        self.pConst_HadZH_BKG_MCFM_JECDn = t.pConst_HadZH_BKG_MCFM_JECDn
        self.pConst_HadWH_BKG_MCFM_JECDn = t.pConst_HadWH_BKG_MCFM_JECDn
        self.pConst_JJQCD_BKG_MCFM_JECDn = t.pConst_JJQCD_BKG_MCFM_JECDn

        self.PFMET_jesDn = t.PFMET_jesDn
        self.nCleanedJetsPt30_jecDn = t.nCleanedJetsPt30_jecDn
        ptDn = sorted(((pt*(1+sigma), i) for i, (pt, sigma) in enumerate(izip(t.JetPt, t.JetSigma))), reverse=True)
        jecDnIndices = [i for pt, i in ptDn if pt>30]

        self.Jet1Pt_jecDn = ptDn[0][0] if ptDn else None
        if self.nCleanedJetsPt30_jecDn < 2:
            self.DiJetMass_jecDn = self.DiJetDEta_jecDn = self.HjjPt_jecDn = -99
        else:
            jets_jecDn = [tlvfromptetaphim(pt, t.JetEta[i], t.JetPhi[i], t.JetMass[i]) for pt, i in ptDn[:2]]
            self.DiJetMass_jecDn = (jets_jecDn[0]+jets_jecDn[1]).M()
            self.DiJetDEta_jecDn = jets_jecDn[0].Eta() - jets_jecDn[1].Eta()
            self.HjjPt_jecDn = (Hvector + jets_jecDn[0] + jets_jecDn[1]).Pt() if nCleanedJets >= 2 else -99

        if self.nCleanedJetsPt30 == self.nCleanedJetsPt30_jecDn:
            self.nCleanedJetsPt30BTagged_bTagSF_jecDn = self.nCleanedJetsPt30BTagged_bTagSF
            self.jetQGLikelihood_jecDn = self.jetQGLikelihood
            self.jetPhi_jecDn = self.jetPhi
            #note that this automatically takes care of the dummyfloatstar case, since if there are no jets
            # nCleanedJetsPt30 will be the same as nCleanedJetsPt30_jecDn
        else:
            self.nCleanedJetsPt30BTagged_bTagSF_jecDn = int(sum(t.JetIsBtaggedWithSF[i] for i in jecDnIndices))
            if jecDnIndices == range(self.nCleanedJetsPt30_jecDn):
                #should happen almost always
                self.jetQGLikelihood_jecDn = self.jetQGLikelihood
                self.jetPhi_jecDn = self.jetPhi
            else:
                self.jetQGLikelihood_jecDn = array("f", [t.JetQGLikelihood[i] for i in jecDnIndices])
                self.jetPhi_jecDn = array("f", [t.JetPhi[i] for i in jecDnIndices])

        self.notdijet_JECDn = self.M2g1_VBF_JECDn <= 0

        return self

    @cache_instancemethod
    def __len__(self):
        if self.isdummy:
            return 0
        elif self.maxevent is None or self.maxevent >= self.tree.GetEntries():
            return self.tree.GetEntries() - self.minevent
        else:
            return self.maxevent - self.minevent + 1

    def Show(self, *args, **kwargs):
        self.tree.Show(*args, **kwargs)

##########
#Category#
##########

    categorizations = []
    for JEC in enums.JECsystematics:
      for btag in enums.btagsystematics:
        if JEC != "Nominal" and btag != "Nominal": continue
        append = [
            categorization.SingleCategorizationgm4l("SM", JEC, btag),
            categorization.SingleCategorizationgm4l("0-", JEC, btag),
            categorization.SingleCategorizationgm4l("a2", JEC, btag),
            categorization.SingleCategorizationgm4l("L1", JEC, btag),
            categorization.SingleCategorizationgm4l("L1Zg", JEC, btag),
        ]
        append += [
            categorization.MultiCategorization("0P_or_0M_or_a2_or_L1_or_L1Zg"+btag.appendname+JEC.appendname, *append)
        ]
        categorizations += append
    categorizations.append(categorization.NoCategorization())
    assert len({_.category_function_name for _ in categorizations}) == len(categorizations)
    del append, btag, JEC

    def p_Gen_ttH_SIG_kappa_1_JHUGen(self):
      h = self.treesample.hffhypothesis
      if h == "Hff0+": return 1
      if h == "fCP0.5": return 0
      assert False
    def p_Gen_ttH_SIG_kappa_tilde_1_JHUGen(self):
      h = self.treesample.hffhypothesis
      if h == "Hff0-": return 1/1.6**2
      if h == "fCP0.5": return 0
      assert False
    def p_Gen_ttH_SIG_kappa_1_kappa_tilde_1_JHUGen(self):
      h = self.treesample.hffhypothesis
      if h == "Hff0+": return 0
      if h == "Hff0-": return 0
      if h == "fCP0.5": return 1/1.6
      assert False

#############
#Init things#
#############

    @classmethod
    def initcategoryfunctions(cls):
        for _ in cls.categorizations:
            setattr(cls, _.category_function_name, _.get_category_function())
            if isinstance(_, categorization.BaseSingleCategorization):
                setattr(cls, _.pVBF_function_name, _.get_pVBF_function())
                setattr(cls, _.pHJJ_function_name, _.get_pHJJ_function())
                setattr(cls, _.pZH_function_name, _.get_pZH_function())
                setattr(cls, _.pWH_function_name, _.get_pWH_function())
                #setattr(cls, _.pHJ_function_name, _.get_pHJ_function())
                #setattr(cls, _.pVBF1j_function_name, _.get_pVBF1j_function())
                #setattr(cls, _.pAux_function_name, _.get_pAux_function())

    @classmethod
    def initsystematics(cls):
        for name, discriminant in inspect.getmembers(cls, predicate=lambda x: isinstance(x, MakeSystematics)):
            nominal = discriminant.getnominal()
            alternates = discriminant.getalternates()

            setattr(cls, nominal.__name__, nominal)
            for alternate in alternates:
                assert not hasattr(cls, alternate.__name__), alternate.__name__
                setattr(cls, alternate.__name__, alternate)

    def initlists(self):
        self.toaddtotree = [
            "D_bkg",
            "D_bkg_ResUp",
            "D_bkg_ResDown",
            "D_bkg_ScaleUp",
            "D_bkg_ScaleDown",
            "D_bkg_VBFdecay_ResUp",
            "D_bkg_VBFdecay_ResDown",
            "D_bkg_VBFdecay_ScaleUp",
            "D_bkg_VBFdecay_ScaleDown",
            "D_bkg_HadVHdecay_ResUp",
            "D_bkg_HadVHdecay_ResDown",
            "D_bkg_HadVHdecay_ScaleUp",
            "D_bkg_HadVHdecay_ScaleDown",
            "D_2jet_0plus",
            "D_2jet_0minus",
            "D_2jet_a2",
            "D_2jet_L1",
            "D_2jet_L1Zg",
            "D_HadWH_0plus",
            "D_HadWH_0minus",
            "D_HadWH_a2",
            "D_HadWH_L1",
            "D_HadWH_L1Zg",
            "D_HadZH_0plus",
            "D_HadZH_0minus",
            "D_HadZH_a2",
            "D_HadZH_L1",
            "D_HadZH_L1Zg",
            "D_0minus_decay",
            "D_CP_decay",
            "D_CP_decay_new",
            "D_0hplus_decay",
            "D_int_decay",
            "D_int_decay_new",
            "D_L1_decay",
            "D_L1int_decay",
            "D_L1int_decay_new",
            "D_L1Zg_decay",
            "D_L1Zgint_decay",
            "D_L1Zgint_decay_new",
            "D_L1L1Zg_decay",
            "D_L1L1Zgint_decay",
            "D_L1L1Zgint_decay_new",
        ]

        self.exceptions = [
            "D_eL_decay",
            "D_eLint_decay",
            "D_eLint_decay_new",
            "D_eR_decay",
            "D_eRint_decay",
            "D_eRint_decay_new",
            "D_eLeR_decay",
            "D_eLeRint_decay",
            "D_eLeRint_decay_new",

            "D_4couplings_general",
            "D_4couplings_general_raw",
            "binning_4couplings_decay",
            "binning_4couplings_VBFdecay",
            "binning_4couplings_VBFdecay_JECUp",
            "binning_4couplings_VBFdecay_JECDn",
            "binning_4couplings_HadVHdecay",
            "binning_4couplings_HadVHdecay_JECUp",
            "binning_4couplings_HadVHdecay_JECDn",
            "foldbins_4couplings_decay",
            "foldbins_4couplings_VBFdecay",
            "foldbins_4couplings_HadVHdecay",

            "D_bkg_kin_VBFdecay",
            "D_bkg_kin_VBFdecay_JECUp",
            "D_bkg_kin_VBFdecay_JECDn",
            "D_bkg_kin_HadVHdecay",
            "D_bkg_kin_HadVHdecay_JECUp",
            "D_bkg_kin_HadVHdecay_JECDn",

            "categorizations",
            "cconstantforDbkg",
            "cconstantforD2jet",
            "cconstantforDHadWH",
            "cconstantforDHadZH",
            "checkfunctions",
            "counters",
            "cutoffs",
            "definitelyexists",
            "exceptions",
            "f",
            "failedtree",
            "filename",
            "GEN",
            "hypothesis",
            "initcategoryfunctions",
            "initlists",
            "initsystematics",
            "isalternate",
            "isbkg",
            "isdata",
            "isdummy",
            "isZX",
            "kfactors",
            "maxevent",
            "minevent",
            "nevents",
            "nevents2e2mu",
            "next",
            "passesblindcut",
            "per_event_scale_factor",
            "printevery",
            "productionmode",
            "toaddtotree",
            "toaddtotree_float",
            "toaddtotree_int",
            "tree",
            "treesample",
            "Show",
            "unblind",
            "xsec",
            "year",
        ]
        self.toaddtotree_float = []
        self.toaddtotree_int = [
            "D_4couplings_decay_raw",
            "D_4couplings_decay",
        ]

        proddiscriminants = [
            "D_bkg_{prod}decay",
            "D_0minus_{prod}",
            "D_CP_{prod}",
            "D_CP_{prod}_new",
            "D_0hplus_{prod}",
            "D_int_{prod}",
            "D_int_{prod}_new",
            "D_L1_{prod}",
            "D_L1int_{prod}",
            "D_L1int_{prod}_new",
            "D_L1Zg_{prod}",
            "D_L1Zgint_{prod}",
            "D_L1Zgint_{prod}_new",
            "D_0minus_{prod}decay",
            "D_0hplus_{prod}decay",
            "D_L1_{prod}decay",
            "D_L1Zg_{prod}decay",
        ]
        STXSdiscriminants = [
            "D_STXS_stage1p1",
            "D_STXS_stage1p1_bTagSFUp",
            "D_STXS_stage1p1_bTagSFDn",
            "D_4couplings_VBFdecay_raw",
            "D_4couplings_VBFdecay",
            "D_4couplings_HadVHdecay_raw",
            "D_4couplings_HadVHdecay",
        ]

        for JEC in "", "_JECUp", "_JECDn":
            for prod in ("VBF", "HadVH"):
                self.toaddtotree += [_.format(prod=prod)+JEC for _ in proddiscriminants]
            self.toaddtotree_int += [_+JEC for _ in STXSdiscriminants if not (JEC and "bTagSF" in _)]

        for _ in self.categorizations:
            if self.GEN and _.issystematic:
                self.exceptions.append(_.category_function_name)
            else:
                self.toaddtotree_int.append(_.category_function_name)
            if isinstance(_, categorization.BaseSingleCategorization):
                for name in (
                    _.pHJJ_function_name,
                    _.pVBF_function_name,
                    _.pZH_function_name,
                    _.pWH_function_name,
                    #_.pHJ_function_name,
                    #_.pVBF1j_function_name,
                    #_.pAux_function_name,
                ):
                    if name not in self.exceptions:
                        self.exceptions.append(name)

        if self.treesample.productionmode == "data":
            self.exceptions.append("MC_weight_nominal")
        else:
            self.toaddtotree.append("MC_weight_nominal")

        if self.GEN:
            for lst in self.toaddtotree, self.toaddtotree_int, self.toaddtotree_float:
                for _ in lst[:]:
                    if "jecup" in _.lower() or "jecdn" in _.lower():
                        lst.remove(_)
                        self.exceptions.append(_)

        self.kfactors = []
        if not self.GEN:
            if self.treesample.productionmode == "qqZZ": self.kfactors += ["KFactor_EW_qqZZ", "KFactor_QCD_qqZZ_M"]
            if self.treesample.productionmode == "ggZZ": self.kfactors += ["KFactor_QCD_ggZZ_Nominal"]

        if self.treesample.productionmode == "ttH" and getattr(self.treesample, "alternategenerator", None) == None:
            if self.treesample.hffhypothesis == "Hff0+":
                self.toaddtotree.append("p_Gen_ttH_SIG_kappa_1_JHUGen")
                self.exceptions.append("p_Gen_ttH_SIG_kappa_tilde_1_JHUGen")
                self.exceptions.append("p_Gen_ttH_SIG_kappa_1_kappa_tilde_1_JHUGen")
            elif self.treesample.hffhypothesis == "Hff0-":
                self.exceptions.append("p_Gen_ttH_SIG_kappa_1_JHUGen")
                self.toaddtotree.append("p_Gen_ttH_SIG_kappa_tilde_1_JHUGen")
                self.exceptions.append("p_Gen_ttH_SIG_kappa_1_kappa_tilde_1_JHUGen")
            elif self.treesample.hffhypothesis == "fCP0.5":
                self.toaddtotree.append("p_Gen_ttH_SIG_kappa_1_JHUGen")
                self.toaddtotree.append("p_Gen_ttH_SIG_kappa_tilde_1_JHUGen")
                self.toaddtotree.append("p_Gen_ttH_SIG_kappa_1_kappa_tilde_1_JHUGen")
            else:
                assert False
        else:
            self.exceptions.append("p_Gen_ttH_SIG_kappa_1_JHUGen")
            self.exceptions.append("p_Gen_ttH_SIG_kappa_tilde_1_JHUGen")
            self.exceptions.append("p_Gen_ttH_SIG_kappa_1_kappa_tilde_1_JHUGen")

    passesblindcut = config.blindcut

def TreeWrapperFactory(treesample, minevent=0, maxevent=None, LSF=Fake_LSF_creating()):
    if treesample.production.LHE:
        from lhewrapper import LHEWrapper
        cls = LHEWrapper
    else:
        cls = TreeWrapper
    return cls(treesample, minevent=minevent, maxevent=maxevent, LSF=LSF)
