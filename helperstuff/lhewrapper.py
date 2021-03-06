#!/usr/bin/env python

from abc import abstractmethod
import inspect
import os
import random
import sys

import ROOT

import CJLSTscripts
import config
from makesystematics import MakeSystematics
from samples import ReweightingSample
from treewrapper import TreeWrapperBase
from utilities import cache, cache_instancemethod, callclassinitfunctions, product, requirecmsenv, tlvfromptetaphim

requirecmsenv(os.path.join(config.repositorydir, "CMSSW_8_0_20"))

from ZZMatrixElement.MELA.mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar

class LHEEvent(object):
  def __init__(self, event, bkg4l=False):
    self.eventstr = event

    #from https://github.com/hroskes/HiggsAnalysis-ZZMatrixElement/blob/b0ad77b/PythonWrapper/python/mela.py#L338-L366
    lines = event.split("\n")
    lines = [line for line in lines if not ("<event>" in line or "</event>" in line or not line.split("#")[0].strip())]
    nparticles, _, weight, _, _, _ = lines[0].split()
    nparticles = int(nparticles)
    self.weight = float(weight)
    if nparticles != len(lines)-1:
      raise ValueError("Wrong number of particles! Should be {}, have {}".replace(nparticles, len(lines)-1))
    self.gendaughters, self.genmothers, self.genassociated = daughters, mothers, associated = [], [], []
    ids = [None]
    mother1s = [None]
    mother2s = [None]
    for line in lines[1:]:
      id, status, mother1, mother2 = (int(_) for _ in line.split()[0:4])
      ids.append(id)
      mother1s.append(mother1)
      mother2s.append(mother2)
      if status == -1:
        mothers.append(SimpleParticle_t(line))
      elif status == 1 and (1 <= abs(id) <= 6 or 11 <= abs(id) <= 16 or abs(id) in (21, 22)):
        while True:
          if 11 <= abs(id) <= 16 and bkg4l:
            daughters.append(SimpleParticle_t(line))
            break
          if mother1 != mother2 or mother1 is None:
            associated.append(SimpleParticle_t(line))
            break
          if ids[mother1] == 25:
            daughters.append(SimpleParticle_t(line))
            break
          mother2 = mother2s[mother1]
          mother1 = mother1s[mother1]
      elif bkg4l and id == 25:
        raise ValueError("Can't have explicit Higgs for bkg4l")

    if bkg4l and len(daughters) != 4:
      raise ValueError("{} leptons instead of 4 for bkg4l!".format(len(daughters)))

  @staticmethod
  def smear(particle):
    id, p = particle
    if abs(id) == 11:
      newid = id
      smearpt = config.LHEsmearptelectron
    elif abs(id) == 13:
      newid = id
      smearpt = config.LHEsmearptmuon
    elif 1 <= abs(id) <= 5 or abs(id) == 21:
      newid = 0
      smearpt = config.LHEsmearptjet
    elif abs(id) == 15:
      newid = id
      smearpt = 0  #taus are dropped anyway in passparticlecuts
    else:
      raise ValueError("Don't know how to smear {}".format(id))

    newp = tlvfromptetaphim(random.gauss(p.Pt(), smearpt), p.Eta(), p.Phi(), p.M())
    return SimpleParticle_t(newid, newp)

  @staticmethod
  def passparticlecuts(particle):
    id, p = particle
    if abs(id) == 11:
      if p.Pt() < 7 or abs(p.Eta()) > 2.5: return False
    elif abs(id) == 13:
      if p.Pt() < 5 or abs(p.Eta()) > 2.4: return False
    elif abs(id) == 15:
      return False
    elif id == 0 or 1 <= abs(id) <= 5 or id == 21:
      if p.Pt() < 30 or abs(p.Eta()) > 4.7: return False
    else:
      raise ValueError("Don't know how to cut {}".format(id))
    return True

  @property
  @cache_instancemethod
  def recodaughters(self):
    result = [self.smear(_) for _ in self.gendaughters]
    result = [_ for _ in result if self.passparticlecuts(_)]
    return result

  @property
  @cache_instancemethod
  def recoassociated(self):
    result = [self.smear(_) for _ in self.genassociated]
    result = [_ for _ in result if self.passparticlecuts(_)]
    return result

  @property
  @cache_instancemethod
  def ZZMass(self):
    return sum((p for id, p in self.recodaughters), ROOT.TLorentzVector()).M()

  @property
  def genmelaargs(self):
    return SimpleParticleCollection_t(self.gendaughters), SimpleParticleCollection_t(self.genassociated), SimpleParticleCollection_t(self.genmothers), True

  @property
  def recomelaargs(self):
    return SimpleParticleCollection_t(self.recodaughters), SimpleParticleCollection_t(self.recoassociated), 0, False

  @property
  @cache_instancemethod
  def passcuts(self):
    if len(self.recodaughters) < 4: return False
    if not config.m4lmin < self.ZZMass < config.m4lmax: return False
    return True

  def __str__(self):
    return self.eventstr

  @property
  @cache_instancemethod
  def ZZFlav(self):
    return product(id for id, p in self.recodaughters)

  @property
  @cache_instancemethod
  def njets(self):
    return sum(id == 0 for id, p in self.recoassociated)

class TreeWrapperMELA(TreeWrapperBase):
  def __init__(self, treesample, minevent=0, maxevent=None):
    self.mela = self.initmela()
    self.dofa3stuff = False#(treesample.productionmode == "qqZZ")
    super(TreeWrapperMELA, self).__init__(treesample, minevent, maxevent)

  @staticmethod
  @cache
  def initmela(*args, **kwargs):
    return Mela(*args, **kwargs)

  def calcrecoMEs(self, *setinputeventargs):
    m = self.mela
    m.setInputEvent(*setinputeventargs)

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
    m.ghg2 = m.ghz1 = 1
    self.M2g1_decay = m.computeP(True)

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
    m.ghg2 = m.ghz2 = 1
    self.M2g2_decay = m.computeP(True)

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
    m.ghg2 = 1; m.ghz1_prime2 = 1e4
    self.M2g1prime2_decay = m.computeP(True)

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
    m.ghg2 = 1; m.ghzgs1_prime2 = 1e4
    self.M2ghzgs1prime2_decay = m.computeP(True)

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
    m.ghg2 = m.ghz1 = 1; m.ghz1_prime2 = 1e4
    self.M2g1g1prime2_decay = m.computeP(True) - self.M2g1_decay - self.M2g1prime2_decay

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
    m.ghg2 = m.ghz1 = 1; m.ghzgs1_prime2 = 1e4
    self.M2g1ghzgs1prime2_decay = m.computeP(True) - self.M2g1_decay - self.M2ghzgs1prime2_decay

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
    m.ghg2 = 1; m.ghz1_prime2 = m.ghzgs1_prime2 = 1e4
    self.M2g1prime2ghzgs1prime2_decay = m.computeP(True) - self.M2g1prime2_decay - self.M2ghzgs1prime2_decay

    ###############
    #contact terms#
    ###############

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
    m.ghg2 = m.ghzzp1 = m.ezp_L_E = m.ezp_L_M = m.ezp_L_T = 1
    self.M2eL_decay = m.computeP(True)

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
    m.ghg2 = m.ghzzp1 = m.ezp_R_E = m.ezp_R_M = m.ezp_R_T = 1
    self.M2eR_decay = m.computeP(True)

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
    m.ghg2 = m.ghz1 = m.ghzzp1 = m.ezp_L_E = m.ezp_L_M = m.ezp_L_T = 1
    self.M2g1eL_decay = m.computeP(True) - self.M2eL_decay - self.M2g1_decay

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
    m.ghg2 = m.ghz1 = m.ghzzp1 = m.ezp_R_E = m.ezp_R_M = m.ezp_R_T = 1
    self.M2g1eR_decay = m.computeP(True) - self.M2eR_decay - self.M2g1_decay

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
    m.ghg2 = m.ghzzp1 = m.ezp_L_E = m.ezp_L_M = m.ezp_L_T = m.ezp_R_E = m.ezp_R_M = m.ezp_R_T = 1
    self.M2eLeR_decay = m.computeP(True) - self.M2eL_decay - self.M2eR_decay

    ############
    #Dbkg stuff#
    ############

    m.setProcess(TVar.HSMHiggs, TVar.JHUGen, TVar.ZZGG)
    self.p_m4l_SIG = m.computePM4l(TVar.SMSyst_None)

    m.setProcess(TVar.bkgZZ, TVar.JHUGen, TVar.ZZGG)
    self.p_m4l_BKG = m.computePM4l(TVar.SMSyst_None)

    m.setProcess(TVar.bkgZZ, TVar.MCFM, TVar.ZZQQB)
    self.M2qqZZ = m.computeP(True)

    assert not config.useQGTagging
    self.cconstantforDbkg = CJLSTscripts.getDbkgConstant(self.ZZFlav(), self.ZZMass())
    self.cconstantforD2jet = CJLSTscripts.getDVBF2jetsConstant_shiftWP(self.ZZMass(), config.useQGTagging, 0.5)
    self.cconstantforDHadWH = CJLSTscripts.getDWHhConstant_shiftWP(self.ZZMass(), config.useQGTagging, 0.5)
    self.cconstantforDHadZH = CJLSTscripts.getDZHhConstant_shiftWP(self.ZZMass(), config.useQGTagging, 0.5)

    self.M2g1prime2_decay /= 1e4**2
    self.M2ghzgs1prime2_decay /= 1e4**2
    self.M2g1g1prime2_decay /= 1e4
    self.M2g1ghzgs1prime2_decay /= 1e4
    self.M2g1prime2ghzgs1prime2_decay /= 1e4**2

    self.decayangles = m.computeDecayAngles()

    if self.dofa3stuff:
      m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
      m.ghg2 = m.ghz4 = 1
      self.M2g4_decay = m.computeP(True)

      m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
      m.ghg2 = m.ghz1 = m.ghz4 = 1
      self.M2g1g4_decay = m.computeP(True) - self.M2g1_decay - self.M2g4_decay

      if self.event.njets >= 2:
        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
        m.ghz1 = 1
        self.M2g1_VBF = m.computeProdP(True)

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
        m.ghz4 = 1
        self.M2g4_VBF = m.computeProdP(True)

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
        m.ghz1 = m.ghz4 = 1
        self.M2g1g4_VBF = m.computeProdP(True) - self.M2g1_VBF - self.M2g4_VBF

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.Had_ZH)
        m.ghz1 = 1
        self.M2g1_HadZH = m.computeProdP(True)

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.Had_ZH)
        m.ghz4 = 1
        self.M2g4_HadZH = m.computeProdP(True)

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.Had_ZH)
        m.ghz1 = m.ghz4 = 1
        self.M2g1g4_HadZH = m.computeProdP(True) - self.M2g1_HadZH - self.M2g4_HadZH

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.Had_WH)
        m.ghz1 = 1
        self.M2g1_HadWH = m.computeProdP(True)

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.Had_WH)
        m.ghz4 = 1
        self.M2g4_HadWH = m.computeProdP(True)

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.Had_WH)
        m.ghz1 = m.ghz4 = 1
        self.M2g1g4_HadWH = m.computeProdP(True) - self.M2g1_HadWH - self.M2g4_HadWH

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJQCD)
        m.ghg2 = 1
        self.M2g2_HJJ = m.computeProdP(True)

        self.notdijet = False
      else:
        self.notdijet = True

    m.resetInputEvent()

  @classmethod
  def initsystematics(cls):
    for name, discriminant in inspect.getmembers(cls, predicate=lambda x: isinstance(x, MakeSystematics)):
      nominal = discriminant.getnominal()
      setattr(cls, discriminant.name, nominal)

  @abstractmethod
  def ZZMass(self): return self.decayangles.qH

  def Z1Mass(self): return min((self.decayangles.m1, self.decayangles.m2), key=lambda x: abs(x-91.1876))
  def Z2Mass(self): return max((self.decayangles.m1, self.decayangles.m2), key=lambda x: abs(x-91.1876))
  def costhetastar(self): return self.decayangles.costhetastar
  def costheta1(self): return self.decayangles.costheta1
  def costheta2(self): return self.decayangles.costheta2
  def Phi1(self): return self.decayangles.Phi1
  def Phi(self): return self.decayangles.Phi

  @abstractmethod
  def ZZFlav(self): pass

  def category_0P_or_0M(self):
    assert False, self
    if self.D_2jet_0plus() > 0.5 or self.D_2jet_0minus() > 0.5:
      return CJLSTscripts.VBF2jTaggedAC19
    elif self.D_HadZH_0plus() > 0.5 or self.D_HadZH_0minus() > 0.5 or self.D_HadWH_0plus() > 0.5 or self.D_HadWH_0minus() > 0.5:
      return CJLSTscripts.VHHadrTaggedAC19
    elif self.ZZPt > 120:
      return CJLSTscripts.BoostedAC19
    return CJLSTscripts.UntaggedAC19

@callclassinitfunctions("initweightfunctions", "initsystematics")
class TreeWrapperPythia(TreeWrapperMELA):
  def __init__(self, treesample, minevent=0, maxevent=None):
    assert minevent == 0 and maxevent is None
    self.event = None
    self.sumofweights = None
    self.bkg4l = True
    super(TreeWrapperPythia, self).__init__(treesample, minevent, maxevent)
    self.preliminaryloop()
    self.tree = ROOT.TChain("demo/GenEvents")
    self.tree.Add(treesample.pythiafile)

  @cache_instancemethod
  def __len__(self):
    return self.tree.GetEntries()

  def preliminaryloop(self):
    self.sumofweights = len(self)

  @classmethod
  def initweightfunctions(cls):
    pass
  @classmethod
  def initsystematics(cls):
    pass

@callclassinitfunctions("initweightfunctions", "initsystematics")
class LHEWrapper(TreeWrapperMELA):
  def __init__(self, treesample, minevent=0, maxevent=None):
    assert minevent == 0 and maxevent is None
    self.event = None
    self.sumofweights = None
    self.bkg4l = (treesample.productionmode == "qqZZ")
    super(LHEWrapper, self).__init__(treesample, minevent, maxevent)
    self.preliminaryloop()

  @cache_instancemethod
  def __len__(self):
    with open(self.treesample.LHEfile) as f:
      return f.read().count("</event>")

  @property
  @cache_instancemethod
  def xsec(self):
    if self.isdata: return None
    with open(self.treesample.LHEfile) as f:
      for line in f:
        if "<init>" in line:
          break
      next(f)
      line = next(f)
      return float(line.split()[0])

  def preliminaryloop(self):
    if self.isdummy: return
    i = 0
    sumofweights = 0
    weightfunction = self.treesample.get_MC_weight_function(reweightingonly=True, LHE=True)
    with open(self.treesample.LHEfile) as f:
      event = ""
      for line in f:
        if "<event>" not in line and not event: continue
        event += line
        if "</event>" in line:
          self.event = event = LHEEvent(event, bkg4l=self.bkg4l)
          sumofweights += weightfunction(self)
          event = ""
          i += 1
          if i % self.printevery == 0 or i == len(self):
            print i, "/", len(self), "(preliminary loop)"
            #break
    self.sumofweights = sumofweights
    self.event = None

  def __iter__(self):
    self.__i = 0
    self.__f = open(self.treesample.LHEfile)
    return super(LHEWrapper, self).__iter__()

  def next(self):
    if self.__i > len(self): assert False, (self.__i, len(self))
    self.__i += 1
    m = self.mela
    m.resetInputEvent()
    event = ""
    if self.__i % self.printevery == 0 or self.__i == len(self):
      print self.__i, "/", len(self)
      #raise StopIteration
    for line in self.__f:
      if "<event>" not in line and not event:
        continue
      event += line
      if "</event>" in line:
        self.event = event = LHEEvent(event, bkg4l=self.bkg4l)
        if not event.passcuts:
          return next(self)

        self.calcrecoMEs(*event.recomelaargs)
        return self
    raise StopIteration

  def ZZMass(self):
    return self.event.ZZMass #has to happen before doing all the MELA stuff, so can't use computeDecayAngles
  def ZZFlav(self):
    return self.event.ZZFlav

  def __del__(self):
    self.mela.resetInputEvent()

  def initlists(self):
    self.toaddtotree = [
      "ZZMass",
      "Z1Mass",
      "Z2Mass",
      "costhetastar",
      "costheta1",
      "costheta2",
      "Phi1",
      "Phi",
      "D_bkg",
      "D_0hplus_decay",
      "D_L1_decay",
      "D_L1int_decay",
      "D_L1Zg_decay",
      "D_L1Zgint_decay",
      "D_L1L1Zg_decay",
      "D_L1L1Zgint_decay",
      "D_eL_decay",
      "D_eR_decay",
      "D_eLint_decay",
      "D_eRint_decay",
      "D_eLeR_decay",
      "D_eLeRint_decay",
    ]
    self.exceptions = [
      "D_bkg_ResUp",
      "D_bkg_ResDown",
      "D_bkg_ScaleUp",
      "D_bkg_ScaleDown",
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
      "D_int_decay",

      "allsamples",
      "bkg4l",
      "cconstantforDbkg",
      "cconstantforD2jet",
      "cconstantforDHadWH",
      "cconstantforDHadZH",
      "checkfunctions",
      "dofa3stuff",
      "event",
      "exceptions",
      "hypothesis",
      "initlists",
      "initmela",
      "initsystematics",
      "initweightfunctions",
      "isalternate",
      "isbkg",
      "isdata",
      "isZX",
      "maxevent",
      "mela",
      "minevent",
      "next",
      "preliminaryloop",
      "printevery",
      "productionmode",
      "Show",
      "sumofweights",
      "toaddtotree",
      "toaddtotree_float",
      "toaddtotree_int",
      "treesample",
      "unblind",
      "xsec",
    ]
    proddiscriminants = [
      "D_0minus_{prod}",
      "D_CP_{prod}",
      "D_0hplus_{prod}",
      "D_int_{prod}",
      "D_L1_{prod}",
      "D_L1int_{prod}",
      "D_L1Zg_{prod}",
      "D_L1Zgint_{prod}",
      "D_0minus_{prod}decay",
      "D_0hplus_{prod}decay",
      "D_L1_{prod}decay",
      "D_L1Zg_{prod}decay",
    ]
    for prod in ("VBF", "HadVH"):
      self.exceptions += [_.format(prod=prod) for _ in proddiscriminants]

    self.toaddtotree_int = [
      "ZZFlav",
    ]

    self.toaddtotree_float = []

    for sample in self.allsamples:
      if sample == self.treesample.reweightingsample:
        self.toaddtotree.append(sample.weightname())
      else:
        self.exceptions.append(sample.weightname())
    if self.treesample.reweightingsample not in self.allsamples and not self.isdata:
      raise ValueError("{} not in allsamples!".format(self.treesample))

    fa3stuff = [
      "D_0minus_VBFdecay",
      "D_0minus_HadVHdecay",
      "D_0minus_decay",
      "D_CP_VBF",
      "D_CP_HadVH",
      "D_CP_decay",
    ]

    if self.dofa3stuff:
      for thing in fa3stuff:
        self.exceptions.remove(thing)
        self.toaddtotree_float.append(thing)
      self.toaddtotree_int.append("category_0P_or_0M")
    else:
      self.exceptions.append("category_0P_or_0M")

  def Show(self):
    print self.event

  allsamples = [
    ReweightingSample("ggH", "0+"),
    ReweightingSample("ggH", "L1"),
    ReweightingSample("ggH", "fL10.5"),
    ReweightingSample("ggH", "L1Zg"),
    ReweightingSample("ggH", "fL1Zg0.5"),
    ReweightingSample("ggH", "fL10.5fL1Zg0.5"),
    ReweightingSample("qqZZ"),
  ]

if __name__ == "__main__":
  from samples import Sample
  LHEWrapper(Sample("ggH", "L1Zg", "LHE_170509"))
