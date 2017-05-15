#!/usr/bin/env python

import os
import random
import sys

import ROOT

import CJLSTscripts
import config
from treewrapper import TreeWrapperBase
from utilities import cache, cache_instancemethod, product, requirecmsenv, tlvfromptetaphim

requirecmsenv(os.path.join(config.repositorydir, "CMSSW_8_0_20"))

from ZZMatrixElement.PythonWrapper.mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar

class LHEEvent(object):
  def __init__(self, event):
    self.eventstr = event

    #from https://github.com/hroskes/HiggsAnalysis-ZZMatrixElement/blob/b0ad77b/PythonWrapper/python/mela.py#L338-L366
    lines = event.split("\n")
    lines = [line for line in lines if not ("<event>" in line or "</event>" in line or not line.split("#")[0].strip())]
    nparticles, _, _, _, _, _ = lines[0].split()
    nparticles = int(nparticles)
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
          if mother1 != mother2 or mother1 is None:
            associated.append(SimpleParticle_t(line))
            break
          if ids[mother1] == 25:
            daughters.append(SimpleParticle_t(line))
            break
          mother2 = mother2s[mother1]
          mother1 = mother1s[mother1]

  @staticmethod
  def smear(particle):
    id, p = particle
    if abs(id) == 11:
      newid = id
      smearpt = config.smearptelectron
    elif abs(id) == 13:
      newid = id
      smearpt = config.smearptmuon
    elif 1 <= abs(id) <= 5 or abs(id) == 21:
      newid = 0
      smearpt = config.smearptjet
    else:
      assert False

    newp = tlvfromptetaphim(random.gauss(p.Pt(), smearpt), p.Eta(), p.Phi(), p.M())
    return SimpleParticle_t(newid, newp)

  @staticmethod
  def passparticlecuts(particle):
    id, p = particle
    if abs(id) == 11:
      if p.Pt() < 7 or abs(p.Eta()) > 2.5: return False
    elif abs(id) == 13:
      if p.Pt() < 5 or abs(p.Eta()) > 2.4: return False
    elif 1 <= abs(id) <= 5:
      if p.Pt() < 30 or abs(p.Eta()) > 4.7: return False
    return True

  @property
  @cache_instancemethod
  def recodaughters(self):
    return [self.smear(_) for _ in self.gendaughters if self.passparticlecuts(_)]

  @property
  @cache_instancemethod
  def recoassociated(self):
    return [self.smear(_) for _ in self.genassociated if self.passparticlecuts(_)]

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

class LHEWrapper(TreeWrapperBase):
  def __init__(self, treesample, minevent=0, maxevent=None):
    assert minevent == 0 and maxevent is None
    self.mela = self.initmela()
    self.event = None
    super(LHEWrapper, self).__init__(treesample, minevent, maxevent)

  @staticmethod
  @cache
  def initmela(*args, **kwargs):
    return Mela(*args, **kwargs)

  @cache_instancemethod
  def __len__(self):
    with open(self.treesample.LHEfile) as f:
      return f.read().count("</event>")

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
    for line in self.__f:
      if "<event>" not in line and not event:
        continue
      event += line
      if "</event>" in line:
        self.event = event = LHEEvent(event)
        if not event.passcuts:
          return next(self)

        m.setInputEvent(*event.recomelaargs)
        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
        m.ghg2 = m.ghz1 = 1
        self.M2g1_decay = m.computeP(True)

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
        m.ghg2 = 1; m.ghz1_prime2 = 1e4
        self.M2g1prime2_decay = m.computeP(True)

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
        m.ghg2 = 1; m.ghzgs1_prime2 = 1e4
        self.M2ghzgs1prime2_decay = m.computeP(True)

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
        m.ghg2 = m.g1 = 1; m.ghz1_prime2 = 1e4
        self.M2g1g1prime2_decay = m.computeP(True) - self.M2g1_decay - self.M2g1prime2_decay

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
        m.ghg2 = m.g1 = 1; m.ghzgs1_prime2 = 1e4
        self.M2g1ghzgs1prime2_decay = m.computeP(True) - self.M2g1_decay - self.M2ghzgs1prime2_decay

        m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZGG)
        m.ghg2 = 1; m.ghz1_prime2 = m.ghzgs1_prime2 = 1e4
        self.M2g1prime2ghzgs1prime2_decay = m.computeP(True) - self.M2g1prime2_decay - self.M2ghzgs1prime2_decay

        m.setProcess(TVar.HSMHiggs, TVar.JHUGen, TVar.ZZGG)
        self.p_m4l_SIG = m.computePM4l(TVar.SMSyst_None)

        m.setProcess(TVar.bkgZZ, TVar.JHUGen, TVar.ZZGG)
        self.p_m4l_BKG = m.computePM4l(TVar.SMSyst_None)

        m.setProcess(TVar.bkgZZ, TVar.MCFM, TVar.ZZQQB)
        self.M2qqZZ = m.computeP(True)

        self.cconstantforDbkg = CJLSTscripts.getDbkgConstant(self.ZZFlav(), self.ZZMass())

        self.M2g1prime2_decay /= 1e4**2
        self.M2ghzgs1prime2_decay /= 1e4**2
        self.M2g1g1prime2_decay /= 1e4
        self.M2g1ghzgs1prime2_decay /= 1e4
        self.M2g1prime2ghzgs1prime2_decay /= 1e4**2

        return self
    raise StopIteration

  def ZZMass(self):
    return self.event.ZZMass

  def ZZFlav(self):
    return self.event.ZZFlav

  def __del__(self):
    self.mela.resetInputEvent()

  def initlists(self):
    self.toaddtotree = [
      "ZZMass",
      "D_bkg",
      "D_L1_decay",
      "D_L1int_decay",
      "D_L1Zg_decay",
      "D_L1Zgint_decay",
      "D_L1L1Zg_decay",
      "D_L1L1Zgint_decay",
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
      "D_0hplus_decay",
      "D_int_decay",

      "cconstantforDbkg",
      "cconstantforD2jet",
      "cconstantforDHadWH",
      "cconstantforDHadZH",
      "checkfunctions",
      "event",
      "exceptions",
      "hypothesis",
      "initlists",
      "initmela",
      "isalternate",
      "isbkg",
      "isdata",
      "isZX",
      "maxevent",
      "mela",
      "minevent",
      "next",
      "printevery",
      "productionmode",
      "toaddtotree",
      "toaddtotree_int",
      "treesample",
      "Show",
      "unblind",
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

  def Show(self):
    print self.event

if __name__ == "__main__":
  from samples import Sample
  LHEWrapper(Sample("ggH", "L1Zg", "LHE_170509"))
