#!/usr/bin/env python
from abc import ABCMeta, abstractproperty
from collections import Counter
import itertools
from math import copysign, sqrt
import numpy
import os

import ROOT

import config
import constants
from enums import AlternateGenerator, analyses, Analysis, Flavor, flavors, purehypotheses, HffHypothesis, hffhypotheses, Hypothesis, MultiEnum, MultiEnumABCMeta, Production, ProductionMode, productions, PythiaSystematic, pythiasystematics
from utilities import cache, product
from weightshelper import WeightsHelper


class SampleBase(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def g1(self):
        pass
    @abstractproperty
    def g2(self):
        pass
    @abstractproperty
    def g4(self):
        pass
    @abstractproperty
    def g1prime2(self):
        pass
    @abstractproperty
    def ghzgs1prime2(self):
        pass
    @abstractproperty
    def ghg2(self):
        pass
    @abstractproperty
    def ghg4(self):
        pass
    @abstractproperty
    def kappa(self):
        pass
    @abstractproperty
    def kappa_tilde(self):
        pass

    def ValueError(self, functionname):
        return ValueError("invalid sample {} in function {}".format(self, functionname))

    @property
    def JHUxsec(self):
        if self.photoncut:
            JHUXSggH2L2la1 = constants.JHUXSggH2L2la1_photoncut
            JHUXSVBFa1 = constants.JHUXSVBFa1_photoncut
            JHUXSZHa1 = constants.JHUXSZHa1_photoncut
            JHUXSWHa1 = constants.JHUXSWHa1_photoncut
            assert self.g2 == self.g4 == self.g1prime2 == 0, self
        else:
            JHUXSggH2L2la1 = constants.JHUXSggH2L2la1
            JHUXSVBFa1 = constants.JHUXSVBFa1
            JHUXSZHa1 = constants.JHUXSZHa1
            JHUXSWHa1 = constants.JHUXSWHa1
            assert self.ghzgs1prime2 == 0
        if self.productionmode == "ggH":
            return (
                                JHUXSggH2L2la1*self.g1**2
                    + constants.JHUXSggH2L2la2*self.g2**2
                    + constants.JHUXSggH2L2la3*self.g4**2
                    + constants.JHUXSggH2L2lL1*self.g1prime2**2
                    + constants.JHUXSggH2L2lL1Zg*self.ghzgs1prime2**2
                    + constants.JHUXSggH2L2la1a2 * self.g1*self.g2 / constants.g2decay
                    + constants.JHUXSggH2L2la1a3 * self.g1*self.g4 / constants.g4decay
                    + constants.JHUXSggH2L2la1L1 * self.g1*self.g1prime2 / constants.g1prime2decay_gen
                    + constants.JHUXSggH2L2la1L1Zg * self.g1*self.ghzgs1prime2 / constants.ghzgs1prime2decay_gen
                   )

        if self.productionmode == "VBF":
            return (
                                JHUXSVBFa1 * self.g1**2
                    + constants.JHUXSVBFa2 * self.g2**2
                    + constants.JHUXSVBFa3 * self.g4**2
                    + constants.JHUXSVBFL1 * self.g1prime2**2
                    + constants.JHUXSVBFL1Zg * self.ghzgs1prime2**2
                    + constants.JHUXSVBFa1a2 * self.g1*self.g2 / constants.g2VBF
                    + constants.JHUXSVBFa1a3 * self.g1*self.g4 / constants.g4VBF
                    + constants.JHUXSVBFa1L1 * self.g1*self.g1prime2 / constants.g1prime2VBF_gen
                    + constants.JHUXSVBFa1L1Zg * self.g1*self.ghzgs1prime2 / constants.ghzgs1prime2VBF_gen
                   ) * (
                                JHUXSggH2L2la1 * self.g1**2
                    + constants.JHUXSggH2L2la2 * self.g2**2
                    + constants.JHUXSggH2L2la3 * self.g4**2
                    + constants.JHUXSggH2L2lL1 * self.g1prime2**2
                    + constants.JHUXSggH2L2lL1Zg * self.ghzgs1prime2**2
                    + constants.JHUXSggH2L2la1a2 * self.g1*self.g2 / constants.g2decay
                    + constants.JHUXSggH2L2la1a3 * self.g1*self.g4 / constants.g4decay
                    + constants.JHUXSggH2L2la1L1 * self.g1*self.g1prime2 / constants.g1prime2decay_gen
                    + constants.JHUXSggH2L2la1L1Zg * self.g1*self.ghzgs1prime2 / constants.ghzgs1prime2decay_gen
                   )

        if self.productionmode == "ZH":
            return (
                                JHUXSZHa1 * self.g1**2
                    + constants.JHUXSZHa2 * self.g2**2
                    + constants.JHUXSZHa3 * self.g4**2
                    + constants.JHUXSZHL1 * self.g1prime2**2
                    + constants.JHUXSZHL1Zg * self.ghzgs1prime2**2
                    + constants.JHUXSZHa1a2 * self.g1*self.g2 / constants.g2ZH
                    + constants.JHUXSZHa1a3 * self.g1*self.g4 / constants.g4ZH
                    + constants.JHUXSZHa1L1 * self.g1*self.g1prime2 / constants.g1prime2ZH_gen
                    + constants.JHUXSZHa1L1Zg * self.g1*self.ghzgs1prime2 / constants.ghzgs1prime2ZH_gen
                   ) * (
                                JHUXSggH2L2la1 * self.g1**2
                    + constants.JHUXSggH2L2la2 * self.g2**2
                    + constants.JHUXSggH2L2la3 * self.g4**2
                    + constants.JHUXSggH2L2lL1 * self.g1prime2**2
                    + constants.JHUXSggH2L2lL1Zg * self.ghzgs1prime2**2
                    + constants.JHUXSggH2L2la1a2 * self.g1*self.g2 / constants.g2decay
                    + constants.JHUXSggH2L2la1a3 * self.g1*self.g4 / constants.g4decay
                    + constants.JHUXSggH2L2la1L1 * self.g1*self.g1prime2 / constants.g1prime2decay_gen
                    + constants.JHUXSggH2L2la1L1Zg * self.g1*self.ghzgs1prime2 / constants.ghzgs1prime2decay_gen
                   )

        if self.productionmode == "WH":
            return (
                                JHUXSWHa1 * self.g1**2
                    + constants.JHUXSWHa2 * self.g2**2
                    + constants.JHUXSWHa3 * self.g4**2
                    + constants.JHUXSWHL1 * self.g1prime2**2
                    + constants.JHUXSWHL1Zg * self.ghzgs1prime2**2
                    + constants.JHUXSWHa1a2 * self.g1*self.g2 / constants.g2WH
                    + constants.JHUXSWHa1a3 * self.g1*self.g4 / constants.g4WH
                    + constants.JHUXSWHa1L1 * self.g1*self.g1prime2 / constants.g1prime2WH_gen
                    + constants.JHUXSWHa1L1Zg * self.g1*self.ghzgs1prime2 / constants.ghzgs1prime2WH_gen
                   ) * (
                                JHUXSggH2L2la1 * self.g1**2
                    + constants.JHUXSggH2L2la2 * self.g2**2
                    + constants.JHUXSggH2L2la3 * self.g4**2
                    + constants.JHUXSggH2L2lL1 * self.g1prime2**2
                    + constants.JHUXSggH2L2lL1Zg * self.ghzgs1prime2**2
                    + constants.JHUXSggH2L2la1a2 * self.g1*self.g2 / constants.g2decay
                    + constants.JHUXSggH2L2la1a3 * self.g1*self.g4 / constants.g4decay
                    + constants.JHUXSggH2L2la1L1 * self.g1*self.g1prime2 / constants.g1prime2decay_gen
                    + constants.JHUXSggH2L2la1L1Zg * self.g1*self.ghzgs1prime2 / constants.ghzgs1prime2decay_gen
                   )
        if self.productionmode == "HJJ":
            return (
                      constants.JHUXSHJJa2 * self.ghg2**2
                    + constants.JHUXSHJJa3 * self.ghg4**2
                    + constants.JHUXSHJJa2a3 * self.ghg2*self.ghg4 / constants.ghg4HJJ
                   ) * (
                                JHUXSggH2L2la1 * self.g1**2
                    + constants.JHUXSggH2L2la2 * self.g2**2
                    + constants.JHUXSggH2L2la3 * self.g4**2
                    + constants.JHUXSggH2L2lL1 * self.g1prime2**2
                    + constants.JHUXSggH2L2lL1Zg * self.ghzgs1prime2**2
                    + constants.JHUXSggH2L2la1a2 * self.g1*self.g2 / constants.g2decay
                    + constants.JHUXSggH2L2la1a3 * self.g1*self.g4 / constants.g4decay
                    + constants.JHUXSggH2L2la1L1 * self.g1*self.g1prime2 / constants.g1prime2decay_gen
                    + constants.JHUXSggH2L2la1L1Zg * self.g1*self.ghzgs1prime2 / constants.ghzgs1prime2decay_gen
                   )
        if self.productionmode == "ttH":
            return (
                      constants.JHUXSttHkappa * self.kappa**2
                    + constants.JHUXSttHkappatilde * self.kappa_tilde**2
                    + constants.JHUXSttHkappakappatilde * self.kappa*self.kappa_tilde / constants.kappa_tilde_ttH
                   ) * (
                                JHUXSggH2L2la1 * self.g1**2
                    + constants.JHUXSggH2L2la2 * self.g2**2
                    + constants.JHUXSggH2L2la3 * self.g4**2
                    + constants.JHUXSggH2L2lL1 * self.g1prime2**2
                    + constants.JHUXSggH2L2lL1Zg * self.ghzgs1prime2**2
                    + constants.JHUXSggH2L2la1a2 * self.g1*self.g2 / constants.g2decay
                    + constants.JHUXSggH2L2la1a3 * self.g1*self.g4 / constants.g4decay
                    + constants.JHUXSggH2L2la1L1 * self.g1*self.g1prime2 / constants.g1prime2decay_gen
                    + constants.JHUXSggH2L2la1L1Zg * self.g1*self.ghzgs1prime2 / constants.ghzgs1prime2decay_gen
                   )
        assert False

    @property
    def xsec(self):
        if self.productionmode == "ggH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSggH2L2l
            return constants.SMXSggH2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "VBF":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSVBF2L2l
            return constants.SMXSVBF2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "ZH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSZH2L2l
            return constants.SMXSZH2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "WH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSWH2L2l
            return constants.SMXSWH2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "WplusH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSWpH2L2l
        if self.productionmode == "WminusH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSWmH2L2l
        if self.productionmode == "HJJ":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return 1
            return constants.SMXSHJJ2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+", "Hff0+").JHUxsec
        if self.productionmode == "ttH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSttH2L2l
            return constants.SMXSttH2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+", "Hff0+").JHUxsec
        raise self.ValueError("xsec")

    @property
    def SMxsec(self):
        hffhypothesis = "Hff0+" if self.hffhypothesis else None
        return ReweightingSample(self.productionmode, "SM", hffhypothesis).xsec

    def __eq__(self, other):
        return all((
                   self.productionmode == other.productionmode,
                   self.photoncut == other.photoncut,
                   self.g1 == other.g1,
                   self.g2 == other.g2,
                   self.g4 == other.g4,
                   self.g1prime2 == other.g1prime2,
                   self.ghzgs1prime2 == other.ghzgs1prime2,
               )) and (
                   self.productionmode != "HJJ" or all((
                       self.ghg2 == other.ghg2,
                       self.ghg4 == other.ghg4,
                   ))
               ) and (
                   self.productionmode != "ttH" or all((
                       self.kappa == other.kappa,
                       self.kappa_tilde == other.kappa_tilde,
                   ))
               )
    def __ne__(self, other):
        return not self == other

    def linearcombinationofreweightingsamples(self, hypotheses):
        g1 = self.g1
        if   self.g2 == self.g1prime2 == self.ghzgs1prime2 == 0: analysis = Analysis("fa3");   gi = self.g4
        elif self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0: analysis = Analysis("fa2");   gi = self.g2
        elif self.g2 == self.g4       == self.ghzgs1prime2 == 0: analysis = Analysis("fL1");   gi = self.g1prime2
        elif self.g2 == self.g4       == self.g1prime2     == 0: analysis = Analysis("fL1Zg"); gi = self.ghzgs1prime2
        else: assert False

        if hypotheses == "templates":
            from templates import TemplatesFile
            templatesfile = TemplatesFile(analysis, str(self.productionmode).lower(), "2e2mu", "Untagged", config.productionsforcombine[0])
            hypotheses = [_.hypothesis for _ in templatesfile.templates()]

        basis = SampleBasis(hypotheses, self.productionmode, analysis)
        vectorofreweightingsamples = [ReweightingSample(self.productionmode, hypothesis) for hypothesis in hypotheses]
        invertedmatrix = basis.invertedmatrix
        maxpower = invertedmatrix.shape[0]-1
        vectorTofg1xgiy = numpy.matrix([[g1**(maxpower-i)*gi**i for i in range(maxpower+1)]])
        vectorTmatrix = vectorTofg1xgiy*invertedmatrix
        vectorTmatrixaslist = numpy.array(vectorTmatrix).flatten().tolist()
        return {
                reweightingsample: factor
                 for reweightingsample, factor in zip(vectorofreweightingsamples, vectorTmatrixaslist)
                 if factor
               }

    @property
    def MC_weight_terms(self):
        factors = []

        weightshelper = WeightsHelper(self)
        if weightshelper.weightprodstring is not None:
            prod = Counter()
            SMcoupling  = getattr(self, weightshelper.prodSMcoupling )
            BSMcoupling = getattr(self, weightshelper.prodBSMcoupling) / float(weightshelper.prodBSMcouplingvalue)
            prod[weightshelper.prodweightSM]  += SMcoupling**2            - SMcoupling * BSMcoupling
            prod[weightshelper.prodweightmix] += SMcoupling * BSMcoupling
            prod[weightshelper.prodweightBSM] += BSMcoupling**2           - SMcoupling * BSMcoupling
            factors.append([
                            (weightname, couplingsq)
                                  for weightname, couplingsq in prod.iteritems()
                                   if couplingsq
                                         and weightname is not None
                           ])

        if weightshelper.weightdecaystring is not None:
            decay = Counter()
            SMcoupling  = getattr(self, weightshelper.decaySMcoupling )
            BSMcoupling = getattr(self, weightshelper.decayBSMcoupling) / float(weightshelper.decayBSMcouplingvalue)
            decay[weightshelper.decayweightSM]  += SMcoupling**2            - SMcoupling * BSMcoupling
            decay[weightshelper.decayweightmix] += SMcoupling * BSMcoupling
            decay[weightshelper.decayweightBSM] += BSMcoupling**2           - SMcoupling * BSMcoupling
            factors.append([
                            (weightname, couplingsq)
                                  for weightname, couplingsq in decay.iteritems()
                                   if couplingsq
                                         and weightname is not None
                           ])

        return factors

    @property
    def MC_weight(self):
        if self.photoncut: raise NotImplementedError
        result = "*".join(
                          "("+
                          "+".join(
                                   "({}*{})".format(weightname, couplingsq)
                                         for weightname, couplingsq in factor
                                  )
                          +")" for factor in self.MC_weight_terms
                         )
        return result

    def get_MC_weight_function(self_sample, functionname=None, reweightingonly=False):
        if functionname is None:
            functionname = self_sample.weightname()

        strsample = str(self_sample)

        if self_sample.productionmode == "ggZZ":
            def MC_weight_function(self_tree):
                KFactor = self_tree.KFactor_QCD_ggZZ_Nominal
                return self_tree.overallEventWeight * self_tree.xsec * KFactor / self_tree.nevents

        elif self_sample.productionmode == "qqZZ":
            def MC_weight_function(self_tree):
                if hasattr(self_tree, "KFactor_EW_qqZZ"): #qqZZ->itself
                    KFactor = self_tree.KFactor_EW_qqZZ * self_tree.KFactor_QCD_qqZZ_M
                    return self_tree.overallEventWeight * self_tree.xsec * KFactor #/ self_tree.nevents
                                                                                   #Do NOT divide by nevents.  This will happen
                                                                                   #in creating the templates when we divide by
                                                                                   #sum(effectiveentries[==nevents])
                else:                                     #VBFbkg->qqZZ
                    result = self_tree.p_Gen_JJQCD_BKG_MCFM
                    if not reweightingonly and result != 0:
                        cutoff = self_tree.cutoffs[strsample]
                        if result > cutoff:
                            result = cutoff**2 / result
                        result *= (
                                     self_tree.overallEventWeight
                                   * self_tree.multiplyweight[strsample]
                                  )
                    return result


        elif self_sample.productionmode == "ZX":
            import ZX
            def MC_weight_function(self_tree):
                LepPt, LepEta, LepLepId = self_tree.tree.LepPt, self_tree.tree.LepEta, self_tree.tree.LepLepId
                return ZX.fakeRate13TeV(LepPt[2],LepEta[2],LepLepId[2]) * ZX.fakeRate13TeV(LepPt[3],LepEta[3],LepLepId[3])

        elif (self_sample.productionmode == "VBF bkg"
                   or hasattr(self_sample, "alternategenerator") and self_sample.alternategenerator in ("POWHEG", "MINLO", "NNLOPS")
                   or hasattr(self_sample, "pythiasystematic") and self_sample.pythiasystematic is not None):
            def MC_weight_function(self_tree):
                return self_tree.overallEventWeight * self_tree.xsec / self_tree.nevents

        elif self_sample.issignal:

            factors = self_sample.MC_weight_terms

            photoncut_decay = self_sample.photoncut
            photoncut_ZH = self_sample.photoncut and self_sample.productionmode == "ZH"
            photoncut_VBF = self_sample.photoncut and self_sample.productionmode == "VBF"

            from utilities import tlvfromptetaphim
            from itertools import product as cartesianproduct

            def MC_weight_function(self_tree):
                if photoncut_decay:
                    leptons = [(id, tlvfromptetaphim(pt, eta, phi, m)) for id, pt, eta, phi, m in zip(self_tree.LHEDaughterId, self_tree.LHEDaughterPt, self_tree.LHEDaughterEta, self_tree.LHEDaughterPhi, self_tree.LHEDaughterMass)]
                    for (id1, p1), (id2, p2) in cartesianproduct(leptons, leptons):
                        if id1 == -id2 and (p1+p2).M() < 4: return 0
                if photoncut_ZH:
                    Zdecay = [tlvfromptetaphim(pt, eta, phi, m) for pt, eta, phi, m in zip(self_tree.LHEAssociatedParticlePt, self_tree.LHEAssociatedParticleEta, self_tree.LHEAssociatedParticlePhi, self_tree.LHEAssociatedParticleMass)]
                    assert len(Zdecay) == 2
                    if (Zdecay[0]+Zdecay[1]).M() < 4: return 0
                if photoncut_VBF:
                    jets = [tlvfromptetaphim(pt, eta, phi, m) for pt, eta, phi, m in zip(self_tree.LHEAssociatedParticlePt, self_tree.LHEAssociatedParticleEta, self_tree.LHEAssociatedParticlePhi, self_tree.LHEAssociatedParticleMass)]
                    assert len(jets) == 2
                    if any(jet.Pt() < 15 for jet in jets): return 0
                result = product(
                                 sum(
                                     getattr(self_tree, weightname) * couplingsq
                                           for weightname, couplingsq in factor
                                    ) for factor in factors
                                )

                if not reweightingonly and result != 0:
                    cutoff = self_tree.cutoffs[strsample]
                    if result > cutoff:
                        result = cutoff**2 / result
                    result *= (
                                 self_tree.overallEventWeight
                               * self_tree.multiplyweight[strsample]
                              )
                return result

        else:
            raise ValueError("No MC weight function defined for {}".format(self_sample))

        MC_weight_function.__name__ = functionname
        return MC_weight_function

    def fai(self, productionmode, analysis, withdecay=False):
        from combinehelpers import mixturesign
        analysis = Analysis(analysis)
        hypothesis = analysis.purehypotheses[1]
        mixhypothesis = analysis.mixdecayhypothesis
        if productionmode == "VH":
            productionmodes = [ProductionMode("ZH"), ProductionMode("WH")]
        else:
            productionmodes = [ProductionMode(productionmode)]

        if productionmode == "ggH":
            withdecay = True

        power = 4 if productionmode != "ggH" and withdecay else 2
        numerator = None
        denominator = 0
        for h in purehypotheses:
            if not getattr(self, h.couplingname) and h != hypothesis: continue
            if analysis.photoncut and h == "0+": continue
            if not analysis.photoncut and h == "0+_photoncut": continue
            kwargs = {_.couplingname: 0 for _ in purehypotheses}
            kwargs[h.couplingname] = 1
            kwargs["photoncut"] = analysis.photoncut
            term = 0
            for productionmode in productionmodes:
                term += getattr(self, h.couplingname)**power * ArbitraryCouplingsSample(productionmode, **kwargs).xsec
            if not withdecay:
                term /= ArbitraryCouplingsSample("ggH", **kwargs).xsec
            denominator += term
            if h == hypothesis:
                numerator = term
        return copysign(
                        numerator / denominator,
                        #mixturesign: multiply by -1 for fL1
                        mixturesign(analysis)*getattr(self, hypothesis.couplingname)
                       )

class ArbitraryCouplingsSample(SampleBase):
    def __init__(self, productionmode, g1, g2, g4, g1prime2, ghzgs1prime2, ghg2=None, ghg4=None, kappa=None, kappa_tilde=None, photoncut=False):
        self.productionmode = ProductionMode(productionmode)
        self.__g1, self.__g2, self.__g4, self.__g1prime2, self.__ghzgs1prime2 = g1, g2, g4, g1prime2, ghzgs1prime2
        self.photoncut = photoncut or (ghzgs1prime2 != 0)
        if self.productionmode not in ("ggH", "VBF", "ZH", "WH", "HJJ", "ttH"):
            raise ValueError("Bad productionmode {}".format(self.productionmode))
        if sum(bool(g) for g in (g2, g4, g1prime2, ghzgs1prime2)) > 1:
            raise ValueError("Can only set at most one of g2, g4, g1prime2, ghzgs1prime2")

        if self.productionmode == "HJJ":
            self.__ghg2, self.__ghg4 = ghg2, ghg4
            if ghg2 is None or ghg4 is None:
                raise ValueError("Have to set ghg2 and ghg4 for HJJ")
        elif self.productionmode == "ggH":
            if ghg2 is None:
                ghg2 = 1
            if ghg4 is None:
                ghg4 = 0
            if ghg2 != 1:
                raise ValueError("ghg2 has to = 1 for {}".format(self.productionmode))
            if ghg4 != 0:
                raise ValueError("ghg4 has to = 0 for {}".format(self.productionmode))
        else:
            if ghg2 is not None or ghg4 is not None:
                raise ValueError("Can't set ghg2 or ghg4 for {}".format(self.productionmode))

        if self.productionmode == "ttH":
            self.__kappa, self.__kappa_tilde = kappa, kappa_tilde
            if kappa is None or kappa_tilde is None:
                raise ValueError("Have to set kappa and kappa_tilde for ttH")
        else:
            if kappa is not None or kappa_tilde is not None:
                raise ValueError("Can't set kappa or kappa_tilde for {}".format(self.productionmode))

    @property
    def g1(self):
        return self.__g1
    @property
    def g2(self):
        return self.__g2
    @property
    def g4(self):
        return self.__g4
    @property
    def g1prime2(self):
        return self.__g1prime2
    @property
    def ghzgs1prime2(self):
        return self.__ghzgs1prime2
    @property
    def ghg2(self):
        if self.__ghg2 is None: raise ValueError("{} doesn't have ghg2!".format(self))
        return self.__ghg2
    @property
    def ghg4(self):
        if self.__ghg4 is None: raise ValueError("{} doesn't have ghg4!".format(self))
        return self.__ghg4
    @property
    def kappa(self):
        return self.__kappa
    @property
    def kappa_tilde(self):
        return self.__kappa_tilde

    def __repr__(self):
        kwargs = ()
        if self.productionmode == "HJJ": kwargs += ("ghg2", "ghg4")
        if self.productionmode == "ttH": kwargs += ("kappa", "kappa_tilde")
        kwargs += ("g1", "g2", "g4", "g1prime2", "ghzgs1prime2", "photoncut")
        return "{}({}, {})".format(
                                   type(self).__name__,
                                   repr(self.productionmode.item.name),
                                   ", ".join(
                                             "{}={}".format(kwarg, getattr(self, kwarg))
                                             for kwarg in kwargs
                                            ),
                                  )

def samplewithfai(productionmode, analysis, fai, withdecay=False, productionmodeforfai=None):
    from combinehelpers import mixturesign, sigmaioversigma1
    analysis = Analysis(analysis)
    productionmode = ProductionMode(productionmode)

    if isinstance(withdecay, (basestring, ProductionMode)) and productionmodeforfai is None:
        raise TypeError(
                        "withdecay is a {thetype}!  Maybe you meant\n"
                        "  samplewithfai({productionmode!r}, {analysis!r}, {fai!r}, productionmodeforfai={withdecay!r})\n"
                        "or\n"
                        "  samplewithfai({productionmode!r}, {analysis!r}, {fai!r}, bool({withdecay!r}))"
                        .format(thetype=type(withdecay).__name__, productionmode=productionmode, analysis=analysis, fai=fai, withdecay=withdecay)
                       )

    if productionmodeforfai is None:
        productionmodeforfai = productionmode

    if productionmodeforfai == "ggH":
        withdecay = True

    if str(productionmodeforfai) == "VH":
        pass
    else:
        productionmodeforfai = ProductionMode(productionmodeforfai)

    kwargs = {coupling: 0 for coupling in ("g1", "g2", "g4", "g1prime2", "ghzgs1prime2")}
    kwargs["photoncut"] = analysis.photoncut
    power = (.25 if productionmodeforfai != "ggH" and withdecay else .5)
    kwargs["g1"] = (1-abs(fai))**power
    xsecratio = sigmaioversigma1(analysis, productionmodeforfai)
    if not withdecay:
        xsecratio /= sigmaioversigma1(analysis, "ggH")
    mixhypothesis = analysis.mixdecayhypothesis
    kwargs[analysis.couplingname] = copysign(
                                             (abs(fai) / xsecratio)**power,
                                             #mixturesign: multiply by -1 for fL1
                                             mixturesign(analysis)*fai
                                            )
    return ArbitraryCouplingsSample(productionmode, **kwargs)

class ReweightingSample(MultiEnum, SampleBase):
    __metaclass__ = MultiEnumABCMeta
    enumname = "reweightingsample"
    enums = [ProductionMode, Hypothesis, HffHypothesis, Flavor]

    def __eq__(self, other):
        try:
            return MultiEnum.__eq__(self, other)
        except BaseException as e1:
            try:
                return SampleBase.__eq__(self, other)
            except BaseException as e2:
                raise ValueError(
                                 "Can't compare {} == {}!  Trying to compare as MultiEnums:\n"
                                 "  {}\n"
                                 "Trying to compare couplings:\n"
                                 "  {}"
                                 .format(self, other, e1, e2)
                                )

    def check(self, *args):
        if self.productionmode is None:
            raise ValueError("No option provided for productionmode\n{}".format(args))
        elif self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "HJJ", "ttH"):
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hypothesis not in self.productionmode.validhypotheses:
                raise ValueError("{} hypothesis can't be {}\n{}".format(self.productionmode, self.hypothesis, args))
            if self.flavor is not None:
                raise ValueError("Flavor provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.productionmode in ("HJJ", "ttH"):
                if self.hffhypothesis is None:
                    raise ValueError("Hff hypothesis not provided for {} productionmode\n{}".format(self.productionmode, args))
            else:
                if self.hffhypothesis is not None:
                    raise ValueError("Hff hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
        elif self.productionmode in ("ggZZ", "VBF bkg"):
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hffhypothesis is not None:
                raise ValueError("Hff hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.flavor is None:
                raise ValueError("No flavor provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.productionmode == "VBF bkg" and self.flavor.hastaus:
                raise ValueError("No {} samples with taus\n{}".format(self.productionmode, args))
        elif self.productionmode == "qqZZ":
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hffhypothesis is not None:
                raise ValueError("Hff hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.flavor is not None and self.flavor.hastaus:
                raise ValueError("No {} samples with taus\n{}".format(self.productionmode, args))
        elif self.productionmode in ("ZX", "data"):
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hffhypothesis is not None:
                raise ValueError("Hff hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.flavor is not None:
                raise ValueError("Flavor provided for {} productionmode\n{}".format(self.productionmode, args))
        else:
            raise ValueError("Bad productionmode {}\n{}".format(self.productionmode, args))

    @property
    def photoncut(self):
        if self.hypothesis is None: return False
        return self.hypothesis.photoncut

    def reweightingsamples(self):
        if self.productionmode == "qqZZ" and self.flavor is not None:
            return [self]
        if self.productionmode in ("ggZZ", "qqZZ", "ZX") or self.alternategenerator in ("POWHEG", "MINLO", "NNLOPS") or self.pythiasystematic is not None:
            return [self]
        if self.productionmode == "VBF bkg":
            return [self, ReweightingSample("qqZZ", self.flavor)]
        if self.productionmode == "ttH":  #ttH reweighting doesn't work unfortunately
            return [
                    ReweightingSample(self.productionmode, hypothesis, hffhypothesis)
                            for hypothesis in self.productionmode.validhypotheses
                            for hffhypothesis in (hffhypotheses if hypothesis in ("0+", "0-", "fa30.5") else ("Hff0+",))
                            if hffhypothesis == self.hffhypothesis
                   ]
        if self.productionmode in ("ttH", "HJJ"):
            return [
                    ReweightingSample(self.productionmode, hypothesis, hffhypothesis)
                            for hypothesis in self.productionmode.validhypotheses
                            for hffhypothesis in (hffhypotheses if hypothesis in ("0+", "0-", "fa30.5") else ("Hff0+",))
                   ]
        elif self.productionmode in ("ggH", "VBF", "ZH", "WH"):
            return [ReweightingSample(self.productionmode, hypothesis) for hypothesis in self.productionmode.validhypotheses]
        elif self.productionmode == "data":
            return []
        raise self.ValueError("reweightingsamples")

    @property
    def isbkg(self):
        return self.productionmode.isbkg
    @property
    def issignal(self):
        return not self.isbkg

    def isZX(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBF bkg", "data", "VBF", "ZH", "WH", "ttH", "HJJ", "WplusH", "WminusH"):
            return False
        elif self.productionmode == "ZX":
            return True
        raise self.ValueError("isZX")

    def isdata(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBF bkg", "ZX", "VBF", "ZH", "WH", "ttH", "HJJ", "WplusH", "WminusH"):
            return False
        elif self.productionmode == "data":
            return True
        raise self.ValueError("isdata")

    def weightname(self):
        if self.productionmode == "ggH":
            if self.hypothesis == "0+":
                return "MC_weight_ggH_g1"
            elif self.hypothesis == "0+_photoncut":
                return "MC_weight_ggH_g1_photoncut"
            elif self.hypothesis == "a2":
                return "MC_weight_ggH_g2"
            elif self.hypothesis == "0-":
                return "MC_weight_ggH_g4"
            elif self.hypothesis == "L1":
                return "MC_weight_ggH_g1prime2"
            elif self.hypothesis == "L1Zg":
                return "MC_weight_ggH_ghzgs1prime2"
            elif self.hypothesis == "fa20.5":
                return "MC_weight_ggH_g1g2"
            elif self.hypothesis == "fa2dec-0.9":
                return "MC_weight_ggH_g1g2_minuspoint9"
            elif self.hypothesis == "fa30.5":
                return "MC_weight_ggH_g1g4"
            elif self.hypothesis == "fL10.5":
                return "MC_weight_ggH_g1g1prime2"
            elif self.hypothesis == "fL1Zg0.5":
                return "MC_weight_ggH_g1ghzgs1prime2"
            for a in analyses:
                if self.hypothesis == "{}-0.5".format(a):
                    return ReweightingSample(self.productionmode, "{}0.5".format(a)).weightname()+"_pi"
        elif self.productionmode in ("VBF", "ZH", "WH", "WplusH", "WminusH"):
            if self.hypothesis == "0+":
                return "MC_weight_{}_g1".format(self.productionmode)
            elif self.hypothesis == "0+_photoncut":
                return "MC_weight_{}_g1_photoncut".format(self.productionmode)
            elif self.hypothesis == "a2":
                return "MC_weight_{}_g2".format(self.productionmode)
            elif self.hypothesis == "0-":
                return "MC_weight_{}_g4".format(self.productionmode)
            elif self.hypothesis == "L1":
                return "MC_weight_{}_g1prime2".format(self.productionmode)
            elif self.hypothesis == "L1Zg":
                return "MC_weight_{}_ghzgs1prime2".format(self.productionmode)

            elif self.hypothesis == "fa2dec0.5":
                return "MC_weight_{}_g1g2_dec".format(self.productionmode)
            elif self.hypothesis == "fa3dec0.5":
                return "MC_weight_{}_g1g4_dec".format(self.productionmode)
            elif self.hypothesis == "fL1dec0.5":
                return "MC_weight_{}_g1g1prime2_dec".format(self.productionmode)
            elif self.hypothesis == "fL1Zgdec0.5":
                return "MC_weight_{}_g1ghzgs1prime2_dec".format(self.productionmode)

            elif self.hypothesis == "fa2prod0.5":
                return "MC_weight_{}_g1g2_prod".format(self.productionmode)
            elif self.hypothesis == "fa3prod0.5":
                return "MC_weight_{}_g1g4_prod".format(self.productionmode)
            elif self.hypothesis == "fL1prod0.5":
                return "MC_weight_{}_g1g1prime2_prod".format(self.productionmode)
            elif self.hypothesis == "fL1Zgprod0.5":
                return "MC_weight_{}_g1ghzgs1prime2_prod".format(self.productionmode)

            elif self.hypothesis == "fa2proddec0.5":
                return "MC_weight_{}_g1g2_proddec".format(self.productionmode)
            elif self.hypothesis == "fa3proddec0.5":
                return "MC_weight_{}_g1g4_proddec".format(self.productionmode)
            elif self.hypothesis == "fL1proddec0.5":
                return "MC_weight_{}_g1g1prime2_proddec".format(self.productionmode)
            elif self.hypothesis == "fL1Zgproddec0.5":
                return "MC_weight_{}_g1ghzgs1prime2_proddec".format(self.productionmode)

            elif self.hypothesis == "fa2dec-0.9":
                return "MC_weight_{}_g1g2_dec_minuspoint9".format(self.productionmode)

            for a in analyses:
                for b in "prod", "dec", "proddec":
                    if self.hypothesis == "{}{}-0.5".format(a, b):
                        return ReweightingSample(self.productionmode, "{}{}0.5".format(a, b)).weightname()+"_pi"

        elif self.productionmode in ("HJJ", "ttH"):
            if self.productionmode == "HJJ":
                if self.hffhypothesis == "Hff0+":
                    result = "MC_weight_HJJ_ghg2"
                elif self.hffhypothesis == "Hff0-":
                    result = "MC_weight_HJJ_ghg4"
                elif self.hffhypothesis == "fCP0.5":
                    result = "MC_weight_HJJ_ghg2ghg4"
            elif self.productionmode == "ttH":
                if self.hffhypothesis == "Hff0+":
                    result = "MC_weight_ttH_kappa"
                elif self.hffhypothesis == "Hff0-":
                    result = "MC_weight_ttH_kappatilde"
                elif self.hffhypothesis == "fCP0.5":
                    result = "MC_weight_ttH_kappakappatilde"
            result += "_"
            if self.hypothesis == "0+":
                return result + "g1"
            elif self.hypothesis == "0+_photoncut":
                return result + "g1_photoncut"
            elif self.hypothesis == "a2":
                return result + "g2"
            elif self.hypothesis == "0-":
                return result + "g4"
            elif self.hypothesis == "L1":
                return result + "g1prime2"
            elif self.hypothesis == "L1Zg":
                return result + "ghzgs1prime2"
            elif self.hypothesis == "fa20.5":
                return result + "g1g2"
            elif self.hypothesis == "fa30.5":
                return result + "g1g4"
            elif self.hypothesis == "fL10.5":
                return result + "g1g1prime2"
            elif self.hypothesis == "fL1Zg0.5":
                return result + "g1ghzgs1prime2"
            elif self.hypothesis == "fa2dec-0.9":
                return result + "g1g2_minuspoint9"
            for a in analyses:
                if self.hypothesis == "{}-0.5".format(a):
                    return ReweightingSample(self.productionmode, self.hffhypothesis, "{}0.5".format(a)).weightname()+"_pi"
            
        elif self.productionmode == "ggZZ":
            return "MC_weight_ggZZ"
        elif self.productionmode == "qqZZ":
            return "MC_weight_qqZZ"
        elif self.productionmode == "VBF bkg":
            return "MC_weight_VBFbkg"
        elif self.productionmode == "ZX":
            return "MC_weight_ZX"
        raise self.ValueError("weightname")

    @property
    def MC_weight(self):
        return self.weightname()

    def TDirectoryname(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBFbkg", "data", "VBF", "ZH", "WH", "ttH", "HJJ", "WplusH", "WminusH") or self.productionmode == "ZX" and not config.usedata:
            return "ZZTree"
        if self.productionmode == "ZX":
            return "CRZLLTree"
        raise self.ValueError("TDirectoryname")

    def color(self):
        if self.productionmode == "ggH":
            if self.hypothesis == "0+":
                return 1
            elif self.hypothesis == "0-":
                return 2
            elif self.hypothesis == "fa30.5":
                return ROOT.kRed+3
            elif self.hypothesis == "a2":
                return ROOT.kCyan
            elif self.hypothesis == "fa20.5":
                return 4
            elif self.hypothesis == "L1":
                return 3
            elif self.hypothesis == "fL10.5":
                return ROOT.kGreen+3
        raise self.ValueError("color")

    def onlyweights(self):
        """True if this sample is not needed for making templates,
           and only the weight and ZZMass should be recorded in the tree"""
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "ggZZ", "qqZZ", "VBF bkg", "ZX", "data", "HJJ", "ttH", "WplusH", "WminusH"):
            return False
        if self.productionmode in ():
            return True
        raise self.ValueError("onlyweights")

    def weightingredients(self):
        """only needs to be defined if self.onlyweights() is True.
           gives a list of variable names that are used to calculate
           the weight"""
        raise self.ValueError("weightingredients")

    @property
    def g1(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "ttH", "HJJ"):
            if self.hypothesis in (
                                   ["0+", "0+_photoncut"]
                                   + [
                                      "{}{}{}0.5".format(a, b, c) for a in ("fa2", "fa3", "fL1", "fL1Zg")
                                                                  for b in ("dec", "prod", "proddec")
                                                                  for c in ("+", "-")
                                     ] + ["fa2dec-0.9"]
                                  ):
                return 1
            if self.hypothesis in ("a2", "0-", "L1", "L1Zg"):
                return 0
        raise self.ValueError("g1")

    @property
    def g2(self):
        if self.hypothesis in (
                               ["0+", "0+_photoncut", "0-", "L1", "L1Zg"]
                             + ["{}{}{}0.5".format(a, b, c) for a in ("fa3", "fL1", "fL1Zg") for b in ("dec", "prod", "proddec") for c in ("+", "-")]
                              ):
            return 0
        if self.hypothesis == "a2":
            if self.productionmode in ("ggH", "ttH"):
                return constants.g2decay
            else:
                return 1

        if self.hypothesis == "fa2dec0.5":
            return constants.g2decay

        if self.hypothesis == "fa2dec-0.9":
            return -constants.g2decay * 3

        if self.productionmode == "VBF":
            if self.hypothesis == "fa2prod0.5":
                return constants.g2VBF
            if self.hypothesis == "fa2proddec0.5":
                return sqrt(constants.g2VBF*constants.g2decay)

        if self.productionmode == "ZH":
            if self.hypothesis == "fa2prod0.5":
                return constants.g2ZH
            if self.hypothesis == "fa2proddec0.5":
                return sqrt(constants.g2ZH*constants.g2decay)

        if self.productionmode == "WH":
            if self.hypothesis == "fa2prod0.5":
                return constants.g2WH
            if self.hypothesis == "fa2proddec0.5":
                return sqrt(constants.g2WH*constants.g2decay)

        if self.hypothesis == "fa2prod-0.5": return -ReweightingSample(self.productionmode, "fa2prod0.5", self.hffhypothesis).g2
        if self.hypothesis == "fa2dec-0.5": return -ReweightingSample(self.productionmode, "fa2dec0.5", self.hffhypothesis).g2
        if self.hypothesis == "fa2proddec-0.5": return -ReweightingSample(self.productionmode, "fa2proddec0.5", self.hffhypothesis).g2

        raise self.ValueError("g2")

    @property
    def g4(self):
        if self.hypothesis in (
                               ["0+", "0+_photoncut", "a2", "L1", "L1Zg"]
                             + ["{}{}{}0.5".format(a, b, c) for a in ("fa2", "fL1", "fL1Zg") for b in ("dec", "prod", "proddec") for c in "+-"]
                             + ["fa2dec-0.9"]
                              ):
            return 0
        if self.hypothesis == "0-":
            if self.productionmode in ("ggH", "ttH"):
                return constants.g4decay
            else:
                return 1

        if self.hypothesis == "fa3dec0.5":
            return constants.g4decay

        if self.productionmode == "VBF":
            if self.hypothesis == "fa3prod0.5":
                return constants.g4VBF
            if self.hypothesis == "fa3proddec0.5":
                return sqrt(constants.g4VBF*constants.g4decay)

        if self.productionmode == "ZH":
            if self.hypothesis == "fa3prod0.5":
                return constants.g4ZH
            if self.hypothesis == "fa3proddec0.5":
                return sqrt(constants.g4ZH*constants.g4decay)

        if self.productionmode == "WH":
            if self.hypothesis == "fa3prod0.5":
                return constants.g4WH
            if self.hypothesis == "fa3proddec0.5":
                return sqrt(constants.g4WH*constants.g4decay)

        if self.hypothesis == "fa3prod-0.5": return -ReweightingSample(self.productionmode, "fa3prod0.5", self.hffhypothesis).g4
        if self.hypothesis == "fa3dec-0.5": return -ReweightingSample(self.productionmode, "fa3dec0.5", self.hffhypothesis).g4
        if self.hypothesis == "fa3proddec-0.5": return -ReweightingSample(self.productionmode, "fa3proddec0.5", self.hffhypothesis).g4

        raise self.ValueError("g4")

    @property
    def g1prime2(self):
        if self.hypothesis in (
                               ["0+", "0+_photoncut", "a2", "0-", "L1Zg"]
                             + ["{}{}{}0.5".format(a, b, c) for a in ("fa2", "fa3", "fL1Zg") for b in ("dec", "prod", "proddec") for c in "+-"]
                             + ["fa2dec-0.9"]
                              ):
            return 0

        if self.hypothesis == "L1":
            if self.productionmode in ("ggH", "ttH"):
                return constants.g1prime2decay_gen
            else:
                return 1

        if self.hypothesis == "fL1dec0.5":
            return constants.g1prime2decay_gen

        if self.productionmode == "VBF":
            if self.hypothesis == "fL1prod0.5":
                return constants.g1prime2VBF_gen
            if self.hypothesis == "fL1proddec0.5":
                return -sqrt(constants.g1prime2VBF_gen*constants.g1prime2decay_gen)

        if self.productionmode == "ZH":
            if self.hypothesis == "fL1prod0.5":
                return constants.g1prime2ZH_gen
            if self.hypothesis == "fL1proddec0.5":
                return -sqrt(constants.g1prime2ZH_gen*constants.g1prime2decay_gen)

        if self.productionmode == "WH":
            if self.hypothesis == "fL1prod0.5":
                return constants.g1prime2WH_gen
            if self.hypothesis == "fL1proddec0.5":
                return -sqrt(constants.g1prime2WH_gen*constants.g1prime2decay_gen)

        if self.hypothesis == "fL1prod-0.5": return -ReweightingSample(self.productionmode, "fL1prod0.5", self.hffhypothesis).g1prime2
        if self.hypothesis == "fL1dec-0.5": return -ReweightingSample(self.productionmode, "fL1dec0.5", self.hffhypothesis).g1prime2
        if self.hypothesis == "fL1proddec-0.5": return -ReweightingSample(self.productionmode, "fL1proddec0.5", self.hffhypothesis).g1prime2

        raise self.ValueError("g1prime2")

    @property
    def ghzgs1prime2(self):
        if self.hypothesis in (
                               ["0+", "0+_photoncut", "a2", "0-", "L1"]
                             + ["{}{}{}0.5".format(a, b, c) for a in ("fa2", "fa3", "fL1") for b in ("dec", "prod", "proddec") for c in "+-"]
                             + ["fa2dec-0.9"]
                              ):
            return 0

        if self.hypothesis == "L1Zg":
            if self.productionmode in ("ggH", "ttH"):
                return constants.ghzgs1prime2decay_gen
            else:
                return 1

        if self.hypothesis == "fL1Zgdec0.5":
            return constants.ghzgs1prime2decay_gen

        if self.productionmode == "VBF":
            if self.hypothesis == "fL1Zgprod0.5":
                return constants.ghzgs1prime2VBF_gen
            if self.hypothesis == "fL1Zgproddec0.5":
                return -sqrt(constants.ghzgs1prime2VBF_gen*constants.ghzgs1prime2decay_gen)

        if self.productionmode == "ZH":
            if self.hypothesis == "fL1Zgprod0.5":
                return constants.ghzgs1prime2ZH_gen
            if self.hypothesis == "fL1Zgproddec0.5":
                return -sqrt(constants.ghzgs1prime2ZH_gen*constants.ghzgs1prime2decay_gen)

        if self.productionmode == "WH":        #use the ZH numbers.  Shouldn't matter whatever we use.
            if self.hypothesis == "fL1Zgprod0.5":
                return constants.ghzgs1prime2ZH_gen
            if self.hypothesis == "fL1Zgproddec0.5":
                return -sqrt(constants.ghzgs1prime2ZH_gen*constants.ghzgs1prime2decay_gen)

        if self.hypothesis == "fL1Zgprod-0.5": return -ReweightingSample(self.productionmode, "fL1Zgprod0.5", self.hffhypothesis).ghzgs1prime2
        if self.hypothesis == "fL1Zgdec-0.5": return -ReweightingSample(self.productionmode, "fL1Zgdec0.5", self.hffhypothesis).ghzgs1prime2
        if self.hypothesis == "fL1Zgproddec-0.5": return -ReweightingSample(self.productionmode, "fL1Zgproddec0.5", self.hffhypothesis).ghzgs1prime2

        raise self.ValueError("ghzgs1prime2")

    @property
    def ghg2(self):
        if self.productionmode == "ggH":
            return 1
        if self.productionmode == "HJJ":
            if self.hffhypothesis in ("Hff0+", "fCP0.5"):
                return 1
            if self.hffhypothesis == "Hff0-":
                return 0
        raise self.ValueError("ghg2")

    @property
    def ghg4(self):
        if self.productionmode == "ggH":
            return 0
        if self.productionmode == "HJJ":
            if self.hffhypothesis == "Hff0+":
                return 0
            if self.hffhypothesis == "Hff0-":
                return 1
            if self.hffhypothesis == "fCP0.5":
                return constants.ghg4HJJ
        raise self.ValueError("ghg4")

    @property
    def kappa(self):
        if self.productionmode == "ttH":
            if self.hffhypothesis in ("Hff0+", "fCP0.5"):
                return 1
            if self.hffhypothesis == "Hff0-":
                return 0
        raise self.ValueError("kappa")

    @property
    def kappa_tilde(self):
        if self.productionmode == "ttH":
            if self.hffhypothesis == "Hff0+":
                return 0
            if self.hffhypothesis == "Hff0-":
                return 1
            if self.hffhypothesis == "fCP0.5":
                return constants.kappa_tilde_ttH

class ReweightingSamplePlus(ReweightingSample):
    enums = [ReweightingSample, AlternateGenerator, PythiaSystematic]
    enumname = "reweightingsampleplus"

    def check(self, *args):
        if (
            self.alternategenerator == "POWHEG"
            and (
                 self.hypothesis != "0+"
                 or self.productionmode not in ("VBF", "ZH", "WplusH", "WminusH", "ttH")
                )
           ) or (
            self.alternategenerator in ("MINLO", "NNLOPS")
            and (
                 self.hypothesis != "0+"
                 or self.productionmode != "ggH"
                )
           ):
            raise ValueError("No {} {} sample produced with {}\n{}".format(self.productionmode, self.hypothesis, self.alternategenerator, args))

        if (
            self.alternategenerator == "POWHEG"
            and (
                 self.hffhypothesis != "Hff0+"
                 and self.hffhypothesis is not None
                )
           ):
            raise ValueError("No {} {} sample produced with {}\n{}".format(self.productionmode, self.hffhypothesis, self.alternategenerator, args))

        if self.pythiasystematic is not None:
            if (
                self.hypothesis != "0+"
                or self.productionmode not in ("ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH")
               ):
                raise ValueError("No {} {} sample produced with pythia {}\n{}".format(self.productionmode, self.hypothesis, self.pythiasystematic, args))
            if self.productionmode != "ggH" and self.alternategenerator != "POWHEG":
                raise ValueError("{} sample with pythia {} is produced with POWHEG\n{}".format(self.productionmode, self.pythiasystematic, args))

        super(ReweightingSamplePlus, self).check(*args)

    def weightname(self):
        result = super(ReweightingSamplePlus, self).weightname()
        if self.alternategenerator:
            result += "_" + str(self.alternategenerator)
        if self.pythiasystematic:
            result += "_" + str(self.pythiasystematic)
        return result

class Sample(ReweightingSamplePlus):
    enums = [ReweightingSamplePlus, Production]

    def check(self, *args):
        if self.production is None:
            raise ValueError("No option provided for production\n{}".format(args))

        if self.hypothesis is not None and self.hypothesis not in self.productionmode.generatedhypotheses:
            raise ValueError("No {} sample produced with hypothesis {}!\n{}".format(self.productionmode, self.hypothesis, args))

        if self.productionmode == "qqZZ" and self.flavor is not None:
            raise ValueError("No {} sample produced with separate flavors!\n{}".format(self.productionmode, args))

        if self.productionmode in ("WplusH", "WminusH") and self.alternategenerator != "POWHEG":
            raise ValueError("Separate {} sample is produced with POWHEG.  Maybe you meant to specify POWHEG, or WH?\n{}".format(self.productionmode, args))

        if self.pythiasystematic is None and self.alternategenerator == "MINLO":
            raise ValueError("MINLO nominal sample is not ready yet :(\n{}".format(args))

        if self.pythiasystematic is not None and self.alternategenerator == "NNLOPS":
            raise ValueError("No NNLOPS samples with systematics!\n{}".format(args))

        super(Sample, self).check(*args)

    def CJLSTmaindir(self):
        if self.productionmode == "ggH":
            return self.production.CJLSTdir_anomalous()
        if self.productionmode == "VBF":
            return self.production.CJLSTdir_anomalous_VBF()
        if self.productionmode in ("ZH", "WH"):
            return self.production.CJLSTdir_anomalous_VH()
        if self.productionmode in ("data", "ZX"):
            return self.production.CJLSTdir_data()
        return self.production.CJLSTdir()

    def CJLSTdirname(self):
        if self.alternategenerator is not None:
            if self.alternategenerator == "POWHEG":
                if self.productionmode in ("VBF", "ZH", "WplusH", "WminusH", "ttH"):
                    s = str(self.productionmode)
                    if self.productionmode == "VBF": s = "VBFH"
                    if self.hypothesis == "0+": result = "{}125".format(s)
                if self.pythiasystematic is not None:
                    result += self.pythiasystematic.appendname
                return result
        if self.alternategenerator in ("MINLO", "NNLOPS") or self.pythiasystematic is not None:
            if self.productionmode == "ggH":
                result = "ggH125"
                if self.pythiasystematic is not None:
                    result += self.pythiasystematic.appendname
                if self.alternategenerator == "MINLO":
                    result += "_minloHJJ"
                return result
            raise self.ValueError("CJLSTdirname")
        if self.productionmode in ("ggH", "VBF", "ZH", "WH"):
            s = str(self.productionmode)
            if self.productionmode == "VBF": s = "VBFH"
            if self.hypothesis == "0+": result = "{}0PM_M125".format(s)
            if self.hypothesis == "a2": result = "{}0PH_M125".format(s)
            if self.hypothesis == "0-": result = "{}0M_M125".format(s)
            if self.hypothesis == "L1": result = "{}0L1_M125".format(s)
            if self.hypothesis in ("fa20.5", "fa2prod0.5"): result = "{}0PHf05ph0_M125".format(s)
            if self.hypothesis in ("fa30.5", "fa3prod0.5"): result = "{}0Mf05ph0_M125".format(s)
            if self.hypothesis in ("fL10.5", "fL1prod0.5"): result = "{}0L1f05ph0_M125".format(s)
            return result
        if self.productionmode in ("HJJ", "ttH"):
            s = str(self.productionmode)
            if self.hffhypothesis == "Hff0+": return "{}0PM_M125".format(s)
            if self.hffhypothesis == "Hff0-": return "{}0M_M125".format(s)
            if self.hffhypothesis == "fCP0.5": return "{}0Mf05ph0_M125".format(s)
        if self.productionmode == "ggZZ":
            return "ggTo{}_Contin_MCFM701".format(self.flavor)
        if self.productionmode == "VBF bkg":
            return "VBFTo{}JJ_Contin_phantom128".format(self.flavor)
        if self.productionmode == "qqZZ":
            return "ZZTo4l"
        if self.productionmode == "ZX" or self.productionmode == "data":
            return "AllData"
        raise self.ValueError("CJLSTdirname")

    def CJLSTfile(self):
        return os.path.join(self.CJLSTmaindir(), self.CJLSTdirname(), "ZZ4lAnalysis.root")

    def withdiscriminantsfile(self):
        result = os.path.join(config.repositorydir, "step3_withdiscriminants", "{}.root".format(self).replace(" ", ""))
        return result
        raise self.ValueError("withdiscriminantsfile")

    @staticmethod
    def effectiveentries(reweightfrom, reweightto):
        from utilities import tfiles
        if reweightto.productionmode in ("ggZZ", "VBF bkg", "ZX", "data"):
            assert reweightfrom.reweightingsample == reweightto
            return 1
        f = tfiles[reweightfrom.withdiscriminantsfile()]
        t = f.effectiveentries
        assert t.GetEntries() == 1, "{}???".format(t.GetEntries())
        t.GetEntry(0)
        try:
            return getattr(t, reweightto.weightname())
        except:
            t.Show()
            raise

    @property
    def xsec(self):
        try:
            return super(Sample, self).xsec
        except ValueError:
            t = ROOT.TChain("candTree")
            t.Add(self.withdiscriminantsfile())
            t.GetEntry(1)
            result = t.xsec
            if t.xsec <= 0:
                raise ValueError("Can't get xsec for {!r}".format(self))
            return result

class SampleBasis(MultiEnum):
    enums = [ProductionMode, Analysis]
    def __init__(self, hypotheses, *args):
        self.hypotheses = [Hypothesis(_) for _ in hypotheses]
        super(SampleBasis, self).__init__(*args)

    def check(self, *args):
        args = (self.hypotheses,)+args
        if self.productionmode in ("ggH", "ttH"):
            dimension = 3
        elif self.productionmode in ("VBF", "WH", "ZH"):
            dimension = 5
        else:
            raise ValueError("Bad productionmode {}\n{}".format(self.productionmode, args))

        if len(self.hypotheses) != dimension:
            raise ValueError("Wrong number of hypotheses ({}, should be {})\n{}".format(len(self.hypotheses), dimension, args))
        if len(set(self.hypotheses)) != len(self.hypotheses):
            raise ValueError("Duplicate hypothesis\n{}".format(args))

    @property
    @cache
    def matrix(self):
        hffhypothesis = None
        if self.productionmode == "ttH": hffhypothesis = "Hff0+"
        dimension = len(self.hypotheses)
        maxpower = dimension-1
        samples = [ReweightingSample(self.productionmode, _, hffhypothesis) for _ in self.hypotheses]
        return numpy.matrix(
                            [
                             [
                              sample.g1**(maxpower-i) * getattr(sample, self.analysis.couplingname)**i
                                  for i in range(dimension)
                             ]
                                 for sample in samples
                            ]
                           )
    @property
    @cache
    def invertedmatrix(self):
        return self.matrix.I

def allsamples():
    for production in productions:
        for productionmode in "ggH", "VBF", "ZH", "WH":
            for hypothesis in ProductionMode(productionmode).generatedhypotheses:
                yield Sample(productionmode, hypothesis, production)
        for hypothesis in hffhypotheses:
            yield Sample("HJJ", hypothesis, "0+", production)
            yield Sample("ttH", hypothesis, "0+", production)
        for flavor in flavors:
            yield Sample("ggZZ", flavor, production)
            if not flavor.hastaus:
                yield Sample("VBF bkg", flavor, production)
        yield Sample("ggH", "0+", "NNLOPS", production)
        for productionmode in "VBF", "ZH", "WplusH", "WminusH":
            yield Sample(productionmode, "0+", "POWHEG", production)
        yield Sample("ttH", "Hff0+", "0+", "POWHEG", production)
        for systematic in pythiasystematics:
            for productionmode in "VBF", "ZH", "WplusH", "WminusH":
                yield Sample(productionmode, "0+", "POWHEG", production, systematic)
            yield Sample("ttH", "Hff0+", "0+", "POWHEG", production, systematic)
            yield Sample("ggH", "0+", production, systematic)
            yield Sample("ggH", "0+", "MINLO", production, systematic)
        yield Sample("qqZZ", production)
        yield Sample("ZX", production)
        yield Sample("data", production)

def xcheck():
    CJLST = {_: _.CJLSTfile() for _ in allsamples()}
    withdiscs = {_: _.withdiscriminantsfile() for _ in allsamples()}

    for (k1, v1), (k2, v2) in itertools.product(CJLST.iteritems(), CJLST.iteritems()):
        if k1.productionmode == "ZX" and k2.productionmode == "data" or k2.productionmode == "ZX" and k1.productionmode == "data": continue
        if k1 != k2 and v1 == v2:
            raise ValueError("CJLST files for {} and {} are the same:\n{}".format(k1, k2, v1))

    for (k1, v1), (k2, v2) in itertools.product(withdiscs.iteritems(), withdiscs.iteritems()):
        if k1 != k2 and v1 == v2:
            raise ValueError("with discriminants files for {} and {} are the same:\n{}".format(k1, k2, v1))

xcheck()
del xcheck
