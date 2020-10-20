#!/usr/bin/env python
import os, tempfile
from abc import ABCMeta, abstractproperty
from collections import Counter
from itertools import combinations, permutations, product as cartesianproduct
from math import copysign, sqrt

import numpy

import ROOT

import config
import constants
import eft
from enums import AlternateGenerator, AlternateWeight, analyses, Analysis, Extension, Flavor, flavors, ggZZoffshellproductionmodes, purehypotheses, HffHypothesis, hffhypotheses, Hypothesis, MultiEnum, MultiEnumABCMeta, Production, ProductionMode, productions, PythiaSystematic, pythiasystematics
from extendedcounter import ExtendedCounter
from utilities import cache, cache_file, deprecate, generatortolist, product, TFile, tlvfromptetaphim, WriteOnceDict
from weightshelper import WeightsHelper

class SumOfSamplesBase(object):
  __metaclass__ = ABCMeta
  @abstractproperty
  def samplesandfactors(self): pass
  def __add__(self, other):
    return SumOfSamples(self.samplesandfactors + other.samplesandfactors)
  def __neg__(self):
    return SumOfSamples(-self.samplesandfactors)
  def __sub__(self, other):
    return SumOfSamples(self.samplesandfactors - other.samplesandfactors)
  def __mul__(self, scalar):
    return SumOfSamples(self.samplesandfactors * scalar)
  def __rmul__(self, scalar):
    return self * scalar
  def __div__(self, scalar):
    return SumOfSamples(self.samplesandfactors / scalar)
  @property
  def MC_weight_terms_expanded(self):
    terms = [(weight, weightfactor*factor) for s, factor in self.samplesandfactors.iteritems() for weight, weightfactor in s.MC_weight_terms_expanded]
    c = ExtendedCounter()

    for weight, factor in terms:
      c[weight] += factor

    maxfactor = max(abs(factor) for factor in c.values())
    if maxfactor:
      for weight, factor in c.items():
        if abs(factor/maxfactor) < 1e-10: del c[weight]

    return [(weight, factor) for weight, factor in c.iteritems()]

  @property
  def MC_weight_terms(self):
    return [self.MC_weight_terms_expanded]  #for compatibility with SampleBase below

  @property
  def MC_weight(self):
    result = "*".join(
                      "("+
                      "+".join(
                               "({}*{})".format(weightname, couplingsq)
                                   for weightname, couplingsq in sorted(factor)
                              )
                      +")" if factor else "0"
                      for factor in self.MC_weight_terms
                     )
    return result


class SumOfSamples(SumOfSamplesBase):
  def __init__(self, samplesandfactors=None):
    if samplesandfactors is None: samplesandfactors = ExtendedCounter()
    self.__samplesandfactors = samplesandfactors
  @property
  def samplesandfactors(self): return self.__samplesandfactors

class SampleBase(SumOfSamplesBase):
    @property
    def samplesandfactors(self): return ExtendedCounter({self: 1})
    @abstractproperty
    def ghz1(self):
        pass
    @abstractproperty
    def ghz2(self):
        pass
    @abstractproperty
    def ghz4(self):
        pass
    @abstractproperty
    def ghz1prime2(self):
        pass

    @abstractproperty
    def ghw1(self):
        pass
    @abstractproperty
    def ghw2(self):
        pass
    @abstractproperty
    def ghw4(self):
        pass
    @abstractproperty
    def ghw1prime2(self):
        pass

    @abstractproperty
    def ghzgs1prime2(self):
        pass

    @property
    def g1(self):
      if self.ghz1 != self.ghw1: raise ValueError("ghz1 = {self.ghz1}, ghw1 = {self.ghw1}".format(self=self))
      return self.ghz1
    @property
    def g2(self):
      if self.ghz2 != self.ghw2: raise ValueError("ghz2 = {self.ghz2}, ghw2 = {self.ghw2}".format(self=self))
      return self.ghz2
    @property
    def g4(self):
      if self.ghz4 != self.ghw4: raise ValueError("ghz4 = {self.ghz4}, ghw4 = {self.ghw4}".format(self=self))
      return self.ghz4
    @property
    def g1prime2(self):
      if self.ghz1prime2 != self.ghw1prime2: raise ValueError("ghz1prime2 = {self.ghz1prime2}, ghw1prime2 = {self.ghw1prime2}".format(self=self))
      return self.ghz1prime2

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
    @abstractproperty
    def pdf(self):
        pass

    def ValueError(self, functionname):
        return ValueError("invalid sample {} in function {}".format(self, functionname))

    @property
    def nominalJHUxsec(self):
        result = constants.JHUXSHZZ2e2mua1.nominal_value
        if self.g2 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
            result *= self.g1**2
        elif self.g1 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
            result *= self.g2**2 / constants.g2HZZ**2
        elif self.g1 == self.g2 == self.g1prime2 == self.ghzgs1prime2 == 0:
            result *= self.g4**2 / constants.g4HZZ**2
        elif self.g1 == self.g2 == self.g4 == self.ghzgs1prime2 == 0:
            result *= self.g1prime2**2 / constants.g1prime2HZZ**2
        elif self.g1 == self.g2 == self.g4 == self.g1prime2 == 0:
            result *= self.ghzgs1prime2**2 / constants.ghzgs1prime2HZZ**2
        else:
            raise ValueError("Nominal xsec is only defined for samples with no interference")

        if self.productionmode == "ggH":
            if self.ghg4 == 0:
                result *= self.ghg2**2
            elif self.ghg2 == 0:
                result *= self.ghg4**2 / constants.ghg4HJJ**2
            else:
                raise ValueError("Nominal xsec is only defined for samples with no interference")
        elif self.productionmode in ("bbH", "tqH"):
            pass
        elif self.productionmode == "VBF":
            result *= constants.NNPDF30_lo_as_0130.JHUXSVBFa1.nominal_value
            if self.g2 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g1**2
            elif self.g1 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g2**2 / constants.g2VBF**2
            elif self.g1 == self.g2 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g4**2 / constants.g4VBF**2
            elif self.g1 == self.g2 == self.g4 == self.ghzgs1prime2 == 0:
                result *= self.g1prime2**2 / constants.g1prime2VBF**2
            elif self.g1 == self.g2 == self.g4 == self.g1prime2 == 0:
                result *= self.ghzgs1prime2**2 / constants.ghzgs1prime2VBF**2
            else:
                raise ValueError("Nominal xsec is only defined for samples with no interference")
        elif self.productionmode == "ZH":
            result *= constants.NNPDF30_lo_as_0130.JHUXSZHa1.nominal_value
            if self.g2 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g1**2
            elif self.g1 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g2**2 / constants.g2ZH**2
            elif self.g1 == self.g2 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g4**2 / constants.g4ZH**2
            elif self.g1 == self.g2 == self.g4 == self.ghzgs1prime2 == 0:
                result *= self.g1prime2**2 / constants.g1prime2ZH**2
            elif self.g1 == self.g2 == self.g4 == self.g1prime2 == 0:
                result *= self.ghzgs1prime2**2 / constants.ghzgs1prime2ZH**2
            else:
                raise ValueError("Nominal xsec is only defined for samples with no interference")
        elif self.productionmode == "WH":
            result *= constants.NNPDF30_lo_as_0130.JHUXSWHa1.nominal_value
            if self.g2 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g1**2
            elif self.g1 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g2**2 / constants.g2WH**2
            elif self.g1 == self.g2 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g4**2 / constants.g4WH**2
            elif self.g1 == self.g2 == self.g4 == self.ghzgs1prime2 == 0:
                result *= self.g1prime2**2 / constants.g1prime2WH**2
            elif self.g1 == self.g2 == self.g4 == self.g1prime2 == 0:
                result *= 0
            else:
                raise ValueError("Nominal xsec is only defined for samples with no interference")
        elif self.productionmode == "ttH":
            result *= constants.NNPDF30_lo_as_0130.JHUXSttHkappa.nominal_value
            if self.kappa_tilde == 0:
                result *= self.kappa**2
            elif self.kappa == 0:
                result *= self.kappa_tilde**2 / constants.kappa_tilde_ttH**2
            else:
                raise ValueError("Nominal xsec is only defined for samples with no interference")
        else:
            raise self.ValueError("nominalJHUxsec")

        if not (self.productionmode == "WH" and self.g1 == self.g2 == self.g4 == self.g1prime2 == 0) and result == 0 != self.JHUxsec:
            raise ValueError("Something is wrong")

        return result

    @property
    def JHUxsec(self):
        try:
            return self.JHUxsec_pdf(self.pdf)
        except ValueError as e:
            if "with no pdf" in str(e):
                raise ValueError("{} has no pdf, so can't get JHUxsec".format(self))
            raise

    def JHUxsec_pdf(self, pdf):
        from constants import JHUXSHZZ2e2mua1, JHUXSHZZ2e2mua2, JHUXSHZZ2e2mua3, JHUXSHZZ2e2muL1, JHUXSHZZ2e2muL1Zg,\
                              JHUXSHZZ2e2mua1a2, JHUXSHZZ2e2mua1a3, JHUXSHZZ2e2mua1L1, JHUXSHZZ2e2mua1L1Zg,         \
                              JHUXSHZZ2e2mua2a3, JHUXSHZZ2e2mua2L1, JHUXSHZZ2e2mua2L1Zg, JHUXSHZZ2e2mua3L1,         \
                              JHUXSHZZ2e2mua3L1Zg, JHUXSHZZ2e2muL1L1Zg

        if pdf is not None:
            exec 'from helperstuff.constants.'+pdf+' import ' + """                                            \
                              JHUXSVBFa1, JHUXSVBFa2, JHUXSVBFa3, JHUXSVBFL1, JHUXSVBFL1Zg,                    \
                              JHUXSVBFa1a2, JHUXSVBFa1a3, JHUXSVBFa1L1, JHUXSVBFa1L1Zg,                        \
                              JHUXSVBFa2a3, JHUXSVBFa2L1, JHUXSVBFa2L1Zg, JHUXSVBFa3L1,                        \
                              JHUXSVBFa3L1Zg, JHUXSVBFL1L1Zg,                                                  \
                              JHUXSZHa1, JHUXSZHa2, JHUXSZHa3, JHUXSZHL1, JHUXSZHL1Zg,                         \
                              JHUXSZHa1a2, JHUXSZHa1a3, JHUXSZHa1L1, JHUXSZHa1L1Zg,                            \
                              JHUXSZHa2a3, JHUXSZHa2L1, JHUXSZHa2L1Zg, JHUXSZHa3L1,                            \
                              JHUXSZHa3L1Zg, JHUXSZHL1L1Zg,                                                    \
                              JHUXSWHa1, JHUXSWHa2, JHUXSWHa3, JHUXSWHL1, JHUXSWHL1Zg,                         \
                              JHUXSWHa1a2, JHUXSWHa1a3, JHUXSWHa1L1, JHUXSWHa1L1Zg,                            \
                              JHUXSWHa2a3, JHUXSWHa2L1, JHUXSWHa2L1Zg, JHUXSWHa3L1,                            \
                              JHUXSWHa3L1Zg, JHUXSWHL1L1Zg,                                                    \
                              JHUXSHJJa2, JHUXSHJJa3, JHUXSHJJa2a3,                                            \
                              JHUXSttHkappa, JHUXSttHkappatilde, JHUXSttHkappakappatilde
            """

        if self.productionmode == "bbH" or self.productionmode == "ggH" and self.ghg2 == 1 and self.ghg4 == 0:
            return (
                      JHUXSHZZ2e2mua1*self.g1**2
                    + JHUXSHZZ2e2mua2*self.g2**2
                    + JHUXSHZZ2e2mua3*self.g4**2
                    + JHUXSHZZ2e2muL1*self.g1prime2**2
                    + JHUXSHZZ2e2muL1Zg*self.ghzgs1prime2**2

                    + JHUXSHZZ2e2mua1a2 * self.g1*self.g2
                    + JHUXSHZZ2e2mua1a3 * self.g1*self.g4
                    + JHUXSHZZ2e2mua1L1 * self.g1*self.g1prime2
                    + JHUXSHZZ2e2mua1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua2a3 * self.g4*self.g2
                    + JHUXSHZZ2e2mua2L1 * self.g2*self.g1prime2
                    + JHUXSHZZ2e2mua2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua3L1 * self.g4*self.g1prime2
                    + JHUXSHZZ2e2mua3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSHZZ2e2muL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   )

        if pdf is None:
            raise ValueError("Can't get xsec for {} with no pdf".format(self.productionmode))


        if self.productionmode == "ggH":
            return (
                      JHUXSHJJa2 * self.ghg2**2
                    + JHUXSHJJa3 * self.ghg4**2
                    + JHUXSHJJa2a3 * self.ghg2*self.ghg4
                   ) / JHUXSHJJa2 * (
                      JHUXSHZZ2e2mua1*self.g1**2
                    + JHUXSHZZ2e2mua2*self.g2**2
                    + JHUXSHZZ2e2mua3*self.g4**2
                    + JHUXSHZZ2e2muL1*self.g1prime2**2
                    + JHUXSHZZ2e2muL1Zg*self.ghzgs1prime2**2

                    + JHUXSHZZ2e2mua1a2 * self.g1*self.g2
                    + JHUXSHZZ2e2mua1a3 * self.g1*self.g4
                    + JHUXSHZZ2e2mua1L1 * self.g1*self.g1prime2
                    + JHUXSHZZ2e2mua1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua2a3 * self.g4*self.g2
                    + JHUXSHZZ2e2mua2L1 * self.g2*self.g1prime2
                    + JHUXSHZZ2e2mua2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua3L1 * self.g4*self.g1prime2
                    + JHUXSHZZ2e2mua3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSHZZ2e2muL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   )

        if self.productionmode == "VBF":
            return (
                      JHUXSVBFa1*self.g1**2
                    + JHUXSVBFa2*self.g2**2
                    + JHUXSVBFa3*self.g4**2
                    + JHUXSVBFL1*self.g1prime2**2
                    + JHUXSVBFL1Zg*self.ghzgs1prime2**2

                    + JHUXSVBFa1a2 * self.g1*self.g2
                    + JHUXSVBFa1a3 * self.g1*self.g4
                    + JHUXSVBFa1L1 * self.g1*self.g1prime2
                    + JHUXSVBFa1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSVBFa2a3 * self.g4*self.g2
                    + JHUXSVBFa2L1 * self.g2*self.g1prime2
                    + JHUXSVBFa2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSVBFa3L1 * self.g4*self.g1prime2
                    + JHUXSVBFa3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSVBFL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   ) * (
                      JHUXSHZZ2e2mua1*self.g1**2
                    + JHUXSHZZ2e2mua2*self.g2**2
                    + JHUXSHZZ2e2mua3*self.g4**2
                    + JHUXSHZZ2e2muL1*self.g1prime2**2
                    + JHUXSHZZ2e2muL1Zg*self.ghzgs1prime2**2

                    + JHUXSHZZ2e2mua1a2 * self.g1*self.g2
                    + JHUXSHZZ2e2mua1a3 * self.g1*self.g4
                    + JHUXSHZZ2e2mua1L1 * self.g1*self.g1prime2
                    + JHUXSHZZ2e2mua1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua2a3 * self.g4*self.g2
                    + JHUXSHZZ2e2mua2L1 * self.g2*self.g1prime2
                    + JHUXSHZZ2e2mua2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua3L1 * self.g4*self.g1prime2
                    + JHUXSHZZ2e2mua3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSHZZ2e2muL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   )

        if self.productionmode == "ZH":
            return (
                      JHUXSZHa1*self.g1**2
                    + JHUXSZHa2*self.g2**2
                    + JHUXSZHa3*self.g4**2
                    + JHUXSZHL1*self.g1prime2**2
                    + JHUXSZHL1Zg*self.ghzgs1prime2**2

                    + JHUXSZHa1a2 * self.g1*self.g2
                    + JHUXSZHa1a3 * self.g1*self.g4
                    + JHUXSZHa1L1 * self.g1*self.g1prime2
                    + JHUXSZHa1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSZHa2a3 * self.g4*self.g2
                    + JHUXSZHa2L1 * self.g2*self.g1prime2
                    + JHUXSZHa2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSZHa3L1 * self.g4*self.g1prime2
                    + JHUXSZHa3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSZHL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   ) * (
                      JHUXSHZZ2e2mua1*self.g1**2
                    + JHUXSHZZ2e2mua2*self.g2**2
                    + JHUXSHZZ2e2mua3*self.g4**2
                    + JHUXSHZZ2e2muL1*self.g1prime2**2
                    + JHUXSHZZ2e2muL1Zg*self.ghzgs1prime2**2

                    + JHUXSHZZ2e2mua1a2 * self.g1*self.g2
                    + JHUXSHZZ2e2mua1a3 * self.g1*self.g4
                    + JHUXSHZZ2e2mua1L1 * self.g1*self.g1prime2
                    + JHUXSHZZ2e2mua1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua2a3 * self.g4*self.g2
                    + JHUXSHZZ2e2mua2L1 * self.g2*self.g1prime2
                    + JHUXSHZZ2e2mua2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua3L1 * self.g4*self.g1prime2
                    + JHUXSHZZ2e2mua3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSHZZ2e2muL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   )

        if self.productionmode == "WH":
            return (
                      JHUXSWHa1*self.g1**2
                    + JHUXSWHa2*self.g2**2
                    + JHUXSWHa3*self.g4**2
                    + JHUXSWHL1*self.g1prime2**2
                    + JHUXSWHL1Zg*self.ghzgs1prime2**2

                    + JHUXSWHa1a2 * self.g1*self.g2
                    + JHUXSWHa1a3 * self.g1*self.g4
                    + JHUXSWHa1L1 * self.g1*self.g1prime2
                    + JHUXSWHa1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSWHa2a3 * self.g4*self.g2
                    + JHUXSWHa2L1 * self.g2*self.g1prime2
                    + JHUXSWHa2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSWHa3L1 * self.g4*self.g1prime2
                    + JHUXSWHa3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSWHL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   ) * (
                      JHUXSHZZ2e2mua1*self.g1**2
                    + JHUXSHZZ2e2mua2*self.g2**2
                    + JHUXSHZZ2e2mua3*self.g4**2
                    + JHUXSHZZ2e2muL1*self.g1prime2**2
                    + JHUXSHZZ2e2muL1Zg*self.ghzgs1prime2**2

                    + JHUXSHZZ2e2mua1a2 * self.g1*self.g2
                    + JHUXSHZZ2e2mua1a3 * self.g1*self.g4
                    + JHUXSHZZ2e2mua1L1 * self.g1*self.g1prime2
                    + JHUXSHZZ2e2mua1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua2a3 * self.g4*self.g2
                    + JHUXSHZZ2e2mua2L1 * self.g2*self.g1prime2
                    + JHUXSHZZ2e2mua2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua3L1 * self.g4*self.g1prime2
                    + JHUXSHZZ2e2mua3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSHZZ2e2muL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   )
        if self.productionmode == "ttH":
            return (
                      JHUXSttHkappa * self.kappa**2
                    + JHUXSttHkappatilde * self.kappa_tilde**2
                    + JHUXSttHkappakappatilde * self.kappa*self.kappa_tilde
                   ) * (
                      JHUXSHZZ2e2mua1*self.g1**2
                    + JHUXSHZZ2e2mua2*self.g2**2
                    + JHUXSHZZ2e2mua3*self.g4**2
                    + JHUXSHZZ2e2muL1*self.g1prime2**2
                    + JHUXSHZZ2e2muL1Zg*self.ghzgs1prime2**2

                    + JHUXSHZZ2e2mua1a2 * self.g1*self.g2
                    + JHUXSHZZ2e2mua1a3 * self.g1*self.g4
                    + JHUXSHZZ2e2mua1L1 * self.g1*self.g1prime2
                    + JHUXSHZZ2e2mua1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua2a3 * self.g4*self.g2
                    + JHUXSHZZ2e2mua2L1 * self.g2*self.g1prime2
                    + JHUXSHZZ2e2mua2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua3L1 * self.g4*self.g1prime2
                    + JHUXSHZZ2e2mua3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSHZZ2e2muL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   )

        if self.productionmode == "tqH":
            assert self.kappa_tilde == self.g4 == self.g2 == self.g1prime2 == self.ghzgs1prime2 == 0, self
            return 1

        assert False

    @property
    def xsec(self):
        if self.productionmode == "ggH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSggH2e2mu
            return constants.SMXSggH2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+", "Hff0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "VBF":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSVBF2e2mu
            return constants.SMXSVBF2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "ZH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSZH2e2mu
            return constants.SMXSZH2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "WH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSWH2e2mu
            return constants.SMXSWH2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "VH":
            if isinstance(self, ReweightingSampleWithPdf):
                return ReweightingSampleWithPdf("ZH", self.hypothesis, self.production).xsec + ReweightingSampleWithPdf("WH", self.hypothesis, self.production).xsec
        if self.productionmode == "WplusH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSWpH2e2mu
        if self.productionmode == "WminusH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSWmH2e2mu
        if self.productionmode == "ttH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM" and self.hffhypothesis == "Hff0+": return constants.SMXSttH2e2mu
            return constants.SMXSttH2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+", "Hff0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "bbH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSbbH2e2mu
            return constants.SMXSbbH2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "tqH":
            if hasattr(self, hypothesis) and self.hypothesis == "SM": return constants.SMXStqH2e2mu
            return constants.SMXStqH2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+", "Hff0+").JHUxsec_pdf(self.pdf)
        raise self.ValueError("xsec")

    @property
    def SMxsec(self):
        hffhypothesis = "Hff0+" if self.hffhypothesis else None
        return ReweightingSample(self.productionmode, "SM", hffhypothesis).xsec

    def __eq__(self, other):
        return all((
                   self.productionmode == other.productionmode,
                   self.g1 == other.g1,
                   self.g2 == other.g2,
                   self.g4 == other.g4,
                   self.g1prime2 == other.g1prime2,
                   self.ghzgs1prime2 == other.ghzgs1prime2,
               )) and (
                   self.productionmode != "ggH" or all((
                       self.ghg2 == other.ghg2,
                       self.ghg4 == other.ghg4,
                   ))
               ) and (
                   self.productionmode != "ttH" or all((
                       self.kappa == other.kappa,
                       self.kappa_tilde == other.kappa_tilde,
                   ))
               ) and (
                   self.pdf is None or other.pdf is None or self.pdf == other.pdf
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
        if self.productionmode.isbkg: return [[("1", 1)]]

        factors = []

        weightshelper = WeightsHelper(self)

        for prodordec in "prod", "dec":
          if weightshelper.useproddec(prodordec):
            counter = Counter()
            for couplingname, couplingvalue, weight in weightshelper.couplingsandweights(prodordec, mix=False):
              coupling = getattr(self, couplingname) / float(couplingvalue)
              counter[weight] += coupling**2
            for ((coupling1name, coupling1value, weight1),
                 (coupling2name, coupling2value, weight2),
                 weightint) in weightshelper.couplingsandweights(prodordec, mix=True):
              coupling1 = getattr(self, coupling1name ) / float(coupling1value)
              coupling2 = getattr(self, coupling2name ) / float(coupling2value)
              counter[weight1]   -= coupling1*coupling2
              counter[weight2]   -= coupling1*coupling2
              counter[weightint] += coupling1*coupling2

            factors.append([
                            (weightname, couplingsq)
                                  for weightname, couplingsq in counter.iteritems()
                                   if couplingsq
                                         and weightname is not None
                           ])

        assert factors
        return factors

    @property
    def MC_weight_terms_expanded(self):
      factors = self.MC_weight_terms
      result = []
      for individualfactors in cartesianproduct(*factors):
        result.append(("*".join(weightname for weightname, multiplier in individualfactors), product(multiplier for weightname, multiplier in individualfactors)))
      return result

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
            if h.isEFT: continue
            if not getattr(self, h.couplingname) and h != hypothesis: continue
            kwargs = {_.couplingname: 0 for _ in purehypotheses if not _.isEFT}
            kwargs["pdf"] = "NNPDF30_lo_as_0130"
            kwargs[h.couplingname] = 1
            term = 0
            for productionmode in productionmodes:
                term += getattr(self, h.couplingname)**power * ArbitraryCouplingsSample(productionmode, **kwargs).nominalJHUxsec \
                         * (constants.nominal_normalize_WH_to_ZH if productionmode == "WH" else 1)
            if not withdecay:
                term /= ArbitraryCouplingsSample("ggH", **kwargs).nominalJHUxsec
            denominator += term
            if h == hypothesis:
                numerator = term
        return copysign(
                        numerator / denominator,
                        #mixturesign: multiply by -1 for fL1
                        mixturesign(analysis)*getattr(self, hypothesis.couplingname)*self.g1
                       )

    @property
    def a1eLeR(self):
        assert self.g2 == self.g4 == 0
        ghz1 = self.g1 + 2 * self.g1prime2 * (constants.M_Z**2 - 1j*constants.M_Z*constants.Ga_Z)/constants.L1**2
        ghz1 = ghz1.real  #<-- !!!
        ghzzp1 = constants.M_Z**2/constants.L1**2
        ezp_L = constants.aL * self.g1prime2 + constants.e * self.ghzgs1prime2
        ezp_R = constants.aR * self.g1prime2 + constants.e * self.ghzgs1prime2
        return ghz1, ghzzp1*ezp_L, ghzzp1*ezp_R

    @property
    def feL(self):
        a1, eL, eR = self.a1eLeR
        return (1 if a1*eL>0 else -1) * eL**2 * constants.JHUXSHZZ2e2mueL / (eL**2 * constants.JHUXSHZZ2e2mueL + eR**2 * constants.JHUXSHZZ2e2mueR + abs(a1)**2 * constants.JHUXSHZZ2e2mua1)

    @property
    def feR(self):
        a1, eL, eR = self.a1eLeR
        return (1 if a1*eR>0 else -1) * eR**2 * constants.JHUXSHZZ2e2mueR / (eL**2 * constants.JHUXSHZZ2e2mueL + eR**2 * constants.JHUXSHZZ2e2mueR + abs(a1)**2 * constants.JHUXSHZZ2e2mua1)

    @property
    def hasZZ(self):
        return not self.ghz1 == self.ghz2 == self.ghz4 == self.ghz1prime2 == 0
    @property
    def hasZg(self):
        return not self.ghzgs1prime2 == 0
    @property
    def hasgg(self):
        return False

class ArbitraryCouplingsSample(SampleBase):
    def __init__(self, productionmode, **kwargs):
        g1 = kwargs.pop("g1", None)
        g2 = kwargs.pop("g2", None)
        g4 = kwargs.pop("g4", None)
        g1prime2 = kwargs.pop("g1prime2", None)

        ghz1 = kwargs.pop("ghz1", None)
        ghz2 = kwargs.pop("ghz2", None)
        ghz4 = kwargs.pop("ghz4", None)
        ghz1prime2 = kwargs.pop("ghz1prime2", None)

        ghw1 = kwargs.pop("ghw1", None)
        ghw2 = kwargs.pop("ghw2", None)
        ghw4 = kwargs.pop("ghw4", None)
        ghw1prime2 = kwargs.pop("ghw1prime2", None)

        ghzgs1prime2 = kwargs.pop("ghzgs1prime2", None)

        ghg2 = kwargs.pop("ghg2", None)
        ghg4 = kwargs.pop("ghg4", None)

        kappa = kwargs.pop("kappa", None)
        kappa_tilde = kwargs.pop("kappa_tilde", None)

        pdf = kwargs.pop("pdf", None)

        if kwargs: raise TypeError("Extra kwargs: "+str(kwargs))

        self.productionmode = ProductionMode(productionmode)
        if self.productionmode not in ("ggH", "VBF", "ZH", "WH", "ttH", "bbH"):
            raise ValueError("Bad productionmode {}".format(self.productionmode))

        if self.productionmode in ("ggH", "ZH", "ttH", "bbH"):
            if not (ghw1 is ghw2 is ghw4 is ghw1prime2 is None):
                raise ValueError("Don't provide WW couplings for {}".format(self.productionmode))
            if g1 is None: g1, ghz1 = ghz1, None
            if g2 is None: g2, ghz2 = ghz2, None
            if g4 is None: g4, ghz4 = ghz4, None
            if g1prime2 is None: g1prime2, ghz1prime2 = ghz1prime2, None

        if ghz1 is ghw1 is None: ghz1 = ghw1 = g1; g1 = None
        if ghz2 is ghw2 is None: ghz2 = ghw2 = g2; g2 = None
        if ghz4 is ghw4 is None: ghz4 = ghw4 = g4; g4 = None
        if ghz1prime2 is ghw1prime2 is None: ghz1prime2 = ghw1prime2 = g1prime2; g1prime2 = None

        if g1 is not None: raise TypeError("Provided both g1 and ghz1/ghw1")
        if g2 is not None: raise TypeError("Provided both g2 and ghz2/ghw2")
        if g4 is not None: raise TypeError("Provided both g4 and ghz4/ghw4")
        if g1prime2 is not None: raise TypeError("Provided both g1prime2 and ghz1prime2/ghw1prime2")

        if ghz1 is None: raise TypeError("Have to provide ghz1 or g1")
        if ghz2 is None: raise TypeError("Have to provide ghz2 or g2")
        if ghz4 is None: raise TypeError("Have to provide ghz4 or g4")
        if ghz1prime2 is None: raise TypeError("Have to provide ghz1prime2 or g1prime2")

        if ghw1 is None: raise TypeError("Have to provide ghw1 or g1")
        if ghw2 is None: raise TypeError("Have to provide ghw2 or g2")
        if ghw4 is None: raise TypeError("Have to provide ghw4 or g4")
        if ghw1prime2 is None: raise TypeError("Have to provide ghw1prime2 or g1prime2")

        if ghzgs1prime2 is None: raise TypeError("Have to provide ghzgs1prime2")

        self.__ghz1 = ghz1
        self.__ghz2 = ghz2
        self.__ghz4 = ghz4
        self.__ghz1prime2 = ghz1prime2

        self.__ghw1 = ghw1
        self.__ghw2 = ghw2
        self.__ghw4 = ghw4
        self.__ghw1prime2 = ghw1prime2

        self.__ghzgs1prime2 = ghzgs1prime2

        if self.productionmode == "ggH":
            if ghg2 is ghg4 is None:
                ghg2 = 1
                ghg4 = 0
            if ghg2 is None or ghg4 is None:
                raise ValueError("Have to set ghg2 and ghg4 for ggH, or leave both to the defaults of 1 and 0")
        else:
            if ghg2 is not None or ghg4 is not None:
                raise ValueError("Can't set ghg2 or ghg4 for {}".format(self.productionmode))

        self.__ghg2 = ghg2
        self.__ghg4 = ghg4

        if self.productionmode == "ttH":
            self.__kappa, self.__kappa_tilde = kappa, kappa_tilde
            if kappa is None or kappa_tilde is None:
                raise ValueError("Have to set kappa and kappa_tilde for ttH")
        else:
            if kappa is not None or kappa_tilde is not None:
                raise ValueError("Can't set kappa or kappa_tilde for {}".format(self.productionmode))

        self.__kappa = kappa
        self.__kappa_tilde = kappa_tilde

        self.__pdf = pdf
        if pdf is not None:
            try:
                exec "from constants import "+pdf in globals(), locals()
            except ImportError:
                raise ValueError("Unknown pdf "+pdf)

    @property
    def ghz1(self):
        return self.__ghz1
    @property
    def ghz2(self):
        return self.__ghz2
    @property
    def ghz4(self):
        return self.__ghz4
    @property
    def ghz1prime2(self):
        return self.__ghz1prime2
    @property
    def ghw1(self):
        return self.__ghw1
    @property
    def ghw2(self):
        return self.__ghw2
    @property
    def ghw4(self):
        return self.__ghw4
    @property
    def ghw1prime2(self):
        return self.__ghw1prime2
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
    @property
    def pdf(self):
        return self.__pdf

    def __repr__(self):
        kwargs = ()
        if self.productionmode == "ggH": kwargs += ("ghg2", "ghg4")
        if self.productionmode == "ttH": kwargs += ("kappa", "kappa_tilde")
        if self.pdf is not None: kwargs += ("pdf",)
        kwargs += ("ghz1", "ghz2", "ghz4", "ghz1prime2", "ghzgs1prime2")
        if self.productionmode in ("VBF", "WH"):
            kwargs += ("ghw1", "ghw2", "ghw4", "ghw1prime2")
        return "{}({}, {})".format(
                                   type(self).__name__,
                                   repr(self.productionmode.item.name),
                                   ", ".join(
                                             "{}={}".format(kwarg, getattr(self, kwarg))
                                             for kwarg in kwargs
                                            ),
                                  )

def samplewithfai(productionmode, analysis, fai, withdecay=False, productionmodeforfai=None, pdf=None):
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
    kwargs["pdf"] = pdf
    return ArbitraryCouplingsSample(productionmode, **kwargs)

def samplewitha1eLeR(a1, eL, eR):
    kwargs = {
      "g1": a1 - 2/(constants.aL - constants.aR) * (eL - eR),
      "g1prime2": constants.L1**2 / (constants.M_Z**2 * (constants.aL - constants.aR)) * (eL-eR),
      "ghzgs1prime2": -constants.L1**2 / (constants.M_Z**2 * constants.e * (constants.aL - constants.aR)) * (constants.aR*eL-constants.aL*eR),
      "g2": 0,
      "g4": 0,
    }
    result = ArbitraryCouplingsSample("ggH", **kwargs)
    return result

def samplewithfeLfeR(feL, feR):
    fa1 = 1-abs(feL)-abs(feR)
    a1 = sqrt(fa1)
    eL = (1 if feL>0 else -1) * sqrt(abs(feL) * constants.JHUXSHZZ2e2mua1 / constants.JHUXSHZZ2e2mueL)
    eR = (1 if feR>0 else -1) * sqrt(abs(feR) * constants.JHUXSHZZ2e2mua1 / constants.JHUXSHZZ2e2mueR)
    return samplewitha1eLeR(a1, eL, eR)

def samplewithfL1fL1Zg(fL1, fL1Zg):
    fa1 = 1-abs(fL1)-abs(fL1Zg)
    a1 = sqrt(fa1)
    L1 = -(1 if fL1>0 else -1) * sqrt(abs(fL1) * constants.JHUXSHZZ2e2mua1 / constants.JHUXSHZZ2e2muL1)
    L1Zg = -(1 if fL1Zg>0 else -1) * sqrt(abs(fL1Zg) * constants.JHUXSHZZ2e2mua1 / constants.JHUXSHZZ2e2muL1Zg)
    return ArbitraryCouplingsSample("ggH", g1=a1, g1prime2=L1, ghzgs1prime2=L1Zg, g2=0, g4=0)

class ReweightingSample(MultiEnum, SampleBase):
    __metaclass__ = MultiEnumABCMeta
    enumname = "reweightingsample"
    enums = [ProductionMode, Hypothesis, HffHypothesis]

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
        elif self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "ttH", "bbH", "tqH", "VH"):
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hypothesis not in self.productionmode.validhypotheses:
                raise ValueError("{} hypothesis can't be {}\n{}".format(self.productionmode, self.hypothesis, args))
            if self.productionmode in ("ggH", "ttH", "tqH"):
                if self.hffhypothesis is None:
                    if self.productionmode == "ggH":
                        self.hffhypothesis = HffHypothesis("Hff0+")
                    else:
                        raise ValueError("Hff hypothesis not provided for {} productionmode\n{}".format(self.productionmode, args))
                if self.productionmode == "tqH" and self.hffhypothesis != "Hff0+":
                    raise ValueError("Bad hffhypothesis {} for {} productionmode\n{}".format(self.hffhypothesis, self.productionmode, args))
            else:
                if self.hffhypothesis is not None:
                    raise ValueError("Hff hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
        elif self.productionmode in ("ggZZ", "VBF bkg"):
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hffhypothesis is not None:
                raise ValueError("Hff hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
        elif self.productionmode in ("qqZZ", "TTZZ", "ZZZ", "WZZ", "WWZ", "TTWW", "TTZJets_M10_MLM", "TTZToLLNuNu_M10", "TTZToLL_M1to10_MLM", "EW") + tuple(ggZZoffshellproductionmodes):
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hffhypothesis is not None:
                raise ValueError("Hff hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
        elif self.productionmode in ("ZX", "data"):
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hffhypothesis is not None:
                raise ValueError("Hff hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
        else:
            raise ValueError("Bad productionmode {}\n{}".format(self.productionmode, args))

    @property
    def isbkg(self):
        return self.productionmode.isbkg
    @property
    def issignal(self):
        return not self.isbkg

    def isZX(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBF bkg", "data", "VBF", "ZH", "WH", "ttH", "WplusH", "WminusH", "bbH", "tqH", "TTZZ", "ZZZ", "WZZ", "WWZ", "TTWW", "TTZJets_M10_MLM", "TTZToLLNuNu_M10", "TTZToLL_M1to10_MLM") + tuple(ggZZoffshellproductionmodes):
            return False
        elif self.productionmode == "ZX":
            return True
        raise self.ValueError("isZX")

    def isdata(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBF bkg", "ZX", "VBF", "ZH", "WH", "ttH", "WplusH", "WminusH", "bbH", "tqH", "TTZZ", "ZZZ", "WZZ", "WWZ", "TTWW", "TTZJets_M10_MLM", "TTZToLLNuNu_M10", "TTZToLL_M1to10_MLM") + tuple(ggZZoffshellproductionmodes):
            return False
        elif self.productionmode == "data":
            return True
        raise self.ValueError("isdata")

    def TDirectoryname(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBFbkg", "data", "VBF", "ZH", "WH", "ttH", "WplusH", "WminusH", "bbH", "tqH", "TTZZ", "ZZZ", "WZZ", "WWZ", "TTWW", "TTZJets_M10_MLM", "TTZToLLNuNu_M10", "TTZToLL_M1to10_MLM") + tuple(ggZZoffshellproductionmodes) or self.productionmode == "ZX" and not config.usedata:
            return "ZZTree"
        if self.productionmode == "ZX":
            return "CRZLLTree"
        raise self.ValueError("TDirectoryname")

    @property
    @cache
    def ghz1(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "ttH", "bbH", "tqH"):
            if self.hypothesis == "0+": return 1
            if self.hypothesis in ("0-", "a2", "L1", "L1Zg", "a3EFT", "a2EFT", "L1EFT"): return 0
            for EFT, couplings in ("", ("fa3", "fa2", "fL1", "fL1Zg")), ("EFT", ("fa3", "fa2", "fL1")):
                for proddec in "prod", "dec", "proddec":
                    minus = "-" if proddec == "proddec" else ""
                    doplushalf = (not EFT) or proddec != "proddec"
                    dominushalf = (not EFT) or proddec == "proddec"
                    for a in couplings:
                        if doplushalf and self.hypothesis == a+EFT+proddec+"0.5":
                            return 1
                        if dominushalf and self.hypothesis == a+EFT+proddec+"-0.5":
                            return 1
                    for a, b in combinations(couplings, 2):
                        if self.hypothesis == a+EFT+proddec+"0.5"+b+EFT+proddec+minus+"0.5":
                            return 0
                        if self.hypothesis == a+EFT+proddec+"0.33"+b+EFT+proddec+minus+"0.33":
                            return 1
                    for a, b, c in combinations(couplings, 3):
                        if self.hypothesis == a+EFT+proddec+"0.33"+b+EFT+proddec+"0.33"+c+EFT+proddec+minus+"0.33":
                            return 0
                        if proddec != "proddec": continue
                        if self.hypothesis == a+EFT+proddec+"0.25"+b+EFT+proddec+"0.25"+c+EFT+proddec+"0.25":
                            return 1
                    if proddec != "proddec": continue
                    for a, b, c, d in combinations(couplings, 4):
                        if self.hypothesis == a+EFT+proddec+"0.25"+b+EFT+proddec+"0.25"+c+EFT+proddec+"0.25"+d+EFT+proddec+"0.25":
                            return 0

            if self.hypothesis == "fa2dec-0.9": return 1
            if self.hypothesis == "fa3VBF0.5": return 1
            if self.hypothesis == "fa3VH0.5": return 1
            if self.hypothesis == "fa2VBF0.5": return 1
            if self.hypothesis == "fa2VH0.5": return 1
            if self.hypothesis == "fL1VBF0.5": return 1
            if self.hypothesis == "fL1VH0.5": return 1
            if self.hypothesis == "BestFit19009": return 0.504778 ** 0.5

        raise self.ValueError("ghz1")

    @property
    @cache
    def ghz2(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "ttH", "bbH"):
            if self.hypothesis in ("a2", "a2EFT"): return 1
            if self.hypothesis in ("0+", "0-", "L1", "L1Zg", "a3EFT", "L1EFT"): return 0

            ghz2 = {"dec": constants.g2HZZ}
            if self.productionmode == "VBF": ghz2["prod"] = constants.g2VBF
            elif self.productionmode == "ZH":  ghz2["prod"] = constants.g2ZH
            elif self.productionmode == "WH":  ghz2["prod"] = constants.g2WH
            else: ghz2["prod"] = None
            ghz2["proddec"] = None if ghz2["prod"] is None else sqrt(ghz2["prod"]*ghz2["dec"])

            for EFT, couplings in ("", ("fa3", "fa2", "fL1", "fL1Zg")), ("EFT", ("fa3", "fa2", "fL1")):
                for proddec in "prod", "dec", "proddec":
                    doplushalf = (not EFT) or proddec != "proddec"
                    dominushalf = (not EFT) or proddec == "proddec"
                    for a in couplings:
                        if doplushalf and self.hypothesis == a+EFT+proddec+"0.5":
                            if a == "fa2":
                                return ghz2[proddec]
                            return 0
                        if dominushalf and self.hypothesis == a+EFT+proddec+"-0.5":
                            if a == "fa2":
                                return -ghz2[proddec]
                            return 0

                    minus = "-" if proddec == "proddec" else ""
                    for a, b in combinations(couplings, 2):
                        if self.hypothesis in (a+EFT+proddec+"0.5"+b+EFT+proddec+minus+"0.5",
                                               a+EFT+proddec+"0.33"+b+EFT+proddec+minus+"0.33"):
                            if a == "fa2":
                                return ghz2[proddec]
                            if b == "fa2":
                                return float(minus+"1") * ghz2[proddec]
                            return 0
                    for a, b, c in combinations(couplings, 3):
                        if (self.hypothesis == a+EFT+proddec+"0.33"+b+EFT+proddec+"0.33"+c+EFT+proddec+minus+"0.33"
                              or proddec=="proddec"
                              and self.hypothesis == a+EFT+proddec+"0.25"+b+EFT+proddec+"0.25"+c+EFT+proddec+"0.25"
                        ):
                            if a == "fa2" or b == "fa2" or "0.25" in str(self.hypothesis) and c == "fa2":
                                return ghz2[proddec]
                            if c == "fa2":
                                return float(minus+"1") * ghz2[proddec]
                            return 0

                    if proddec != "proddec": continue
                    for a, b, c, d in combinations(couplings, 4):
                        if self.hypothesis == a+EFT+proddec+"0.25"+b+EFT+proddec+"0.25"+c+EFT+proddec+"0.25"+d+EFT+proddec+"0.25":
                            if a == "fa2" or b == "fa2" or c == "fa2" or d == "fa2":
                                return ghz2[proddec]
                            return 0

            if self.hypothesis == "fa2dec-0.9": return -constants.g2HZZ * 3
            if self.hypothesis == "fa3VBF0.5": return 0
            if self.hypothesis == "fa3VH0.5": return 0
            if self.hypothesis == "fa2VBF0.5": return constants.g2VBF
            if self.hypothesis == "fa2VH0.5": return constants.g2VH
            if self.hypothesis == "fL1VBF0.5": return 0
            if self.hypothesis == "fL1VH0.5": return 0
            if self.hypothesis == "BestFit19009": return constants.g2HZZ * -0.289103**0.5

        raise self.ValueError("ghz2")

    @property
    @cache
    def ghz4(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "ttH", "bbH"):
            if self.hypothesis in ("0-", "a3EFT"): return 1
            if self.hypothesis in ("0+", "a2", "L1", "L1Zg", "a2EFT", "L1EFT"): return 0

            ghz4 = {"dec": constants.g4HZZ}
            if self.productionmode == "VBF": ghz4["prod"] = constants.g4VBF
            elif self.productionmode == "ZH":  ghz4["prod"] = constants.g4ZH
            elif self.productionmode == "WH":  ghz4["prod"] = constants.g4WH
            else: ghz4["prod"] = None
            ghz4["proddec"] = None if ghz4["prod"] is None else sqrt(ghz4["prod"]*ghz4["dec"])

            for EFT, couplings in ("", ("fa3", "fa2", "fL1", "fL1Zg")), ("EFT", ("fa3", "fa2", "fL1")):
                for proddec in "prod", "dec", "proddec":
                    doplushalf = (not EFT) or proddec != "proddec"
                    dominushalf = (not EFT) or proddec == "proddec"
                    for a in couplings:
                        if doplushalf and self.hypothesis == a+EFT+proddec+"0.5":
                            if a == "fa3":
                                return ghz4[proddec]
                            return 0
                        if dominushalf and self.hypothesis == a+EFT+proddec+"-0.5":
                            if a == "fa3":
                                return -ghz4[proddec]
                            return 0

                    minus = "-" if proddec == "proddec" else ""
                    for a, b in combinations(couplings, 2):
                        if self.hypothesis in (a+EFT+proddec+"0.5"+b+EFT+proddec+minus+"0.5",
                                               a+EFT+proddec+"0.33"+b+EFT+proddec+minus+"0.33"):
                            if a == "fa3":
                                return ghz4[proddec]
                            if b == "fa3":
                                return float(minus+"1") * ghz4[proddec]
                            return 0
                    for a, b, c in combinations(couplings, 3):
                        if (self.hypothesis == a+EFT+proddec+"0.33"+b+EFT+proddec+"0.33"+c+EFT+proddec+minus+"0.33"
                              or proddec=="proddec"
                              and self.hypothesis == a+EFT+proddec+"0.25"+b+EFT+proddec+"0.25"+c+EFT+proddec+"0.25"
                        ):
                            if a == "fa3" or b == "fa3" or "0.25" in str(self.hypothesis) and c == "fa3":
                                return ghz4[proddec]
                            if c == "fa3":
                                return float(minus+"1") * ghz4[proddec]
                            return 0

                    if proddec != "proddec": continue
                    for a, b, c, d in combinations(couplings, 4):
                        if self.hypothesis == a+EFT+proddec+"0.25"+b+EFT+proddec+"0.25"+c+EFT+proddec+"0.25"+d+EFT+proddec+"0.25":
                            if a == "fa3" or b == "fa3" or c == "fa3" or d == "fa3":
                                return ghz4[proddec]
                            return 0

            if self.hypothesis == "fa2dec-0.9": return 0

            if self.hypothesis == "fa3VBF0.5": return constants.g4VBF
            if self.hypothesis == "fa3VH0.5": return constants.g4VH
            if self.hypothesis == "fa2VBF0.5": return 0
            if self.hypothesis == "fa2VH0.5": return 0
            if self.hypothesis == "fL1VBF0.5": return 0
            if self.hypothesis == "fL1VH0.5": return 0
            if self.hypothesis == "BestFit19009": return constants.g4HZZ * -0.01**0.5

        raise self.ValueError("ghz4")

    @property
    @cache
    def ghz1prime2(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "ttH", "bbH"):
            if self.hypothesis in ("L1", "L1EFT"): return 1e4
            if self.hypothesis in ("0+", "a2", "0-", "L1Zg", "a2EFT", "a3EFT"): return 0

            ghz1prime2 = {"dec": constants.g1prime2HZZ}
            if self.productionmode == "VBF": ghz1prime2["prod"] = constants.g1prime2VBF
            elif self.productionmode == "ZH":  ghz1prime2["prod"] = constants.g1prime2ZH
            elif self.productionmode == "WH":  ghz1prime2["prod"] = constants.g1prime2WH
            else: ghz1prime2["prod"] = None
            ghz1prime2["proddec"] = None if ghz1prime2["prod"] is None else -sqrt(ghz1prime2["prod"]*ghz1prime2["dec"])

            for EFT, couplings in ("", ("fa3", "fa2", "fL1", "fL1Zg")), ("EFT", ("fa3", "fa2", "fL1")):
                for proddec in "prod", "dec", "proddec":
                    doplushalf = (not EFT) or proddec != "proddec"
                    dominushalf = (not EFT) or proddec == "proddec"
                    for a in couplings:
                        if doplushalf and self.hypothesis == a+EFT+proddec+"0.5":
                            if a == "fL1":
                                return ghz1prime2[proddec]
                            return 0
                        if dominushalf and self.hypothesis == a+EFT+proddec+"-0.5":
                            if a == "fL1":
                                return -ghz1prime2[proddec]
                            return 0

                    minus = "-" if proddec == "proddec" else ""
                    for a, b in combinations(couplings, 2):
                        if self.hypothesis in (a+EFT+proddec+"0.5"+b+EFT+proddec+minus+"0.5",
                                               a+EFT+proddec+"0.33"+b+EFT+proddec+minus+"0.33"):
                            if a == "fL1":
                                return ghz1prime2[proddec]
                            if b == "fL1":
                                return float(minus+"1") * ghz1prime2[proddec]
                            return 0
                    for a, b, c in combinations(couplings, 3):
                        if (self.hypothesis == a+EFT+proddec+"0.33"+b+EFT+proddec+"0.33"+c+EFT+proddec+minus+"0.33"
                              or proddec=="proddec"
                              and self.hypothesis == a+EFT+proddec+"0.25"+b+EFT+proddec+"0.25"+c+EFT+proddec+"0.25"
                        ):
                            if a == "fL1" or b == "fL1" or "0.25" in str(self.hypothesis) and c == "fL1":
                                return ghz1prime2[proddec]
                            if c == "fL1":
                                return float(minus+"1") * ghz1prime2[proddec]
                            return 0

                    if proddec != "proddec": continue
                    for a, b, c, d in combinations(couplings, 4):
                        if self.hypothesis == a+EFT+proddec+"0.25"+b+EFT+proddec+"0.25"+c+EFT+proddec+"0.25"+d+EFT+proddec+"0.25":
                            if a == "fL1" or b == "fL1" or c == "fL1" or d == "fL1":
                                return ghz1prime2[proddec]
                            return 0

            if self.hypothesis == "fa2dec-0.9": return 0

            if self.hypothesis == "fa3VBF0.5": return 0
            if self.hypothesis == "fa3VH0.5": return 0
            if self.hypothesis == "fa2VBF0.5": return 0
            if self.hypothesis == "fa2VH0.5": return 0
            if self.hypothesis == "fL1VBF0.5": return constants.g1prime2VBF
            if self.hypothesis == "fL1VH0.5": return constants.g1prime2VH
            if self.hypothesis == "BestFit19009": return constants.g1prime2HZZ * 0.133879**0.5

        raise self.ValueError("ghz1prime2")

    @property
    def ghw1(self): return self.ghz1
    @property
    def ghw2(self):
        if self.hypothesis.isEFT: return eft.ghw2(ghz1=self.ghz1, ghz2=self.ghz2, ghz4=self.ghz4, ghz1prime2=self.ghz1prime2)
        return self.ghz2
    @property
    def ghw4(self):
        if self.hypothesis.isEFT: return eft.ghw4(ghz1=self.ghz1, ghz2=self.ghz2, ghz4=self.ghz4, ghz1prime2=self.ghz1prime2)
        return self.ghz4
    @property
    def ghw1prime2(self):
        if self.hypothesis.isEFT: return eft.ghw1prime2(ghz1=self.ghz1, ghz2=self.ghz2, ghz4=self.ghz4, ghz1prime2=self.ghz1prime2)
        return self.ghz1prime2

    @property
    @cache
    def ghzgs1prime2(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "ttH", "bbH"):
            if self.hypothesis.isEFT: return eft.ghzgs1prime2(ghz1=self.ghz1, ghz2=self.ghz2, ghz4=self.ghz4, ghz1prime2=self.ghz1prime2)
            if self.hypothesis == "L1Zg": return 1e4
            if self.hypothesis in ("0+", "a2", "0-", "L1"): return 0

            ghzgs1prime2 = {"dec": constants.ghzgs1prime2HZZ}
            if self.productionmode == "VBF": ghzgs1prime2["prod"] = constants.ghzgs1prime2VBF
            elif self.productionmode == "ZH":  ghzgs1prime2["prod"] = constants.ghzgs1prime2ZH
            elif self.productionmode == "WH":  ghzgs1prime2["prod"] = -1000  #dummy value, should remove this entirely
            else: ghzgs1prime2["prod"] = None
            ghzgs1prime2["proddec"] = None if ghzgs1prime2["prod"] is None else -sqrt(ghzgs1prime2["prod"]*ghzgs1prime2["dec"])

            couplings = "fa3", "fa2", "fL1", "fL1Zg"
            for proddec in "prod", "dec", "proddec":
                for a in couplings:
                    if self.hypothesis == a+proddec+"0.5":
                        if a == "fL1Zg":
                            return ghzgs1prime2[proddec]
                        return 0
                    if self.hypothesis == a+proddec+"-0.5":
                        if a == "fL1Zg":
                            return -ghzgs1prime2[proddec]
                        return 0

                minus = "-" if proddec == "proddec" else ""
                for a, b in combinations(couplings, 2):
                    if self.hypothesis in (a+proddec+"0.5"+b+proddec+minus+"0.5",
                                           a+proddec+"0.33"+b+proddec+minus+"0.33"):
                        if a == "fL1Zg":
                            return ghzgs1prime2[proddec]
                        if b == "fL1Zg":
                            return float(minus+"1") * ghzgs1prime2[proddec]
                        return 0
                for a, b, c in combinations(couplings, 3):
                    if (self.hypothesis == a+proddec+"0.33"+b+proddec+"0.33"+c+proddec+minus+"0.33"
                          or proddec=="proddec"
                          and self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"
                    ):
                        if a == "fL1Zg" or b == "fL1Zg" or "0.25" in str(self.hypothesis) and c == "fL1Zg":
                            return ghzgs1prime2[proddec]
                        if c == "fL1Zg":
                            return float(minus+"1") * ghzgs1prime2[proddec]
                        return 0

                if proddec != "proddec": continue
                for a, b, c, d in combinations(couplings, 4):
                    if self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"+d+proddec+"0.25":
                        if a == "fL1Zg" or b == "fL1Zg" or c == "fL1Zg" or d == "fL1Zg":
                            return ghzgs1prime2[proddec]
                        return 0

            if self.hypothesis == "fa2dec-0.9": return 0

            if self.hypothesis == "fa3VBF0.5": return 0
            if self.hypothesis == "fa3VH0.5": return 0
            if self.hypothesis == "fa2VBF0.5": return 0
            if self.hypothesis == "fa2VH0.5": return 0
            if self.hypothesis == "fL1VBF0.5": return 0
            if self.hypothesis == "fL1VH0.5": return 0
            if self.hypothesis == "BestFit19009": return constants.ghzgs1prime2HZZ * -0.0622403**0.5

        raise self.ValueError("ghzgs1prime2")

    @property
    def ghg2(self):
        if self.productionmode == "ggH":
            if self.hffhypothesis in ("Hff0+", "fCP0.5"):
                return 1
            if self.hffhypothesis == "Hff0-":
                return 0
        raise self.ValueError("ghg2")

    @property
    def ghg4(self):
        if self.productionmode == "ggH":
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
        if self.productionmode == "tqH":
            assert self.hffhypothesis == "Hff0+"
            return 1
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
        if self.productionmode == "tqH":
            assert self.hffhypothesis == "Hff0+"
            return 0
    @property
    def pdf(self): return None

class ReweightingSampleWithPdf(ReweightingSample):
    enums = [ReweightingSample, Production]
    @property
    def pdf(self): return self.production.pdf

class ReweightingSampleWithFlavor(ReweightingSample):
    enums = [ReweightingSample, Flavor]
    def check(self, *args):
        if self.productionmode in ("ggZZ",):
            if self.flavor is None:
                raise ValueError("No flavor provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.productionmode == "VBF bkg" and self.flavor.hastaus:
                raise ValueError("No {} samples with taus\n{}".format(self.productionmode, args))
        else:
            if self.flavor is not None:
                raise ValueError("Flavor provided for {} productionmode\n{}".format(self.productionmode, args))

class ReweightingSamplePlus(ReweightingSample):
    enums = [ReweightingSample, AlternateGenerator, PythiaSystematic, Extension]
    enumname = "reweightingsampleplus"

    def check(self, *args):
        super(ReweightingSamplePlus, self).check(*args)

        if self.hffhypothesis is None and self.productionmode == "ggH": self.hffhypothesis = HffHypothesis("Hff0+")

        if self.pythiasystematic is not None:
            if (
                self.hypothesis != "0+"
                or self.productionmode not in ("ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH")
               ):
                raise ValueError("No {} {} sample produced with pythia {}\n{}".format(self.productionmode, self.hypothesis, self.pythiasystematic, args))
            if self.alternategenerator not in ("POWHEG", "MINLO"):
                raise ValueError("{} sample with pythia {} is produced with POWHEG{}\n{}".format(self.productionmode, self.pythiasystematic, " or MINLO" if self.productionmode == "ggH" else "", args))

        if self.pythiasystematic is not None:
            if (
                self.hffhypothesis != "Hff0+"
                 and self.hffhypothesis is not None
               ):
                raise ValueError("No {} {} sample produced with pythia {}\n{}".format(self.productionmode, self.hffhypothesis, self.pythiasystematic, args))

class ReweightingSamplePlusWithFlavor(ReweightingSamplePlus):
    enumname = "reweightingsamplepluswithflavor"
    enums = [ReweightingSamplePlus, Flavor]

    def check(self, *args):
        self.reweightingsamplewithflavor = ReweightingSampleWithFlavor(self.reweightingsample, self.flavor)
        super(ReweightingSamplePlusWithFlavor, self).check(*args)

class Sample(ReweightingSamplePlusWithFlavor):
    enums = [ReweightingSamplePlusWithFlavor, Production]

    def check(self, *args):
        if self.production is None:
            raise ValueError("No option provided for production\n{}".format(args))

        if self.hypothesis is not None and self.hypothesis not in self.productionmode.generatedhypotheses(self.production):
            raise ValueError("No {} sample produced with hypothesis {}!\n{}".format(self.productionmode, self.hypothesis, args))

        if self.productionmode in ("WplusH", "WminusH") and self.alternategenerator != "POWHEG":
            raise ValueError("Separate {} sample is produced with POWHEG.  Maybe you meant to specify POWHEG, or WH?\n{}".format(self.productionmode, args))

        self.reweightingsamplewithpdf = ReweightingSampleWithPdf(self.reweightingsample, self.production)

        if self.pythiasystematic is not None and self.alternategenerator == "NNLOPS":
            raise ValueError("No NNLOPS samples with systematics!\n{}".format(args))

        super(Sample, self).check(*args)

        if {
          "POWHEG": (
             self.hypothesis != "0+"
             or self.productionmode not in ("ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH")
             or None is not self.hffhypothesis != "Hff0+"
           ),
           "MINLO": (
             self.hypothesis != "0+"
             or self.productionmode != "ggH"
             or self.hffhypothesis != "Hff0+"
           ),
           "NNLOPS": (
             self.hypothesis != "0+"
             or self.productionmode != "ggH"
             or self.hffhypothesis != "Hff0+"
           ),
           "MCatNLO": (
             self.hypothesis != "0+"
             or self.productionmode != "ggH"
           ),
           "JHUGen": (
             self.hypothesis != "0+"
             or self.productionmode != "ggH"
           ),
           "None": (
             self.productionmode == "ggH" and self.hffhypothesis != "Hff0+"
           )
        }[str(self.alternategenerator)]:
            raise ValueError("No {} sample produced with {}\n{}".format(self.reweightingsample, self.alternategenerator, args))

    def CJLSTmaindir(self):
      if self.alternategenerator is None or self.alternategenerator in ("MCatNLO", "JHUGen"):
        if self.productionmode in ("ggH", "ttH"):
          return self.production.CJLSTdir_anomalous()
        if self.productionmode == "VBF":
          return self.production.CJLSTdir_anomalous_VBF()
        if self.productionmode in ("ZH", "WH"):
          return self.production.CJLSTdir_anomalous_VH()
        if self.productionmode in ("data", "ZX"):
          return self.production.CJLSTdir_data()
      if self.productionmode == "ggH" and self.alternategenerator == "MINLO" and self.pythiasystematic is None:
        return self.production.CJLSTdir_MINLO()
      if self.productionmode == "VBF bkg":
        return self.production.CJLSTdir_VBFbkg()
      if self.productionmode in ggZZoffshellproductionmodes:
        return self.production.CJLSTdir_ggZZoffshell()
      return self.production.CJLSTdir()

    def CJLSTdirname(self):
        if self.alternategenerator == "POWHEG":
            if self.productionmode in ("ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH"):
                s = str(self.productionmode)
                if self.productionmode == "VBF": s = "VBFH"
                if self.hypothesis == "0+": result = "{}125".format(s)
            if self.pythiasystematic is not None:
                result += self.pythiasystematic.appendname
            if self.extension is not None: result += str(self.extension)
            return result
        if self.alternategenerator in ("MINLO", "NNLOPS"):
            if self.productionmode == "ggH":
                result = "ggH125"
                if self.pythiasystematic is not None:
                    result += self.pythiasystematic.appendname
                if self.alternategenerator == "MINLO":
                    result += "_minloHJJ"
                if self.alternategenerator == "NNLOPS":
                    result += "_NNLOPS"
                return result
            raise self.ValueError("CJLSTdirname")
        if self.alternategenerator in ("JHUGen", "MCatNLO") and self.productionmode == "ggH":
            if self.hffhypothesis == "Hff0+": result = "HJJ0PM_M125"
            if self.hffhypothesis == "Hff0-": result = "HJJ0M_M125"
            if self.hffhypothesis == "fCP0.5": result = "HJJ0Mf05ph0_M125"
            if self.alternategenerator == "MCatNLO":
                result += "_mcanlo"
                result = result.replace("ph0", "")
            if self.extension is not None: result += "_" + str(self.extension)
            return result
        if self.productionmode in ("ggH", "VBF", "ZH", "WH"):
            s = str(self.productionmode)
            if self.productionmode == "VBF": s = "VBFH"
            if self.hypothesis == "0+": result = "{}0PM_M125".format(s)
            if self.hypothesis == "a2": result = "{}0PH_M125".format(s)
            if self.hypothesis == "0-": result = "{}0M_M125".format(s)
            if self.hypothesis == "L1": result = "{}0L1_M125".format(s)
            if self.hypothesis == "L1Zg": result = "{}0L1Zg_M125".format(s)
            if self.hypothesis in ("fa20.5", "fa2prod0.5"): result = "{}0PHf05ph0_M125".format(s)
            if self.hypothesis in ("fa30.5", "fa3prod0.5"): result = "{}0Mf05ph0_M125".format(s)
            if self.hypothesis in ("fL10.5", "fL1prod0.5"): result = "{}0L1f05ph0_M125".format(s)
            if self.hypothesis in ("fL1Zg0.5", "fL1Zgprod0.5"): result = "{}0L1Zgf05ph0_M125".format(s)
            if self.extension is not None: result += "_" + str(self.extension)
            return result
        if self.productionmode == "ttH":
            if self.hffhypothesis == "Hff0+": result = "ttH0PM_M125"
            if self.hffhypothesis == "Hff0-": result = "ttH0M_M125"
            if self.hffhypothesis == "fCP0.5": result = "ttH0Mf05ph0_M125"
            if self.extension is not None: result += "_" + str(self.extension)
            return result
        if self.productionmode == "bbH": return "bbH125"
        if self.productionmode == "tqH": return "tqH125"
        if self.productionmode == "ggZZ":
            return "ggTo{}_Contin_MCFM701".format(self.flavor)
        if self.productionmode == "VBF bkg":
            return "{}_mc_VBFToContinToZZTo4l_ecdaf558".format(self.production.year)
        if self.productionmode == "qqZZ":
            result = "ZZTo4l"
            if self.extension is not None: result += str(self.extension)
            return result
        if self.productionmode in ("TTZZ", "ZZZ", "WZZ", "WWZ", "TTWW", "TTZJets_M10_MLM", "TTZToLLNuNu_M10", "TTZToLL_M1to10_MLM") + tuple(ggZZoffshellproductionmodes):
            result = str(self.productionmode)
            if self.extension is not None: result += str(self.extension)
            return result
        if self.productionmode == "ZX" or self.productionmode == "data":
            return "AllData"
        raise self.ValueError("CJLSTdirname")

    def CJLSTfile(self):
        if self.production.LHE: raise ValueError("Can't get the CJLST file when in LHE mode!")
        return os.path.join(self.CJLSTmaindir(), self.CJLSTdirname(), "ZZ4lAnalysis.root")

    @property
    def LHEfile(self):
        if not self.production.LHE: raise ValueError("Can't get the lhe file when not in LHE mode!")
        raise self.ValueError("LHEfile")

    def withdiscriminantsfile(self):
        if self.copyfromothersample: return self.copyfromothersample.withdiscriminantsfile()
        result = os.path.join(config.repositorydir, "step3_withdiscriminants", str(self.production), "{}.root".format(self).replace(" ", ""))
        return result
        raise self.ValueError("withdiscriminantsfile")

    @staticmethod
    @cache
    def effectiveentries(reweightfrom, reweightto):
        from utilities import tfiles
        if (reweightto.productionmode in ("ggZZ", "VBF bkg", "ZX", "data")
                or reweightfrom.alternategenerator is not None
                or reweightfrom.production.LHE):
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

    @property
    def copyfromothersample(self):
        otherproduction = None

        kwargs = {enum.enumname: getattr(self, enum.enumname) for enum in type(self).needenums}

        if self.production == "190821_2016" and self.pythiasystematic in ("ScaleUp", "ScaleDown"):
            del kwargs["pythiasystematic"]
            kwargs["extension"] = "ext"
            otherproduction = "190821_2017"
        if self == Sample("190821_2017", "WplusH", "0+", "POWHEG", "ext"):
            otherproduction = "190821_2018"
            del kwargs["extension"]
        if self == Sample("190821_2017", "ggH", "0+", "POWHEG", "ext"):
            otherproduction = "190821_2018"
            del kwargs["extension"]
        if self == Sample("190821_2017", "ggH", "0+", "POWHEG", "TuneUp"):
            otherproduction = "190821_2016"
        if self == Sample("190821_2017", "ggH", "0+", "POWHEG", "TuneDown"):
            otherproduction = "190821_2016"
        if self == Sample("190821_2018", "ggH", "0+", "POWHEG", "TuneUp"):
            otherproduction = "190821_2017"
        if self == Sample("190821_2018", "ggH", "0+", "POWHEG", "TuneDown"):
            otherproduction = "190821_2017"
        if self == Sample("190821_2017", "ZH", "0+", "POWHEG", "TuneUp"):
            otherproduction = "190821_2018"
        if self == Sample("190821_2017", "ZH", "0+", "POWHEG", "TuneDown"):
            otherproduction = "190821_2018"
        if self == Sample("190821_2017", "WplusH", "0+", "POWHEG", "TuneUp"):
            otherproduction = "190821_2018"
        if self.productionmode == "ggZZ" and self.production == "190821_2018" and self.flavor == "2mu2tau":
            otherproduction = "190821_2017"

        if self.production == "200205_2016" and self.pythiasystematic in ("ScaleUp", "ScaleDown"):
            del kwargs["pythiasystematic"]
            kwargs["extension"] = "ext"
            otherproduction = "200205_2017"

        if otherproduction is None: return None

        kwargs["production"] = otherproduction
        return Sample(*kwargs.values())

    @property
    def pdf(self): return self.reweightingsamplewithpdf.pdf

    @property
    @generatortolist
    def alternateweights(self):
      for alternateweight in AlternateWeight.items():
        if alternateweight != "1" and self.production.GEN: continue
        if alternateweight != "1" and self.productionmode == "ZX": continue
        if alternateweight.isTHUggH and self.productionmode != "ggH": continue
        if alternateweight in ("EWcorrUp", "EWcorrDn") and self.productionmode != "qqZZ": continue
        if (
          alternateweight == "PythiaScaleUp" or alternateweight == "PythiaScaleDown"
        ) and (
          self.production.year == 2016
          or self.production.year == 2017 and self.extension is None
        ): continue
        yield alternateweight

@cache_file(os.path.join(config.repositorydir, "data", "invertedmatrices.pkl"), lambda x: bytes(x.data))
def invertnumpymatrix(numpymatrix):
    return numpymatrix.I

class SampleBasis(MultiEnum):
    enums = [ProductionMode, Analysis]
    def __init__(self, hypotheses, *args):
        self.hypotheses = [Hypothesis(_) for _ in hypotheses]
        super(SampleBasis, self).__init__(*args)

    def __hash__(self):
        return hash((super(SampleBasis, self).__hash__(), tuple(self.hypotheses)))
    def __eq__(self, other):
        return super(SampleBasis, self).__eq__(other) and self.hypotheses == other.hypotheses
    #__ne__ is already defined as not == in MultiEnum

    def check(self, *args):
        args = (self.hypotheses,)+args
        if self.productionmode in ("ggH", "ttH", "bbH"):
            if self.analysis.dimensions == 4:
                dimension = 15
            elif self.analysis.dimensions == 3:
                dimension = 10
            elif self.analysis.dimensions == 2:
                dimension = 6
            else:
                dimension = 3
        elif self.productionmode in ("VBF", "WH", "ZH"):
            if self.analysis.dimensions == 4:
                dimension = 70
            elif self.analysis.dimensions == 3:
                dimension = 35
            elif self.analysis.dimensions == 2:
                assert False  #have to figure this out
            else:
                dimension = 5
        else:
            raise ValueError("Bad productionmode {}\n{}".format(self.productionmode, args))

        if len(self.hypotheses) != dimension:
            raise ValueError("Wrong number of hypotheses ({}, should be {})\n{}".format(len(self.hypotheses), dimension, args))
        if len(set(self.hypotheses)) != len(self.hypotheses):
            raise ValueError("Duplicate hypothesis\n{}".format(args))

    def scaleby(self, i):
        if self.productionmode in ("VBF", "ZH", "WH", "ggH", "ttH", "bbH"):
            return {
                Analysis("fa3"): 1,
                Analysis("fa2"): 1,
                Analysis("fL1"): 1e-4,
                Analysis("fL1Zg"): 1e-4,
            }[self.analysis.fais[i]]
        assert False

    @property
    @cache
    def matrix(self):
        hffhypothesis = None
        if self.productionmode == "ttH": hffhypothesis = "Hff0+"
        dimension = len(self.hypotheses)
        samples = [ReweightingSample(self.productionmode, _, hffhypothesis) for _ in self.hypotheses]
        if self.analysis.dimensions == 4:
            if dimension == 15:
                maxpower = 2
            elif dimension == 70:
                maxpower = 4
            else:
                assert False, dimension
            assert len(self.analysis.couplingnames) == 4
            return numpy.matrix(
                                [
                                 [
                                  sample.g1**(maxpower-i-j-k-l)
                                  * (getattr(sample, self.analysis.couplingnames[0]) * self.scaleby(0))**i
                                  * (getattr(sample, self.analysis.couplingnames[1]) * self.scaleby(1))**j
                                  * (getattr(sample, self.analysis.couplingnames[2]) * self.scaleby(2))**k
                                  * (getattr(sample, self.analysis.couplingnames[3]) * self.scaleby(3))**l
                                     for l in range(maxpower+1)
                                     for k in range(maxpower+1-l)
                                     for j in range(maxpower+1-k-l)
                                     for i in range(maxpower+1-j-k-l)
                                 ]
                                    for sample in samples
                                ]
                               )

        elif self.analysis.dimensions == 3:
            if dimension == 10:
                maxpower = 2
            elif dimension == 35:
                maxpower = 4
            else:
                assert False, dimension
            assert len(self.analysis.couplingnames) == 3
            assert self.analysis.isEFT
            return numpy.matrix(
                                [
                                 [
                                  sample.g1**(maxpower-i-j-k)
                                  * (getattr(sample, self.analysis.couplingnames[0].replace("g", "ghz")) * self.scaleby(0))**i
                                  * (getattr(sample, self.analysis.couplingnames[1].replace("g", "ghz")) * self.scaleby(1))**j
                                  * (getattr(sample, self.analysis.couplingnames[2].replace("g", "ghz")) * self.scaleby(2))**k
                                     for k in range(maxpower+1)
                                     for j in range(maxpower+1-k)
                                     for i in range(maxpower+1-j-k)
                                 ]
                                    for sample in samples
                                ]
                               )

        elif self.analysis.dimensions == 2:
            assert dimension == 6  #not VBF or VH
            maxpower = 2
            return numpy.matrix(
                                [
                                 [
                                  sample.g1**(maxpower-i-j)
                                  * (getattr(sample, self.analysis.couplingnames[0]) * self.scaleby(0))**i
                                  * (getattr(sample, self.analysis.couplingnames[1]) * self.scaleby(1))**j
                                     for j in range(maxpower+1)
                                     for i in range(maxpower+1-j)
                                 ]
                                    for sample in samples
                                ]
                               )
        else:
            maxpower = dimension-1
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
        self.matrix.flags.writeable = False
        return invertnumpymatrix(self.matrix)

def allsamples(doxcheck=True):
    if doxcheck: __xcheck(*allsamples(doxcheck=False))

    for production in productions:
        if production.GEN:
            for productionmode in "ggH", "VBF", "ZH", "WH", "bbH":
                for hypothesis in ProductionMode(productionmode).generatedhypotheses(production):
                    yield Sample(productionmode, hypothesis, production)
            for hypothesis in hffhypotheses:
                yield Sample("ggH", "JHUGen", hypothesis, "0+", production)
                yield Sample("ttH", hypothesis, "0+", production)
            yield Sample("tqH", "0+", "Hff0+", production)
            yield Sample("ggH", "0+", "MINLO", production)
            for productionmode in "ggH", "VBF", "ZH", "WplusH", "WminusH":
                yield Sample(productionmode, "0+", "POWHEG", production)
                yield Sample(productionmode, "0+", "POWHEG", "ext", production)
            yield Sample("ttH", "Hff0+", "0+", "POWHEG", production)
            yield Sample("ttH", "Hff0+", "0+", "POWHEG", "ext", production)
            yield Sample("qqZZ", production)
            yield Sample("qqZZ", production, "ext2")
            yield Sample("data", production)  #(not real data, just an empty dummy file)
            continue

        if production.LHE:
            for hypothesis in "0+", "L1Zg", "fL1Zg0.5", "fL10.5fL1Zg0.5":
                yield Sample("ggH", hypothesis, production)
            yield Sample("qqZZ", production)
            yield Sample("data", production)  #(not real data, just an empty dummy file)
            continue

        for productionmode in "ggH", "VBF", "ZH", "WH", "bbH":
            for hypothesis in ProductionMode(productionmode).generatedhypotheses(production):
                yield Sample(productionmode, hypothesis, production)
                if production.year == 2017 and (
                   productionmode == "ggH" and hypothesis in (
                       "fL1Zg0.5", "fa30.5",
                   ) or productionmode == "VBF" and hypothesis in (
                       "fL1prod0.5", "a3", "fa3prod0.5", "a2", "a1",
                   ) or productionmode == "ZH" and hypothesis in (
                       "a1", "a2", "a3", "L1", "L1Zg", "fa3prod0.5", "fL1Zgprod0.5",
                   ) or productionmode == "WH"
                ):
                    yield Sample(productionmode, hypothesis, production, "ext1")
        for hypothesis in hffhypotheses:
            yield Sample("ggH", "JHUGen", hypothesis, "0+", production)
            yield Sample("ggH", "JHUGen", hypothesis, "0+", production, "ext1")
            if production.year == 2017: yield Sample("ggH", "JHUGen", hypothesis, "0+", production, "ext2")
            yield Sample("ggH", hypothesis, "0+", production, "MCatNLO")
            yield Sample("ttH", hypothesis, "0+", production)
            if production.year == 2017: yield Sample("ttH", hypothesis, "0+", production, "ext1")
        yield Sample("tqH", "0+", "Hff0+", production)
        for flavor in flavors:
            yield Sample("ggZZ", flavor, production)
        yield Sample("VBF bkg", production)
        yield Sample("ggH", "0+", "NNLOPS", production)
        yield Sample("ggH", "0+", "MINLO", production)
        for productionmode in "ggH", "VBF", "ZH", "WplusH", "WminusH":
            yield Sample(productionmode, "0+", "POWHEG", production)
            if production.year == 2017: yield Sample(productionmode, "0+", "POWHEG", "ext", production)
        yield Sample("ttH", "Hff0+", "0+", "POWHEG", production)
        if production.year == 2017: yield Sample("ttH", "Hff0+", "0+", "POWHEG", "ext", production)
        for systematic in pythiasystematics:
            if systematic in ("ScaleUp", "ScaleDown") and production.year >= 2017: continue
            for productionmode in "ggH", "VBF", "ZH", "WplusH", "WminusH":
                yield Sample(productionmode, "0+", "POWHEG", production, systematic)
            yield Sample("ttH", "Hff0+", "0+", "POWHEG", production, systematic)
            yield Sample("ggH", "0+", "MINLO", production, systematic)
        yield Sample("qqZZ", production)
        if production.year == 2018:
            yield Sample("qqZZ", "ext1", production)
            yield Sample("qqZZ", "ext2", production)
        else:
            yield Sample("qqZZ", "ext", production)
        for _ in "TTZZ", "ZZZ", "WZZ", "WWZ", "TTWW", "TTZToLL_M1to10_MLM":
            yield Sample(_, production)
        for _ in "TTZJets_M10_MLM", "TTZToLLNuNu_M10":
            yield Sample(_, production)
            yield Sample(_, production, "ext1")
        yield Sample("ZX", production)
        yield Sample("data", production)

        for _ in ggZZoffshellproductionmodes:
            yield Sample(_, production)

@cache_file(os.path.join(config.repositorydir, "data", "samples_xcheck.pkl"))
def __xcheck(*samples):
    if os.path.exists(os.path.join(config.repositorydir, "data", "samples_xcheck.pkl")) and not os.path.exists(os.path.join(config.repositorydir, "data", "samples_xcheck.pkl.tmp")):
        os.remove(os.path.join(config.repositorydir, "data", "samples_xcheck.pkl")) #no need to have >1 list in the cache

    def CJLSTfile(sample):
        if sample.production.LHE:
             return sample.LHEfile
        else:
             return sample.CJLSTfile()

    CJ2wd = WriteOnceDict("Multiple withdiscriminantsfiles:\n  {oldvalue}\n  {newvalue}\nfrom the same CJLST file:\n  {key}")
    wd2CJ = WriteOnceDict("Multiple CJLST files:\n  {oldvalue}\n  {newvalue}\nthat give the same withdiscriminantsfile:\n  {key}")

    for sample in samples:
        if sample.productionmode == "data": continue
        if sample.copyfromothersample: continue
        print sample
        try:
            CJ2wd[CJLSTfile(sample)] = sample.withdiscriminantsfile()
            wd2CJ[sample.withdiscriminantsfile()] = CJLSTfile(sample)
        except KeyError as e:
            print e
            raise

if __name__ == "__main__":
    for s in (
      ReweightingSample("ggH", "fa2dec0.5fa3dec0.5"),
      ReweightingSample("VBF", "fa2prod0.5fa3prod0.5"),
      ReweightingSample("ZH",  "fa2proddec-0.33fa3proddec0.33"),
      ReweightingSample("VBF", "fa2proddec0.25fa3proddec0.25fL1proddec0.25"),
    ):
        print s, s.g1, s.g2, s.g4, s.g1prime2, s.ghzgs1prime2
