from abc import ABCMeta, abstractmethod, abstractproperty
from CJLSTscripts import categoryIchep16, VBF2jTaggedIchep16, VHHadrTaggedIchep16
import config
from enums import Category, categories, HffHypothesis, Hypothesis, JECSystematic
from samples import ArbitraryCouplingsSample
import re

class BaseCategorization(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def category_function_name(self): pass
    @abstractmethod
    def get_category_function(self_categorization): pass

    def __hash__(self): return hash(self.category_function_name)
    def __eq__(self, other): return self.category_function_name == other.category_function_name
    def __ne__(self, other): return not self == other

class BaseSingleCategorization(BaseCategorization):
    def __init__(self, JEC):
        self.JEC = JECSystematic(JEC)
    @abstractproperty
    def g1(self): pass
    @abstractproperty
    def g2(self): pass
    @abstractproperty
    def g4(self): pass
    @abstractproperty
    def g1prime2(self): pass
    @abstractproperty
    def ghzgs1prime2(self): pass
    @abstractproperty
    def ghg2(self): pass
    @abstractproperty
    def ghg4(self): pass
    @abstractproperty
    def pHJJ_function_name(self): pass
    @abstractproperty
    def pVBF_function_name(self): pass
    @abstractproperty
    def pZH_function_name(self): pass
    @abstractproperty
    def pWH_function_name(self): pass

    @property
    def pHJ_variable_name(self): return "p_JQCD_SIG_ghg2_1_JHUGen_{}".format(self.JEC)
    @property
    def pVBF1j_variable_name(self): return "p_JVBF_SIG_ghv1_1_JHUGen_{}".format(self.JEC)
    @property
    def pAux_variable_name(self): return "pAux_JVBF_SIG_ghv1_1_JHUGen_{}".format(self.JEC)

    @property
    def njets_variable_name(self): return "nCleanedJetsPt30" + self.JEC.appendname.replace("JEC", "jec")
    @property
    def nbtagged_variable_name(self): return "nCleanedJetsPt30BTagged_bTagSF" + self.JEC.appendname.replace("JEC", "jec")
    @property
    def phi_variable_name(self): return "jetPhi" + self.JEC.appendname.replace("JEC", "jec")
    @property
    def QGL_variable_name(self): return "jetQGLikelihood" + self.JEC.appendname.replace("JEC", "jec")

    @staticmethod
    def get_p_function(terms, multiplier, name):
        def function(self_tree):
            return sum(getattr(self_tree, k)*v for k, v in terms) * multiplier
        function.__name__ = name
        return function

    def get_pHJJ_function(self):
        terms = {
                 "p_JJQCD_SIG_ghg2_1_JHUGen_{}".format(self.JEC): self.ghg2**2,
                 "p_JJQCD_SIG_ghg4_1_JHUGen_{}".format(self.JEC): self.ghg4**2,
                 "p_JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_{}".format(self.JEC): self.ghg2*self.ghg4,
                }
        terms = tuple((k, v) for k, v in terms.iteritems() if v)
        multiplier = ArbitraryCouplingsSample("HJJ", 1, 0, 0, 0, 0, ghg2=1, ghg4=0).JHUxsec / ArbitraryCouplingsSample("HJJ", 1, 0, 0, 0, 0, ghg2=self.ghg2, ghg4=self.ghg4).JHUxsec
        return self.get_p_function(terms, multiplier, self.pHJJ_function_name)

    def get_pVBF_function(self):
        terms = {
                 "p_JJVBF_SIG_ghv1_1_JHUGen_{}".format(self.JEC): self.g1**2,
                 "p_JJVBF_SIG_ghv2_1_JHUGen_{}".format(self.JEC): self.g2**2,
                 "p_JJVBF_SIG_ghv4_1_JHUGen_{}".format(self.JEC): self.g4**2,
                 "p_JJVBF_SIG_ghv1prime2_1_JHUGen_{}".format(self.JEC): self.g1prime2**2,
                 "p_JJVBF_SIG_ghzgs1prime2_1_JHUGen_{}".format(self.JEC): self.ghzgs1prime2**2,
                 "p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_{}".format(self.JEC): self.g1 * self.g2,
                 "p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_{}".format(self.JEC): self.g1 * self.g4,
                 "p_JJVBF_SIG_ghv1_1_ghv1prime2_1_JHUGen_{}".format(self.JEC): self.g1 * self.g1prime2,
                 "p_JJVBF_SIG_ghv1_1_ghzgs1prime2_1_JHUGen_{}".format(self.JEC): self.g1 * self.ghzgs1prime2,
                }
        terms = tuple((k, v) for k, v in terms.iteritems() if v)
        multiplier = (
                      (ArbitraryCouplingsSample("VBF", 1, 0, 0, 0, 0).xsec / ArbitraryCouplingsSample("VBF", self.g1, self.g2, self.g4, self.g1prime2, self.ghzgs1prime2).xsec)
                     /
                      (ArbitraryCouplingsSample("ggH", 1, 0, 0, 0, 0).xsec / ArbitraryCouplingsSample("ggH", self.g1, self.g2, self.g4, self.g1prime2, self.ghzgs1prime2).xsec)
                     )
        return self.get_p_function(terms, multiplier, self.pVBF_function_name)

    def get_pZH_function(self):
        terms = {
                 "p_HadZH_SIG_ghz1_1_JHUGen_{}".format(self.JEC): self.g1**2,
                 "p_HadZH_SIG_ghz2_1_JHUGen_{}".format(self.JEC): self.g2**2,
                 "p_HadZH_SIG_ghz4_1_JHUGen_{}".format(self.JEC): self.g4**2,
                 "p_HadZH_SIG_ghz1prime2_1_JHUGen_{}".format(self.JEC): self.g1prime2**2,
                 "p_HadZH_SIG_ghzgs1prime2_1_JHUGen_{}".format(self.JEC): self.ghzgs1prime2**2,
                 "p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_{}".format(self.JEC): self.g1 * self.g2,
                 "p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_{}".format(self.JEC): self.g1 * self.g4,
                 "p_HadZH_SIG_ghz1_1_ghz1prime2_1_JHUGen_{}".format(self.JEC): self.g1 * self.g1prime2,
                 "p_HadZH_SIG_ghz1_1_ghzgs1prime2_1_JHUGen_{}".format(self.JEC): self.g1 * self.ghzgs1prime2,
                }
        terms = tuple((k, v) for k, v in terms.iteritems() if v)
        multiplier = (
                      (ArbitraryCouplingsSample("ZH", 1, 0, 0, 0, 0).xsec / ArbitraryCouplingsSample("ZH", self.g1, self.g2, self.g4, self.g1prime2, self.ghzgs1prime2).xsec)
                     /
                      (ArbitraryCouplingsSample("ggH", 1, 0, 0, 0, 0).xsec / ArbitraryCouplingsSample("ggH", self.g1, self.g2, self.g4, self.g1prime2, self.ghzgs1prime2).xsec)
                     )
        return self.get_p_function(terms, multiplier, self.pZH_function_name)

    def get_pWH_function(self):
        if self.ghzgs1prime2 != 0:
            def result(self_tree): return 0
            result.__name__ = self.pWH_function_name
            return result
        terms = {
                 "p_HadWH_SIG_ghw1_1_JHUGen_{}".format(self.JEC): self.g1**2,
                 "p_HadWH_SIG_ghw2_1_JHUGen_{}".format(self.JEC): self.g2**2,
                 "p_HadWH_SIG_ghw4_1_JHUGen_{}".format(self.JEC): self.g4**2,
                 "p_HadWH_SIG_ghw1prime2_1_JHUGen_{}".format(self.JEC): self.g1prime2**2,
                 "p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_{}".format(self.JEC): self.g1 * self.g2,
                 "p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_{}".format(self.JEC): self.g1 * self.g4,
                 "p_HadWH_SIG_ghw1_1_ghw1prime2_1_JHUGen_{}".format(self.JEC): self.g1 * self.g1prime2,
                }
        terms = tuple((k, v) for k, v in terms.iteritems() if v)
        multiplier = (
                      (ArbitraryCouplingsSample("WH", 1, 0, 0, 0, 0).xsec / ArbitraryCouplingsSample("WH", self.g1, self.g2, self.g4, self.g1prime2, self.ghzgs1prime2).xsec)
                     /
                      (ArbitraryCouplingsSample("ggH", 1, 0, 0, 0, 0).xsec / ArbitraryCouplingsSample("ggH", self.g1, self.g2, self.g4, self.g1prime2, self.ghzgs1prime2).xsec)
                     )
        return self.get_p_function(terms, multiplier, self.pWH_function_name)

    def get_category_function(self_categorization):
        pHJJ_function_name = self_categorization.pHJJ_function_name
        pVBF_function_name = self_categorization.pVBF_function_name
        pZH_function_name = self_categorization.pZH_function_name
        pWH_function_name = self_categorization.pWH_function_name

        pHJ_variable_name = self_categorization.pHJ_variable_name
        pVBF1j_variable_name = self_categorization.pVBF1j_variable_name
        pAux_variable_name = self_categorization.pAux_variable_name

        njets_variable_name = self_categorization.njets_variable_name
        nbtagged_variable_name = self_categorization.nbtagged_variable_name
        phi_variable_name = self_categorization.phi_variable_name
        QGL_variable_name = self_categorization.QGL_variable_name

        def function(self_tree):
            result = self_categorization.lastvalue = categoryIchep16(
                self_tree.nExtraLep,
                self_tree.nExtraZ,
                  getattr(self_tree, njets_variable_name),
                  getattr(self_tree, nbtagged_variable_name),
                  getattr(self_tree, QGL_variable_name),
                  getattr(self_tree, pHJJ_function_name)(),
                  getattr(self_tree, pHJ_variable_name),
                  getattr(self_tree, pVBF_function_name)(),
                  getattr(self_tree, pVBF1j_variable_name),
                  getattr(self_tree, pAux_variable_name),
                  getattr(self_tree, pWH_function_name)(),
                  getattr(self_tree, pZH_function_name)(),
                  getattr(self_tree, phi_variable_name),
                self_tree.ZZMass,
                config.useQGTagging,
            )
            return result
        function.__name__ = self_categorization.category_function_name
        return function

class ArbitraryCouplingsSingleCategorization(BaseSingleCategorization):
    def __init__(self, g1, g2, g4, g1prime2, ghg2=1, ghg4=0, JEC="Nominal"):
        self.__g1 = g1
        self.__g2 = g2
        self.__g4 = g4
        self.__g1prime2 = g1prime2
        self.__ghzgs1prime2 = ghzgs1prime2
        if len(x for x in (g2, g4, g1prime2, ghzgs1prime2) if x) > w:
            raise ValueError("Can't use more than 1 of g2, g4, g1prime2, ghzgs1prime2")
        self.__ghg2 = ghg2
        self.__ghg4 = ghg4
        super(ArbitraryCouplingsSingleCategorization, self).__init__(JEC)

    @property
    def g1(self): return self.__g1
    @property
    def g2(self): return self.__g2
    @property
    def g4(self): return self.__g4
    @property
    def g1prime2(self): return self.__g1prime2
    @property
    def g1prime2(self): return self.__ghzgs1prime2
    @property
    def ghg2(self): return self.__ghg2
    @property
    def ghg4(self): return self.__ghg4

    @property
    def pHJJ_function_name(self): return "category_p_JJQCD_SIG_{}_JHUGen_{}".format("_".join("{}_{}".format(name, getattr(self, name)) for name in ("ghg2", "ghg4")), self.JEC)
    @property
    def pVBF_function_name(self): return "category_p_JJVBF_SIG_{}_JHUGen_{}".format("_".join("{}_{}".format(name.replace("g", "ghv"), getattr(self, name)) for name in ("g1", "g2", "g4", "g1prime2")), self.JEC)
    @property
    def pZH_function_name(self): return "category_p_HadZH_SIG_{}_JHUGen_{}".format("_".join("{}_{}".format(name.replace("g", "ghz"), getattr(self, name)) for name in ("g1", "g2", "g4", "g1prime2")), self.JEC)
    @property
    def pWH_function_name(self): return "category_p_HadWH_SIG_{}_JHUGen_{}".format("_".join("{}_{}".format(name.replace("g", "ghw"), getattr(self, name)) for name in ("g1", "g2", "g4", "g1prime2")), self.JEC)

    @property
    def category_function_name(self): return "category_{}".format("_".join("{}_{}".format(name, getattr(self, name)) for name in ("ghg2", "ghg4", "g1", "g2", "g4", "g1prime2"))) + self.JEC.appendname

class SingleCategorizationFromSample(BaseSingleCategorization):
    def __init__(self, sample, JEC="Nominal"):
        self.sample = sample
        super(SingleCategorizationFromSample, self).__init__(JEC)

    @property
    def g1(self): return self.sample.g1
    @property
    def g2(self): return self.sample.g2
    @property
    def g4(self): return self.sample.g4
    @property
    def g1prime2(self): return self.sample.g1prime2
    @property
    def ghzgs1prime2(self): return self.sample.ghzgs1prime2
    @property
    def ghg2(self):
        try:
            return self.sample.ghg2
        except ValueError:
            return 1
    @property
    def ghg4(self):
        try:
            return self.sample.ghg4
        except ValueError:
            return 0

    @property
    def hypothesisname(self):
        return self.nicename(self.sample.hypothesis)
    @property
    def hffhypothesisname(self):
        return self.nicename(self.sample.hffhypothesis if self.sample.hffhypothesis is not None else HffHypothesis("Hff0+"))
    def nicename(self, s):
        result = str(s).replace("prod", str(self.sample.productionmode)).replace("+", "P").replace("-", "M").replace(".", "p")
        assert re.match("[\w]+", result)
        return result

    @property
    def pHJJ_function_name(self): return "category_p_JJQCD_SIG_{}_JHUGen_{}".format(self.hffhypothesisname, self.JEC)
    @property
    def pVBF_function_name(self): return "category_p_JJVBF_SIG_{}_JHUGen_{}".format(self.hypothesisname, self.JEC)
    @property
    def pZH_function_name(self): return "category_p_HadZH_SIG_{}_JHUGen_{}".format(self.hypothesisname, self.JEC)
    @property
    def pWH_function_name(self): return "category_p_HadWH_SIG_{}_JHUGen_{}".format(self.hypothesisname, self.JEC)

    @property
    def category_function_name(self):
        result = "category_{}".format(self.hypothesisname) + self.JEC.appendname
        if self.sample.hffhypothesis:
            result += "_" + self.hffhypothesisname
        return result

    def __str__(self):
        assert self.sample.hffhypothesis is None or self.sample.hffhypothesis == "Hff0+"
        return self.hypothesisname

class MultiCategorization(BaseCategorization):
    def __init__(self, name, *singles):
        self.name = name
        self.singles = frozenset(singles)
#        for single in self.singles:
#            single.nmulticategorizations += 1
    def __hash__(self):
        return hash(self.singles)
    def __eq__(self, other):
        return self.singles == other.singles
    def __ne__(self, other):
        return self.singles != other.singles
    def __gt__(self, other):
        return self.singles > other.singles
    def __ge__(self, other):
        return self.singles >= other.singles
    def __lt__(self, other):
        return self.singles < other.singles
    def __le__(self, other):
        return self.singles <= other.singles

    def __str__(self):
        return self.name.replace("_", " ")

    def get_category_function(self_categorization):
        singles = self_categorization.singles
        def function(self_tree):
            lastvalues = {single.lastvalue for single in singles}
            if any(_ == VBF2jTaggedIchep16 for _ in lastvalues):
                return VBF2jTaggedIchep16
            if any(_ == VHHadrTaggedIchep16 for _ in lastvalues):
                return VHHadrTaggedIchep16
            assert len(lastvalues) == 1
            return lastvalues.pop()
        function.__name__ = self_categorization.category_function_name
        return function

    @property
    def category_function_name(self): return "category_{}".format(self.name)
