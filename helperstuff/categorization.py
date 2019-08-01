from abc import ABCMeta, abstractmethod, abstractproperty
from CJLSTscripts import categoryMor18, categoryAC19, UntaggedAC19, VBF2jTaggedAC19, VHHadrTaggedAC19, BoostedAC19
import config
from enums import BTagSystematic, Category, categories, HffHypothesis, Hypothesis, JECSystematic
from samples import ArbitraryCouplingsSample
from utilities import setname
import re

class BaseCategorization(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def category_function_name(self): pass
    @abstractmethod
    def get_category_function(self_categorization): pass
    @abstractproperty
    def issystematic(self): pass

    def __hash__(self): return hash((type(self), self.category_function_name))
    def __eq__(self, other): return type(self) == type(other) and self.category_function_name == other.category_function_name
    def __ne__(self, other): return not self == other

class NoCategorization(BaseCategorization):
    @property
    def category_function_name(self): return "category_nocategorization"
    def get_category_function(self_categorization):
        @setname(self_categorization.category_function_name)
        def function(self_tree): return UntaggedAC19
        return function
    @property
    def issystematic(self): return False

class BaseSingleCategorization(BaseCategorization):
    def __init__(self, JEC, btag):
        self.JEC = JECSystematic(JEC)
        self.btag = BTagSystematic(btag)
        if self.btag != "Nominal" and self.JEC != "Nominal":
            raise ValueError("Can't have systematics on both btag and JEC at the same time! {}, {}".format(self.btag, self.JEC))
    @property
    def issystematic(self): return self.btag != "Nominal" or self.JEC != "Nominal"

    @property
    def pHJ_variable_name(self): return "p_JQCD_SIG_ghg2_1_JHUGen_{}".format(self.JEC)
    @property
    def pVBF1j_variable_name(self): return "p_JVBF_SIG_ghv1_1_JHUGen_{}".format(self.JEC)
    @property
    def pAux_variable_name(self): return "pAux_JVBF_SIG_ghv1_1_JHUGen_{}".format(self.JEC)
    @property
    def p_HadWH_mavjj_variable_name(self): return "p_HadWH_mavjj_{}".format(self.JEC)
    @property
    def p_HadWH_mavjj_true_variable_name(self): return "p_HadWH_mavjj_true_{}".format(self.JEC)
    @property
    def p_HadZH_mavjj_variable_name(self): return "p_HadZH_mavjj_{}".format(self.JEC)
    @property
    def p_HadZH_mavjj_true_variable_name(self): return "p_HadZH_mavjj_true_{}".format(self.JEC)

    @property
    def njets_variable_name(self): return "nCleanedJetsPt30" + self.JEC.njetsappendname
    @property
    def nbtagged_variable_name(self): return "nCleanedJetsPt30BTagged" + self.btag.njetsappendname + self.JEC.njetsappendname
    @property
    def phi_variable_name(self): return "jetPhi" + self.JEC.njetsappendname
    @property
    def QGL_variable_name(self): return "jetQGLikelihood" + self.JEC.njetsappendname

    @abstractmethod
    def get_pHJJ_function(self): pass
    @abstractmethod
    def get_pVBF_function(self): pass
    @abstractmethod
    def get_pZH_function(self): pass
    @abstractmethod
    def get_pWH_function(self): pass

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

        p_HadWH_mavjj_variable_name = self_categorization.p_HadWH_mavjj_variable_name
        p_HadWH_mavjj_true_variable_name = self_categorization.p_HadWH_mavjj_true_variable_name
        p_HadZH_mavjj_variable_name = self_categorization.p_HadZH_mavjj_variable_name
        p_HadZH_mavjj_true_variable_name = self_categorization.p_HadZH_mavjj_true_variable_name

        @setname(self_categorization.category_function_name)
        def function(self_tree):
            result = self_categorization.lastvalue = categoryAC19(
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
                  getattr(self_tree, p_HadWH_mavjj_variable_name),
                  getattr(self_tree, p_HadWH_mavjj_true_variable_name),
                  getattr(self_tree, p_HadZH_mavjj_variable_name),
                  getattr(self_tree, p_HadZH_mavjj_true_variable_name),
                  getattr(self_tree, phi_variable_name),
                self_tree.ZZMass,
                self_tree.ZZPt,
                self_tree.PFMET,
                config.useVHMETTagged,
                config.useQGTagging,
            )
            return result

        return function

class BaseSingleCategorizationCouplings(BaseSingleCategorization):
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

    @staticmethod
    def get_p_function(terms, multiplier, name):
        @setname(name)
        def function(self_tree):
            return sum(getattr(self_tree, k)*v for k, v in terms) * multiplier
        return function

    def get_pHJJ_function(self):
        terms = {
                 "p_JJQCD_SIG_ghg2_1_JHUGen_{}".format(self.JEC): self.ghg2**2,
                 "p_JJQCD_SIG_ghg4_1_JHUGen_{}".format(self.JEC): self.ghg4**2,
                 "p_JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_{}".format(self.JEC): self.ghg2*self.ghg4,
                }
        terms = tuple((k, v) for k, v in terms.iteritems() if v)
        multiplier = ArbitraryCouplingsSample("HJJ", 1, 0, 0, 0, 0, ghg2=1, ghg4=0, pdf="NNPDF30_lo_as_0130").JHUxsec / ArbitraryCouplingsSample("HJJ", 1, 0, 0, 0, 0, ghg2=self.ghg2, ghg4=self.ghg4, pdf="NNPDF30_lo_as_0130").JHUxsec
        return self.get_p_function(terms, multiplier.nominal_value, self.pHJJ_function_name)


    def get_pVBF_function(self):
        terms = {
                 "p_JJVBF_SIG_ghv1_1_JHUGen_{}".format(self.JEC): self.g1**2,
                 "p_JJVBF_SIG_ghv2_1_JHUGen_{}".format(self.JEC): self.g2**2,
                 "p_JJVBF_SIG_ghv4_1_JHUGen_{}".format(self.JEC): self.g4**2,
                 "p_JJVBF_SIG_ghv1prime2_1_JHUGen_{}".format(self.JEC): self.g1prime2**2,
                 "p_JJVBF_SIG_ghza1prime2_1_JHUGen_{}".format(self.JEC): self.ghzgs1prime2**2,
                 "p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_{}".format(self.JEC): self.g1 * self.g2,
                 "p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_{}".format(self.JEC): self.g1 * self.g4,
                 "p_JJVBF_SIG_ghv1_1_ghv1prime2_1_JHUGen_{}".format(self.JEC): self.g1 * self.g1prime2,
                 "p_JJVBF_SIG_ghv1_1_ghza1prime2_1_JHUGen_{}".format(self.JEC): self.g1 * self.ghzgs1prime2,
                }
        terms = tuple((k, v) for k, v in terms.iteritems() if v)
        multiplier = (
                      (ArbitraryCouplingsSample("VBF", 1, 0, 0, 0, 0, pdf="NNPDF30_lo_as_0130").xsec / ArbitraryCouplingsSample("VBF", self.g1, self.g2, self.g4, self.g1prime2, self.ghzgs1prime2, pdf="NNPDF30_lo_as_0130").xsec)
                     /
                      (ArbitraryCouplingsSample("ggH", 1, 0, 0, 0, 0, pdf="NNPDF30_lo_as_0130").xsec / ArbitraryCouplingsSample("ggH", self.g1, self.g2, self.g4, self.g1prime2, self.ghzgs1prime2, pdf="NNPDF30_lo_as_0130").xsec)
                     )
        return self.get_p_function(terms, multiplier.nominal_value, self.pVBF_function_name)

    def get_pZH_function(self):
        terms = {
                 "p_HadZH_SIG_ghz1_1_JHUGen_{}".format(self.JEC): self.g1**2,
                 "p_HadZH_SIG_ghz2_1_JHUGen_{}".format(self.JEC): self.g2**2,
                 "p_HadZH_SIG_ghz4_1_JHUGen_{}".format(self.JEC): self.g4**2,
                 "p_HadZH_SIG_ghz1prime2_1_JHUGen_{}".format(self.JEC): self.g1prime2**2,
                 "p_HadZH_SIG_ghza1prime2_1_JHUGen_{}".format(self.JEC): self.ghzgs1prime2**2,
                 "p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_{}".format(self.JEC): self.g1 * self.g2,
                 "p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_{}".format(self.JEC): self.g1 * self.g4,
                 "p_HadZH_SIG_ghz1_1_ghz1prime2_1_JHUGen_{}".format(self.JEC): self.g1 * self.g1prime2,
                 "p_HadZH_SIG_ghz1_1_ghza1prime2_1_JHUGen_{}".format(self.JEC): self.g1 * self.ghzgs1prime2,
                }
        terms = tuple((k, v) for k, v in terms.iteritems() if v)
        multiplier = (
                      (ArbitraryCouplingsSample("ZH", 1, 0, 0, 0, 0, pdf="NNPDF30_lo_as_0130").xsec / ArbitraryCouplingsSample("ZH", self.g1, self.g2, self.g4, self.g1prime2, self.ghzgs1prime2, pdf="NNPDF30_lo_as_0130").xsec)
                     /
                      (ArbitraryCouplingsSample("ggH", 1, 0, 0, 0, 0, pdf="NNPDF30_lo_as_0130").xsec / ArbitraryCouplingsSample("ggH", self.g1, self.g2, self.g4, self.g1prime2, self.ghzgs1prime2, pdf="NNPDF30_lo_as_0130").xsec)
                     )
        return self.get_p_function(terms, multiplier.nominal_value, self.pZH_function_name)

    def get_pWH_function(self):
        if self.ghzgs1prime2 != 0:
            @setname(self.pWH_function_name)
            def result(self_tree): return 0
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
                      (ArbitraryCouplingsSample("WH", 1, 0, 0, 0, 0, pdf="NNPDF30_lo_as_0130").xsec / ArbitraryCouplingsSample("WH", self.g1, self.g2, self.g4, self.g1prime2, self.ghzgs1prime2, pdf="NNPDF30_lo_as_0130").xsec)
                     /
                      (ArbitraryCouplingsSample("ggH", 1, 0, 0, 0, 0, pdf="NNPDF30_lo_as_0130").xsec / ArbitraryCouplingsSample("ggH", self.g1, self.g2, self.g4, self.g1prime2, self.ghzgs1prime2, pdf="NNPDF30_lo_as_0130").xsec)
                     )
        return self.get_p_function(terms, multiplier.nominal_value, self.pWH_function_name)

class ArbitraryCouplingsSingleCategorization(BaseSingleCategorizationCouplings):
    def __init__(self, g1, g2, g4, g1prime2, ghg2=1, ghg4=0, JEC="Nominal", btag="Nominal"):
        self.__g1 = g1
        self.__g2 = g2
        self.__g4 = g4
        self.__g1prime2 = g1prime2
        self.__ghzgs1prime2 = ghzgs1prime2
        if len(x for x in (g2, g4, g1prime2, ghzgs1prime2) if x) > w:
            raise ValueError("Can't use more than 1 of g2, g4, g1prime2, ghzgs1prime2")
        self.__ghg2 = ghg2
        self.__ghg4 = ghg4
        super(ArbitraryCouplingsSingleCategorization, self).__init__(JEC=JEC, btag=btag)

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
    def category_function_name(self): return "category_{}".format("_".join("{}_{}".format(name, getattr(self, name)) for name in ("ghg2", "ghg4", "g1", "g2", "g4", "g1prime2"))) + self.btag.appendname + self.JEC.appendname

class SingleCategorizationFromSample(BaseSingleCategorizationCouplings):
    def __init__(self, sample, JEC="Nominal", btag="Nominal"):
        self.sample = sample
        super(SingleCategorizationFromSample, self).__init__(JEC=JEC, btag=btag)

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
        result = "category_{}".format(self.hypothesisname) + self.btag.appendname + self.JEC.appendname
        if self.sample.hffhypothesis:
            result += "_" + self.hffhypothesisname
        return result

    def __str__(self):
        assert self.sample.hffhypothesis is None or self.sample.hffhypothesis == "Hff0+"
        return self.hypothesisname

    def __repr__(self):
        return "{}({!r}, {!r}, {!r})".format(type(self).__name__, self.sample, self.JEC, self.btag)

class SingleCategorizationgm4l(BaseSingleCategorization):
    def __init__(self, hypothesis, JEC="Nominal", btag="Nominal"):
        self.hypothesis = Hypothesis(hypothesis)
        if not self.hypothesis.ispure: raise ValueError("Invalid hypothesis {}".format(hypothesis))
        super(SingleCategorizationgm4l, self).__init__(JEC=JEC, btag=btag)

    @property
    def pHJJ_function_name(self): return "category_p_JJQCD_SIG_Hff0P_JHUGen_{}_gm4l".format(self.JEC)
    @property
    def pVBF_function_name(self): return "category_p_JJVBF_SIG_{}_JHUGen_{}_gm4l".format(self.hypothesisname, self.JEC)
    @property
    def pZH_function_name(self): return "category_p_HadZH_SIG_{}_JHUGen_{}_gm4l".format(self.hypothesisname, self.JEC)
    @property
    def pWH_function_name(self): return "category_p_HadWH_SIG_{}_JHUGen_{}_gm4l".format(self.hypothesisname, self.JEC)

    @staticmethod
    def get_p_function(probability, gconstant, name):
        usegconstant = not (gconstant == 1)
        @setname(name)
        def function(self_tree):
            result = getattr(self_tree, probability)
            if usegconstant:
              result *= getattr(self_tree, gconstant)**2
            return result
        return function

    def get_pHJJ_function(self):
        probability, gconstant = "p_JJQCD_SIG_ghg2_1_JHUGen_{}".format(self.JEC), 1
        return self.get_p_function(probability, gconstant, self.pHJJ_function_name)


    def get_pVBF_function(self):
        if self.hypothesis == "0+": probability, gconstant = "p_JJVBF_SIG_ghv1_1_JHUGen_{}".format(self.JEC), 1
        if self.hypothesis == "a2": probability, gconstant = "p_JJVBF_SIG_ghv2_1_JHUGen_{}".format(self.JEC), "g2VBF_m4l"
        if self.hypothesis == "a3": probability, gconstant = "p_JJVBF_SIG_ghv4_1_JHUGen_{}".format(self.JEC), "g4VBF_m4l"
        if self.hypothesis == "L1": probability, gconstant = "p_JJVBF_SIG_ghv1prime2_1_JHUGen_{}".format(self.JEC), "g1prime2VBF_m4l"
        if self.hypothesis == "L1Zg": probability, gconstant = "p_JJVBF_SIG_ghza1prime2_1_JHUGen_{}".format(self.JEC), "ghzgs1prime2VBF_m4l"
        return self.get_p_function(probability, gconstant, self.pVBF_function_name)

    def get_pZH_function(self):
        if self.hypothesis == "0+": probability, gconstant = "p_HadZH_SIG_ghz1_1_JHUGen_{}".format(self.JEC), 1
        if self.hypothesis == "a2": probability, gconstant = "p_HadZH_SIG_ghz2_1_JHUGen_{}".format(self.JEC), "g2ZH_m4l"
        if self.hypothesis == "a3": probability, gconstant = "p_HadZH_SIG_ghz4_1_JHUGen_{}".format(self.JEC), "g4ZH_m4l"
        if self.hypothesis == "L1": probability, gconstant = "p_HadZH_SIG_ghz1prime2_1_JHUGen_{}".format(self.JEC), "g1prime2ZH_m4l"
        if self.hypothesis == "L1Zg": probability, gconstant = "p_HadZH_SIG_ghza1prime2_1_JHUGen_{}".format(self.JEC), "ghzgs1prime2ZH_m4l"
        return self.get_p_function(probability, gconstant, self.pZH_function_name)

    def get_pWH_function(self):
        if self.hypothesis == "L1Zg":
            @setname(self.pWH_function_name)
            def result(self_tree): return 0
            return result
        if self.hypothesis == "0+": probability, gconstant = "p_HadWH_SIG_ghw1_1_JHUGen_{}".format(self.JEC), 1
        if self.hypothesis == "a2": probability, gconstant = "p_HadWH_SIG_ghw2_1_JHUGen_{}".format(self.JEC), "g2WH_m4l"
        if self.hypothesis == "a3": probability, gconstant = "p_HadWH_SIG_ghw4_1_JHUGen_{}".format(self.JEC), "g4WH_m4l"
        if self.hypothesis == "L1": probability, gconstant = "p_HadWH_SIG_ghw1prime2_1_JHUGen_{}".format(self.JEC), "g1prime2WH_m4l"
        return self.get_p_function(probability, gconstant, self.pWH_function_name)

    def nicename(self, h):
        result = str(h).replace("+", "P").replace("-", "M")
        assert re.match("[\w]+", result)
        return result

    @property
    def hypothesisname(self):
        return self.nicename(self.hypothesis)

    @property
    def category_function_name(self):
        result = "category_" + self.hypothesisname + self.btag.appendname + self.JEC.appendname
        return result

class MultiCategorization(BaseCategorization):
    def __init__(self, name, *singles):
        self.name = name
        self.singles = frozenset(singles)
    @property
    def issystematic(self):
        result = {_.issystematic for _ in self.singles}
        assert len(result) == 1, result
        return result.pop()
    def __hash__(self):
        return hash(self.singles)
    def __eq__(self, other):
        try:
            return self.singles == other.singles
        except AttributeError:
            return self.singles == {other}
    def __ne__(self, other):
        return not self == other
    def __gt__(self, other):
        try:
            return self.singles > other.singles
        except AttributeError:
            return isinstance(other, NoCategorization) or other in self.singles
    def __lt__(self, other):
        try:
            return self.singles > other.singles
        except AttributeError:
            return False
    def __ge__(self, other):
        return self == other or self > other
    def __le__(self, other):
        return self == other or self < other

    def __str__(self):
        return self.name.replace("_", " ")

    def get_category_function(self_categorization):
        singles = self_categorization.singles
        @setname(self_categorization.category_function_name)
        def function(self_tree):
            lastvalues = {single.lastvalue for single in singles}
            if any(_ == VBF2jTaggedAC19 for _ in lastvalues):
                return VBF2jTaggedAC19
            if any(_ == VHHadrTaggedAC19 for _ in lastvalues):
                return VHHadrTaggedAC19
            assert len(lastvalues) == 1
            return lastvalues.pop()
        return function

    @property
    def category_function_name(self): return "category_{}".format(self.name)

    def __repr__(self):
        #return ("{}(" + ", ".join(["{!r}"]*(len(self.singles)+1))+")").format(type(self).__name__, self.name, *self.singles)
        return self.name
