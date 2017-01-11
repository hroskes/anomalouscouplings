#!/usr/bin/env python
from collections import Counter
from helperstuff import config
from helperstuff.CJLSTscripts import categoryIchep16
from helperstuff.enums import Hypothesis, purehypotheses
from helperstuff.samples import Sample
from helperstuff.treewrapper import dummyfloatstar
from helperstuff.utilities import tfiles

def category(tree, hypothesis):
    hypothesis = Hypothesis(hypothesis)
    nExtraLep = tree.nExtraLep
    nExtraZ = tree.nExtraZ
    nCleanedJetsPt30 = tree.nCleanedJetsPt30
    nCleanedJetsPt30BTagged = tree.nCleanedJetsPt30BTagged
    jetQGLikelihood = tree.JetQGLikelihood.data()
    jetPhi = tree.JetPhi.data()
    ZZMass = tree.ZZMass
    p_JQCD_SIG_ghg2_1_JHUGen_JECNominal = tree.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal
    p_JVBF_SIG_ghv1_1_JHUGen_JECNominal = tree.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal
    pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal = tree.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal


    if nCleanedJetsPt30 == 0:
        jetQGLikelihood = jetPhi = dummyfloatstar

    if hypothesis == "0+": coupling = "gh{}1_1"
    elif hypothesis == "0-": coupling = "gh{}4_1"
    elif hypothesis == "a2": coupling = "gh{}2_1"
    elif hypothesis == "L1": coupling = "gh{}1prime2_1E4"
    else: assert False

    p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal = getattr(tree, "p_JJVBF_SIG_{}_JHUGen_JECNominal".format(coupling).format("v"))
    p_HadWH_SIG_ghw1_1_JHUGen_JECNominal = getattr(tree, "p_HadWH_SIG_{}_JHUGen_JECNominal".format(coupling).format("w"))
    p_HadZH_SIG_ghz1_1_JHUGen_JECNominal = getattr(tree, "p_HadZH_SIG_{}_JHUGen_JECNominal".format(coupling).format("z"))

    if hypothesis == "0-":
        p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal = tree.p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal
    else:
        p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal = tree.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal
    

    result = categoryIchep16(
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
          p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
          p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
          jetPhi,
          ZZMass,
          config.useQGTagging,
         )

    if hypothesis == "0+":
        assert result == tree.category

    return result

def count(*sample):
    sample = Sample(*sample)
    t = tfiles[sample.withdiscriminantsfile()].candTree
    counters = {hypothesis: Counter({_: 0 for _ in range(6)}) for hypothesis in purehypotheses}
    total = 0
    weightname = sample.weightname()
    countersitems = counters.items()
    for entry in t:
        weight = getattr(t, weightname)
        for hypothesis, counter in countersitems:
            counter[category(t, hypothesis)] += weight
        total += weight

    for counter in counters.values():
        for k, v in counter.iteritems():
            counter[k] = v/total

    return counters

if __name__ == "__main__":
    counters = count("ggH", "0+", "161221")
    for hypothesis in purehypotheses:
        print hypothesis
        for category in range(6):
            print "{}    {}%".format(category, counters[hypothesis][category]*100)
        print
