#!/usr/bin/env python
from collections import Counter
from helperstuff import config
from helperstuff import constants
from helperstuff.CJLSTscripts import categoryIchep16, VBF2jTaggedIchep16, VHHadrTaggedIchep16
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

    if hypothesis == "0+":
        p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal = tree.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal
        p_HadWH_SIG_ghw1_1_JHUGen_JECNominal = tree.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal
        p_HadZH_SIG_ghz1_1_JHUGen_JECNominal = tree.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal
    elif hypothesis == "0-":
        p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal = tree.p_JJVBF_SIG_ghv4_1_JHUGen_JECNominal * constants.g4VBF**2
        p_HadWH_SIG_ghw1_1_JHUGen_JECNominal = tree.p_HadWH_SIG_ghw4_1_JHUGen_JECNominal * constants.g4WH**2
        p_HadZH_SIG_ghz1_1_JHUGen_JECNominal = tree.p_HadZH_SIG_ghz4_1_JHUGen_JECNominal * constants.g4ZH**2
    elif hypothesis == "a2":
        p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal = tree.p_JJVBF_SIG_ghv2_1_JHUGen_JECNominal * constants.g2VBF**2
        p_HadWH_SIG_ghw1_1_JHUGen_JECNominal = tree.p_HadWH_SIG_ghw2_1_JHUGen_JECNominal * constants.g2WH**2
        p_HadZH_SIG_ghz1_1_JHUGen_JECNominal = tree.p_HadZH_SIG_ghz2_1_JHUGen_JECNominal * constants.g2ZH**2
    elif hypothesis == "L1":
        p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal = tree.p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal * constants.g1prime2VBF_reco**2 / 1e8
        p_HadWH_SIG_ghw1_1_JHUGen_JECNominal = tree.p_HadWH_SIG_ghw1prime2_1E4_JHUGen_JECNominal * constants.g1prime2WH_reco**2 / 1e8
        p_HadZH_SIG_ghz1_1_JHUGen_JECNominal = tree.p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECNominal * constants.g1prime2ZH_reco**2 / 1e8
    else: assert False

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

def count(desiredcategory, *sample):
    print sample
    sample = Sample(*sample)
    t = tfiles[sample.withdiscriminantsfile()].candTree
    counters = {hypothesis: Counter({_: 0 for _ in range(-1, 2)}) for hypothesis in purehypotheses}
    total = 0
    weightname = sample.weightname()
    countersitems = counters.items()
    length = t.GetEntries()
    for i, entry in enumerate(t, start=1):
        weight = getattr(t, weightname)
        for hypothesis, counter in countersitems:
            if not (t.nExtraLep==0 and (((t.nCleanedJetsPt30==2 or t.nCleanedJetsPt30==3) and t.nCleanedJetsPt30BTagged<=1) or (t.nCleanedJetsPt30>=4 and t.nCleanedJetsPt30BTagged==0))):
                counter[-1] += weight
            else:
                counter[category(t, hypothesis)==desiredcategory] += weight
        total += weight
        if i % 10000 == 0 or i == length:
            print i, "/", length

    for counter in counters.values():
        for k, v in counter.iteritems():
            counter[k] = v/total

    return counters

if __name__ == "__main__":
    counters = count(VBF2jTaggedIchep16, "ggH", "a2", "161221")
    for hypothesis in purehypotheses:
        print hypothesis
        for category in range(-1, 2):
            print "{}    {}%".format(category, counters[hypothesis][category]*100)
        print
