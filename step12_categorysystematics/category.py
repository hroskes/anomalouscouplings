#!/usr/bin/env python
from collections import Counter
from helperstuff import config
from helperstuff import constants
from helperstuff.CJLSTscripts import categoryIchep16, VBF2jTaggedIchep16, VHHadrTaggedIchep16
from helperstuff.enums import Hypothesis, purehypotheses
from helperstuff.samples import ReweightingSample, Sample
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

    if hypothesis == "0+": assert result == tree.category_0P
    if hypothesis == "0-": assert result == tree.category_0M
    if hypothesis == "a2": assert result == tree.category_a2, "{} {} {}".format(result, tree.category_a2, tree.Show())
    if hypothesis == "L1": assert result == tree.category_L1

    return result

def count(desiredcategory, sampleforfile, *tohypotheses):
    tohypotheses = [Hypothesis(_) for _ in tohypotheses]
    assert len(tohypotheses) == len(set(tohypotheses))
    productionmode = sampleforfile.productionmode
    samples = [ReweightingSample(productionmode, hypothesis) for hypothesis in tohypotheses]
    t = tfiles[sampleforfile.withdiscriminantsfile()].candTree
    SM = Hypothesis("0+")

    categoryhypotheses = sorted(purehypotheses, key=lambda x: x != "0+")

    counters = {(usehypothesis, hypothesis): Counter({_: 0 for _ in range(-1, 2)}) for hypothesis in categoryhypotheses for usehypothesis in tohypotheses}
    counters_or = {(usehypothesis, (SM, hypothesis)): Counter({_: 0 for _ in range(-1, 2)}) for hypothesis in categoryhypotheses for usehypothesis in tohypotheses if hypothesis != SM}
    total = Counter()
    weightnames = {hypothesis: sample.weightname() for hypothesis, sample in zip(tohypotheses, samples)}.items()
    length = t.GetEntries()
    tmpresult = {}
    for i, entry in enumerate(t, start=1):
        for hypothesis in categoryhypotheses:
            if not (t.nExtraLep==0 and (((t.nCleanedJetsPt30==2 or t.nCleanedJetsPt30==3) and t.nCleanedJetsPt30BTagged<=1) or (t.nCleanedJetsPt30>=4 and t.nCleanedJetsPt30BTagged==0))):
                tmpresult[hypothesis] = -1
            else:
                tmpresult[hypothesis] = category(t, hypothesis)==desiredcategory

        for usehypothesis, weightname in weightnames:
            weight = getattr(t, weightname)
            total[usehypothesis] += weight

            for hypothesis in categoryhypotheses:
                counters[usehypothesis, hypothesis][tmpresult[hypothesis]] += weight
                if hypothesis != SM:
                    counters_or[usehypothesis, (SM, hypothesis)][max(tmpresult[SM], tmpresult[hypothesis])] += weight

        if i % 10000 == 0 or i == length:
            print i, "/", length
            #break

    counters.update(counters_or)
    for usehypothesis in tohypotheses:
        for hypothesis in categoryhypotheses:
            counter = counters[usehypothesis, hypothesis]
            for k, v in counter.items():
                counter[k] = v/total[usehypothesis]
            if hypothesis != SM:
                counter = counters[usehypothesis, (SM, hypothesis)]
                for k, v in counter.items():
                    counter[k] = v/total[usehypothesis]

    return counters

if __name__ == "__main__":
    SM = Hypothesis("0+")

    tohypotheses = [Hypothesis(_) for _ in ["0+", "0-", "a2", "L1", "fa2prod-0.5", "fL1prod0.5"]]

    counters = None

    fromhypotheses = "0+", #"0-", "a2", "L1", "fa3prod0.5", "fa2prod0.5", "fL1prod0.5"

    for fromhypothesis in fromhypotheses:
        theupdate = count(VBF2jTaggedIchep16, Sample("VBF", fromhypothesis, config.productionsforcombine[0]), *tohypotheses)
        for counter in theupdate.values():
            for k in counter.keys():
                counter[k] /= len(fromhypotheses)

        if counters is None:
            counters = theupdate
        else:
            assert counters.keys() == theupdate.keys()
            for key in counters.keys():
                counters[key] += theupdate[key]

    fmt = "{:11} {:5.1%}"
    print
    print "Not enough jets:"
    for tohypothesis in tohypotheses:
        print fmt.format(tohypothesis, counters[tohypothesis, SM][-1])
    print

    fmt = " ".join(["{:>11}"]*8)
    print fmt.format(*[""] + [_ for _ in purehypotheses] + ["0+ or {}".format(_) for _ in purehypotheses if _ != "0+"])

    fmt = " ".join(["{:11}"] + ["{:11.1%}"]*7)
    for tohypothesis in tohypotheses:
        print fmt.format(*[tohypothesis]
                        + [counters[tohypothesis, hypothesis][1] for hypothesis in purehypotheses]
                        + [counters[tohypothesis, (SM, hypothesis)][1] for hypothesis in purehypotheses if hypothesis != SM]
                        )
