#!/usr/bin/env python

from itertools import product

from helperstuff import config

from helperstuff.enums import categories, Category, JECSystematic, ProductionMode
from helperstuff.samples import ReweightingSample, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import MultiplyCounter, tfiles

def count(fromsample, tosamples, categorizations):
    print fromsample
    f = tfiles[fromsample.withdiscriminantsfile()]
    t = f.candTree
    length = t.GetEntries()
    result = MultiplyCounter()
    t.SetBranchStatus("*", 0)
    t.SetBranchStatus("category_*", 1)
    t.SetBranchStatus("MC_weight_*", 1)

    theproduct = list(product(tosamples, categorizations))

    weight, category = {}, {}
    for i, entry in enumerate(t, start=1):
        for tosample in tosamples:
            weight[tosample] = getattr(t, tosample.weightname())
            result[tosample] += weight[tosample]
        for categorization in categorizations:
            category[categorization] = getattr(t, categorization.category_function_name)

        for tosample, categorization in theproduct:
            result[tosample, categorization, category[categorization]] += weight[tosample]
        if i % 10000 == 0 or i == length:
            print i, "/", length
            #break
    return result

def findsystematic(categorizations, categorization, JEC):
    name = categorization.category_function_name.replace("_JECUp", "").replace("_JECDn", "")
    name += JECSystematic(JEC).appendname
    result = {_ for _ in categorizations if _.category_function_name == name}
    assert len(result) == 1, result
    return result.pop()

def maketable(productionmode):
    assert len(config.productionsforcombine) == 1
    productionmode = ProductionMode(productionmode)
    fromsamples = [Sample(config.productionsforcombine[0], productionmode, _) for _ in productionmode.generatedhypotheses]
    tosamples = [ReweightingSample(productionmode, _) for _ in productionmode.validhypotheses]
    categorizations = [_ for _ in TreeWrapper.categorizations if "category_0P_or_0M" in _.category_function_name]

    result = sum((count(fromsample, tosamples, categorizations) for fromsample in fromsamples), MultiplyCounter()) / len(fromsamples)

    for key in result:
        if isinstance(key, ReweightingSample): continue
        elif isinstance(key, tuple):
            try:
                result[key] /= result[key[0]]
            except ZeroDivisionError:
                pass
        else: assert False

    rowfmt = "| {:15} | " + " | ".join("{:12.1%} {:+10.1%} {:+10.1%}" for category in categories) + " |"
    headerfmt = rowfmt.replace(".1%", "").replace("+", "").replace(":10", ":^10").replace(":12", ":^12")
    line = "-"*len(headerfmt.format(*[""]*headerfmt.count("{")))

    tosample = ReweightingSample(productionmode, "SM")
    categorization = {_ for _ in categorizations if _.category_function_name == "category_0P_or_0M"}
    assert len(categorization) == 1, categorization
    categorization = categorization.pop()

    print line
    print headerfmt.format(productionmode, *sum(([category, "JEC up", "JEC down"] for category in categories), []))
    print line

    for tosample in tosamples:
        fmtargs = [tosample.hypothesis]
    #or replace last 2 lines with
    #for categorization in categorizations:
    #   if categorization.JEC != "JECNominal": continue
    #   fmtargs = [categorization]
        for category in categories:
            nominal = sum(result[tosample, categorization, _] for _ in category.idnumbers)
            JECUp   = sum(result[tosample, findsystematic(categorizations, categorization, "JECUp"), _] for _ in category.idnumbers)
            JECDn   = sum(result[tosample, findsystematic(categorizations, categorization, "JECDn"), _] for _ in category.idnumbers)
            JECUp = (JECUp - nominal)/nominal if nominal != 0 else float("nan")
            JECDn = (JECDn - nominal)/nominal if nominal != 0 else float("nan")
            fmtargs += [nominal, JECUp, JECDn]
        print rowfmt.format(*fmtargs)

    print line

if __name__ == "__main__":
    maketable("ZH")
