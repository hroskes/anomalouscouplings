#!/usr/bin/env python

from itertools import product
import os

import ROOT

from helperstuff import config, utilities

from helperstuff.enums import analyses, Analysis, BTagSystematic, categories, Category, JECSystematic, ProductionMode
from helperstuff.samples import ReweightingSample, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import MultiplyCounter

def count(fromsamples, tosamples, categorizations):
    t = ROOT.TChain("candTree")
    for fromsample in fromsamples:
        t.Add(fromsample.withdiscriminantsfile())
    length = t.GetEntries()
    result = MultiplyCounter()
    t.SetBranchStatus("*", 0)
    t.SetBranchStatus("category_*", 1)
    t.SetBranchStatus("MC_weight_*", 1)

    c = ROOT.TCanvas()
    for tosample, categorization in product(tosamples, categorizations):
        t.Draw(categorization.category_function_name, tosample.weightname())
        h = c.GetListOfPrimitives()[0]
        for i in range(6):
            result[tosample, categorization, Category.fromid(i)] += h.GetBinContent(i+1)

        if tosample not in result:
            result[tosample] = sum(h.GetBinContent(i+1) for i in range(6))

    return result

def findsystematic(categorizations, categorization, JEC, btag):
    name = categorization.category_function_name.replace("_JECUp", "").replace("_JECDn", "").replace("_btagSFUp", "").replace("_btagSFDn", "")
    name += BTagSystematic(btag).appendname + JECSystematic(JEC).appendname
    result = {_ for _ in categorizations if _.category_function_name == name}
    assert len(result) == 1, result
    return result.pop()

def maketable(analysis):
    analysis = Analysis(analysis)
    assert len(config.productionsforcombine) == 1
    production = config.productionsforcombine[0]
    categorizations = [_ for _ in TreeWrapper.categorizations if analysis.categoryname in _.category_function_name]

    BSM = analysis.purehypotheses[1]

    tosamples = [
                 ReweightingSample("ggH", "0+"),
                 ReweightingSample("VBF", "0+"),
                 ReweightingSample("ZH", "0+"),
                 ReweightingSample("WH", "0+"),
                 ReweightingSample("ttH", "0+", "Hff0+"),
                 ReweightingSample("ggH", BSM),
                 ReweightingSample("VBF", BSM),
                 ReweightingSample("ZH", BSM),
                 ReweightingSample("WH", BSM),
                 ReweightingSample("ttH", BSM, "Hff0+"),
                 ReweightingSample("qqZZ"),
                 ReweightingSample("ggZZ", "2e2mu"),
                 ReweightingSample("VBF bkg", "2e2mu"),
                 ReweightingSample("ZX"),
                ]
    if ReweightingSample("WH", "L1Zg") in tosamples: tosamples.remove(ReweightingSample("WH", "L1Zg"))

    result = MultiplyCounter()

    for tosample in tosamples:
        print tosample
        if tosample in result: continue
        result += count(tosample.productionmode.allsamples(production), [_ for _ in tosamples if _.productionmode==tosample.productionmode], categorizations)

    for key in result:
        if isinstance(key, ReweightingSample): continue
        elif isinstance(key, tuple):
            try:
                result[key] /= result[key[0]]
            except ZeroDivisionError:
                pass
        else: assert False

    rowfmt = "| {:15} | " + " | ".join("{:12.1%} {:+10.1%} {:+10.1%} {:+10.1%} {:+10.1%}" for category in categories) + " |"
    headerfmt = rowfmt.replace(".1%", "").replace("+", "").replace(":10", ":^10").replace(":12", ":^12")
    line = "-"*len(headerfmt.format(*[""]*headerfmt.count("{")))

    categorization = {_ for _ in categorizations if _.category_function_name == "category_"+analysis.categoryname}
    assert len(categorization) == 1, categorization
    categorization = categorization.pop()

    tablerows = []

    tablerows.append(line)
    tablerows.append(headerfmt.format(analysis, *sum(([category, "JEC up", "JEC down", "btag SF up", "btag SF dn"] for category in categories), [])))
    tablerows.append(line)

    for tosample in tosamples:
        fmtargs = [str(tosample).replace(" 2e2mu", "")]
    #or replace last 2 lines with
    #for categorization in categorizations:
    #   if categorization.JEC != "JECNominal": continue
    #   fmtargs = [categorization]
        for category in categories:
            nominal = result[tosample, categorization, category]
            JECUp   = result[tosample, findsystematic(categorizations, categorization, "JECUp", "Nominal"), category]
            JECDn   = result[tosample, findsystematic(categorizations, categorization, "JECDn", "Nominal"), category]
            btSFUp  = result[tosample, findsystematic(categorizations, categorization, "Nominal", "bTagSFUp"), category]
            btSFDn  = result[tosample, findsystematic(categorizations, categorization, "Nominal", "bTagSFDn"), category]
            JECUp   = (JECUp - nominal)/nominal if nominal != 0 else float("nan")
            JECDn   = (JECDn - nominal)/nominal if nominal != 0 else float("nan")
            btSFUp  = (btSFUp - nominal)/nominal if nominal != 0 else float("nan")
            btSFDn  = (btSFDn - nominal)/nominal if nominal != 0 else float("nan")
            fmtargs += [nominal, JECUp, JECDn, btSFUp, btSFDn]
        tablerows.append(rowfmt.format(*fmtargs))

    tablerows.append(line)

    table = "\n".join(tablerows)

    print table

    outdir = os.path.join(config.plotsbasedir, "categorization", str(analysis))
    utilities.mkdir_p(outdir)
    with open(os.path.join(outdir, "table.txt"), "w") as f:
        f.write(table)

if __name__ == "__main__":
    for analysis in analyses:
        maketable(analysis)
