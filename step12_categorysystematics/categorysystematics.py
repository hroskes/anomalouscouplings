#!/usr/bin/env python

from itertools import product
import os

import ROOT

from helperstuff import config, utilities

from helperstuff.enums import analyses, Analysis, BTagSystematic, categories, Category, JECSystematic, pythiasystematics
from helperstuff.samples import ReweightingSample, ReweightingSamplePlus, ReweightingSampleWithFlavor, Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import MultiplyCounter
from helperstuff.yields import count

def findsystematic(categorizations, categorization, JEC, btag):
    name = categorization.category_function_name.replace("_JECUp", "").replace("_JECDn", "").replace("_btagSFUp", "").replace("_btagSFDn", "")
    name += BTagSystematic(btag).appendname + JECSystematic(JEC).appendname
    result = {_ for _ in categorizations if _.category_function_name == name}
    assert len(result) == 1, result
    return result.pop()

def maketable(analysis, production):
    analysis = Analysis(analysis)
    categorizations = [_ for _ in TreeWrapper.categorizations if analysis.categoryname in _.category_function_name]

    BSM = analysis.purehypotheses[1]

    tosamples = [
                     ReweightingSample("ggH", "0+"),
                 ReweightingSamplePlus("ggH", "0+", "MINLO"),
                     ReweightingSample("HJJ", "0+", "Hff0+"),
                     ReweightingSample("VBF", "0+"),
                 ReweightingSamplePlus("VBF", "0+", "POWHEG"),
                     ReweightingSample("ZH", "0+"),
                 ReweightingSamplePlus("ZH", "0+", "POWHEG"),
                     ReweightingSample("WH", "0+"),
                 ReweightingSamplePlus("WplusH", "0+", "POWHEG"),
                 ReweightingSamplePlus("WminusH", "0+", "POWHEG"),
                     ReweightingSample("ttH", "0+", "Hff0+"),
                 ReweightingSamplePlus("ttH", "0+", "Hff0+", "POWHEG"),
                     ReweightingSample("ggH", BSM),
                     ReweightingSample("VBF", BSM),
                     ReweightingSample("ZH", BSM),
                     ReweightingSample("WH", BSM),
                     ReweightingSample("ttH", BSM, "Hff0+"),
                     ReweightingSample("qqZZ"),
                     ReweightingSampleWithFlavor("ggZZ", "2e2mu"),
                     ReweightingSampleWithFlavor("VBF bkg", "2e2mu"),
                     ReweightingSample("ZX"),
                ]

    othertosamples = [
                      ReweightingSamplePlus(_, systematic) for _ in tosamples for systematic in pythiasystematics
                                                           if isinstance(_, ReweightingSamplePlus)
                     ] + [
                      ReweightingSamplePlus(_, systematic, "POWHEG") for _ in tosamples for systematic in pythiasystematics
                                                           if _.productionmode == "ggH" and _.hypothesis == "0+" and not isinstance(_, ReweightingSamplePlus)
                     ]

    result = MultiplyCounter()

    for tosample in tosamples+othertosamples:
        print tosample
        if tosample in result: continue

        ####################################
        if isinstance(tosample, ReweightingSamplePlus) and tosample.alternategenerator == "MINLO" and tosample.pythiasystematic is None:
            try:
                Sample(tosample, production)
            except ValueError:
                pass
            else:
                assert False, "Delete this section!"
            continue
        ####################################

        if isinstance(tosample, ReweightingSamplePlus):
            fromsamples = [Sample(tosample, production)]
            usetosamples = [tosample]
        else:
            fromsamples = tosample.productionmode.allsamples(production)
            usetosamples = [_ for _ in tosamples if _.productionmode == tosample.productionmode and not isinstance(_, ReweightingSamplePlus)]

        result += count(fromsamples, usetosamples, categorizations)

    for key in result:
        if isinstance(key, ReweightingSample): continue
        elif isinstance(key, tuple):
            try:
                result[key] /= result[key[0]]
            except ZeroDivisionError:
                pass
        else: assert False

    rowfmt = "| {:20} | " + " | ".join("{:12.1%} {:+10.1%} {:+10.1%} {:+10.1%} {:+10.1%} {:+10.1%} {:+10.1%} {:+10.1%} {:+10.1%}" for category in categories) + " |"
    headerfmt = rowfmt.replace(".1%", "").replace("+", "").replace(":10", ":^10").replace(":12", ":^12")
    line = "-"*len(headerfmt.format(*[""]*headerfmt.count("{")))

    categorization = {_ for _ in categorizations if _.category_function_name == "category_"+analysis.categoryname}
    assert len(categorization) == 1, categorization
    categorization = categorization.pop()

    tablerows = []

    tablerows.append(line)
    tablerows.append(headerfmt.format(analysis, *sum(([category, "JEC up", "JEC down", "btag SF up", "btag SF dn", "scale up", "scale dn", "tune up", "tune dn"] for category in categories), [])))

    lasthypothesis = None

    for tosample in tosamples:
        if lasthypothesis != tosample.hypothesis:
            tablerows.append(line)
        lasthypothesis = tosample.hypothesis
        ####################################
        if isinstance(tosample, ReweightingSamplePlus) and tosample.alternategenerator == "MINLO" and tosample.pythiasystematic is None:
            try:
                Sample(tosample, production)
            except ValueError:
                pass
            else:
                assert False, "Delete this section!"
            for ctgrztn, category in product(categorizations, categories):
                result[tosample,ctgrztn,category] = .25 * sum(
                    result[ReweightingSamplePlus(tosample, _), ctgrztn, category]
                       for _ in pythiasystematics
                )
        ####################################

        if tosample == ReweightingSample("WH", "L1Zg"):
            for ctgrztn, category in product(categorizations, categories):
                result[tosample, ctgrztn, category] = float("nan")

        thisrowfmt = rowfmt
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
            btSFUp  = (btSFUp - nominal)/nominal if nominal != 0 and tosample.productionmode != "ZX" else float("nan")
            btSFDn  = (btSFDn - nominal)/nominal if nominal != 0 and tosample.productionmode != "ZX" else float("nan")

            if isinstance(tosample, ReweightingSamplePlus) or tosample.productionmode == "ggH" and tosample.hypothesis == "0+":
                ScaleUp = result[ReweightingSamplePlus(tosample, "ScaleUp"), categorization, category]
                ScaleDn = result[ReweightingSamplePlus(tosample, "ScaleDn"), categorization, category]
                TuneUp = result[ReweightingSamplePlus(tosample, "TuneUp"), categorization, category]
                TuneDn = result[ReweightingSamplePlus(tosample, "TuneDn"), categorization, category]
                ScaleUp   = (ScaleUp - nominal)/nominal if nominal != 0 else float("nan")
                ScaleDn   = (ScaleDn - nominal)/nominal if nominal != 0 else float("nan")
                TuneUp  = (TuneUp - nominal)/nominal if nominal != 0 else float("nan")
                TuneDn  = (TuneDn - nominal)/nominal if nominal != 0 else float("nan")
            else:
                ScaleUp = ScaleDn = TuneUp = TuneDn = float("nan")

            fmtargs += [nominal, JECUp, JECDn, btSFUp, btSFDn, ScaleUp, ScaleDn, TuneUp, TuneDn]
        tablerows.append(thisrowfmt.format(*fmtargs)
                                                    .replace("+nan%",  "     ")
                                                    .replace("-nan%",  "     ")
                                                    .replace( "nan%",   "    ")
                                                                               )

    tablerows.append(line)

    table = "\n".join(tablerows)

    print table

    outdir = os.path.join(config.plotsbasedir, "categorization", str(analysis))
    utilities.mkdir_p(outdir)
    with open(os.path.join(outdir, "table.txt"), "w") as f:
        f.write(table)

if __name__ == "__main__":
  for analysis in analyses:
    for production in config.productionsforcombine:
      maketable(analysis, production)
