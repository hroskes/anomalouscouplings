#!/usr/bin/env python

import operator, glob, re

def prod(iterable):
  return reduce(operator.mul, iterable, 1)

def genxsec(filename):
  import ROOT
  f = ROOT.TFile(filename)
  t = f.candTree
  t.GetEntry(0)
  return t.genxsec * t.genBR

def printgenxsecs(productionmode):
  filenames = [_ for _ in sorted(glob.iglob(productionmode+"*.root")) if "POWHEG" not in _ and "MINLO" not in _]
  maxlen = max(len(_) for _ in filenames)
  fmt = "{:%d} {:.4e}"%(maxlen+5)

  for _ in filenames: print fmt.format(_, genxsec(_))

def sumofweights(filename, weightformula, cutformula="1"):
  import numpy as np
  import ROOT, rootoverloads
  from uncertainties import ufloat

  f = ROOT.TFile(filename)
  t = f.candTree
  t.SetBranchStatus("*", 0)
  for formula in weightformula, cutformula:
    for branch in re.findall(r"\b[a-zA-Z_][a-zA-Z0-9_]+\b", formula):
      t.SetBranchStatus(branch, 1)


  weight = ROOT.TTreeFormula("weight", weightformula, t)
  cut = ROOT.TTreeFormula("cut", cutformula, t)

  a = []

  for entry in t:
    if cut.EvalInstance():
      wt = weight.EvalInstance()
      a.append(wt)

  a = np.array(a)

  return ufloat(sum(a), sum(a**2)**.5)

def run(productionmode, weightformula, cutformula="1"):
  filenames = [_ for _ in sorted(glob.iglob(productionmode+"*.root*")) if "POWHEG" not in _ and "MINLO" not in _]

  assert len(filenames) == 7 + 2*(productionmode != "WH")

  filenames = [_ for _ in filenames if ".tmp" not in _]
  if "ghv" in weightformula:
    filenames = [_ for _ in filenames if "VBFL1Zg" not in _]

  sums = []

  maxlen = max(len(_) for _ in filenames)
  fmt = "{:%d} {:.4f}"%(maxlen+5)

  for filename in filenames:
    sums.append(sumofweights(filename, weightformula, cutformula))
    print fmt.format(filename, sums[-1])

  from uncertainties import ufloat

  print
  average = ufloat(
    sum(_.nominal_value / _.std_dev**2 for _ in sums) / sum(1 / _.std_dev**2 for _ in sums),
    sum(1 / _.std_dev**2 for _ in sums) ** -0.5,
  )
  print
  print fmt.format("Average:", average)
  print
  for filename, _ in zip(filenames, sums):
    print fmt.format(filename, _ - average)

  print
  print

#run("VBF", "MC_weight_nominal * p_Gen_Dec_SIG_ghz1_1_JHUGen * p_Gen_VBF_SIG_ghv1_1_JHUGen")
#run("ggH", "MC_weight_nominal * "
#           "(p_Gen_GG_SIG_ghg2_1_ghz1prime2_1E4_ghz2_1_JHUGen - p_Gen_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen - p_Gen_GG_SIG_ghg2_1_ghz2_1_JHUGen)"
#)
#run("VBF", "MC_weight_nominal * "
#           "(p_Gen_Dec_SIG_ghz1prime2_1E4_ghz2_1_JHUGen - p_Gen_Dec_SIG_ghz1prime2_1E4_JHUGen - p_Gen_Dec_SIG_ghz2_1_JHUGen) * "
#           "(p_Gen_VBF_SIG_ghv1_1_ghza1prime2_1E4_JHUGen - p_Gen_VBF_SIG_ghv1_1_JHUGen - p_Gen_VBF_SIG_ghza1prime2_1E4_JHUGen)"
#)
#run("WH", "MC_weight_nominal * "
#           "(p_Gen_Dec_SIG_ghz1prime2_1E4_ghz2_1_JHUGen - p_Gen_Dec_SIG_ghz1prime2_1E4_JHUGen - p_Gen_Dec_SIG_ghz2_1_JHUGen) * "
#           "p_Gen_WH_SIG_ghw1prime2_1E4_JHUGen"
#)
#run("ZH", "MC_weight_nominal * "
#           "(p_Gen_Dec_SIG_ghz2_1_ghza1prime2_1E4_JHUGen - p_Gen_Dec_SIG_ghza1prime2_1E4_JHUGen - p_Gen_Dec_SIG_ghz2_1_JHUGen) * "
#           "p_Gen_ZH_SIG_ghza1prime2_1E4_JHUGen"
#           " * (D_4couplings_HadVHdecay == 0)"
#)
#run("ZH", "MC_weight_nominal * "
#           "p_Gen_Dec_SIG_ghza1prime2_1E4_JHUGen * "
#           "p_Gen_ZH_SIG_ghza1prime2_1E4_JHUGen"
#)
#run("WH", "MC_weight_nominal * "
#           "p_Gen_Dec_SIG_ghza1prime2_1E4_JHUGen * "
#           "p_Gen_WH_SIG_ghw1prime2_1E4_ghw2_1_JHUGen"
#)
run("ZH", "MC_weight_nominal * "
           "p_Gen_Dec_SIG_ghz1_1_JHUGen * "
           "p_Gen_ZH_SIG_ghz1prime2_1E4_ghz2_1_JHUGen"
)
#run("VBF", "MC_weight_nominal * "
#            "p_Gen_Dec_SIG_ghza1prime2_1E4_JHUGen * "
#            "p_Gen_VBF_SIG_ghza1prime2_1E4_JHUGen")
#run("VBF", "MC_weight_nominal * "
#           "(p_Gen_Dec_SIG_ghz1prime2_1E4_ghz2_1_JHUGen - p_Gen_Dec_SIG_ghz1prime2_1E4_JHUGen - p_Gen_Dec_SIG_ghz2_1_JHUGen) * "
#           "p_Gen_VBF_SIG_ghza1prime2_1E4_JHUGen"
#)

#printgenxsecs("ZH")
