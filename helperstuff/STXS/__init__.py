import os

import utilities

thisfolder = os.path.dirname(os.path.abspath(__file__))
utilities.LoadMacro(os.path.join(thisfolder, "stage1.cc+"))

import ROOT

from ..CJLSTscripts import DVBF1j_ME as __DVBF1j_ME

def stage1_reco_sync(Njets, pTj1, mjj, deta_jj, H_pt, category, pt_hjj, p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass):
  D1jet = __DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass)
  return ROOT.stage1_reco_sync(Njets, pTj1, mjj, deta_jj, H_pt, category, D1jet, pt_hjj)
