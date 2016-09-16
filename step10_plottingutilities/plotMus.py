from ROOT import *

from helperstuff import config, style
from collections import namedtuple
import os

gStyle.SetPadTickY(0)
gStyle.SetTickLength(0.0, "Y");
gStyle.SetEndErrorSize(10)

from array import array
from limits import printlimits

a_mus_up = {}; a_mus_dn = {}; v_mus_up = {}; v_mus_dn = {}; g_mus = {};
a_y = {}; v_y = {};

a_zero = array('d',[0])
v_zero = TVectorD(len(a_zero),a_zero)

Lumi = namedtuple("Lumi", "lumi foldername fa2plotname fa3plotname fL1plotname")
lumis = [
         Lumi("300fb",      "September12_forECFA_300",
                            "limit_100,-0.05,0.05", "limit_100,-0.2,0.2", "limit_50,-0.2,0.2"),
         Lumi("300fb_sc2",  "September12_forECFA_300_scenario2",
                            "limit_100,-0.05,0.05", "limit_100,-0.2,0.2", "limit_50,-0.2,0.2"),
         Lumi("3000fb",     "September12_forECFA",
                            "limit_80,-0.004,0.004", "limit_60,-0.12,0.12", "limit_100,-0.02,0.02"),
         Lumi("3000fb_sc2", "September12_forECFA_scenario2",
                            "limit_80,-0.004,0.004", "limit_60,-0.12,0.12", "limit_100,-0.02,0.02"),
        ]

for lumi, foldername, fa2plotname, fa3plotname, fL1plotname in lumis:
  for i, (mu, plotname) in enumerate(zip(['fL1','fa2','fa3'], (fL1plotname, fa2plotname, fa3plotname)), start=1):

    print lumi, mu
    minimum, results = printlimits(mu, foldername, plotname=plotname)
    assert len(results[3.84]) == 1

    a_y[lumi+'_'+mu] = array('d',[i])
    v_y[lumi+'_'+mu] = TVectorD(len(a_y[lumi+'_'+mu]),a_y[lumi+'_'+mu])

    a_mus_up[lumi+'_'+mu] = array('d',[ results[3.84][0][1]])
    a_mus_dn[lumi+'_'+mu] = array('d',[-results[3.84][0][0]])
    v_mus_up[lumi+'_'+mu] = TVectorD(len(a_mus_up[lumi+'_'+mu]),a_mus_up[lumi+'_'+mu])
    v_mus_dn[lumi+'_'+mu] = TVectorD(len(a_mus_dn[lumi+'_'+mu]),a_mus_dn[lumi+'_'+mu])

    g_mus[lumi+'_'+mu] = TGraphAsymmErrors(v_zero,v_y[lumi+'_'+mu],v_mus_dn[lumi+'_'+mu],v_mus_up[lumi+'_'+mu],v_zero,v_zero)

    g_mus[lumi+'_'+mu].SetLineWidth(3)

    if ('sc2' in lumi): g_mus[lumi+'_'+mu].SetLineColor(kRed)
    else: g_mus[lumi+'_'+mu].SetLineColor(kGreen+2)

for lumi in ['300fb','3000fb']:

  c1 = TCanvas("c1","c1",1000,800)
  c1.SetRightMargin(0.05)
  c1.SetLeftMargin(0.05)

  if (lumi=='300fb'): dummy = TH1D("dummy","dummy",1,-0.2,0.2)
  if (lumi=='3000fb'): dummy = TH1D("dummy","dummy",1,-0.02,0.02)

  dummy.SetMinimum(0.0)
  dummy.SetMaximum(7.0)
  dummy.SetLineColor(0)
  dummy.SetMarkerColor(0)
  dummy.SetLineWidth(0)
  dummy.SetMarkerSize(0)
  dummy.GetYaxis().SetTitle("")
  dummy.GetYaxis().SetLabelSize(0)
  dummy.GetXaxis().SetTitle("expected uncertainty")
  dummy.GetXaxis().SetLabelSize(0.04)
  dummy.SetNdivisions(510, "X");
  dummy.Draw()

  latex2 = TLatex()
  latex2.SetNDC()
  latex2.SetTextSize(0.8*c1.GetTopMargin())
  latex2.SetTextFont(42)
  latex2.SetTextAlign(11) # align right
  latex2.DrawLatex(0.05, 0.95, "CMS Projection")

  latex2.SetTextSize(0.55*c1.GetTopMargin())
  latex2.DrawLatex(0.08,0.85, "Expected uncertainties on")
  latex2.DrawLatex(0.08,0.81, "Higgs Boson anomalous couplings")

  latex2.SetTextSize(0.55*c1.GetTopMargin())
  latex2.DrawLatex(0.65,0.66, "H #rightarrow ZZ* #rightarrow 4#font[12]{l}")
  latex2.DrawLatex(0.65,0.60, "m_{H} = 125 GeV")

  latex2.SetTextSize(0.55*c1.GetTopMargin())

  for mu in ["fa3", "fa2", "fL1"]:
    g_mus[lumi+'_'+mu].Draw("|same")
    g_mus[lumi+'_sc2_'+mu].Draw("|same")
    latex2.DrawLatex(0.08,a_y[lumi+'_'+mu][0]/8.7+0.13, mu.replace("f", "f_{").replace("L", "#Lambda")+"}^{ZZ}")

  line =  TLine(0.0,0.0, 0.0, 5.5)
  line.SetLineWidth(2)
  line.SetLineColor(1)
  line.Draw("same")

  legend = TLegend(.55,.75,.94,.90)
  legend.SetBorderSize(0)
  legend.AddEntry(g_mus[lumi+'_'+mu], lumi.replace("fb","")+" fb^{-1} at #sqrt{s}=13 TeV Scenario 1", "l")
  legend.AddEntry(g_mus[lumi+'_sc2_'+mu], lumi.replace("fb","")+" fb^{-1} at #sqrt{s}=13 TeV Scenario 2", "l")
  legend.Draw("same")

  saveasdir = os.path.join(config.plotsbasedir, "ECFAsummary")
  try: os.makedirs(saveasdir)
  except OSError: pass

  for ext in "png eps root pdf".split():
    c1.SaveAs(os.path.join(saveasdir, "faisummary_{}.{}".format(lumi, ext)))




c1 = TCanvas("c1","c1",1000,800)
c1.SetRightMargin(0.05)
c1.SetLeftMargin(0.05)

dummy = TH1D("dummy","dummy",1,-0.2,0.2)

dummy.SetMinimum(0.0)
dummy.SetMaximum(7.0)
dummy.SetLineColor(0)
dummy.SetMarkerColor(0)
dummy.SetLineWidth(0)
dummy.SetMarkerSize(0)
dummy.GetYaxis().SetTitle("")
dummy.GetYaxis().SetLabelSize(0)
dummy.GetXaxis().SetTitle("expected uncertainty")
dummy.GetXaxis().SetLabelSize(0.04)
dummy.SetNdivisions(510, "X");
dummy.Draw()

latex2 = TLatex()
latex2.SetNDC()
latex2.SetTextSize(0.8*c1.GetTopMargin())
latex2.SetTextFont(42)
latex2.SetTextAlign(11) # align right
latex2.DrawLatex(0.05, 0.95, "CMS Projection")

latex2.SetTextSize(0.55*c1.GetTopMargin())
latex2.DrawLatex(0.08,0.85, "Expected uncertainties on")
latex2.DrawLatex(0.08,0.81, "Higgs Boson anomalous couplings")

latex2.SetTextSize(0.55*c1.GetTopMargin())
latex2.DrawLatex(0.65,0.66, "H #rightarrow ZZ* #rightarrow 4#font[12]{l}")
latex2.DrawLatex(0.65,0.60, "m_{H} = 125 GeV")

latex2.SetTextSize(0.55*c1.GetTopMargin())

for mu in ["fa3", "fa2", "fL1"]:
  g_mus['300fb_'+mu].SetLineColor(4)
  g_mus['300fb_sc2_'+mu].SetLineColor(6)
  for lumi in "300fb", "3000fb":
    g_mus[lumi+'_'+mu].Draw("|same")
    #g_mus[lumi+'_sc2_'+mu].Draw("|same")

  latex2.DrawLatex(0.08,a_y[lumi+'_'+mu][0]/8.7+0.13, mu.replace("f", "f_{").replace("L", "#Lambda")+"}^{ZZ}")

line =TLine(0.0,0.0, 0.0, 5.5)
line.SetLineWidth(2)
line.SetLineColor(1)
line.Draw("same")

legend = TLegend(.55,.75,.94,.90)
legend.SetBorderSize(0)
for lumi in "300fb", "3000fb":
  legend.AddEntry(g_mus[lumi+'_'+mu], lumi.replace("fb","")+" fb^{-1} at #sqrt{s}=13 TeV Scenario 1", "l")
  #legend.AddEntry(g_mus[lumi+'_sc2_'+mu], lumi.replace("fb","")+" fb^{-1} at #sqrt{s}=13 TeV Scenario 2", "l")
legend.Draw("same")

saveasdir = os.path.join(config.plotsbasedir, "ECFAsummary")
try: os.makedirs(saveasdir)
except OSError: pass

for ext in "png eps root pdf".split():
  c1.SaveAs(os.path.join(saveasdir, "faisummary_both.{}".format(ext)))


