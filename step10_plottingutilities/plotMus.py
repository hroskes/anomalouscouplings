from ROOT import *

from helperstuff import style

gStyle.SetPadTickY(0)
gStyle.SetTickLength(0.0, "Y");
gStyle.SetEndErrorSize(10)

from array import array

from results import results

a_mus = {}; v_mus = {}; g_mus = {};
a_y = {}; v_y = {};

a_zero = array('d',[0])
v_zero = TVectorD(len(a_zero),a_zero)

for lumi in ['300fb','300fb_sc2','3000fb','3000fb_sc2']:
  i=5
  for mu in ['r_tot','r_ggH','r_VBF','r_VH','r_ttH']:

    a_y[lumi+'_'+mu] = array('d',[i])
    v_y[lumi+'_'+mu] = TVectorD(len(a_y[lumi+'_'+mu]),a_y[lumi+'_'+mu])

    a_mus[lumi+'_'+mu] = array('d',[(results[lumi+'_'+mu]['up']/results[lumi+'_'+mu]['central']+
                                    results[lumi+'_'+mu]['dn']/results[lumi+'_'+mu]['central'])/2] )
    v_mus[lumi+'_'+mu] = TVectorD(len(a_mus[lumi+'_'+mu]),a_mus[lumi+'_'+mu])

    g_mus[lumi+'_'+mu] = TGraphAsymmErrors(v_zero,v_y[lumi+'_'+mu],v_zero,v_mus[lumi+'_'+mu],v_zero,v_zero)

    g_mus[lumi+'_'+mu].SetLineWidth(3)

    if ('sc2' in lumi): g_mus[lumi+'_'+mu].SetLineColor(kRed)
    else: g_mus[lumi+'_'+mu].SetLineColor(kGreen+2)

    i-=1

for lumi in ['300fb','3000fb']:

  c1 = TCanvas("c1","c1",1000,800)
  c1.SetRightMargin(0.05)
  c1.SetLeftMargin(0.05)

  if (lumi=='300fb'): dummy = TH1D("dummy","dummy",1,-0.159,0.879)
  if (lumi=='3000fb'): dummy = TH1D("dummy","dummy",1,-0.049,0.299)

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
  latex2.DrawLatex(0.08,0.81, "Higgs Boson signal strengths")
  
  latex2.SetTextSize(0.55*c1.GetTopMargin())
  latex2.DrawLatex(0.65,0.66, "H #rightarrow ZZ* #rightarrow 4#font[12]{l}")
  latex2.DrawLatex(0.65,0.60, "m_{H} = 125.09 GeV")
  
  latex2.SetTextSize(0.55*c1.GetTopMargin())
  
  for mu in ['r_tot','r_ggH','r_VBF','r_VH','r_ttH']:
    g_mus[lumi+'_'+mu].Draw("|same")
    g_mus[lumi+'_sc2_'+mu].Draw("|same")

    latex2.DrawLatex(0.08,a_y[lumi+'_'+mu][0]/8.7+0.13, "#mu_{"+mu.replace("r_","")+"}^{ZZ}")

  line =  TLine(0.0,0.0, 0.0, 5.5)
  line.SetLineWidth(2)
  line.SetLineColor(1)
  line.Draw("same")

  legend = TLegend(.55,.75,.94,.90)
  legend.SetBorderSize(0)
  legend.AddEntry(g_mus[lumi+'_'+mu], lumi.replace("fb","")+" fb^{-1} at #sqrt{s}=13 TeV Scenario 1", "l")
  legend.AddEntry(g_mus[lumi+'_sc2_'+mu], lumi.replace("fb","")+" fb^{-1} at #sqrt{s}=13 TeV Scenario 2", "l")
  legend.Draw("same")
  
  c1.SaveAs("muPerProc_"+lumi+".pdf")
  

