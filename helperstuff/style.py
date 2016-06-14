import ROOT

from ROOT import gStyle, gPad, gROOT, kWhite

def fixOverlay():
  gPad.RedrawAxis()


def setTDRStyle(force):

  gStyle.SetCanvasBorderMode(0)
  gStyle.SetCanvasColor(kWhite)
  gStyle.SetCanvasDefH(600)
  gStyle.SetCanvasDefW(600)
  gStyle.SetCanvasDefX(0)
  gStyle.SetCanvasDefY(0)

  gStyle.SetPadBorderMode(0)
  gStyle.SetPadColor(kWhite)
  gStyle.SetPadGridX(False)
  gStyle.SetPadGridY(False)
  gStyle.SetGridColor(0)
  gStyle.SetGridStyle(3)
  gStyle.SetGridWidth(1)

  gStyle.SetFrameBorderMode(0)
  gStyle.SetFrameBorderSize(1)
  gStyle.SetFrameFillColor(0)
  gStyle.SetFrameFillStyle(0)
  gStyle.SetFrameLineColor(1)
  gStyle.SetFrameLineStyle(1)
  gStyle.SetFrameLineWidth(1)

  if force:
      gStyle.SetHistLineColor(1)
      gStyle.SetHistLineStyle(0)
      gStyle.SetHistLineWidth(1)


  gStyle.SetEndErrorSize(2)
  gStyle.SetErrorX(0.)

  gStyle.SetMarkerStyle(20)

  gStyle.SetOptFit(1)
  gStyle.SetFitFormat("5.4g")
  gStyle.SetFuncColor(2)
  gStyle.SetFuncStyle(1)
  gStyle.SetFuncWidth(1)

  gStyle.SetOptDate(0)

  gStyle.SetOptFile(0)
  gStyle.SetOptStat(0)
  gStyle.SetStatColor(kWhite)
  gStyle.SetStatFont(42)
  gStyle.SetStatFontSize(0.04)
  gStyle.SetStatTextColor(1)
  gStyle.SetStatFormat("6.4g")
  gStyle.SetStatBorderSize(1)
  gStyle.SetStatH(0.1)
  gStyle.SetStatW(0.2)


  gStyle.SetPadTopMargin(0.05)
  gStyle.SetPadBottomMargin(0.13)
  gStyle.SetPadLeftMargin(0.16)
  gStyle.SetPadRightMargin(0.04)


  gStyle.SetOptTitle(0)
  gStyle.SetTitleFont(42)
  gStyle.SetTitleColor(1)
  gStyle.SetTitleTextColor(1)
  gStyle.SetTitleFillColor(10)
  gStyle.SetTitleFontSize(0.05)


  gStyle.SetTitleColor(1, "XYZ")
  gStyle.SetTitleFont(42, "XYZ")
  gStyle.SetTitleSize(0.06, "XYZ")
  gStyle.SetTitleXOffset(0.9)
  gStyle.SetTitleYOffset(1.25)


  gStyle.SetLabelColor(1, "XYZ")
  gStyle.SetLabelFont(42, "XYZ")
  gStyle.SetLabelOffset(0.007, "XYZ")
  gStyle.SetLabelSize(0.05, "XYZ")


  gStyle.SetAxisColor(1, "XYZ")
  gStyle.SetStripDecimals(True)
  gStyle.SetTickLength(0.03, "XYZ")
  gStyle.SetNdivisions(510, "XYZ")
  gStyle.SetPadTickX(1)
  gStyle.SetPadTickY(1)

  gStyle.SetOptLogx(0)
  gStyle.SetOptLogy(0)
  gStyle.SetOptLogz(0)

  gStyle.SetPaperSize(20.,20.)


  gROOT.ForceStyle()



def tdrstyle(force=True):
    setTDRStyle(force)


tdrstyle()









def drawlines():
    line68 = ROOT.TLine()
    line68.SetLineStyle(9)
    line68.DrawLine(-1,1,1,1)
    line95 = ROOT.TLine()
    line95.SetLineStyle(9)
    line95.SetLineColor(4)
    line95.DrawLine(-1,3.84,1,3.84)

    oneSig = ROOT.TPaveText(0.18,0.16,0.28,0.19,"NDC")
    oneSig.SetFillColor(0)
    oneSig.SetFillStyle(0)
    oneSig.SetTextFont(42)
    oneSig.SetBorderSize(0)
    oneSig.AddText("68\% CL")
    oneSig.Draw()

    twoSig = ROOT.TPaveText(0.18,0.24,0.28,0.27,"NDC")
    twoSig.SetFillColor(0)
    twoSig.SetFillStyle(0)
    twoSig.SetTextFont(42)
    twoSig.SetBorderSize(0)
    twoSig.AddText("95\% CL")
    twoSig.Draw()
