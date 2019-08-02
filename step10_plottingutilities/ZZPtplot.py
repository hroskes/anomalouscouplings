#!/usr/bin/env python

from discriminantplots import *

if __name__ == "__main__":
  with PlotCopier() as pc:
    plots = [
#      Plot(
#        name="ZZPt",
#        xtitle="p_{T}^{4l}",
#        ytitle="Events / bin",
#        hypothesislines=purehypothesislines + [a3mixVBF, a2mixVBF],
#        xformula="min(ZZPt, 300)",
#        cutformula=untaggedenrichcut,
#        binning = np.array([0, 60., 120, 200, 400]),
#        legendargs=(.2, .5, .9, .9),
#        categorylabel="Untagged",
#        legendcolumns=2,
#        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
#        ymax=50,
#        plotcopier=pc,
#        CMStext="",
#      ),
#      Plot(
#        name="ZZPt_zoom",
#        xtitle="p_{T}^{4l}",
#        ytitle="Events / bin",
#        hypothesislines=purehypothesislines + [a3mixVBF, a2mixVBF],
#        xformula="min(ZZPt, 650)",
#        cutformula=untaggedenrichcut,
#        binning = np.array([120, 200., 300, 400, 500, 600, 700]),
#        legendargs=(.2, .5, .9, .9),
#        categorylabel="Untagged",
#        legendcolumns=2,
#        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
#        ymax=50,
#        plotcopier=pc,
#        CMStext="",
#      ),
      Plot(
        name="D_STXS_stage1p1",
        xtitle="D_{STXS}^{1.1}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a3mixVBF, a2mixVBF],
        xformula="D_STXS_stage1p1",
        cutformula=masscut + " && (category_0P == 2 ? D_bkg_VBFdecay > 0.2 : (category_0P == 4 ? D_bkg_HadVHdecay > 0.2 : D_bkg > 0.7))",
        binning = np.arange(0., 23., 1.),
        legendargs=(.2, .5, .9, .9),
        categorylabel=None,
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="",
      ),
    ]
    for plot in plots: plot.makeplot()

