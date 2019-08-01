#!/usr/bin/env python

from discriminantplots import *

if __name__ == "__main__":
  with PlotCopier() as pc:
    Plot(
      name="ZZPt",
      xtitle="ZZPt",
      ytitle="Events / bin",
      hypothesislines=purehypothesislines + [a3mixVBF, a2mixVBF],
      xformula="min(ZZPt, 300)",
      cutformula=untaggedenrichcut,
      binning = np.array([0, 60., 120, 200, 400]),
      legendargs=(.2, .5, .9, .9),
      categorylabel="Untagged",
      legendcolumns=2,
      saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
      ymax=50,
      plotcopier=pc,
      CMStext="",
    ).makeplot()

