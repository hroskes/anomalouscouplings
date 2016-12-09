#!/usr/bin/env python
from __future__ import print_function
from helperstuff import config
from helperstuff.combinehelpers import getrate
from helperstuff.enums import categories, channels
from helperstuff.filemanager import tfiles
from helperstuff.optimizesmoothing import ControlPlot
from helperstuff.templates import templatesfiles
from itertools import product
import os

print("neffectiveentries", "%contribution", "channel", "templategroup", "analysis", "production", "category", "productionmode", "hypothesis", "link")

totalrates = {}
for category, channel in product(categories, channels):
    totalrates[category,channel] = sum(
                                       getrate(productionmode, category, channel, "forexpectedscan", config.productionsforcombine[0])
                                          for productionmode in ("ggH", "VBF", "ZH", "WH", "ttH", "ggZZ", "qqZZ", "VBF bkg", "ZX")
                                      )

for tf in templatesfiles:
    for t in tf.templates():
        print(
              ControlPlot(t.discriminants[0], t).GetEffectiveEntries("raw"),
              str(getrate(t.productionmode, t.category, t.channel, "forexpectedscan", t.production) / totalrates[t.category, t.channel]*100)+"%",
              str(t).replace("D_int_prod ", "").replace("0-", "'0-"),
              "'" if t.hypothesis is None else "",
              '=HYPERLINK("https://hroskes.web.cern.ch/hroskes/anomalouscouplings_production/templateprojections/controlplots/{}/?match={}","link")'.format(os.path.basename(tf.templatesfile()).replace(".root", ""), t.templatename(final=False))
             )
