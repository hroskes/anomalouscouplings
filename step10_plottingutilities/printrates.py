#!/usr/bin/env python

from collections import Counter
import sys
analysis = sys.argv[1]

from helperstuff import config
from helperstuff.combinehelpers import getdatatree, getrate
from helperstuff.enums import analyses, categories, channels

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

header = "{:10} {:10} {:10} {:10}"
row = "{:10} {:10.2f} {:10.2f} {:10.2f}"
obs = "{:10} {:10.0f} {:10.0f} {:10.0f}"

print "="*len(str(analysis))
print analysis
print "="*len(str(analysis))

print header.format("", *categories)

total = Counter({c: 0 for c in categories})
for p in "ggH", "VBF", "ZH", "WH", "ttH", "qqZZ", "ggZZ", "VBF bkg", "ZX":
    result = Counter({category: sum(getrate(p, channel, category, analysis, production, "fordata") for channel in channels) for category in categories})
    total += result
    print row.format(p, *(result[c] for c in categories))

print

print row.format("Expected", *(total[c] for c in categories))
print obs.format("Observed", *(sum(getdatatree(channel, category, analysis, production).GetEntries() for channel in channels) for category in categories))
