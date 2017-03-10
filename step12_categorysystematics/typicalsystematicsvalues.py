#!/usr/bin/env python

from helperstuff.enums import productionmodes, categories, analyses, channels
from helperstuff.yields import YieldSystematic, YieldSystematicValue

for productionmode in "ggH", "qqH", "ZH", "WH", "ttH", "qqZZ", "ggZZ", "VBFbkg", "ZX":
  for category in categories:
    if category == "Untagged": continue
    maxseen = 0
    for syst in YieldSystematic.items(lambda _: "cat" in str(_) or _ in ("JES", "bTagSF")):
      for analysis in analyses:
        for channel in channels:
          value = YieldSystematicValue(analysis, category, productionmode, syst, channel).value
          try:
            len(value)
          except:
            value = [value]
          for _ in value:
            if _ is None: _ = 1
            maxseen = max(maxseen, abs(_-1))
    print "{:10} {:20} {:.2%}".format(productionmode, category, maxseen)
