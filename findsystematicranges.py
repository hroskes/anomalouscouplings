#!/usr/bin/env python

import glob, pprint, collections

badlines = collections.Counter()

for filename in glob.iglob("slurm-*.out"):
  with open(filename) as f:
    for line in f:
      if "Number of events is negative or error @ params" in line:
        line = line.split("(")[1].split(")")[0]
        varsandvalues = line.split(",")
        varsandvalues = {varandvalue.split()[0]: float(varandvalue.split()[2]) for varandvalue in varsandvalues}

        if any(abs(varsandvalues[_]) > 2 for _ in varsandvalues if "CMS_scale_j" in _):
          badlines[0] += 1
        elif any(abs(varsandvalues[_]) > 2 for _ in varsandvalues):
          badlines[1] += 1
        elif any(abs(varsandvalues[_]) > 1 for _ in varsandvalues if "CMS_scale_j" in _):
          badlines[2] += 1
        elif any(abs(varsandvalues[_]) > 2 for _ in varsandvalues if "CMS_scale" in _ or "CMS_res" in _ or "THU" in _):
          badlines[3] += 1
        else:
          badlines[4] += 1

  pprint.pprint(badlines)
