from helperstuff.enums import *
from helperstuff.samples import *
from helperstuff.utilities import *
from collections import Counter

samples = [Sample("160928", "VBF", hypothesis) for hypothesis in prodonlyhypotheses]

othersamples = [ReweightingSample("VBF", h) for h in "fa2prod-0.5", "fa2dec-0.5", "fa2proddec0.5"]

for s in samples:
    print s
    print
    f = tfiles[s.withdiscriminantsfile()]
    t = f.candTree
    c = Counter()
    for entry in t:
        for sample in othersamples:
            if getattr(t, sample.weightname()) < 0:
                c[sample] += 1
    for sample in othersamples:
        print sample, c[sample]
