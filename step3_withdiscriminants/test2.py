from helperstuff.constants import SMXSZH2L2l
from helperstuff.filemanager import tfiles
from helperstuff.samples import Sample

for h in "0+", "fL1prod0.5":
    f = tfiles[Sample(h, "ZH", "160909").withdiscriminantsfile()]
    CJLST = tfiles[Sample(h, "ZH", "160909").CJLSTfile()]
    hc = CJLST.ZZTree.Get("Counters_reweighted")
    n2L2l = hc.GetBinContent(4,1)+hc.GetBinContent(8,1)
    total = 0
    t = f.candTree
    length = t.GetEntries()
    for i, entry in enumerate(t, start=1):
        total += t.MC_weight_ZH_g1
        if i % 10000 == 0 or i == length:
            print i, "/", length
    print h, total, total / SMXSZH2L2l
