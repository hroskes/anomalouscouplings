rates = """ggH 2.03655774801
qqH 0.117870293266
WH 0.050917903693
ZH 0.045285083273
ttH 0.0215524981897
ggH 0.0511681752353
qqH 0.0826491672701
WH 0.0012070057929
ZH 0.00102389826213
ttH 0.00162274800869
qqZZ 2.44504272266
qqZZ 0.0125085590152
ggZZ 0.181148443157
ggZZ 0.00478279942071
zjets 0.687494279508
zjets 0.0198287038378
ggH 0.863880340333
qqH 0.0503571687183
WH 0.0224938993483
ZH 0.0184693428675
ttH 0.00956001176684
ggH 0.0220349203476
qqH 0.0358187997828
WH 0.000541927226647
ZH 0.000417196958726
ttH 0.0007427108979
qqZZ 0.912062997828
qqZZ 0.00204765966691
ggZZ 0.109306444605
ggZZ 0.00288597682839
zjets 0.465475162925
zjets 0.0134251991311
ggH 1.5417487328
qqH 0.0885156589428
WH 0.038228756336
ZH 0.0329028511948
ttH 0.0155734793628
ggH 0.0378175506879
qqH 0.0611885409124
WH 0.000896106444605
ZH 0.000783670528602
ttH 0.00112499275887
qqZZ 2.00359160029
qqZZ 0.0079220333092
ggZZ 0.211314844316
ggZZ 0.00557928892107
zjets 0.518448081101
zjets 0.0149530340333"""


rates = [line.split() for line in rates.split("\n")]

print "VBF", sum(float(rate) for name, rate in rates if name == "qqH")
print "VH", sum(float(rate) for name, rate in rates if name in ("ZH", "WH"))
print "other sig", sum(float(rate) for name, rate in rates if name in ("ggH", "ttH"))
print "EW", sum(float(rate) for name, rate in rates if name in ("qqZZ", "ggZZ"))
print "ZX", sum(float(rate) for name, rate in rates if name == "zjets")
print "total", sum(float(rate) for name, rate in rates if True)
