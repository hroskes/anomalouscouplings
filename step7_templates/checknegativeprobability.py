#!/usr/bin/env python

import json, itertools, numpy as np

from helperstuff.templates import *

tf = TemplatesFile("vh", "190703_2017", "4mu", "VBFtagged", "fa3fa2fL1fL1Zg")
with open(tf.jsonfile()) as f:
  templatenames = json.load(f)["constraints"][0]["templates"]

g1 = 0.7834214446898936
g2 = -0.08089243633861883
g4 = 0.10581969077638284
g1prime2 = -0.7459199848980587
ghzgs1prime2 = 0.035266732443872774

f = ROOT.TFile(tf.templatesfile())
templates = [getattr(f, name) for name in templatenames]
t = templates[0]

for xyz in itertools.product(xrange(1, t.GetNbinsX()+1), xrange(1, t.GetNbinsY()+1), xrange(1, t.GetNbinsZ()+1)):
  p = sum(
    a1*a2*a3*a4*t.GetBinContent(*xyz)
    for (a1, a2, a3, a4), t in itertools.izip_longest(
      itertools.combinations_with_replacement((g1, g4, g2, g1prime2, ghzgs1prime2), 4),
      templates
    )
  )
  if p<0:
    print xyz, p
    print np.array([t.GetBinContent(*xyz) for t in templates])
