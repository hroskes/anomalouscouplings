#!/usr/bin/env python
"""
example script to copy smoothing parameters from one set of templates to another
"""

from templates import *

for channel in channels:
  for hypothesis in "0+_photoncut", "L1_photoncut", "L1Zg", "fL10.5_photoncut", "fL1Zg0.5", "fL10.5fL1Zg0.5":
    Template(channel, "ggH", "fL1fL1Zg", "170825", "Untagged", hypothesis).reweightrebin = [0, 1, 2]
  Template(channel, "qqZZ", "fL1fL1Zg", "170825", "Untagged").reweightrebin = [0, 1, 2]
  Template(channel, "ggZZ", "fL1fL1Zg", "170825", "Untagged").reweightrebin = [0, 1, 2]
  Template(channel, "VBF bkg", "fL1fL1Zg", "170825", "Untagged").reweightrebin = [0, 1, 2]

Template.writedict()
