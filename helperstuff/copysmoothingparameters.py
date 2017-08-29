#!/usr/bin/env python
"""
example script to copy smoothing parameters from one set of templates to another
"""

from templates import *

templates = []

for channel in channels:
  for hypothesis in "0+_photoncut", "L1_photoncut", "L1Zg", "fL10.5_photoncut", "fL1Zg0.5", "fL10.5fL1Zg0.5":
    templates.append(Template(channel, "ggH", "fL1fL1Zg", "170825", "Untagged", hypothesis))
  templates.append(Template(channel, "qqZZ", "fL1fL1Zg", "170825", "Untagged"))
  templates.append(Template(channel, "ggZZ", "fL1fL1Zg", "170825", "Untagged"))
  templates.append(Template(channel, "VBF bkg", "fL1fL1Zg", "170825", "Untagged"))


for template in templates:
  template.reweightrebin = None
  template.reweightaxes = [0, 1, 2]

Template.writedict()
