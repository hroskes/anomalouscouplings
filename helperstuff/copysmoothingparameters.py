#!/usr/bin/env python
"""
example script to copy smoothing parameters from one set of templates to another
"""

from templates import *

for hypothesis in "0+_photoncut", "L1_photoncut", "L1Zg", "fL10.5_photoncut", "fL1Zg0.5", "fL10.5fL1Zg0.5":
  Template('2e2mu', 'ggH', "fL1fL1Zg_DeR_DeLint", 'LHE_170509', 'Untagged', hypothesis).smoothingparameters = \
  Template('2e2mu', 'ggH', "fL1fL1Zg_DL1_DL1L1Zgint", 'LHE_170509', 'Untagged', hypothesis).smoothingparameters
Template('2e2mu', 'qqZZ', "fL1fL1Zg_DeR_DeLint", 'LHE_170509', 'Untagged').smoothingparameters = \
Template('2e2mu', 'qqZZ', "fL1fL1Zg_DL1_DL1L1Zgint", 'LHE_170509', 'Untagged').smoothingparameters

Template.writedict()
