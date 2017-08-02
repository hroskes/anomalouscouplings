#!/usr/bin/env python
"""
example script to copy smoothing parameters from one set of templates to another
"""

from templates import *

for analysis in analyses:
  for hypothesis in "0+_photoncut", "L1_photoncut", "L1Zg", "fL10.5_photoncut", "fL1Zg0.5", "fL10.5fL1Zg0.5":
    Template('4e',    'ggH', analysis, 'LHE_170509', 'Untagged', hypothesis).smoothingparameters = \
    Template('2e2mu', 'ggH', analysis, 'LHE_170509', 'Untagged', hypothesis).smoothingparameters
  Template('4e',    'qqZZ', analysis, 'LHE_170509', 'Untagged').smoothingparameters = \
  Template('2e2mu', 'qqZZ', analysis, 'LHE_170509', 'Untagged').smoothingparameters

Template.writedict()
