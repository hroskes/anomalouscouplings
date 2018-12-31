import os

import utilities

thisfolder = os.path.dirname(os.path.abspath(__file__))
utilities.LoadMacro(os.path.join(thisfolder, "stage1.cc+"))

from ROOT import stage1_reco_stage1
