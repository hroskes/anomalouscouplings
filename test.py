#!/usr/bin/env python

def getrels(*couplings):
  remaining = 1
  for coupling in couplings:
    if remaining == 0:
      assert coupling == 0
      rel = coupling
    else:
      rel = coupling / remaining
    assert -1<=rel<=1
    print rel
    remaining -= abs(coupling)

import argparse
p = argparse.ArgumentParser()
p.add_argument("couplings", nargs="+", type=float)
args = p.parse_args()
getrels(*args.couplings)
