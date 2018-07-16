#!/usr/bin/env python

from helperstuff.constants import *
gp4 = {
  "fa3": g4VBF**2*g4HZZ**2,
  "fa2": g2VBF**2*g2HZZ**2,
  "fL1": g1prime2VBF**2*g1prime2HZZ**2,
  "fL1Zg": ghzgs1prime2VBF**2*ghzgs1prime2HZZ**2,
}

def f(analysis, production, channel, category):
    from helperstuff.templates import TemplatesFile
    tf = TemplatesFile(analysis, production, channel, category, "vbf")
    print "  Heshy  ", tf.templates()[1].gettemplate().Integral() * gp4[analysis] / tf.templates()[0].gettemplate().Integral()
    tf = TemplatesFile(analysis, str(production)+"_Ulascan", channel, category, "vbf")
    print "  Ulascan", tf.templates()[1].gettemplate().Integral() * gp4[analysis] / tf.templates()[0].gettemplate().Integral()


def ff(analysis, channel, production):
    from helperstuff.enums import categories
    for category in categories:
        print category
        f(analysis, production, channel, category)

import argparse
p = argparse.ArgumentParser()
p.add_argument("analysis", choices="fa3 fa2 fL1 fL1Zg".split())
p.add_argument("channel", choices="2e2mu 4e 4mu".split())
p.add_argument("production", choices="180530 180531".split())
args = p.parse_args()
ff(args.analysis, args.channel, args.production)
