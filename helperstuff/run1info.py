from collections import namedtuple
from enums import Analysis

class AllowedInterval(namedtuple("AllowedInterval", "low hi")):
    def __init__(self, *args, **kwargs):
        super(AllowedInterval, self).__init__(*args, **kwargs)
        if self.low >= self.hi:
            raise ValueError("low >= hi! ({}, {})".format(self.low, self.hi))
    def __contains__(self, other):
        return self.low <= other <= self.hi
    

#http://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-14-018/CMS-HIG-14-018_Table_013.pdf
def allowedintervals2sigma(analysis):
    analysis = Analysis(analysis)
    if analysis == "fL1":
        return [AllowedInterval(-.25, .37)]
    if analysis == "fa2":
        return [AllowedInterval(-1, -.999), #this one isn't in the table, included so that -1 is allowed like +1
                AllowedInterval(-.66, -.57), AllowedInterval(-.15, 1)]
    if analysis == "fa3":
        return [AllowedInterval(-.4, .43)]

def isallowed2sigma(analysis, fai):
    if not -1 <= fai <= 1:
        raise ValueError("fai={} doesn't make sense!".format(fai))
    return any(fai in _ for _ in allowedintervals2sigma(analysis))
