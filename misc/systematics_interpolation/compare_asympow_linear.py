#!/usr/bin/env python

import numpy as np
import scipy.interpolate

class AsymPow(object):
  def __init__(self, downnominalup):
    down, nominal, up = downnominalup
    self.__nominal = nominal
    self.__kappaLow = down / nominal
    self.__kappaHigh = up / nominal

  def __call__(self, x):
    return self.__nominal * np.exp(self.logKappaForX(x) * x)

  def logKappaForX(self, x):
    if abs(x) >= 0.5: return np.log(self.__kappaHigh) if x >= 0 else -np.log(self.__kappaLow)
    # interpolate between log(kappaHigh) and -log(kappaLow)
    #    logKappa(x) = avg + halfdiff * h(2x)
    # where h(x) is the 3th order polynomial
    #    h(x) = (3 x^5 - 10 x^3 + 15 x)/8;
    # chosen so that h(x) satisfies the following:
    #      h (+/-1) = +/-1
    #      h'(+/-1) = 0
    #      h"(+/-1) = 0
    logKhi =  np.log(self.__kappaHigh)
    logKlo = -np.log(self.__kappaLow)
    avg = 0.5*(logKhi + logKlo)
    halfdiff = 0.5*(logKhi - logKlo)
    twox = x+x
    twox2 = twox*twox
    alpha = 0.125 * twox * (twox2 * (3*twox2 - 10.) + 15.)
    ret = avg + alpha*halfdiff;
    assert alpha >= -1 and alpha <= 1, "Something is wrong in the interpolation"
    return ret

x = np.array([-1., 0, 1])
y = np.array([9., 10, 12])

linear = scipy.interpolate.interp1d(x, y, fill_value="extrapolate", kind="linear")
lnN = AsymPow(y)

for i in np.arange(-7, 7.5, .5):
  print i, linear(1.*i+.5) / lnN(1.*i+.5)
