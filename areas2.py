from __future__ import division
from math import exp, expm1, pow
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt


#Parameters

rec = 0.2

#define subjective probabilities over two periods taking into account of correlation


def mUU(c, h):
    out = c*h*(1 - h) + pow(h,2)
    return out
def mUD(c,h):
    out = 2*h*(h - 1)*(c - 1)
    return out
def mDD(c,h):
    out = pow((1 - h),2) + c*(h -pow(h,2))
    return out


#area division for period 0
def payoffs0(p,c,h):
    ps0 = p[0]; pa0 = p[1]; pjDD = p[2]; pjUD = p[3]; psDD =p[4]; paUD = p[5]; paDD = p[6]
    cash = 1
    junior = (mUU(c,h)*1 + mUD(c,h)*pjUD + mDD(c,h)*pjDD)/(2*pa0 - ps0)
    asset = (mUU(c,h)* 1 + mUD(c,h)*paUD + mDD(c,h)*paDD - paDD)/(pa0 - paDD)
    senior =(mUU(c,h)* 1 + mUD(c,h)* 1 + mDD(c,h)*psDD - psDD)/(pa0 - psDD)
    payoffs = np.array([cash,junior,asset,senior])
    return payoffs

def Holder0(p,c,h):
    payoffs = -payoffs0(p,c,h)
    index=payoffs.argmin()
    holder = [0,0,0,0]
    holder [index] = 1
    return holder

#area division for period 1, UD, in which senior tranche is the same as cash

def payoffsUD(p,c,h):
    ps0 = p[0]; pa0 = p[1]; pjDD = p[2]; pjUD = p[3]; psDD =p[4]; paUD = p[5]; paDD = p[6]
    cash=1
    junior = (((mUU(c,h)) + mUD(c,h)/2)*1 + (mUD(c,h)/2 + mDD(c,h))*rec - rec)/ (pjUD - rec)
    asset = (((mUU(c,h)) + mUD(c,h)/2)*1 + (mUD(c,h)/2 + mDD(c,h))*((1+rec)/2) - ((1+rec)/2))/(paUD - ((1+rec)/2))
    payoffs = np.array([cash,junior,asset])
    return payoffs

def HolderUD(p,c,h):
    payoffs = -payoffsUD(p,c,h)
    index=payoffs.argmin()
    holder = [0,0,0]
    holder [index] = 1
    return holder

#area division for period 1, DD

def payoffsDD(p,c,h):
    ps0 = p[0]; pa0 = p[1]; pjDD = p[2]; pjUD = p[3]; psDD =p[4]; paUD = p[5]; paDD = p[6]
    cash=1
    junior = (mUU(c, h)*1+ mUD(c, h)*rec + mDD(c, h)*0)/pjDD
    asset = (mUU(c,h)*1 + mUD(c,h)*((1+rec)/2) + mDD(c,h)*rec - rec)/(paDD - rec)
    senior = (mUU(c,h)*1 + mUD(c,h)*1 + mDD(c,h)*2*rec - 2*rec)/(psDD - 2*rec)
    payoffs = np.array([cash,junior,asset,senior])
    return payoffs

def HolderDD(p,c,h):
    payoffs = -payoffsUD(p,c,h)
    index=payoffs.argmin()
    holder = [0,0,0,0]
    holder [index] = 1
    return holder




p = [0.9,0.7,0.25,0.45,0.5,0.81,0.3]

for i in range(0,100):
    print mUU(i/100,i/100), mUD(i/100,i/100), mDD(i/100,i/100)
    print payoffs0(p,i/100,i/100)
    print Holder0(p,i/100,i/100)
    print HolderUD(p,i/100,i/100)
    print HolderDD(p,i/100,i/100)