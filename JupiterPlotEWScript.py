# -*- coding: utf-8 -*-
"""
Created on Mon Apr 01 13:30:40 2019

@author: Steven Hill
"""

import JupiterPlotEW as JPEW


CLRBandIDs=["543CH4","619CH4","705CH4","725CH4"]
NIRBandIDs=["790CH4","842CH4","862CH4","889CH4"]

CLRDates=["20150123UT","20150209UT","20150210UT"]
NIRDates=["20150123UT","20150209UT","20150210UT","20150318UT","20150322UT",
          "20150331UT"]

for date in CLRDates:
    for band in CLRBandIDs:
        JPEW.JupiterPlotEW(date,band,"100lpm-550CLR")
        
for date in NIRDates:
    for band in NIRBandIDs:
        JPEW.JupiterPlotEW(date,band,"100lpm-742NIR")
        