# -*- coding: utf-8 -*-
"""
Created on Fri Mar 01 13:04:54 2019

@author: Steven Hill
"""
import sys
sys.path.append('f:\\Astronomy\Python Play')
sys.path.append('f:\\Astronomy\Python Play\Utils')
sys.path.append('f:\\Astronomy\Python Play\Spectrophotometry\Spectroscopy')
import matplotlib.pyplot as pl
import numpy as np
import scipy
from copy import deepcopy
import EquivWidthUtils as EWU
import PlotUtils as PU
import ConfigFiles as CF
import GeneralSpecUtils as GSU
from PyAstronomy import pyasl #This is where the best smoothing algorithm is!
import datetime

#Retrieve Target Parameters and create data paths
J=CF.Target_Parameters("f:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy/Target_Parameters.txt")
J.loadtargetparams("Jupiter")
JupPath=CF.built_path(J)
JupPath.spectra("20150123UT")

JupiterAlb=PU.PlotSetup("JupiterSpecPlotConfig.txt")
JupiterAlb.loadplotparams("f:","JupiterAlbedo","Spectra")

SPEC1 = scipy.fromfile(file=JupPath.One_D_path+"Jupiter20150123082046UT_100lpm-742NIR_Albedo.txt", dtype=float, count=-1, sep='\t')    
SPEC1=scipy.reshape(SPEC1,[SPEC1.size/2,2])
SPEC2 = scipy.fromfile(file=JupPath.One_D_path+"Jupiter20150123115324UT_100lpm-742NIR_Albedo.txt", dtype=float, count=-1, sep='\t')    
SPEC2=scipy.reshape(SPEC2,[SPEC2.size/2,2])
Ratio=GSU.SpectrumMath(SPEC2,SPEC1,"Divide")

JupiterAlb.Y1=1.2
JupiterAlb.Setup_Plot()
#pl.title(" Ratio ",fontsize=9)
pl.plot(SPEC1[:,0],SPEC1[:,1]*.56,label='Low Air Mass')
pl.plot(SPEC2[:,0],SPEC2[:,1]*.56,label='High Air Mass')
pl.plot(Ratio[:,0],Ratio[:,1]*1.0,label='Ratio')

pl.legend(loc=0,ncol=4, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.15)

