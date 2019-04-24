# -*- coding: utf-8 -*-
"""
Created on Fri Apr 05 12:46:09 2019

@author: Steven Hill
"""

def JupiterEWStatistics():
    
    import sys
    sys.path.append('f:\Astronomy\Projects\Jupiter\Spectral Data')
    sys.path.append('f:\\Astronomy\Python Play')
    sys.path.append('f:\\Astronomy\Python Play\Utils')
    sys.path.append('f:\\Astronomy\Python Play\Spectrophotometry\Spectroscopy')
    
    import numpy as np
    import EquivWidthUtils as EWU

    #CLRBandIDs=["543CH4","619CH4","705CH4","725CH4"]
    CLRBandIDs=["619CH4"]
    CLRDates=["20150123UT","20150209UT","20150210UT"]
    Grating="100lpm-550CLR"
    AlbedoEW=[]
    RawEW=[]

    for date in CLRDates:
        for band in CLRBandIDs:
            AEWs=EWU.EWObservations("../EWs/"+date+"-"+Grating+"-Albedo-EW.txt")
            AEWs.load_records("Jupiter",band)
            AlbedoEW.extend(AEWs.EW)
            REWs=EWU.EWObservations("../EWs/"+date+"-"+Grating+"-RawFlux-EW.txt")
            REWs.load_records("Jupiter",band)
            RawEW.extend(REWs.EW)
            
    AlbedoEW=np.array(AlbedoEW,dtype=float)
    print AlbedoEW
    RawEW=np.array(RawEW,dtype=float)
    print "Average AlbedoEW, Std. Dev., 95% Conf.= ",np.mean(AlbedoEW),np.std(AlbedoEW),\
            1.96*np.std(AlbedoEW)/np.sqrt(len(AlbedoEW)),len(AlbedoEW)    
    print "Average RawEW, Std. Dev., 95% Conf.= ",np.mean(RawEW),np.std(RawEW),\
            1.96*np.std(RawEW)/np.sqrt(len(RawEW)),len(RawEW)    

