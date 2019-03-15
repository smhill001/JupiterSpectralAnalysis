# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 12:48:36 2019

@author: Steven Hill

Make this callable with Target, DateUT, and Band ID


"""
def JupiterPlotEW(DateUT,BandID,Grating):
    import sys
    sys.path.append('f:\Astronomy\Projects\Jupiter\Spectral Data')
    sys.path.append('f:\\Astronomy\Python Play')
    sys.path.append('f:\\Astronomy\Python Play\Utils')
    sys.path.append('f:\\Astronomy\Python Play\Spectrophotometry\Spectroscopy')
    
    import astropy.units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz
    import numpy as np
    import ephem
    import EquivWidthUtils as EWU
    from PyAstronomy import pyasl
    import matplotlib.pyplot as pl
    import PlotUtils as PU
    
    TimePlot=PU.PlotSetup("JupiterEWSpecPlotConfig.txt")
    TimePlot.loadplotparams("f:","Jupiter_"+BandID,"Time")
    AirmassPlot=PU.PlotSetup("JupiterEWSpecPlotConfig.txt")
    AirmassPlot.loadplotparams("f:","Jupiter_"+BandID,"Airmass")
    CMIIPlot=PU.PlotSetup("JupiterEWSpecPlotConfig.txt")
    CMIIPlot.loadplotparams("f:","Jupiter_"+BandID,"CMII")

    Pollux = SkyCoord.from_name('Pollux')
    datestring='2015-01-23 08:20:00'
    observer = ephem.Observer()
    #Location from Google Maps at 483 S Oneida Way, Denver 80224
    observer.lon = ephem.degrees('-104.907985')
    observer.lat = ephem.degrees('39.708200')
    print float(observer.lat)*180/np.pi
    Oneida = EarthLocation(lat=float(observer.lat)*u.rad, 
                           lon=float(observer.lon)*u.rad, height=1500*u.m)
    utcoffset = 0*u.hour  # Universal Time
    time = Time('2015-01-23 05:54:00') - utcoffset
    
    Polluxaltaz = Pollux.transform_to(AltAz(obstime=time,location=Oneida))
    print Polluxaltaz.alt.deg
    print("Pollux's Altitude = {0.alt:.2}".format(Polluxaltaz))
    print Polluxaltaz.secz
    
    EWs=EWU.EWObservations("../EWs/"+DateUT+"-"+Grating+"-Albedo-EW.txt")
    EWs.load_records("Jupiter",BandID)
    #print EWs.DateTimeUTObs
    #print EWs.EW
    TimePlot.Setup_PlotDate(min(EWs.DateTimeUTObs),max(EWs.DateTimeUTObs),1)    

    pl.plot_date(EWs.DateTimeUTObs,EWs.EW)
    print "Average EW, Std. Dev., 95% Conf.= ",np.mean(EWs.EW),np.std(EWs.EW),\
            1.96*np.std(EWs.EW)/np.sqrt(len(EWs.EW))
    
    
    ######################
    
    airmass=[]
    CMII=[]
    for d in EWs.DateTimeUTObs:
        
        observer.date=d
        
        Planet = ephem.Jupiter(d)
        Planet.compute(observer)
        #print float(Planet.alt)*180./np.pi
        
        airmass.extend([pyasl.airmassPP(90.-Planet.alt*180./np.pi)])
        #print Planet.cmlII
        CMII.append(Planet.cmlII*180./np.pi)
        #print airmass
    
    AirmassPlot.Setup_Plot()
    pl.scatter(np.array(airmass),EWs.EW)
    Coefs=np.polyfit(np.array(airmass), EWs.EW, 1)
    EWFit=Coefs[1]+Coefs[0]*np.array(airmass)
    pl.plot(np.array(airmass),EWFit)
    corr=np.corrcoef([np.array(airmass),EWs.EW])
    print "Coefs= ",Coefs
    print "Correlation Coef.= ",corr[0,1]
    #CoefsUT=np.polyfit(np.array(EWs.DateTimeUTObs), EWs.EW, 1)
    CMIIPlot.Setup_Plot()
    pl.scatter(np.array(CMII),EWs.EW)
    