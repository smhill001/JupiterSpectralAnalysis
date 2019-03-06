# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 12:48:36 2019

@author: Steven Hill

Make this callable with Target, DateUT, and Band ID


"""

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

EWs=EWU.EWObservations("../EWs/20150331UT-100lpm-742NIR-Albedo-EW.txt")
EWs.load_records("Jupiter","889CH4")
print EWs.EW
print EWs.DateTimeUTObs
pl.figure(figsize=(8., 4.), dpi=150,facecolor="white")
pl.plot_date(EWs.DateTimeUTObs,EWs.EW)


######################

airmass=[]
for d in EWs.DateTimeUTObs:
    
    observer.date=d
    
    Planet = ephem.Jupiter(d)
    Planet.compute(observer)
    print Planet.cmlII
    print float(Planet.alt)*180./np.pi
    
    airmass.extend([pyasl.airmassPP(90.-Planet.alt*180./np.pi)])
    print airmass

pl.figure(figsize=(8., 4.), dpi=150,facecolor="white")
pl.scatter(np.array(airmass),EWs.EW)
Coefs=np.polyfit(np.array(airmass), EWs.EW, 1)
print Coefs
EWFit=Coefs[1]+Coefs[0]*np.array(airmass)
pl.plot(np.array(airmass),EWFit)
print np.corrcoef([np.array(airmass),EWs.EW])
