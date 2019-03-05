# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 12:48:36 2019

@author: Steven Hill
"""

import sys
sys.path.append('f:\Astronomy\Projects\Jupiter\Spectral Data')

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np
import ephem
from PyAstronomy import pyasl

datestring='2015-01-23 08:20:00'
observer = ephem.Observer()
#Location from Google Maps at 483 S Oneida Way, Denver 80224
observer.lon = ephem.degrees('-104.907985')
observer.lat = ephem.degrees('39.708200')
print float(observer.lat)*180/np.pi
observer.date=datestring

Planet = ephem.Jupiter(datestring)
Planet.compute(observer)
print Planet.cmlII
print float(Planet.alt)*180./np.pi

airmass=pyasl.airmassPP(90.-Planet.alt*180./np.pi)
print airmass

Pollux = SkyCoord.from_name('Pollux')
Oneida = EarthLocation(lat=float(observer.lat)*u.rad, 
                       lon=float(observer.lon)*u.rad, height=1500*u.m)
utcoffset = 0*u.hour  # Universal Time
time = Time('2015-01-23 05:54:00') - utcoffset

Polluxaltaz = Pollux.transform_to(AltAz(obstime=time,location=Oneida))
print Polluxaltaz.alt.deg
print("Pollux's Altitude = {0.alt:.2}".format(Polluxaltaz))
print Polluxaltaz.secz


