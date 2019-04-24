# -*- coding: utf-8 -*-
"""
Created on Mon Apr 08 06:48:03 2019

@author: Steven Hill
"""

import numpy as np

nbins=2000
vertarray=np.linspace(0.0,1.0,nbins)
horzarray=np.zeros(20)
latarray=np.arcsin(vertarray)*180./np.pi
airarray=1.0/np.cos(latarray*np.pi/180.)
annulus=2.0*np.pi*(np.abs(vertarray))*(1.0/float(nbins)) #2*pi*r*dr
cirarray=airarray*annulus

#print latarray,airarray,annulus
print np.mean(airarray[0:nbins-1])
print np.sum(cirarray[0:nbins-1])/np.pi # pi*r^2 

def AirmassFactor_mu(lat,lon):
    import numpy as np
    
    lat0,lon0=0.0,0.0
    
    print np.cos(lat*np.pi/180.),np.cos(lon*np.pi/180.)
    arc=np.arccos(np.sin(lat0*np.pi/180.)*np.sin(lat*np.pi/180.)+
                  np.cos(lat0*np.pi/180.)*np.cos(lat*np.pi/180.)*
                  np.cos(lon*np.pi/180.))*180/np.pi
    
    print arc
    