# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 07:20:03 2019

This program computes the condensation altitudes of methane, water,
ammonia, and ammonium hydrosulfide in the Jovian atmosphere.

I need to put the references for this model here.

@author: Steven Hill
"""

import sys
sys.path.append('f:\\Astronomy\Python Play')
sys.path.append('f:\\Astronomy\Python Play\Util')
sys.path.append('f:\\Astronomy\Python Play\Spectrophotometry\Spectroscopy')
import matplotlib.pyplot as pl
import numpy as np
from thermo.chemical import Chemical
import PlotUtils as PU
import ConfigFiles as CF
import JupiterModelLibrary as JML
from scipy import interpolate

URL="https://atmos.nmsu.edu/pdsd/archive/data/gp-j-asi-3-entry-v10/gp_0001/data/asi/loweratm.tab"
Lower_Atm=JML.Jupiter_Galileo_Profile(URL)
Lower_Atm.load_records()

R_star=3.75 # J/(g-K)

CH4=Chemical('methane')
NH3=Chemical('ammonia')
NH4SH=Chemical('Ammonium hydrosulfide')
H2O=Chemical('water')

P=np.geomspace(100.,1.0e6,num=50,endpoint=True,dtype=float) #in Pascals
T_NH3=[]
T_CH4=[]
T_NH4SH=[]
T_H2O=[]
print P.size
for i in range(0,P.size):
    T_NH3.append(NH3.VaporPressure.solve_prop(2.0e-4*P[i]))
    T_CH4.append(CH4.VaporPressure.solve_prop(0.02*P[i]))
    T_H2O.append(H2O.VaporPressure.solve_prop(2.0e-4*P[i]))
    #T_NH4SH.append(NH4SH.VaporPressure.solve_prop(P[i]))

print "Back in Main"
print Lower_Atm.Press[0:5], Lower_Atm.Temp[0:5], Lower_Atm.Grav[0:5]
print ' '
print R_star*np.array(Lower_Atm.Temp[0:5])/np.array(Lower_Atm.Grav[0:5])
print "Done with Print"

############ NH3 CLOUD
Interp_T_NH3=interpolate.interp1d(P,T_NH3,kind='linear', 
                            copy=True,bounds_error=False, 
                            fill_value=np.NaN,axis=0)  
NH3_T_on_Galileo_P=Interp_T_NH3(np.array(Lower_Atm.Press)*1.0e5)

print "NH3_T_on_Galileo_P=",NH3_T_on_Galileo_P

diff=(np.abs(NH3_T_on_Galileo_P-np.array(Lower_Atm.Temp)))
finite=np.isfinite(diff)
findiff=diff[np.where(np.isfinite(diff))]
finNH3_T=NH3_T_on_Galileo_P[np.where(np.isfinite(diff))]
finP_T=(np.array(Lower_Atm.Press)*1.0e5)[np.where(np.isfinite(diff))]
print np.min(findiff),np.argmin(findiff)
mindex=np.argmin(findiff)
print finNH3_T[mindex],(np.array(Lower_Atm.Press)*1.0e5)[mindex]

############ H2O CLOUD

Interp_T_H2O=interpolate.interp1d(P,T_H2O,kind='linear', 
                            copy=True,bounds_error=False, 
                            fill_value=np.NaN,axis=0)  
H2O_T_on_Galileo_P=Interp_T_H2O(np.array(Lower_Atm.Press)*1.0e5)

print "H2O_T_on_Galileo_P=",H2O_T_on_Galileo_P

diff1=(np.abs(H2O_T_on_Galileo_P-np.array(Lower_Atm.Temp)))
finite1=np.isfinite(diff1)
findiff1=diff1[np.where(np.isfinite(diff1))]
finH2O_T=H2O_T_on_Galileo_P[np.where(np.isfinite(diff1))]
finP_TH2O=(np.array(Lower_Atm.Press)*1.0e5)[np.where(np.isfinite(diff1))]
print np.min(findiff1),np.argmin(findiff1)
mindex1=np.argmin(findiff1)
print finH2O_T[mindex1],(np.array(Lower_Atm.Press)*1.0e5)[mindex1]

pl.figure(figsize=(6.0, 6.0), dpi=150, facecolor="white")
pl.yscale('log')
pl.plot(T_NH3,P,label='NH3')
pl.plot(T_CH4,P,label='CH4')
pl.plot(T_H2O,P,label='H2O')
pl.ylim(10000.,1000000.)
ax = pl.gca() 
ax.set_ylim(ax.get_ylim()[::-1]) #reverse y-axis
ax.set_xlim(50.,300.) #reverse y-axis
pl.xscale('linear')
pl.grid(which='both')
pl.plot(np.array(Lower_Atm.Temp),np.array(Lower_Atm.Press)*1.0e5)

pl.title("Jupiter Cloud Levels")
pl.ylabel("Pressure (Pa)")
pl.xlabel("Temperature (K)")

pl.legend(loc=0,ncol=4, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.12, right=0.96, top=0.95, bottom=0.1)
pl.savefig("Jupiter Atmosphere Model.png",dpi=300)        



