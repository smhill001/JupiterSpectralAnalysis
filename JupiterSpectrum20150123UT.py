# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 07:22:23 2014

@author: steven.hill
"""
import sys
sys.path.append('d:\\Astronomy\Python Play')
import matplotlib.pyplot as pl
import pylab
import numpy as np
from scipy import interpolate
import scipy
from copy import deepcopy
import ComputeEW1
import datetime

Jupiter_Karkoschka1993 = scipy.fromfile(file="../../Saturn Project 2013/Spectroscopy/Karkoschka/1993.tab.txt", dtype=float, count=-1, sep=" ")    
Jupiter_Karkoschka1993=scipy.reshape(Jupiter_Karkoschka1993,[Jupiter_Karkoschka1993.size/8,8])

Jupiter_KarkRef1993=np.zeros((Jupiter_Karkoschka1993.size/8,2))
Jupiter_KarkRef1993[:,0]=Jupiter_Karkoschka1993[:,0]*10.
Jupiter_KarkRef1993[:,1]=Jupiter_Karkoschka1993[:,3]
print "***Saturn_KarkRef1993***", Jupiter_KarkRef1993


# Read and reshape spectral data files    
CLR = scipy.fromfile(file="Data/20150123UT/JupiterSpectrum-20150122LT-CLR-sum50m00s-Rotated-Cropped-WVCal.dat", dtype=float, count=-1, sep='\t')    
NIR = scipy.fromfile(file="Data/20150123UT/JupiterSpectrum-20150122LT-742-sum3h20m-Rotated-Cropped-WVCal.dat", dtype=float, count=-1, sep='\t')    
NormResponsewithWV= scipy.fromfile(file="PolluxResponse20150123UT.txt", dtype=float, count=-1, sep=" ")
CLR=scipy.reshape(CLR,[CLR.size/2,2])
NativeDispersion=(CLR[(CLR.size/2.-1),0]-CLR[0,0])/(CLR.size/2.-1.)
NIR=scipy.reshape(NIR,[NIR.size/2,2])
NIR[:,0]=NIR[:,0]+16.
NRespWV=scipy.reshape(NormResponsewithWV,[NormResponsewithWV.size/2,2])
MasterDispersion=(NRespWV[(NRespWV.size/2.-1),0]-NRespWV[0,0])/(NRespWV.size/2.-1.)

#Load Reference Spectrum: Average G2v for albedo calculations
Ref = scipy.loadtxt("g2v.dat", dtype=float, skiprows=3,usecols=(0,1))

#Interpolate NIR, Response and Reference spectra onto CLR Wavelengths

CLRInterp=interpolate.interp1d(CLR[:,0],CLR[:,1],kind='linear', copy=True,
                         bounds_error=False, fill_value=0.0)  
CLRonRef=CLRInterp(Ref[:,0])

NIRInterp=interpolate.interp1d(NIR[:,0],NIR[:,1],kind='linear', copy=True,
                         bounds_error=False, fill_value=0.0)  
NIRonRef=NIRInterp(Ref[:,0])

NRespInterp=interpolate.interp1d(NRespWV[:,0],NRespWV[:,1],kind='linear', copy=True,
                         bounds_error=False, fill_value=0.0)  
NResponRef=NRespInterp(Ref[:,0])

#Create Master Observed Spectrum by merging CLR and NIR spectra

MergeStartWV=7500.
MergeEndWV=7900.

OverlapRegionIndices=np.where((Ref[:,0] >MergeStartWV) & \
     (Ref[:,0] < MergeEndWV))

NIRScaling2CLR= scipy.mean(CLRonRef[OverlapRegionIndices]) \
    / scipy.mean(NIRonRef[OverlapRegionIndices])

MASTER=deepcopy(Ref)
MASTER[OverlapRegionIndices,1]= \
    (CLRonRef[OverlapRegionIndices]+NIRonRef[OverlapRegionIndices]*NIRScaling2CLR)/2.

DedicatedNIRIndices=np.where(Ref[:,0] >= 7900.)
DedicatedCLRIndices=np.where(Ref[:,0] <= 7500.)
MASTER[DedicatedNIRIndices,1]= NIRonRef[DedicatedNIRIndices]*NIRScaling2CLR
MASTER[DedicatedCLRIndices,1]= CLRonRef[DedicatedCLRIndices]

#Compute EWs for telluric bands from MASTER
"""EWFN="JupiterEW20150123UT.txt"
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Ca II H&K",3920.,3990.,20.,EWFN,False)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Delta",4090.,4115.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Ca I 'g band'",4215.,4240.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"G band 4300",4270.,4330.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Gamma",4330.,4350.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Fe I 'd band'",4365.,4405.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Beta",4815.,4890.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Mg 5170",5140.,5200.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Na D",5870.,5920.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Alpha 6563 band",6530.,6600.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Ca II 8498",8480.,8510.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Ca II 8542",8520.,8560.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Ca II 8662",8625.,8675.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"O2 6300 band",6260.,6360.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"O2 B band",6840.,6980.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H2O 7200a band",6910.,7090.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H2O 7200b band",7150.,7350.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"O2 A band",7550.,7710.,40.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H2O Z band",7800.,8470.,40.,EWFN,True)

BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 6190",6130.,6210.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 7250",7100.,7400.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 8620",8500.,8730.,20.,EWFN,True)
#BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 8890",8700.,9200.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 8890",8750.,9100.,60.,EWFN,True)
"""
NativeDispersionNM=NativeDispersion/10.
MasterDispersionNM=MasterDispersion/10.

#Compute top of atmosphere spectrum

BoxcarSmoothIndices=np.where(Ref[:,0] >= 8500.)
NResponRefSmooth=deepcopy(NResponRef)
NResponRefSmooth[BoxcarSmoothIndices]=scipy.convolve(NResponRef[BoxcarSmoothIndices],np.ones((29,))/29)[(28):]
ToA=deepcopy(MASTER)
ToA=MASTER[:,1]/NResponRefSmooth

#Compute Albedo

Albedo=ToA/Ref[:,1]
mAlbedo = np.ma.masked_invalid(Albedo)
AlbedoNormRangeIndices=np.where((Ref[:,0] >4000.) & \
     (Ref[:,0] < 7500.))

NormAlbedo=Albedo/mAlbedo[AlbedoNormRangeIndices].max()
print NormAlbedo.max()

NormAlbedowithWV=deepcopy(Ref)
NormAlbedowithWV[:,1]=NormAlbedo
np.savetxt("JupiterAlbedo20150123UT.txt",NormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")
#Begin plotting 


EWFN="JupiterKarkoschkaEW.txt"
Target="Jupiter"
DateTime=datetime.datetime.strptime("1993-01-01T00:00:00",'%Y-%m-%dT%H:%M:%S')
print "DateTime=",DateTime   
BandType="Methane"         

    
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Jupiter_KarkRef1993,Target,DateTime,BandType,"CH4 5430",5350.,5470.,60.,EWFN,False)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Jupiter_KarkRef1993,Target,DateTime,BandType,"CH4 6190",6100.,6280.,100.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Jupiter_KarkRef1993,Target,DateTime,BandType,"CH4 7050",6980.,7090.,40.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Jupiter_KarkRef1993,Target,DateTime,BandType,"CH4 7250",7140.,7450.,80.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Jupiter_KarkRef1993,Target,DateTime,BandType,"CH4/NH3? 7900",7750.,8210.,40.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Jupiter_KarkRef1993,Target,DateTime,BandType,"CH4 8420",8350.,8480.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Jupiter_KarkRef1993,Target,DateTime,BandType,"CH4 8620",8500.,8730.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Jupiter_KarkRef1993,Target,DateTime,BandType,"CH4 8890",8750.,9100.,60.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Jupiter_KarkRef1993,Target,DateTime,BandType,"CH4 8890N",8830.,8970.,20.,EWFN,True)



pl.figure(figsize=(6.5, 2.5), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0=350
x1=1050

xtks=15
y0=1.0e4
y1=1.0e9
#    ytks=9
ExposureCLR = 3000. #seconds
ExposureNIR = 12000. #seconds
Aperture = (0.135/22.)**2. #meters^2
# Set x limits
pl.xlim(x0,x1)
# Set x ticks
pl.xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
pl.ylim(y0,y1)
# Set y ticks
pl.yscale('log')
pl.grid()
pl.tick_params(axis='both', which='major', labelsize=7)
pl.ylabel("Counts Per Second",fontsize=7)
pl.xlabel("Wavelength (A)",fontsize=7)
pl.title("Jupiter Spectrum 20150123UT",fontsize=9)
pl.plot(Ref[:,0]/10.,CLRonRef/(ExposureCLR*Aperture*NativeDispersionNM),label='CLR',linewidth=0.5)
pl.plot(Ref[:,0]/10.,NIRonRef/(ExposureNIR*Aperture*NativeDispersionNM),label='NIRonCLR',linewidth=0.5)
pl.plot(MASTER[:,0]/10.,MASTER[:,1]/(ExposureCLR*Aperture*NativeDispersionNM),label='MASTER',color='k',linewidth=1)
pl.plot(MASTER[:,0]/10.,ToA//(ExposureCLR*Aperture*NativeDispersionNM),label='Top of Atm.')
pl.plot(MASTER[:,0]/10.,Ref[:,1]*1e7,label='Solar Ref. x 1e7')
pl.plot(MASTER[:,0]/10.,NormAlbedo*1e7,label='Norm. Albedo x 1e7')

pl.plot(Jupiter_KarkRef1993[:,0]/10.,Jupiter_KarkRef1993[:,1]*1.8e7,label='Karkoschka, 1993 x 1.8e7',linewidth=1,color='0.5')


pl.legend(loc=0,ncol=3, borderaxespad=0.,prop={'size':6})
pylab.savefig('JupiterSpectrum20150123UT.png',dpi=300)


TempMaster=MASTER
TempMaster[:,0]=MASTER[:,0]/10.
TempMaster[:,1]=MASTER[:,1]/(ExposureCLR*Aperture*NativeDispersionNM)
np.savetxt("JupiterSpectrum20150123UT.txt",MASTER,delimiter=" ",fmt="%10.3F %10.7F")
np.savetxt("JupiterAlbedo20150123UT.txt",NormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")