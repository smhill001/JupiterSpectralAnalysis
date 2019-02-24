# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 07:22:23 2014

@author: steven.hill
"""
import sys
sys.path.append('g:\\Astronomy\Python Play')
import matplotlib.pyplot as pl
import pylab
import numpy as np
from scipy import interpolate
import scipy
from copy import deepcopy
import ComputeEW1 as CEW1
import quantities
from PyAstronomy import pyasl #This is where the best smoothing algorithm is!
import datetime

Jupiter_Karkoschka1993 = scipy.fromfile(file="../../Saturn Project 2013/Spectroscopy/Karkoschka/1993.tab.txt", dtype=float, count=-1, sep=" ")    
Jupiter_Karkoschka1993=scipy.reshape(Jupiter_Karkoschka1993,[Jupiter_Karkoschka1993.size/8,8])

Jupiter_KarkRef1993=np.zeros((Jupiter_Karkoschka1993.size/8,2))
Jupiter_KarkRef1993[:,0]=Jupiter_Karkoschka1993[:,0]*10.
Jupiter_KarkRef1993[:,1]=Jupiter_Karkoschka1993[:,3]
print "***Saturn_KarkRef1993***", Jupiter_KarkRef1993

CFNArray=["JupiterSpectrum-20150123082046UT-742-001thru010-sum600s-RotCrop-WVCal.dat", 
          "JupiterSpectrum-20150123083157UT-742-011thru020-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123084307UT-742-021thru030-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123085417UT-742-031thru040-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123090528UT-742-041thru050-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123091638UT-742-051thru060-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123092749UT-742-061thru070-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123093859UT-742-071thru080-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123095010UT-742-081thru090-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123100121UT-742-091thru100-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123101136UT-742-101thru110-SUM480s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123102358UT-742-111thru120-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123103514UT-742-121thru130-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123104625UT-742-131thru140-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123105734UT-742-141thru150-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123110844UT-742-151thru160-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123111954UT-742-161thru170-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123113104UT-742-171thru180-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123114214UT-742-181thru190-sum600s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123115324UT-742-191thru200-sum600s-RotCrop-WVCal.dat"]

KeyArray=["Jupiter20150123082046UT",
          "Jupiter20150123083157UT",
          "Jupiter20150123084307UT",
          "Jupiter20150123085417UT",
          "Jupiter20150123090528UT",
          "Jupiter20150123091638UT",
          "Jupiter20150123092749UT",
          "Jupiter20150123093859UT",
          "Jupiter20150123095010UT",
          "Jupiter20150123100121UT",
          "Jupiter20150123101136UT",
          "Jupiter20150123102358UT",
          "Jupiter20150123103514UT",
          "Jupiter20150123104625UT",
          "Jupiter20150123105734UT",
          "Jupiter20150123110844UT",
          "Jupiter20150123111954UT",
          "Jupiter20150123113104UT",
          "Jupiter20150123114214UT",
          "Jupiter20150123115324UT"]

for time in range (0,len(KeyArray)):
    CFN=CFNArray[time] #Clear file name array
    Key=KeyArray[time]
    
    # Read and reshape spectral data files    
    #NIR = scipy.fromfile(file="Data/20150209UT/JupiterSpectrum-20150209T034323UT-sum200s-001thru010-NIR-Aligned-RotCrop-LinearWV.dat", dtype=float, count=-1, sep='\t')    
    #NIR = scipy.fromfile(file="Data/20150209UT/JupiterSpectrum-20150209T034409UT-sum400s-001thru010-742-Aligned-RotCrop.dat", dtype=float, count=-1, sep='\t')    
    NIR = scipy.fromfile(file="Data/20150123UT/"+CFN, dtype=float, count=-1, sep='\t')    
    
    NormResponsewithWV= scipy.fromfile(file="PolluxResponse20150123UT.txt", dtype=float, count=-1, sep=" ")
    NIR=scipy.reshape(NIR,[NIR.size/2,2])
    NativeDispersion=(NIR[(NIR.size/2.-1),0]-NIR[0,0])/(NIR.size/2.-1.)
#    NIR[:,0]=NIR[:,0]+16.
    NRespWV=scipy.reshape(NormResponsewithWV,[NormResponsewithWV.size/2,2])
    MasterDispersion=(NRespWV[(NRespWV.size/2.-1),0]-NRespWV[0,0])/(NRespWV.size/2.-1.)
    
    #Load Reference Spectrum: Average G2v for albedo calculations
    Ref = scipy.loadtxt("g2v.dat", dtype=float, skiprows=3,usecols=(0,1))
    temp=Ref[:,1]
    temp=pyasl.smooth(temp,11,'flat')
    Ref[:,1]=temp
    
    #Interpolate NIR, Response and Reference spectra onto NIR Wavelengths
    
    NIRInterp=interpolate.interp1d(NIR[:,0],NIR[:,1],kind='linear', copy=True,
                             bounds_error=False, fill_value=0.0)  
    NIRonRef=NIRInterp(Ref[:,0])
        
    NRespInterp=interpolate.interp1d(NRespWV[:,0],NRespWV[:,1],kind='linear', copy=True,
                             bounds_error=False, fill_value=0.0)  
    NResponRef=NRespInterp(Ref[:,0])
    
    #Create Master Observed Spectrum by merging NIR and NIR spectra
       
    MASTER=deepcopy(Ref)
    MASTER[:,1]= NIRonRef
    
    #Compute EWs for telluric bands from MASTER
    EWFN=Key+"-RawFlux-EW.txt"
    Target="Jupiter"
    print "Key=",Key
    DateTime=datetime.datetime.strptime(Key[7:11]+"-"+Key[11:13]+"-" \
            +Key[13:15]+"T"+Key[15:17]+":"+Key[17:19]+":"+Key[19:21], \
            '%Y-%m-%dT%H:%M:%S')
    print "DateTime=",DateTime            

    BandType="Solar"
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"Ca II 8498",8480.,8510.,20.,EWFN,False)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"Ca II 8542",8520.,8560.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"Ca II 8662",8625.,8675.,20.,EWFN,True)
    BandType="Telluric"
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"O2 A band",7550.,7710.,40.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"H2O Z band",7800.,8470.,40.,EWFN,True)
    BandType="Target"
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"CH4 8420",8350.,8480.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"CH4 8620",8500.,8730.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"CH4 8890",8750.,9100.,60.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"CH4 8890N",8830.,8970.,20.,EWFN,True)
    BandName,BandCenter1,BandCenter2,HalfWidth,EW=CEW1.ComputeRatio1(MASTER,Target,DateTime,BandType,"829/889",8290.,8890.,70.,EWFN,True)    
    #WHAT ABOUT 7900 CH4 BAND????
    
    NativeDispersionNM=NativeDispersion/10.
    MasterDispersionNM=MasterDispersion/10.
    
    #Compute top of atmosphere spectrum
    
    BoxcarSmoothIndices=np.where(Ref[:,0] >= 8500.) #Do we need to shift here!!!???
    NResponRefSmooth=deepcopy(NResponRef)
    NResponRefSmooth[BoxcarSmoothIndices]=pyasl.smooth(NResponRefSmooth[BoxcarSmoothIndices],29,'flat')
    ToA=deepcopy(MASTER)
    ToA=MASTER[:,1]/NResponRefSmooth
    
    #Compute Albedo
    
    Albedo=ToA/Ref[:,1]
    mAlbedo = np.ma.masked_invalid(Albedo)
    AlbedoNormRangeIndices=np.where((Ref[:,0] >7400.) & \
         (Ref[:,0] < 8000.))
    
    NormAlbedo=Albedo/mAlbedo[AlbedoNormRangeIndices].max()
    print NormAlbedo.max()
    
    NormAlbedowithWV=deepcopy(Ref)
    NormAlbedowithWV[:,1]=NormAlbedo
    np.savetxt(Key+"Albedo.txt",NormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")

    EWFN=Key+"-Albedo-EW.txt"
    
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(NormAlbedowithWV,Target,DateTime,BandType,"CH4/NH3? 7900",7750.,8210.,40.,EWFN,False)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(NormAlbedowithWV,Target,DateTime,BandType,"CH4 8420",8350.,8480.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(NormAlbedowithWV,Target,DateTime,BandType,"CH4 8620",8500.,8730.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(NormAlbedowithWV,Target,DateTime,BandType,"CH4 8890",8750.,9100.,60.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(NormAlbedowithWV,Target,DateTime,BandType,"CH4 8890N",8830.,8970.,20.,EWFN,True)

    #Begin plotting 
    
    pl.figure(figsize=(6.5, 2.5), dpi=150, facecolor="white")
    
    pl.subplot(1, 1, 1)
    #Plot Layout Configuration
    x0=350
    x1=1050
    
    xtks=15
    y0=1.0e4
    y1=1.0e9
    #    ytks=9
    ExposureNIR = 600. #seconds
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
    pl.title(Key+" Spectrum",fontsize=9)
    pl.plot(Ref[:,0]/10.,NIRonRef/(ExposureNIR*Aperture*NativeDispersionNM),label='NIR',linewidth=0.5)
    pl.plot(MASTER[:,0]/10.,MASTER[:,1]/(ExposureNIR*Aperture*NativeDispersionNM),label='MASTER',color='k',linewidth=1)
    pl.plot(MASTER[:,0]/10.,ToA//(ExposureNIR*Aperture*NativeDispersionNM),label='Top of Atm.')
    pl.plot(MASTER[:,0]/10.,Ref[:,1]*1e7,label='Solar Ref. x 1e7')
    pl.plot(MASTER[:,0]/10.,NormAlbedo*1e7,label='Norm. Albedo x 1e7')
    
    pl.plot(Jupiter_KarkRef1993[:,0]/10.,Jupiter_KarkRef1993[:,1]*1.8e7,label='Karkoschka, 1993 x 1.8e7',linewidth=1,color='0.5')
    
    
    pl.legend(loc=0,ncol=3, borderaxespad=0.,prop={'size':6})
    pylab.savefig(Key+"Spectrum.png",dpi=300)
    
    
    TempMaster=MASTER
    TempMaster[:,0]=MASTER[:,0]/10.
    TempMaster[:,1]=MASTER[:,1]/(ExposureNIR*Aperture*NativeDispersionNM)
    np.savetxt(Key+"Spectrum.txt",MASTER,delimiter=" ",fmt="%10.3F %10.7F")
    np.savetxt(Key+"Albedo.txt",NormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")