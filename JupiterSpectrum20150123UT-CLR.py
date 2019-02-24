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

CFNArray=["JupiterSpectrum-20150123071053UT-CLR-001thru015-sum300s-RotCrop-WVCal.dat", 
          "JupiterSpectrum-20150123071738UT-CLR-016thru030-sum300s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123072424UT-CLR-031thru045-sum300s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123073111UT-CLR-046thru060-sum300s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123073756UT-CLR-061thru075-sum300s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123074440UT-CLR-076thru090-sum300s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123075125UT-CLR-091thru105-sum300s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123075810UT-CLR-106thru120-sum300s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123080455UT-CLR-121thru135-sum300s-RotCrop-WVCal.dat",
          "JupiterSpectrum-20150123081140UT-CLR-136thru150-sum300s-RotCrop-WVCal.dat"]

KeyArray=["Jupiter20150123071053UT",
          "Jupiter20150123071738UT",
          "Jupiter20150123072424UT",
          "Jupiter20150123073111UT",
          "Jupiter20150123073756UT",
          "Jupiter20150123074440UT",
          "Jupiter20150123075125UT",
          "Jupiter20150123075810UT",
          "Jupiter20150123080455UT",
          "Jupiter20150123081140UT"]

for time in range (0,len(KeyArray)):
    CFN=CFNArray[time] #Clear file name array
    Key=KeyArray[time]
    
    # Read and reshape spectral data files    
    #CLR = scipy.fromfile(file="Data/20150209UT/JupiterSpectrum-20150209T034323UT-sum200s-001thru010-CLR-Aligned-RotCrop-LinearWV.dat", dtype=float, count=-1, sep='\t')    
    #NIR = scipy.fromfile(file="Data/20150209UT/JupiterSpectrum-20150209T034409UT-sum400s-001thru010-742-Aligned-RotCrop.dat", dtype=float, count=-1, sep='\t')    
    CLR = scipy.fromfile(file="Data/20150123UT/"+CFN, dtype=float, count=-1, sep='\t')    
    
    NormResponsewithWV= scipy.fromfile(file="PolluxResponse20150123UT.txt", dtype=float, count=-1, sep=" ")
    CLR=scipy.reshape(CLR,[CLR.size/2,2])
    NativeDispersion=(CLR[(CLR.size/2.-1),0]-CLR[0,0])/(CLR.size/2.-1.)
#    NIR[:,0]=NIR[:,0]+16.
    NRespWV=scipy.reshape(NormResponsewithWV,[NormResponsewithWV.size/2,2])
    MasterDispersion=(NRespWV[(NRespWV.size/2.-1),0]-NRespWV[0,0])/(NRespWV.size/2.-1.)
    
    #Load Reference Spectrum: Average G2v for albedo calculations
    Ref = scipy.loadtxt("g2v.dat", dtype=float, skiprows=3,usecols=(0,1))
    temp=Ref[:,1]
    temp=pyasl.smooth(temp,11,'flat')
    Ref[:,1]=temp
    
    #Interpolate NIR, Response and Reference spectra onto CLR Wavelengths
    
    CLRInterp=interpolate.interp1d(CLR[:,0],CLR[:,1],kind='linear', copy=True,
                             bounds_error=False, fill_value=0.0)  
    CLRonRef=CLRInterp(Ref[:,0])
        
    NRespInterp=interpolate.interp1d(NRespWV[:,0],NRespWV[:,1],kind='linear', copy=True,
                             bounds_error=False, fill_value=0.0)  
    NResponRef=NRespInterp(Ref[:,0])
    
    #Create Master Observed Spectrum by merging CLR and NIR spectra
       
    MASTER=deepcopy(Ref)
    MASTER[:,1]= CLRonRef
    
    #Compute EWs for telluric bands from MASTER
    EWFN=Key+"-RawFlux-EW.txt"
    Target="Jupiter"
    print "Key=",Key
    DateTime=datetime.datetime.strptime(Key[7:11]+"-"+Key[11:13]+"-" \
            +Key[13:15]+"T"+Key[15:17]+":"+Key[17:19]+":"+Key[19:21], \
            '%Y-%m-%dT%H:%M:%S')
    print "DateTime=",DateTime            

    BandType="Solar"
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"Fe 3820",3820.,3860.,20.,EWFN,False)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"Ca II H&K",3920.,3990.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"H Delta",4085.,4113.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"Ca I 'g band'",4215.,4240.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"G band 4300",4270.,4330.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"H Gamma",4330.,4350.,10.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"Fe I 'd band'",4365.,4405.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"H Beta",4815.,4890.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"Mg 5170",5140.,5200.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"Fe I 5270",5240.,5290.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"Na D",5870.,5920.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"H Alpha 6563 band",6530.,6590.,20.,EWFN,True)
    BandType="Telluric"
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"O2 6300 band",6260.,6360.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"O2 B band",6850.,6920.,50.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"H2O 7200a band",6910.,7090.,10.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"H2O 7200b band",7150.,7350.,20.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"O2 A band",7540.,7720.,40.,EWFN,True)
    BandType="Target"
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"CH4 5430",5350.,5470.,60.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"CH4 6190",6100.,6280.,100.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"CH4 7050",6980.,7090.,40.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(MASTER,Target,DateTime,BandType,"CH4 7250",7140.,7450.,80.,EWFN,True)
    
    #WHAT ABOUT 7900 CH4 BAND????
    
    NativeDispersionNM=NativeDispersion/10.
    MasterDispersionNM=MasterDispersion/10.
    
    #Compute top of atmosphere spectrum
    
    BoxcarSmoothIndices=np.where(Ref[:,0] >= 8500.) #Do we need to shift here!!!???
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
    np.savetxt(Key+"Albedo.txt",NormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")

    EWFN=Key+"-Albedo-EW.txt"
    
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(NormAlbedowithWV,Target,DateTime,BandType,"CH4 5430",5350.,5470.,60.,EWFN,False)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(NormAlbedowithWV,Target,DateTime,BandType,"CH4 6190",6100.,6280.,100.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(NormAlbedowithWV,Target,DateTime,BandType,"CH4 7050",6980.,7090.,40.,EWFN,True)
    BandName,BandStart,BandEnd,ContWidth,EW=CEW1.ComputeEW1(NormAlbedowithWV,Target,DateTime,BandType,"CH4 7250",7140.,7450.,80.,EWFN,True)

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
    ExposureCLR = 300. #seconds
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
    pl.plot(Ref[:,0]/10.,CLRonRef/(ExposureCLR*Aperture*NativeDispersionNM),label='CLR',linewidth=0.5)
    pl.plot(MASTER[:,0]/10.,MASTER[:,1]/(ExposureCLR*Aperture*NativeDispersionNM),label='MASTER',color='k',linewidth=1)
    pl.plot(MASTER[:,0]/10.,ToA//(ExposureCLR*Aperture*NativeDispersionNM),label='Top of Atm.')
    pl.plot(MASTER[:,0]/10.,Ref[:,1]*1e7,label='Solar Ref. x 1e7')
    pl.plot(MASTER[:,0]/10.,NormAlbedo*1e7,label='Norm. Albedo x 1e7')
    
    pl.plot(Jupiter_KarkRef1993[:,0]/10.,Jupiter_KarkRef1993[:,1]*1.8e7,label='Karkoschka, 1993 x 1.8e7',linewidth=1,color='0.5')
    
    
    pl.legend(loc=0,ncol=3, borderaxespad=0.,prop={'size':6})
    pylab.savefig(Key+"Spectrum.png",dpi=300)
    
    
    TempMaster=MASTER
    TempMaster[:,0]=MASTER[:,0]/10.
    TempMaster[:,1]=MASTER[:,1]/(ExposureCLR*Aperture*NativeDispersionNM)
    np.savetxt(Key+"Spectrum.txt",MASTER,delimiter=" ",fmt="%10.3F %10.7F")
    np.savetxt(Key+"Albedo.txt",NormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")