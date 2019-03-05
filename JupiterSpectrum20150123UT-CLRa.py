# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 07:22:23 2014

@author: steven.hill
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

DateUT="20150123UT"
Target="Jupiter"

#Retrieve Target Parameters and create data paths
J=CF.Target_Parameters("f:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy/Target_Parameters.txt")
J.loadtargetparams(Target)
JupPath=CF.built_path(J)
JupPath.spectra(DateUT)
    
#Load response calibration and solar reference spectrum
Response = scipy.fromfile(file="../PolluxResponse20150123UT.txt", dtype=float, count=-1, sep=" ")
Response=scipy.reshape(Response,[Response.size/2,2])
Response[:,0]=Response[:,0]/10.
Response[:,1]=pyasl.smooth(Response[:,1],3,'flat')
MasterDispersion=(Response[(Response.size/2-1),0]-Response[0,0])/(Response.size/2-1)

Ref_g2v = scipy.loadtxt(JupPath.reference_path+J.SpecType, dtype=float, skiprows=3,usecols=(0,1))
Ref_g2v[:,0]=Ref_g2v[:,0]/10.
Ref_g2v[:,1]=pyasl.smooth(Ref_g2v[:,1],3,'flat')

#Load comparison albedo spectrum from Karkoschka, 1994 (1993 observations)
Jupiter_Karkoschka1993 = scipy.fromfile(file="f:/Astronomy/Projects/Planets/Saturn/Spectral Data/Karkoschka/1993.tab.txt", dtype=float, count=-1, sep=" ")    
Jupiter_Karkoschka1993=scipy.reshape(Jupiter_Karkoschka1993,[Jupiter_Karkoschka1993.size/8,8])
Jupiter_KarkRef1993=np.zeros((Jupiter_Karkoschka1993.size/8,2))
Jupiter_KarkRef1993[:,0]=Jupiter_Karkoschka1993[:,0]
Jupiter_KarkRef1993[:,1]=Jupiter_Karkoschka1993[:,3]

#Load plot parameters
Jupiter=PU.PlotSetup("JupiterSpecPlotConfig.txt")
Jupiter.loadplotparams("f:","Jupiter","Spectra")

#Load observations files and create Target+DateTime keys
O=CF.measurement_list(Jupiter.DataFile)
O.load_records(MeasTgt=Target,DateUTSelect=DateUT)
F=CF.ObsFileNames(O.FileList[0])
F.GetFileNames()

#Load spectral bands to measure
Bands=EWU.LinesBands_to_Measure("Jupiter_ObsBands_135mm100lpm.txt")
Bands.load_records(WVRange=[400.,750.])

for time in range (0,len(F.FNArray)):

    #Make key, read raw spectrum
    Key,DateTime=CF.MakeKeyDate(F.FNArray[time])   
    CLR = scipy.fromfile(file=JupPath.input_path+F.FNArray[time], dtype=float, count=-1, sep='\t')    
    CLR=scipy.reshape(CLR,[CLR.size/2,2])
    CLR[:,0]=CLR[:,0]/10.
    NativeDispersion=(CLR[(CLR.size/2-1),0]-CLR[0,0])/(CLR.size/2-1)
    wave,sig=GSU.uniform_wave_grid(CLR[:,0],CLR[:,1],Extend=False)    
    CLRonRef=np.transpose(np.array([wave,sig]))
    CLRonRef[:,1]=pyasl.smooth(CLRonRef[:,1],3,'flat')
    
    #Compute solar, telluric, and planetary equivalent widths from raw spectrum
    EWFN=JupPath.EW_path+Key+"-RawFlux-EW.txt"
    flag=False
    for B in range(0,len(Bands.ID)):
        Temp=EWU.ComputeEW1(CLRonRef,Target,DateTime,Bands.Type[B],Bands.ID[B],
                            Bands.WV0[B],Bands.WV1[B],Bands.WVCont[B],EWFN,flag)
        flag=True

    #Compute top of atmosphere spectrum and albedo
    ToA=GSU.SpectrumMath(CLRonRef,Response,"Divide")
    Albedo=GSU.SpectrumMath(ToA,Ref_g2v,"Divide")
    mAlbedo = np.ma.masked_invalid(Albedo)
    AlbedoNormRangeIndices=np.where((Ref_g2v[:,0] >400.) & \
         (Ref_g2v[:,0] < 750.))   
    NormAlbedo=deepcopy(Albedo)
    NormAlbedo[:,1]=Albedo[:,1]/mAlbedo[AlbedoNormRangeIndices].max()
    np.savetxt(JupPath.One_D_path+Key+"Albedo.txt",NormAlbedo,delimiter=" ",fmt="%10.3F %10.7F")
    
    #Compute planetary equivalent widths from albedo spectrum
    EWFN=JupPath.EW_path+Key+"-Albedo-EW.txt"
    AlbedoBands=EWU.LinesBands_to_Measure("Jupiter_ObsBands_135mm100lpm.txt")
    AlbedoBands.load_records(Type="Planetary",WVRange=[400.,750.])
    print AlbedoBands.ID
    flag=False
    for B in range(0,len(AlbedoBands.ID)):
        #print "B=",B
        Temp=EWU.ComputeEW1(NormAlbedo,Target,DateTime,AlbedoBands.Type[B],AlbedoBands.ID[B],
                            AlbedoBands.WV0[B],AlbedoBands.WV1[B],AlbedoBands.WVCont[B],EWFN,flag)
        flag=True

    #Begin plotting    
    ExposureCLR = 300. #seconds
    Aperture = (0.135/22.)**2. #meters^2    
    Jupiter.Setup_Plot()    
    pl.title(Target+" "+datetime.datetime.strftime(DateTime,'%Y-%m-%d %H:%M:%S')+
                " Spectrum",fontsize=9)

    pl.plot(CLRonRef[:,0],CLRonRef[:,1]/(ExposureCLR*Aperture*NativeDispersion),label='CLR',linewidth=0.5)
    pl.plot(CLRonRef[:,0],ToA[:,1]/(ExposureCLR*Aperture*NativeDispersion),label='Top of Atm.')
    #pl.plot(Ref_g2v[:,0],Ref_g2v[:,1]*1e7,label='Solar Ref. x 1e7')
    pl.plot(NormAlbedo[:,0],NormAlbedo[:,1]*1e7,label='Norm. Albedo x 1e7')
    pl.plot(Jupiter_KarkRef1993[:,0],Jupiter_KarkRef1993[:,1]*1.8e7,label='Karkoschka, 1993 x 1.8e7',linewidth=1,color='0.5')   
    
    pl.legend(loc=0,ncol=4, borderaxespad=0.,prop={'size':6})
    pl.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.15)

    #Save spectrum plot, raw spectra text file, and albedo text file
    pl.savefig(JupPath.One_D_path+Key+"Spectrum.png",dpi=300)        
    np.savetxt(JupPath.One_D_path+Key+"Spectrum.txt",CLRonRef,delimiter=" ",fmt="%10.3F %10.7F")
    np.savetxt(JupPath.One_D_path+Key+"Albedo.txt",NormAlbedo,delimiter=" ",fmt="%10.3F %10.7F")