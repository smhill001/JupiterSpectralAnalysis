# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 07:22:23 2014

@author: steven.hill
"""
import sys
sys.path.append('f:\\Astronomy\Python Play')
sys.path.append('f:\\Astronomy\Python Play\Spectrophotometry\Spectroscopy')
sys.path.append('f:\\Astronomy\Python Play\Util')
import matplotlib.pyplot as pl
import pylab
import numpy as np
import scipy
from scipy import interpolate
from PyAstronomy import pyasl
import os, fnmatch
import GeneralSpecUtils as GSU
import PlotUtils as PU


Jupiter_1996UT = scipy.fromfile(file="f:/Astronomy/Projects/Planets/Saturn/PROJECT/JUPS.DAT", dtype=float, count=-1, sep='\t')    
Jupiter_1996UT=scipy.reshape(Jupiter_1996UT,[Jupiter_1996UT.size/2,2])
Jupiter_1996UT[:,1]=Jupiter_1996UT[:,1]/Jupiter_1996UT[:,1].max() #Normalize the data
JupiterSlopeCorrection=np.linspace(1.02,0.93,Jupiter_1996UT.size/2)
Jupiter_1996UT[:,1]=Jupiter_1996UT[:,1]*JupiterSlopeCorrection
Jupiter_1996Smooth=pyasl.smooth(Jupiter_1996UT[:,1],9,'flat')

NIRFiles=[]
listOfFiles = os.listdir('../1D Spectra/')
pattern='*100lpm-742NIR_Albedo.txt'
for entry in listOfFiles:
    if fnmatch.fnmatch(entry, pattern):
        NIRFiles.append("../1D Spectra/"+entry)
NIRspecarray=GSU.SpectrumAggregation('f:',NIRFiles,FileList=True)
NIRspecarray.ComputeAverageandStats()    

CLRFiles=[]
listOfFiles = os.listdir('../1D Spectra/')
pattern='*100lpm-550CLR_Albedo.txt'
for entry in listOfFiles:
    if fnmatch.fnmatch(entry, pattern):
        CLRFiles.append("../1D Spectra/"+entry)
CLRspecarray=GSU.SpectrumAggregation('f:',CLRFiles,FileList=True)
CLRspecarray.ComputeAverageandStats()    


KeyArray=["Jupiter20150209034346UT",
      "Jupiter20150209035901UT",
      "Jupiter20150209041414UT",
      "Jupiter20150209042927UT",
      "Jupiter20150209044439UT",
      "Jupiter20150209045954UT",
      "Jupiter20150209051507UT",
      "Jupiter20150209053020UT",
      "Jupiter20150209054448UT",
      "Jupiter20150209060047UT",
      "Jupiter20150209061559UT",
      "Jupiter20150210044442UT",
      "Jupiter20150210045951UT",
      "Jupiter20150210051502UT",
      #"Jupiter20150210053012UT",  Dropped due to anomalous EW for 889nm
      "Jupiter20150210054437UT",
      "Jupiter20150210060033UT",
      "Jupiter20150210061543UT",
      "Jupiter20150210063044UT",
      "Jupiter20150210064606UT",
      "Jupiter20150210070116UT",
      "Jupiter20150210071627UT",
      "Jupiter20150210073137UT",
      "Jupiter20150210074647UT",
      "Jupiter20150210080158UT",
      "Jupiter20150210081709UT",
      "Jupiter20150210083221UT",
      "Jupiter20150210084732UT",
      "Jupiter20150210090243UT",
      "Jupiter20150210092047UT"]

# Read and reshape spectral data files    
First=True

for keyindex in range(0,len(KeyArray)):
    print keyindex, KeyArray[keyindex]
#Jupiter_20150123UT = scipy.fromfile(file="JupiterSpectrum20150123UT.txt", dtype=float, count=-1, sep='\t')    
    TempInput = scipy.fromfile(file="../"+KeyArray[keyindex]+"Albedo.txt", dtype=float, count=-1, sep='\t')    
    TempInput=np.nan_to_num(TempInput)
    TempReshape=scipy.reshape(TempInput,[TempInput.size/2,2])
    
    print First,TempReshape[:,1].mean()
    if First==True:
        TempSum=TempReshape
        First=False
    elif First==False:
        TempSum[:,1]=TempSum[:,1]+TempReshape[:,1]
        print First,TempSum[1200,1]

TempSum[:,1]=TempSum[:,1]/len(KeyArray)


Jupiter_Karkoschka1993 = scipy.fromfile(file="F:/Astronomy/Projects/Planets/Saturn/Spectral Data/Karkoschka/1993.tab.txt", dtype=float, count=-1, sep=" ")    
Jupiter_Karkoschka1993=scipy.reshape(Jupiter_Karkoschka1993,[Jupiter_Karkoschka1993.size/8,8])

Jupiter_KarkRef1993=np.zeros((Jupiter_Karkoschka1993.size/8,2))
Jupiter_KarkRef1993[:,0]=Jupiter_Karkoschka1993[:,0]*10.
Jupiter_KarkRef1993[:,1]=Jupiter_Karkoschka1993[:,3]
print "***Jupiter_KarkRef1993***", Jupiter_KarkRef1993

#Should include Jupiter Red and other Spectra here from 20130117UT
#Maybe include 20130109 spectra - might require some extra work, and they're low dispersion!
#pl.figure(figsize=(6.5, 1.5), dpi=150, facecolor="white")
pl.figure(figsize=(8.0, 4.0), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0=400
x1=950
xtks=23
y0=0.0
y1=0.6
ytks=13

# Set x limits
pl.xlim(x0,x1)
# Set x ticks
pl.xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
pl.ylim(y0,y1)
pl.yticks(np.linspace(y0,y1,ytks, endpoint=True))
# Set y ticks
pl.grid(linewidth=0.2)
pl.tick_params(axis='both', which='major', labelsize=8)
pl.ylabel(r"$Albedo$",fontsize=8,color="black")
pl.xlabel(r"$Wavelength (nm)$",fontsize=8)



"""
pl.title("Jupiter Albedo - All Data",fontsize=9)
pl.plot(JupiterINT_20130612UT[:,0]/10.,JupiterINT_20130612UT[:,1]*0.60,label='INT_20130612UT',linewidth=1)
pl.plot(JupiterDSK_20130612UT[:,0]/10.,JupiterDSK_20130612UT[:,1]*0.60,label='DSK_20130612UT',linewidth=0.5)
pl.plot(JupiterRNG_20130612UT[:,0]/10.,JupiterRNG_20130612UT[:,1]*0.60,label='RNG_20130612UT',linewidth=0.5)
pl.plot(JupiterINT_20140629UT[:,0]/10.,JupiterINT_20140629UT[:,1]*0.60,label='INT_20140629UT',linewidth=1,color='k')
pl.plot(Jupiter2014Smth[:,0]/10.,Jupiter2014Smth[:,1]*0.59,label='2014 Smooth',linewidth=1,color='r')

pl.plot(Jupiter_KarkRef1993[:,0]/10.,Jupiter_KarkRef1993[:,1],label='Karkoschka, 1993',linewidth=1,color='0.5')
"""
Indices=np.where((TempSum[:,0] < 4000.))
TempSum[Indices,1]=np.nan
pl.title("Jupiter",fontsize=9)
#pl.plot(JupiterINT_20130612UT[:,0]/10.,JupiterINT_20130612UT[:,1]*0.60,label='INT_20130612UT',linewidth=1)
#pl.plot(JupiterDSK_20130612UT[:,0]/10.,JupiterDSK_20130612UT[:,1]*0.60,label='DSK_20130612UT',linewidth=0.5)
#pl.plot(JupiterRNG_20130612UT[:,0]/10.,JupiterRNG_20130612UT[:,1]*0.60,label='RNG_20130612UT',linewidth=0.5)
pl.step(TempSum[:,0]/10.,TempSum[:,1]*0.54,label='This Work (2014) 0.006m',linewidth=1,where='mid')
#pl.plot(Jupiter2014Smth[:,0]/10.,Jupiter2014Smth[:,1]*0.59,label='2014 Smooth',linewidth=1,color='g')
pl.scatter(Jupiter_1996UT[:,0]/10.,Jupiter_1996UT[:,1]*0.55,label='Jupiter_1996UT',marker='.',s=0.1,color='g')
pl.step((Jupiter_1996UT[:,0]-7.)/10.,Jupiter_1996Smooth*0.55,label='Jupiter_1996UT-Smth',color='g',linewidth=1.0,where='mid')

pl.step(Jupiter_KarkRef1993[:,0]/10.,Jupiter_KarkRef1993[:,1],label='Karkoschka, 1994',linewidth=1,color='0.5',where='mid')

print NIRspecarray.MeanSpec.shape
PU.Draw_with_Conf_Level(NIRspecarray.MeanSpec,0.55,'r','Test',step=True)
PU.Draw_with_Conf_Level(CLRspecarray.MeanSpec,0.55,'b','Test',step=True)
#pl.step(NIRspecarray.MeanSpec[:,0],NIRspecarray.MeanSpec[:,1]*1.15)
#pl.step(CLRspecarray.MeanSpec[:,0],CLRspecarray.MeanSpec[:,1]*0.55)

LineWVs=np.array([486.0,543.0,576.0,  #H I Balmer
                  597.0,619.0,668.0,683.0,705.0,725.0,
                  790.0,842.0,862.0,889.0,
                  
                  551.0,646.0,750.0,760.0,825.0,930.0])                    #H I Balmer
                  
                  
LineY=np.array([0.3,0.3,0.3,
                0.3,0.3,0.3,0.3,0.3,0.05,
                0.2,0.2,0.05,0.1,
                
                0.3,0.3,0.3,0.3,0.3,0.2])
                                        #H I Paschen
LineLabels=[r'$CH_4$',r'$CH_4$',r'$CH_4$',
            r'$CH_4$',r'$CH_4$',r'$CH_4$',r'$CH_4$',r'$CH_4$',r'$CH_4$',
            r'$CH_4$',r'$CH_4$',r'$CH_4$',r'$CH_4$',
            
            r'$NH_3$',r'$NH_3$',r'$NH_3$',r'$NH_3$',r'$NH_3$',r'$NH_3$',r'$NH_3$']

for l in range(0,LineWVs.size):                
    pl.text(LineWVs[l],LineY[l],LineWVs[l].astype('|S5')+' '+LineLabels[l],fontsize=8,
            verticalalignment='bottom',horizontalalignment='center',
            rotation='vertical')

pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':7})
pl.subplots_adjust(left=0.08, bottom=0.12, right=0.98, top=0.92,
                wspace=None, hspace=None)
                
pylab.savefig('JupiterAlbedoAllYears.png',dpi=300)