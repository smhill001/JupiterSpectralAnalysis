# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 07:22:23 2014

@author: steven.hill
"""
def JupiterSpectrumVisualization(Band='619CH4',SpecType="Albedo"):
    #Options for SpecType are:
        #Albedo
        #Spectrum
        #Response
        #RefStar
    import sys
    sys.path.append('f:\\Astronomy\Python Play')
    sys.path.append('f:\\Astronomy\Python Play\Spectrophotometry\Spectroscopy')
    sys.path.append('f:\\Astronomy\Python Play\Util')
    import matplotlib.pyplot as pl
    import pylab
    import numpy as np
    import scipy
    import datetime
    from PyAstronomy import pyasl
    import os, fnmatch
    import GeneralSpecUtils as GSU
    import PlotUtils as PU
    import EquivWidthUtils as EWU
    
    ResponseFile= {'20150123UT':'PolluxResponse20150123UT.txt',
               '20150209UT':'PolluxResponse20150210UT.txt',
               '20150210UT':'PolluxResponse20150210UT.txt',
               '20150318UT':'PolluxResponse20150318UT.txt',
               '20150322UT':'PolluxResponse20150322UT.txt',
               '20150331UT':'PolluxResponse20150322UT.txt'
               }

    ResponsePath="f:/Astronomy/Projects/Planets/Jupiter/Spectral Data/"
    Response20150123UT = scipy.fromfile(file="f:/Astronomy/Projects/Planets/Jupiter/Spectral Data/"+ResponseFile['20150123UT'], dtype=float, count=-1, sep=" ")
    Response20150123UT=scipy.reshape(Response20150123UT,[Response20150123UT.size/2,2])
    Response20150123UT[:,0]=Response20150123UT[:,0]/10.
    Response20150123UT[:,1]=pyasl.smooth(Response20150123UT[:,1],3,'flat')
    
    RefStarFiles= {'20150123UT':'PolluxSpectrum20150123UT.txt',
               '20150209UT':'PolluxSpectrum20150210UT.txt',
               '20150210UT':'PolluxSpectrum20150210UT.txt',
               '20150318UT':'PolluxSpectrum20150318UT.txt',
               '20150322UT':'PolluxSpectrum20150322UT.txt',
               '20150331UT':'PolluxSpectrum20150322UT.txt'
               }

    RefStarPath="f:/Astronomy/Projects/Planets/Jupiter/Spectral Data/"
    RefStar20150123UT = scipy.fromfile(file="f:/Astronomy/Projects/Planets/Jupiter/Spectral Data/"+RefStarFiles['20150123UT'], dtype=float, count=-1, sep=" ")
    RefStar20150123UT=scipy.reshape(RefStar20150123UT,[RefStar20150123UT.size/2,2])
    #RefStar20150123UT[:,0]=RefStar20150123UT[:,0]
    #RefStar20150123UT[:,1]=pyasl.smooth(RefStar20150123UT[:,1],5,'flat')
    RefStar20150210UT = scipy.fromfile(file="f:/Astronomy/Projects/Planets/Jupiter/Spectral Data/"+RefStarFiles['20150210UT'], dtype=float, count=-1, sep=" ")
    RefStar20150210UT=scipy.reshape(RefStar20150210UT,[RefStar20150210UT.size/2,2])
    RefPath="f:/Astronomy/Python Play/SPLibraries/SpectralReferenceFiles/ReferenceLibrary/"
    K0III=scipy.loadtxt(RefPath+"k0iii.dat", dtype=float, usecols=(0,1))
    G2V=scipy.loadtxt(RefPath+"g2v.dat", dtype=float, usecols=(0,1))

    #Should include Jupiter Red and other Spectra here from 20130117UT
    #Maybe include 20130109 spectra - might require some extra work, and they're low dispersion!
    #pl.figure(figsize=(6.5, 1.5), dpi=150, facecolor="white")
    pl.figure(figsize=(8.0, 4.0), dpi=150, facecolor="white")
    
    pl.subplot(1, 1, 1)
    #Plot Layout Configuration
    
    if SpecType=='Albedo':
        PlotLims={'All':{'x0':400,'x1':1000,'xtks':25,'y0':0.00,'y1':0.60,'ytks':13},
                  '619CH4':{'x0':600,'x1':640,'xtks':21,'y0':0.40,'y1':0.55,'ytks':16},
                  '725CH4':{'x0':700,'x1':750,'xtks':11,'y0':0.15,'y1':0.55,'ytks':9},
                  '889CH4':{'x0':860,'x1':920,'xtks':13,'y0':0.00,'y1':0.55,'ytks':12}
                  }
    elif SpecType=='Spectrum':
        PlotLims={'All':{'x0':400,'x1':1000,'xtks':25,'y0':0.00,'y1':1.0e6,'ytks':13},
                  '619CH4':{'x0':600,'x1':640,'xtks':21,'y0':2.0e5,'y1':8.0e5,'ytks':13},
                  '725CH4':{'x0':700,'x1':750,'xtks':11,'y0':0.00,'y1':1.0e6,'ytks':9},
                  '889CH4':{'x0':860,'x1':920,'xtks':13,'y0':0.00,'y1':1.0e6,'ytks':12}
                  }
    elif SpecType=='Response':
        PlotLims={'All':{'x0':400,'x1':1000,'xtks':25,'y0':0.00,'y1':1.0,'ytks':13},
                  '619CH4':{'x0':600,'x1':640,'xtks':21,'y0':0.60,'y1':0.9,'ytks':16},
                  '725CH4':{'x0':700,'x1':750,'xtks':11,'y0':0.00,'y1':1.0e6,'ytks':9},
                  '889CH4':{'x0':860,'x1':920,'xtks':13,'y0':0.00,'y1':1.0e6,'ytks':12}
                  }
    elif SpecType=='RefStar':
        PlotLims={'All':{'x0':400,'x1':1000,'xtks':25,'y0':0.00,'y1':1.2e7,'ytks':13},
                  '619CH4':{'x0':600,'x1':640,'xtks':21,'y0':7.0e6,'y1':1.0e7,'ytks':16},
                  '725CH4':{'x0':700,'x1':750,'xtks':11,'y0':0.00,'y1':1.0e6,'ytks':9},
                  '889CH4':{'x0':860,'x1':920,'xtks':13,'y0':0.00,'y1':1.0e6,'ytks':12}
                  }

    if SpecType=='Albedo' or SpecType=='Spectrum':
        Jupiter_1996UT,Jupiter_1996Smooth,NIR20150123specarray,CLR20150123specarray,\
           NIR20150209specarray,CLR20150209specarray,\
           NIR20150210specarray,CLR20150210specarray,\
           Jupiter_KarkRef1993=GetDailyAverageSpectraJupiter(SpecType)
           
        #print Jupiter_1996UT,NIR20150209specarray

    
    # Set x limits
    
    pl.xlim(PlotLims[Band]['x0'],PlotLims[Band]['x1'])
    pl.xticks(np.linspace(PlotLims[Band]['x0'],PlotLims[Band]['x1'],
                          PlotLims[Band]['xtks'], endpoint=True))
    # Set y limits
    pl.ylim(PlotLims[Band]['y0'],PlotLims[Band]['y1'])
    pl.yticks(np.linspace(PlotLims[Band]['y0'],PlotLims[Band]['y1'],
                          PlotLims[Band]['ytks'], endpoint=True))
    pl.grid(linewidth=0.2)
    pl.tick_params(axis='both', which='major', labelsize=8)
    pl.ylabel(SpecType,fontsize=8,color="black")
    pl.xlabel(r"$Wavelength (nm)$",fontsize=8)
    
    pl.title('Jupiter '+Band+' '+SpecType,fontsize=9)
    
    if SpecType=='Albedo' or SpecType=='Spectrum':
        pl.step(Jupiter_KarkRef1993[:,0]/10.,Jupiter_KarkRef1993[:,1],label='Karkoschka, 1994',linewidth=1,color='k',where='mid')
        
        PU.Draw_with_Conf_Level(NIR20150123specarray.MeanSpec,0.55,'C1','NIR0123',step=True)
        PU.Draw_with_Conf_Level(CLR20150123specarray.MeanSpec,0.522,'C2','CLR0123',step=True)
        PU.Draw_with_Conf_Level(NIR20150209specarray.MeanSpec,0.55,'C3','NIR0209',step=True)
        PU.Draw_with_Conf_Level(CLR20150209specarray.MeanSpec,0.522,'C4','CLR0209',step=True)
        PU.Draw_with_Conf_Level(NIR20150210specarray.MeanSpec,0.55,'C5','NIR0210',step=True)
        PU.Draw_with_Conf_Level(CLR20150210specarray.MeanSpec,0.522,'C6','CLR0210',step=True)
        pl.step((Jupiter_1996UT[:,0]-7.)/10.,Jupiter_1996Smooth*0.55,label='Jupiter_1996UT-Smth',color='g',linewidth=1.0,where='mid')
    
    if SpecType=='Response':
        pl.step(Response20150123UT[:,0],Response20150123UT[:,1],c='C7',label='Response20150123UT')
        TempRefK0III=K0III
        TempRefK0III[:,0]=TempRefK0III[:,0]/10.
        CalcResponse=GSU.SpectrumMath(RefStar20150123UT,TempRefK0III,'Divide')
        WVRangeIndices=np.where((CalcResponse[:,0] >450.) & \
                 (CalcResponse[:,0] < 650.))   
        CalcResponse[:,1]=CalcResponse[:,1]/CalcResponse[WVRangeIndices,1].max()
        pl.step(CalcResponse[:,0],CalcResponse[:,1],c='C6',label='CalcResponse')
        pl.step(CalcResponse[:,0],pyasl.smooth(CalcResponse[:,1],5,'flat'),c='C5',label='CalcResponseSmooth')

    if SpecType=='RefStar':
        print RefStar20150123UT
        pl.step(RefStar20150123UT[:,0],pyasl.smooth(RefStar20150123UT[:,1],5,'flat'),c='C7',label='RefStar20150123UT')
        pl.step(RefStar20150210UT[:,0],pyasl.smooth(RefStar20150210UT[:,1]*0.85,5,'flat'),c='C6',label='RefStar20150210UT')
        pl.step(K0III[:,0]/10.,pyasl.smooth(K0III[:,1]*9e6,5,'flat'),c='C0',label='K0III')
        pl.step(G2V[:,0]/10.,pyasl.smooth(G2V[:,1]*9e6,5,'flat'),c='C1',label='G2V')

    
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
    
    pl.legend(loc=3,ncol=2, borderaxespad=0.,prop={'size':7})
    pl.subplots_adjust(left=0.08, bottom=0.12, right=0.98, top=0.92,
                    wspace=None, hspace=None)
                    
    pylab.savefig('JupiterDayAverages_'+Band+'_'+SpecType+'.png',dpi=300)
    
    
    ###############################################################################
    
    Bands=EWU.LinesBands_to_Measure("Jupiter_ObsBands_135mm100lpm.txt")
    Bands.load_records(Type="Planetary",WVRange=[600.,650.])
        
    EWFN="Test-EW.txt"
    flagA=False
    DTarray=['20150113T000000','20150209T000000','20150210T000000']
    
    for B in range(0,len(Bands.ID)):
        print B,Bands.ID[B]
        for D in range(0,len(DTarray)):
            DT=DTarray[D]
            DateTime=datetime.datetime.strptime(DT[0:4]+"-"+DT[4:6]+"-" \
                +DT[6:8]+"T"+DT[9:11]+":"+DT[11:13]+":"+DT[13:15], \
                '%Y-%m-%dT%H:%M:%S')
            if D==0:
                TempSpec=CLR20150123specarray
            elif D==1:
                TempSpec=CLR20150209specarray
            elif D==2:
                TempSpec=CLR20150210specarray
            Temp=EWU.ComputeEW1(TempSpec.MeanSpec,
                                'Jupiter',DateTime,'Planetary','619CH4',
                                605.5,626.5,2.0,EWFN,flagA)
            #                    Bands.WV0[B],Bands.WV1[B],Bands.WVCont[B],EWFN,flagA)
            flagA=True
def GetDailyAverageSpectraJupiter(SpecType):
    import sys
    sys.path.append('f:\\Astronomy\Python Play')
    sys.path.append('f:\\Astronomy\Python Play\Spectrophotometry\Spectroscopy')
    sys.path.append('f:\\Astronomy\Python Play\Util')
    import numpy as np
    import scipy
    from PyAstronomy import pyasl
    import os, fnmatch
    import GeneralSpecUtils as GSU
    
    Jupiter_1996UT = scipy.fromfile(file="f:/Astronomy/Projects/Planets/Saturn/PROJECT/JUPS.DAT", dtype=float, count=-1, sep='\t')    
    Jupiter_1996UT=scipy.reshape(Jupiter_1996UT,[Jupiter_1996UT.size/2,2])
    Jupiter_1996UT[:,1]=Jupiter_1996UT[:,1]/Jupiter_1996UT[:,1].max() #Normalize the data
    JupiterSlopeCorrection=np.linspace(1.02,0.93,Jupiter_1996UT.size/2)
    Jupiter_1996UT[:,1]=Jupiter_1996UT[:,1]*JupiterSlopeCorrection
    Jupiter_1996Smooth=pyasl.smooth(Jupiter_1996UT[:,1],9,'flat')
    
    NIR20150123Files=[]
    listOfFiles = os.listdir('../1D Spectra/')
    pattern='Jupiter20150123*100lpm-742NIR_'+SpecType+'.txt'
    for entry in listOfFiles:
        if fnmatch.fnmatch(entry, pattern):
            NIR20150123Files.append("../1D Spectra/"+entry)
    NIR20150123specarray=GSU.SpectrumAggregation('f:',NIR20150123Files,FileList=True)
    NIR20150123specarray.ComputeAverageandStats()    
    
    CLR20150123Files=[]
    listOfFiles = os.listdir('../1D Spectra/')
    pattern='Jupiter20150123*100lpm-550CLR_'+SpecType+'.txt'
    for entry in listOfFiles:
        if fnmatch.fnmatch(entry, pattern):
            CLR20150123Files.append("../1D Spectra/"+entry)
    CLR20150123specarray=GSU.SpectrumAggregation('f:',CLR20150123Files,FileList=True)
    CLR20150123specarray.ComputeAverageandStats()    
    
    
    NIR20150209Files=[]
    listOfFiles = os.listdir('../1D Spectra/')
    pattern='Jupiter20150209*100lpm-742NIR_'+SpecType+'.txt'
    for entry in listOfFiles:
        if fnmatch.fnmatch(entry, pattern):
            NIR20150209Files.append("../1D Spectra/"+entry)
    NIR20150209specarray=GSU.SpectrumAggregation('f:',NIR20150209Files,FileList=True)
    NIR20150209specarray.ComputeAverageandStats()    
    
    CLR20150209Files=[]
    listOfFiles = os.listdir('../1D Spectra/')
    pattern='Jupiter20150209*100lpm-550CLR_'+SpecType+'.txt'
    for entry in listOfFiles:
        if fnmatch.fnmatch(entry, pattern):
            CLR20150209Files.append("../1D Spectra/"+entry)
            print entry
    CLR20150209specarray=GSU.SpectrumAggregation('f:',CLR20150209Files,FileList=True)
    CLR20150209specarray.ComputeAverageandStats()    
    
    
    NIR20150210Files=[]
    listOfFiles = os.listdir('../1D Spectra/')
    pattern='Jupiter20150210*100lpm-742NIR_'+SpecType+'.txt'
    for entry in listOfFiles:
        if fnmatch.fnmatch(entry, pattern):
            NIR20150210Files.append("../1D Spectra/"+entry)
    NIR20150210specarray=GSU.SpectrumAggregation('f:',NIR20150210Files,FileList=True)
    NIR20150210specarray.ComputeAverageandStats()    
    
    CLR20150210Files=[]
    listOfFiles = os.listdir('../1D Spectra/')
    pattern='Jupiter20150210*100lpm-550CLR_'+SpecType+'.txt'
    for entry in listOfFiles:
        if fnmatch.fnmatch(entry, pattern):
            CLR20150210Files.append("../1D Spectra/"+entry)
    CLR20150210specarray=GSU.SpectrumAggregation('f:',CLR20150210Files,FileList=True)
    CLR20150210specarray.ComputeAverageandStats()    
    
    Jupiter_Karkoschka1993 = scipy.fromfile(file="F:/Astronomy/Projects/Planets/Saturn/Spectral Data/Karkoschka/1993.tab.txt", dtype=float, count=-1, sep=" ")    
    Jupiter_Karkoschka1993=scipy.reshape(Jupiter_Karkoschka1993,[Jupiter_Karkoschka1993.size/8,8])
    
    Jupiter_KarkRef1993=np.zeros((Jupiter_Karkoschka1993.size/8,2))
    Jupiter_KarkRef1993[:,0]=Jupiter_Karkoschka1993[:,0]*10.
    Jupiter_KarkRef1993[:,1]=Jupiter_Karkoschka1993[:,3]
    print "***Jupiter_KarkRef1993***", Jupiter_KarkRef1993
    
    return Jupiter_1996UT,Jupiter_1996Smooth,NIR20150123specarray,CLR20150123specarray,\
           NIR20150209specarray,CLR20150209specarray,\
           NIR20150210specarray,CLR20150210specarray,\
           Jupiter_KarkRef1993
        