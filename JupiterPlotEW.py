# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 12:48:36 2019

@author: Steven Hill

Make this callable with Target, DateUT, and Band ID


"""
def JupiterPlotEW(DateUT,BandID,Grating,FluxType='Albedo'):
    #JupiterPlotEW("20150123UT","619CH4","100lpm-550CLR")
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
    import PlotUtils as PU
    import ConfigFiles as CF
    from scipy import interpolate
    
    TimePlot=PU.PlotSetup("JupiterEWSpecPlotConfig.txt")
    TimePlot.loadplotparams("f:","Jupiter_"+BandID,"Time")
    AirmassPlot=PU.PlotSetup("JupiterEWSpecPlotConfig.txt")
    AirmassPlot.loadplotparams("f:","Jupiter_"+BandID,"Airmass")
    CMIIPlot=PU.PlotSetup("JupiterEWSpecPlotConfig.txt")
    CMIIPlot.loadplotparams("f:","Jupiter_"+BandID,"CMII")
    PWVPlot=PU.PlotSetup("JupiterEWSpecPlotConfig.txt")
    PWVPlot.loadplotparams("f:","Earth_PWV","Time")
    PressPlot=PU.PlotSetup("JupiterEWSpecPlotConfig.txt")
    PressPlot.loadplotparams("f:","Earth_Press","Time")
    TempFPlot=PU.PlotSetup("JupiterEWSpecPlotConfig.txt")
    TempFPlot.loadplotparams("f:","Earth_TempF","Time")
    Earth_Airmass=PU.PlotSetup("JupiterEWSpecPlotConfig.txt")
    Earth_Airmass.loadplotparams("f:","Earth_Airmass","Time")

    
    CalibrationStars={'20150123UT':"Pollux",
                      '20150209UT':"Pollux",
                      '20150210UT':"Pollux",
                      '20150318UT':"Pollux",
                      '20150322UT':"Pollux",
                      '20150331UT':"Pollux"}
    CalibrationDates={'20150123UT':'2015-01-23 05:54:00',
                      '20150209UT':'2015-02-10 05:54:00',
                      '20150210UT':'2015-02-10 05:54:00',
                      '20150318UT':'2015-03-18 05:54:00',
                      '20150322UT':'2015-03-22 05:54:00',
                      '20150331UT':'2015-03-22 05:54:00'}
    ResponseFile= {'20150123UT':'PolluxResponse20150123UT.txt',
                   '20150209UT':'PolluxResponse20150210UT.txt',
                   '20150210UT':'PolluxResponse20150210UT.txt',
                   '20150318UT':'PolluxResponse20150318UT.txt',
                   '20150322UT':'PolluxResponse20150322UT.txt',
                   '20150331UT':'PolluxResponse20150322UT.txt'
                   }

    Pollux = SkyCoord.from_name('Pollux')
    observer = ephem.Observer()
    #Location from Google Maps at 483 S Oneida Way, Denver 80224
    observer.lon = ephem.degrees('-104.907985')
    observer.lat = ephem.degrees('39.708200')
    print float(observer.lat)*180/np.pi
    Oneida = EarthLocation(lat=float(observer.lat)*u.rad, 
                           lon=float(observer.lon)*u.rad, height=1500*u.m)
    utcoffset = 0*u.hour  # Universal Time
    time = Time(CalibrationDates[DateUT]) - utcoffset

    Polluxaltaz = Pollux.transform_to(AltAz(obstime=time,location=Oneida))
    print Polluxaltaz.alt.deg
    print("Pollux's Altitude = {0.alt:.2}".format(Polluxaltaz))
    print Polluxaltaz.secz
    
    AEWs=EWU.EWObservations("../EWs/"+DateUT+"-"+Grating+"-Albedo-EW.txt")
    AEWs.load_records("Jupiter",BandID)
    REWs=EWU.EWObservations("../EWs/"+DateUT+"-"+Grating+"-RawFlux-EW.txt")
    REWs.load_records("Jupiter",BandID)
    #print EWs.DateTimeUTObs
    #print EWs.EW
    ########################################################
    TimePlot.Setup_PlotDate(min(AEWs.DateTimeUTObs),max(AEWs.DateTimeUTObs),1,
                            canvas_size=[8.0,6.0],subplot=[4,1,1])    
    lbl=BandID+" REW "+str(round(np.mean(REWs.EW),3))+"$\pm$"+str(round(np.std(REWs.EW),3))
    pl.plot_date(REWs.DateTimeUTObs,REWs.EW,label=lbl)
    lbl=BandID+" AEW "+str(round(np.mean(AEWs.EW),3))+"$\pm$"+str(round(np.std(AEWs.EW),3))
    pl.plot_date(AEWs.DateTimeUTObs,AEWs.EW,label=lbl)
    pl.ylabel("Equiv. Width (nm)",fontsize=7)
    pl.legend(loc=2,ncol=3,fontsize=6,scatterpoints=1)
    pl.title("Jupiter "+DateUT+" "+BandID)
    print "Average AEW, Std. Dev., 95% Conf.= ",np.mean(AEWs.EW),np.std(AEWs.EW),\
            1.96*np.std(AEWs.EW)/np.sqrt(len(AEWs.EW))    
    print "Average REW, Std. Dev., 95% Conf.= ",np.mean(REWs.EW),np.std(REWs.EW),\
            1.96*np.std(REWs.EW)/np.sqrt(len(REWs.EW))    
    ########################################################   
    print PWVPlot.DataFile
    Earth_atmo=CF.Observing_Conditions(PWVPlot.DataFile)
    print DateUT[0:4]+'-'+DateUT[4:6]+'-'+DateUT[6:8]
    Earth_atmo.load_records(DateUT[0:4]+'-'+DateUT[4:6]+'-'+DateUT[6:8])
    AX=PWVPlot.Setup_PlotDate(min(AEWs.DateTimeUTObs),max(AEWs.DateTimeUTObs),1,
                            new_canvas=False,subplot=[4,1,2])    

    DTarray=[]
    DT1array=[]
    for DT in Earth_atmo.ObsDateUT:
        TmpDT=Time(DT[0:19], format='isot', scale='utc')
        DTarray.append(float(TmpDT.jd))
    for DT1 in AEWs.DateTimeUTObs:
        TmpDT1=Time(DT1, format='iso', scale='utc')
        DT1array.append(float(TmpDT1.jd))
    
    InterpPWV=interpolate.interp1d(DTarray,Earth_atmo.PWV,kind='linear', 
                                copy=True,bounds_error=False, 
                                fill_value=np.NaN,axis=0)  
    PWVonGrid=InterpPWV(np.array(DT1array))
    AX.plot_date(AEWs.DateTimeUTObs,PWVonGrid,label="PWV [mm]")
    AX.set_ylabel("PWV (mm)",fontsize=7)
    ########################################################
    InterpPress=interpolate.interp1d(DTarray,Earth_atmo.Press,kind='linear', 
                                copy=True,bounds_error=False, 
                                fill_value=np.NaN,axis=0)  
    PressonGrid=InterpPress(np.array(DT1array))
    
    ax1=AX.twinx()
    #ax1.set_xlim(x0,x1)
    #ax1.set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
    ax1.set_ylim(800.,850.)
    ax1.set_yticks(np.linspace(800.,850.,11, endpoint=True))
    ax1.tick_params(axis='y', which='major', labelsize=7)
    ax1.set_ylabel("Pressure (mb)",fontsize=7)

    
    #PressPlot.Setup_PlotDate(min(EWs.DateTimeUTObs),max(EWs.DateTimeUTObs),1)
    ax1.plot_date(AEWs.DateTimeUTObs,PressonGrid,'C1',markerfacecolor='C1',
                  marker="o",linewidth=0.0,label="Pressure")
    AX.legend(loc=2,ncol=3,fontsize=6,scatterpoints=1)
    ax1.legend(loc=1,ncol=3,fontsize=6,scatterpoints=1)
    
    ########################################################
    InterpTempC=interpolate.interp1d(DTarray,Earth_atmo.TempC,kind='linear', 
                                copy=True,bounds_error=False, 
                                fill_value=np.NaN,axis=0)  
    TempConGrid=InterpTempC(np.array(DT1array))

    InterpRelHum=interpolate.interp1d(DTarray,Earth_atmo.RelHum,kind='linear', 
                                copy=True,bounds_error=False, 
                                fill_value=np.NaN,axis=0)  
    RelHumonGrid=InterpRelHum(np.array(DT1array))
        
    az=TempFPlot.Setup_PlotDate(min(AEWs.DateTimeUTObs),max(AEWs.DateTimeUTObs),1,
                            new_canvas=False,subplot=[4,1,3])   
    
    az.plot_date(AEWs.DateTimeUTObs,TempConGrid*9.0/5.0+32.0,label="Temp.[F]")
    az.plot_date(AEWs.DateTimeUTObs,RelHumonGrid,label="Rel. Hum. [%]")
    pl.legend(loc=2,ncol=3,fontsize=6,scatterpoints=1)
    ########################################################
    airmass=[]
    CMII=[]
    for d in AEWs.DateTimeUTObs:
        observer.date=d       
        Planet = ephem.Jupiter(d)
        Planet.compute(observer)
        airmass.extend([pyasl.airmassPP(90.-Planet.alt*180./np.pi)])
        CMII.append(Planet.cmlII*180./np.pi)
    ########################################################  
    
    EAM=Earth_Airmass.Setup_PlotDate(min(AEWs.DateTimeUTObs),max(AEWs.DateTimeUTObs),1,
                            new_canvas=False,subplot=[4,1,4])    
    #lbl=BandID+" AEW "+str(round(np.mean(AEWs.EW),3))+"$\pm$"+str(round(np.std(AEWs.EW),3))
    EAM.plot_date(AEWs.DateTimeUTObs,np.array(airmass),label=lbl)
    ref_airmass=np.full(len(airmass),Polluxaltaz.secz)
    EAM.plot_date(AEWs.DateTimeUTObs,ref_airmass,
                  linestyle='solid', marker='None',
                  color='C0')
    #PLOT HORIZONTAL LINE WITH REFERENCE AIR MASS HERE
    EAM.set_ylabel("Air Mass",fontsize=7)
    EAM.text(sorted(AEWs.DateTimeUTObs)[0],2.9,CalibrationStars[DateUT],
             color='#000000',fontsize=8,
             verticalalignment='top',horizontalalignment='left',)
    EAM.text(sorted(AEWs.DateTimeUTObs)[0],2.7,CalibrationDates[DateUT],
             color='#000000',fontsize=8,
             verticalalignment='top',horizontalalignment='left',)
    EAM.text(sorted(AEWs.DateTimeUTObs)[0],2.5,"Air Mass="+str(round(Polluxaltaz.secz,3)),
             color='#000000',fontsize=8,
             verticalalignment='top',horizontalalignment='left',)
    print ' '
    print sorted(AEWs.DateTimeUTObs)[0]
    pl.subplots_adjust(left=0.08, right=0.92, top=0.95, bottom=0.05)
    pl.savefig("Jupiter EW Plot"+DateUT+" "+BandID+".png",dpi=300)
    

    ########################################################
    ########################################################  
    AirmassPlot.Setup_Plot()
    pl.scatter(np.array(airmass),AEWs.EW)
    Coefs=np.polyfit(np.array(airmass), AEWs.EW, 1)
    EWFit=Coefs[1]+Coefs[0]*np.array(airmass)
    pl.plot(np.array(airmass),EWFit)
    corr=np.corrcoef([np.array(airmass),AEWs.EW])
    print "Coefs= ",Coefs
    print "Correlation Coef.= ",corr[0,1]
    #CoefsUT=np.polyfit(np.array(EWs.DateTimeUTObs), EWs.EW, 1)
    ########################################################
    CMIIPlot.Setup_Plot()
    pl.scatter(np.array(CMII),AEWs.EW)
    pl.text(10.,20.,'GOES-16/SUVI  Fe195  YYYY-MM-DD HH:MM:SS', color='#000000',fontsize=8,
       verticalalignment='bottom',horizontalalignment='left',)
    ########################################################
    