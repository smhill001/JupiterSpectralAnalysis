# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 07:20:03 2019

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

class Jupiter_Galileo_Profile(CF.readurllines):
    """
    Read observing conditions from the Suomi-net web site
    """
    pass
    def load_records(self):
        self.Time=[]  #Keyword for star identification
        self.Press=[]           #Target, e.g., component of a multiple star
        self.Temp=[]           #Target, e.g., component of a multiple star
        self.Dens=[]
        self.Alt=[]
        self.Grav=[]
        for recordindex in range(0,self.nrecords-1):
            self.Time.append(str(self.URLLines[recordindex][13:20]))
            self.Press.append(np.float(self.URLLines[recordindex][12:18]))
            self.Temp.append(np.float(self.URLLines[recordindex][24:34]))
            self.Dens.append(str(self.URLLines[recordindex][62:66]))                
            self.Alt.append(str(self.URLLines[recordindex][67:72]))                
            self.Grav.append(np.float(self.URLLines[recordindex][60:66]))                




