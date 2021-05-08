#!/usr/bin/env python
# vasc_res_sim.py
import subprocess
import shutil
import os
import glob
import numpy as np
import time
import pdb
import string
import math
import sys
from os import path
import os
import shutil


if __name__=="__main__":
    t0 = time.time()

    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/AngleVariation_3X3Network/'
    CollectedResults = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/AngleVariation_3X3Network/CollectedShearResults/'
    if path.isdir(CollectedResults)==0:
        os.mkdir(CollectedResults)
   
    Density = ['PI_2.2','PI_6' ]   
    Collapse = ['0.0','1.227','2.248', '3.178', '4.17', '5.124', '6.119', '7.08', '8.059', '9.032', '10.0']
    

    NewNumbering = [ '0','1','2', '3','4','5', '6', '7','8','9','10']

    for Level in Density:
        counter = -1
        NewPath = CollectedResults + '/'+Level
        if path.isdir(NewPath)==0:
            os.mkdir(NewPath)
        for i in Collapse:
            counter = counter+1
            print counter
            # Move the file
            ShearStress = TerminalOutputFolder + Level+'/'+i+'/Results/Extracted/surface-traction_3600.vtu'
            NewShearStress = NewPath+'/surface-traction_'+NewNumbering[counter]+'.vtu'
        
            mv =  'cp ' +ShearStress +' '+NewShearStress
            subprocess.call(mv, shell=True)

    Density = ['PI_3','PI_5' ]   
    Collapse = ['0','0.1227','0.2248', '0.3178', '0.417', '0.5124', '0.6119', '0.708', '0.8059', '0.9032', '1.0']
  
    for Level in Density:
        counter = -1
        NewPath = CollectedResults + '/'+Level
        if path.isdir(NewPath)==0:
            os.mkdir(NewPath)
        for i in Collapse:
            counter = counter+1
            print counter
            # Move the file
            ShearStress = TerminalOutputFolder + Level+'/'+i+'/Results/Extracted/surface-traction_1800.vtu'
            NewShearStress = NewPath+'/surface-traction_'+NewNumbering[counter]+'.vtu'
        
            mv =  'cp ' +ShearStress +' '+NewShearStress
            subprocess.call(mv, shell=True)