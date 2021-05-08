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

    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/PressureVariation/'
    Pressure = ["200","50"]
    Collapse = ['0.0','10.0','20.0','30.0','40.0','50.0','60.0','70.0','80.0','90.0','100.0']
            
    Parallel = 20
    SleepyTime = 200
    NewNumbering = [ '0','1','2', '3','4','5', '6', '7','8','9','10']
  
    
    NewPath = TerminalOutputFolder+'CollectedShearResults/'
    if path.isdir(NewPath)==0:
        os.mkdir(NewPath)
    
    for j in Pressure:
        PressurePath = NewPath+"Pressure_"+j
        if path.isdir(PressurePath)==0:
            os.mkdir(PressurePath)
        counter = -1
        for i in Collapse:
            counter = counter +1
            # Move the file
            ShearStress = TerminalOutputFolder +'Pressure_'+j+'/' + i+'/Results/Extracted/surface-traction_3600.vtu'
            NewShearStress = PressurePath+'/surface-traction_'+ NewNumbering[counter] +'.vtu'
        
            mv =  'cp ' +ShearStress +' '+NewShearStress
            subprocess.call(mv, shell=True)

# ---------

