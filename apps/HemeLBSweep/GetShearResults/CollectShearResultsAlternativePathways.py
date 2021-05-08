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

    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/IncreasingAlternativePathways/Levels_'
    CollectedResults = '/data/vascrem/testoutput/HemeLBSweep/IncreasingAlternativePathways/CollectedShearResults/'
    if path.isdir(CollectedResults)==0:
        os.mkdir(CollectedResults)

    AltPathways = ['4','5','6','7','8']  
    Collapse = ['0','0.1227','0.2248', '0.3178', '0.417', '0.5124', '0.6119', '0.708', '0.8059', '0.9032', '1.0']

    NewNumbering = [ '0','1','2', '3','4','5', '6', '7','8','9','10']
 
 
    for Level in AltPathways:
        counter = -1
        NewPath = CollectedResults + 'Levels_'+Level
        if path.isdir(NewPath)==0:
            os.mkdir(NewPath)
        for i in Collapse:
            counter = counter+1
        
            # mHemeLBDirectory = TerminalOutputFolder+'Levels_'+Level+'/'+i+'/'

            # Move the file
            ShearStress = TerminalOutputFolder  +Level+'/'+i+'/Results/Extracted/surface-traction_3600.vtu'
            NewShearStress = NewPath+'/surface-traction_'+NewNumbering[counter]+'.vtu'

            mv =  'cp ' +ShearStress +' '+NewShearStress
            subprocess.call(mv, shell=True)