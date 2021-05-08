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

    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/UpperBranchFolder/'
    
    NewPath = TerminalOutputFolder+'CollectedShearResults/'
    if path.isdir(NewPath)==0:
        os.mkdir(NewPath)
   
    Collapse = ['0','1','2','3','4','5','6','7','8','9','10']

    # for i in Collapse:
    #     # Move the file
    #     ShearStress = TerminalOutputFolder +i+'/Results/Extracted/surface-traction_3600.vtu'
    #     NewShearStress = NewPath+'/surface-traction_'+ i +'.vtu'
    
    #     mv =  'cp ' +ShearStress +' '+NewShearStress
    #     subprocess.call(mv, shell=True)

# ---------

    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/UpperBranchFolder/ArouUpperBranchAround6'
    
    NewPath = TerminalOutputFolder+'CollectedShearResultsAround6/'
    if path.isdir(NewPath)==0:
        os.mkdir(NewPath)
   
    # Collapse = ['0','1','2','3','4','5','6','7','8','9','10']
    Collapse =[ '0.5815', '0.5824' , '0.5924' , '0.5932' , '0.5941', '0.5951', '0.5961', '0.5971', '0.6098', '0.6117' , '0.6119', '0.6127', '0.6137', '0.6146']
    
    for i in Collapse:
        # Move the file
        ShearStress = TerminalOutputFolder +'/'+i+'/Results/Extracted/surface-traction_3600.vtu'
        NewShearStress = NewPath+'/surface-traction_'+ i +'.vtu'
    
        mv =  'cp ' +ShearStress +' '+NewShearStress
        subprocess.call(mv, shell=True)


