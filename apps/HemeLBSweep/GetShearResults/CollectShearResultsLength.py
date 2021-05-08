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

    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/Length_Variation/HorizontalLength_'
    CollectedResults = '/data/vascrem/testoutput/HemeLBSweep/Length_Variation/CollectedShearResults/'
    if path.isdir(CollectedResults)==0:
        os.mkdir(CollectedResults)
   
    Scalling = ['0.2', '0.4', '0.6', '0.8', '1','1.2', '1.4', '1.6', '1.8', '2', '2.2', '2.4', '2.6', '2.8', '3']   
    # Collapse = ['0.0','1.227','2.248', '3.178', '4.17', '5.124', '6.119', '7.08', '8.059', '9.032', '10.0']
    Collapse = ['0','1','2','3','4','5','6','7','8','9','10']

     


    NewNumbering = [ '0','1','2', '3','4','5', '6', '7','8','9','10']


    for Level in Scalling:
        counter = -1
        NewPath = CollectedResults + '/HorizontalLength_'+Level
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