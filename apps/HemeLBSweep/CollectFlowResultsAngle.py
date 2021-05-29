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
import os.path
from os import path



if __name__=="__main__":
    t0 = time.time()



    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/AngleVariation_3X3Network/'
    Collapse = [ '0', '0.1227', '0.2248', '0.3178', '0.417', '0.5124', '0.6119', '0.708', '0.8059', '0.9032', '1.0']
    # Collapse = [ '1.227', '2.248', '3.178', '4.17', '5.124', '6.119', '7.08', '8.059', '9.032', '10.0']
    
    NewNumbering = [ '0','1','2', '3','4','5', '6', '7','8','9','10']
    Scalling = ['PI_6/']
    NewPathTop = '/data/vascrem/testoutput/HemeLBSweep/AngleVariation_3X3Network/AngleCollected/'
    if path.isdir(NewPathTop)==0:
        os.mkdir(NewPathTop)
             

    for j in Scalling:
        NewPath = NewPathTop + j
        if path.isdir(NewPath)==0:
            os.mkdir(NewPath)
        counter =-1
        for i in Collapse:
            counter = counter+1
            print counter                                             
            Pressure = TerminalOutputFolder +j+i+'/Results/Extracted/surface-pressure_9500.vtu'
            
            # Traction = TerminalOutputFolder +j+i+'/Results/Extracted/surface-traction_20000.vtu'
            # Velocity = TerminalOutputFolder +j+i+'/Results/Extracted/wholegeometry-velocity_20000.vtu'
            Mesh = TerminalOutputFolder +j+i+'/config.stl'
            I = str(float(i)/100)
            
            NewPressure = NewPath+'surface-pressure_'+NewNumbering[counter]+'.vtu'
            # NewVelocity = NewPath+'wholegeometry-velocity_'+NewNumbering[counter]+'.vtu'
            # NewTraction = NewPath+'surface-traction_'+NewNumbering[counter]+'.vtu'
            NewMesh = NewPath+'mesh.'+NewNumbering[counter]+'.vtu'

            mv =  'cp ' +Pressure +' '+NewPressure
            subprocess.call(mv, shell=True)

            # mv =  'cp ' +Velocity +' '+NewVelocity
            # subprocess.call(mv, shell=True)
            # mv =  'cp ' +Traction +' '+NewTraction
            # subprocess.call(mv, shell=True)
            mv =  'cp ' +Mesh +' '+NewMesh
            subprocess.call(mv, shell=True)



