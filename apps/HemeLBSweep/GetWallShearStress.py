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
from xml.etree import ElementTree
import tempfile 
from stl import mesh
from datetime import datetime




if __name__=="__main__":
    t0 = time.time()


    # chmod 700 RunHemeLBSweepBash
    # subprocess.call("chmod 700 RunHemeLBSweepBash", shell=True)
    
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/IncreasingDensity/'
    CollectedResults = '/data/vascrem/testoutput/HemeLBSweep/IncreasingDensity/CollectedShearResults/'
    # if path.isdir(CollectedResults)==0:
    #     os.mkdir(CollectedResults)
  
  
    Density = ['1','2','3']   
    Collapse = ['0','0.1227','0.2248', '0.3178', '0.417', '0.5124', '0.6119', '0.708', '0.8059', '0.9032', '1.0']
            
    Parallel = 20
    SleepyTime = 200
    NewNumbering = [ '0','1','2', '3','4','5', '6', '7','8','9','10']
    AvaliablePaths = range(Parallel)
    for Level in Density:
        counter = -1
        # NewPath = CollectedResults + 'CollectedResults/Density_'+Level
        # if path.isdir(NewPath)==0:
        #     os.mkdir(NewPath)
        for i in Collapse:
            counter = counter+1
            Core = AvaliablePaths[0] 
            
            mHemeLBDirectory = TerminalOutputFolder+'Density_'+Level+'/'+i+'/'

            GetWallShearStress = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results/Extracted/surface-traction.xtr "
            WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
            subprocess.Popen(['./RunHemeLBSweepBash', " ", " ", " ", GetWallShearStress,WaitFileGeneration ])


            # Move the file


                   # Velocity = TerminalOutputFolder +'Density_' +j+'/'+ i+'/Results/Extracted/wholegeometry-velocity_3600.vtu'

            # if Level == '1' or Level == '2':
            #     ShearStress = TerminalOutputFolder +'Density_' +Level+'/'+i+'/Results/Extracted/surface-traction_3600.vtu'
            #     I = str(float(i)/100)
            #     NewShearStress = NewPath+'/surface-traction_'+NewNumbering[counter]+'.vtu'
            # else: 
            #     ShearStress = TerminalOutputFolder +'Density_' +Level+'/'+i+'/Results/Extracted/surface-traction_4500.vtu'
            #     I = str(float(i)/100)
            #     NewShearStress = NewPath+'/surface-traction_'+NewNumbering[counter]+'.vtu'

            # mv =  'cp ' +ShearStress +' '+NewShearStress
            # subprocess.call(mv, shell=True)

             
            AvaliablePaths.remove(Core) 
            # Check if all positions are taken
            while len(AvaliablePaths) ==0:
                time.sleep(SleepyTime)
                # print "Awake and checking for spare cores" 
                print "Sleep Time"
                for P in range(Parallel):
                    OutputFile = TerminalOutputFolder+'WaitFile'+str(P)+'.txt'
                    if path.exists(OutputFile):
                        AvaliablePaths.append(P)
                        os.remove(OutputFile)
                if len(AvaliablePaths) >0:
                    print AvaliablePaths, "Have found a spare core or two :-) " 
                    print time.time() - t0, "seconds time"