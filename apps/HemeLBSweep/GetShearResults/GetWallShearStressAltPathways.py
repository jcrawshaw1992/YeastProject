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


    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/IncreasingAlternativePathways/Levels_'
 
    AltPathways = ['4','5','6','7','8']   
    Collapse = [ '0', '0.1227', '0.2248', '0.3178', '0.417', '0.5124', '0.6119', '0.708', '0.8059', '0.9032', '1.0']

    Parallel = 28
    SleepyTime = 100
    AvaliablePaths = range(Parallel)
    for Level in AltPathways:
        counter = -1
        for i in Collapse:
            counter = counter+1
            Core = AvaliablePaths[0] 
            
            mHemeLBDirectory = TerminalOutputFolder+Level+'/'+i+'/'

            GetWallShearStress = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results/Extracted/surface-traction.xtr "
            WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
            subprocess.Popen(['.././RunHemeLBSweepBash', " ", " ", " ", GetWallShearStress,WaitFileGeneration])
             
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
