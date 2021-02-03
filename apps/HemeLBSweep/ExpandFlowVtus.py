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

    # chmod 700 RunFlowvtuBash
    # Currently this code does not generate pr2 or xml files :S 
    # subprocess.call("chmod 700 RunHemeLBCollapse", shell=True)
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymetricCollapse/'
    Collapse = ['1']
    Parallel = 2
    # SleepyTime = 60
    AvaliablePaths = range(Parallel)
    print AvaliablePaths
    for i in Collapse:
        Core = AvaliablePaths[0]
        mHemeLBDirectory = TerminalOutputFolder+i+'/'

        # GmyUnstructuredGridReader = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml"
        # subprocess.call(GmyUnstructuredGridReader, shell=True)

        # Generate the flow vtus
        GenerateFlowVtus = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results3/Extracted/surface-pressure.xtr " + mHemeLBDirectory + "Results3/Extracted/wholegeometry-velocity.xtr"
        subprocess.call(GenerateFlowVtus, shell=True)

'/data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymmetricGrowth/'
        
        GenerateFlowVtus = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py /data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymmetricGrowth/10/config.vtu /data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymmetricGrowth/10/Results3/Extracted/surface-pressure.xtr /data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymmetricGrowth/10/Results3/Extracted/wholegeometry-velocity.xtr"
        subprocess.call(GenerateFlowVtus, shell=True)


        # Generate waitFile
        # WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
        # subprocess.Popen(['./RunFlowvtuBash', GenerateFlowVtus, WaitFileGeneration ])
        # AvaliablePaths.remove(Core) 
        # # # Check if all positions are taken
        # while len(AvaliablePaths) ==0:
        #     time.sleep(SleepyTime)
        #     # print "Awake and checking for spare cores" 
        #     for P in range(Parallel):
        #         OutputFile = TerminalOutputFolder+'WaitFile'+str(P)+'.txt'
        #         if path.exists(OutputFile):
        #             AvaliablePaths.append(P)
        #             os.remove(OutputFile)
        #     if len(AvaliablePaths) >0:
        #         print AvaliablePaths, "Have found a spare core or two :-) " 
        #         print time.time() - t0, "seconds time"