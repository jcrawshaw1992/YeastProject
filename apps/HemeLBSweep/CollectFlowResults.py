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
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymmetricGrowth/'
    # Collapse = ['9','8','7','6','5','4','3','2','1','0']
    Collapse = ['10','11','12','13','14','15','16','17','18','19','20']

    mkdir = "mkdir /data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymmetricGrowth/VelocityFiles/"
    subprocess.call(mkdir, shell=True)

    for i in Collapse:
        mHemeLBDirectory = TerminalOutputFolder+i+'/'

        # mv =  "cp " +mHemeLBDirectory +"Results3/Extracted/surface-pressure_8000.vtu "+ TerminalOutputFolder + "CollectedResults2/surface-pressure_"+i+".vtu"   
        # subprocess.call(mv, shell=True)

        mv =  "cp " +mHemeLBDirectory +"Results3/Extracted/wholegeometry-velocity_16000.vtu " +TerminalOutputFolder + "VelocityFiles/wholegeometry-velocity_"+i+".vtu"   
        subprocess.call(mv, shell=True)

        mv =  "cp " +mHemeLBDirectory +"config.xml " +TerminalOutputFolder + "VelocityFiles/config_"+i+".xml"  
        subprocess.call(mv, shell=True)

        # mv =  "cp " +mHemeLBDirectory +"config.stl " +"/data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymmetricGrowth/stls/config"+i+".stl"   
        # subprocess.call(mv, shell=True)
