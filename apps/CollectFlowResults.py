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
    
    

    # Collapse = ['1','2','3','4','5','6','7','8','9','10','0','5.8','5.5','5.6','5.9','6.1','6.2','6.3','6.4','6.5']
    
    Collapse = ['5.5','5.6','5.7','5.8','5.9','6.0','6.1','6.2','6.3','6.4','6.5']
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/AngleVariation_3X3Network/PI_2.2/'

    NewDirectory =TerminalOutputFolder+"CollectedResults/"
    mkdir = "mkdir  "+NewDirectory
    
    subprocess.call(mkdir, shell=True)
    for i in Collapse:
 
        mHemeLBDirectory = TerminalOutputFolder+i+'/'
        # print mHemeLBDirectory
        # print TerminalOutputFolder

        mv =  'cp ' +mHemeLBDirectory +'Results/Extracted/surface-pressure_3600.vtu '+ TerminalOutputFolder +'CollectedResults/surface-pressure_'+i+'.vtu'   
        subprocess.call(mv, shell=True)

        mv =  'cp ' +mHemeLBDirectory +'Results/Extracted/wholegeometry-velocity_3600.vtu '+ TerminalOutputFolder + 'CollectedResults/wholegeometry-velocity_'+i+'.vtu'   
        subprocess.call(mv, shell=True)

        # mv =  "cp " + TerminalOutputFolder +i+'/'+"Results/Extracted/wholegeometry-velocity_3600.vtu " +NewDirectory +"wholegeometry-velocity_"+i+".vtu"   
        # subprocess.call(mv, shell=True)

        # mv =  "cp " + TerminalOutputFolder +i+'/'+"Results/Extracted/surface-pressure_3600.vtu " +NewDirectory +"surface-pressure_"+i+".vtu"   
        # subprocess.call(mv, shell=True)


        # mv =  "cp " +mHemeLBDirectory +"config.xml " +TerminalOutputFolder + "VelocityFiles/config_"+i+".xml"  
        # subprocess.call(mv, shell=True)

        # mv =  "cp " +mHemeLBDirectory +"config.stl " +"/data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymmetricGrowth/stls/config"+i+".stl"   
        # subprocess.call(mv, shell=True)
