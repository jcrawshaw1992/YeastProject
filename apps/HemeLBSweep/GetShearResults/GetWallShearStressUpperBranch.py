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

    
    # TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/UpperBranchFolder/'
    # # if path.isdir(CollectedResults)==0:
    # #     os.mkdir(CollectedResults)
    # # Collapse = ['0.0','1.227','2.248', '3.178', '4.17', '5.124', '6.119', '7.08', '8.059', '9.032', '10.0']
    # Collapse = ['0','1','2','3','4','5','6','7','8','9','10']
            
    # Parallel = 20
    # SleepyTime = 200
    # NewNumbering = [ '0','1','2', '3','4','5', '6', '7','8','9','10']
    # AvaliablePaths = range(Parallel)
    # for i in Collapse:
        
    #     mHemeLBDirectory = TerminalOutputFolder+i+'/'

    #     GetWallShearStress = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results/Extracted/surface-traction.xtr "
    #     WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(1)+'.txt'
    #     subprocess.Popen(['.././RunHemeLBSweepBash', " ", " ", " ", GetWallShearStress,WaitFileGeneration])
            
 
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/UpperBranchFolder/ArouUpperBranchAround6/'
    # if path.isdir(CollectedResults)==0:
    #     os.mkdir(CollectedResults)
    # Collapse = ['0.0','1.227','2.248', '3.178', '4.17', '5.124', '6.119', '7.08', '8.059', '9.032', '10.0']
    Collapse = ['0','1','2','3','4','5','6','7','8','9','10']
    Collapse =[ '0.5815', '0.5824' , '0.5924' , '0.5932' , '0.5941', '0.5951', '0.5961', '0.5971', '0.6098', '0.6117' , '0.6119', '0.6127', '0.6137', '0.6146']
            
    Parallel = 20
    SleepyTime = 200
    NewNumbering = [ '0','1','2', '3','4','5', '6', '7','8','9','10']
    AvaliablePaths = range(Parallel)
    for i in Collapse:
        
        mHemeLBDirectory = TerminalOutputFolder+i+'/'

        GetWallShearStress = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results/Extracted/surface-traction.xtr "
        WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(1)+'.txt'
        subprocess.Popen(['.././RunHemeLBSweepBash', " ", " ", " ", GetWallShearStress,WaitFileGeneration])
  