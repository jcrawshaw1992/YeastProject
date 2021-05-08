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

    
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/MiddleBranchFolder/'
    # if path.isdir(CollectedResults)==0:
    #     os.mkdir(CollectedResults)
    # Collapse = ['0.0','1.227','2.248', '3.178', '4.17', '5.124', '6.119', '7.08', '8.059', '9.032', '10.0']
    Collapse = ['0','1']#'2','3','4','5','6','7','8','9','10']
            

    for i in Collapse:
        
        mHemeLBDirectory = TerminalOutputFolder+i+'/'

        GetWallShearStress = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results/Extracted/surface-traction.xtr "
        WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(1)+'.txt'
        subprocess.Popen(['.././RunHemeLBSweepBash', " ", " ", " ", GetWallShearStress,WaitFileGeneration])
            
 