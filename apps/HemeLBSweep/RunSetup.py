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
from argparse import ArgumentParser




if __name__=="__main__":

    print 'Hi;'
    hemelb_setup_exe = 'env PYTHONPATH=/home/vascrem/hemelb-dev/Tools:/home/vascrem/hemelb-dev/Tools/setuptool:/home/vascrem/hemelb-dev/Tools/hemeTools/converters:$PYTHONPATH /home/vascrem/hemelb-dev/Tools/setuptool/scripts/hemelb-setup-nogui'

    # Define arguments to be parsed
    parser = ArgumentParser(description='Run a vascular remodelling simulation')
    parser.add_argument('-working_directory', dest='working_directory',type=str, default="/data/vascrem/testoutput/ParameterSweep/Cylinder/SweepTerminalOutputs/", help='Set the archieve files for chaste to continue from ')
    parser.add_argument('-WaitFileGeneration', dest='WaitFileGeneration',type=str, default="/data/vascrem/testoutput/ParameterSweep/Cylinder/SweepTerminalOutputs/", help='Set the archieve files for chaste to continue from ')
       
    # Parse arguments (this will create args.flow_to_vtu etc. variables)
    args = parser.parse_args()
    working_directory = args.working_directory
    WaitFileGeneration = args.WaitFileGeneration
 

    heme_profile_filename = working_directory + 'config.pr2' 
    command = hemelb_setup_exe + ' ' + heme_profile_filename 
    subprocess.call(command, shell=True)


    f = open(WaitFileGeneration, "a")
    f.write("Setup Complete")
    f.close()



 