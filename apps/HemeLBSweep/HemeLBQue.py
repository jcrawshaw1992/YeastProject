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

    # command = 'python RunHemeLBSweepAngleVariation2.py'
    # subprocess.call(command, shell=True)

    # command = 'python RunHemeLBPressureSweep.py'
    # subprocess.call(command, shell=True)

    # command = 'python RunHemeLBSweep3By3Network.py'
    # subprocess.call(command, shell=True)

    command = 'python RunHemeLBSweepAdditionalAlternatPathways.py'
    subprocess.call(command, shell=True)


    # command = 'python RunHemeLBSweepVaryingDensityScalledRadius.py'
    # subprocess.call(command, shell=True)

    command = 'python RunHemeLBSweepVaryingDensityScalledRadius2.py'
    subprocess.call(command, shell=True)


    command = 'python RunHemeLBSweepAdditionalAlternatPathways.py'
    subprocess.call(command, shell=True)




    



#  3319  sudo du -h | sort -rh | head -5
#
