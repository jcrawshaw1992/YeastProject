#!/usr/bin/env python
# vasc_res_sim.py
import subprocess
import vtk
import shutil
import os
from xml.etree import ElementTree
import glob
from argparse import ArgumentParser
import numpy as np
import time
import matplotlib.pyplot as plt
import csv
import pdb
import string
import math
from subprocess import Popen, PIPE



if __name__=="__main__":

    # Idea 1 -- This works when I run this code as a python file, but if I call this python code from chaste, this line wont work 
    subprocess.call(' mpirun -np 1 hemelb -in /Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/config.xml', shell=True)
    
    # Idea 2 -- 
    # This line will open a new terminal and run hemelB when I run this python script directly. But when I call this python code through Chaste,
    # I get the error below   -- I obviously need a display, which is the next thing I go about doing 
    subprocess.call(['xterm', '-e', 'mpirun -np 1 hemelb -in /Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/config.xml'])

    # // I get this error 
    # xterm: Xt error: Can't open display: 
    #  xterm: DISPLAY is not set


    # Idea 3 -- 
    # this sucessffully opens a display and can run code but not sure how to exicute mpirun -np 1 hemelb -in /Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/config.xml from here
    xterm = subprocess.Popen('xterm -display :0', shell=True)
    xterm.communicate('-e emacs')
    # This didnt work 
    # xterm.communicate(['xterm', '-e','mpirun -np 1 hemelb -in /Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/config.xml'])

   print "----- finished? -----\n\n"



#    Hi All, 
#    I am trying to incorportate Miguel's HemeLB 'directly' into Chaste by calling a sequence of terminal commands using std::system() from within a force class. It was 
#    going well until the primary command which uses mpirun.  Chaste seems to skip happily past this command, not crashing, but also not running the command. I am assuming 
#    mpirun isn't working because its clashing with something inhertent to Chaste. Does anyone have any experience or advice with this? So far I am looking 
#    into automatically opening another terminal and running HemeLB from that new terminal and then closing it. Thanks
   
   
   
   , has anyone ever tired to run mpirun using std::system from within Chaste 