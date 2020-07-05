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

# Pope



if __name__=="__main__":

    # Define arguments to be parsed
    hemelb_setup_exe = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'

    parser = ArgumentParser(description='Update the necessary xml file for running HemeLB')
    parser.add_argument('-directory', dest='directory', default='/Users/jcrawshaw/Documents/testoutput/', help='Need to provide a destination to find the xml file to edit ')
   
    # Only the final time output is needed, so this should be the period --- this is because it is time costly
    # to write out files, aaaand until the last time point, I cant be sure 
    # that the simulation is at equalibirum --- prior to equi I just have rubbish
   
    args = parser.parse_args()
    # xml_config_file = args.directory + 'config.xml'

#     xml_config_file = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/config.xml'
# # 
#     # subprocess.call(' mpirun -np 1 hemelb -in /Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/config.xml')
    
    

 
#     # print "**************   hemelb step   **************   "
#     # shutil.rmtree(args.directory + 'results/', ignore_errors=True)
#     command = 'mpirun -np 1 hemelb -in ' + xml_config_file
#     subprocess.call( command, shell=True)

    # subprocess.call(['open', '-W', '-a', 'Terminal.app', '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/apps/update_xml_file.py', '--args', 'update_xml_file.py'])
# python /Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/apps/update_xml_file.py
    # subprocess.call('python bb.py', creationflags=subprocess.CREATE_NEW_CONSOLE)



    # Progress summary-- I need to have a display to run from cpp as I do below, making the display works, but I for some reason cant 
    # run the mpirun, using this, even just from phyton. But the above code works if running straight from python

    # subprocess.call(['xterm', '-e', 'mpirun -np 1 hemelb -in /Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/config.xml'])
    # Popen.communicate(['xterm', '-e', 'tail -f %s' % PIPE_PATH])


    xterm = subprocess.Popen('xterm -display :0', shell=True)
    xterm.communicate('-e emacs')
    # print "dpe"
    # xterm.communicate(['xterm', '-e','mpirun -np 1 hemelb -in /Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/config.xml'])


    # fd = open("output", "w")
    # # Open xterm window with logs from that file
    # p = subprocess.Popen(["xterm", "-e", "tail", "-f", "output"])
    # # Do logging
    # fd.write("Hello, world!")
    # fd.flush()



    # https://stackoverflow.com/questions/34800280/python-send-a-command-to-xterm

    # command = 'mpirun -np 1 hemelb -in /Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/config.xml'
    # subprocess.Popen('xterm -hold -e "%s"' % command)

    # Popen.communicate('mpirun -np 1 hemelb -in /Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/config.xml')
    # subprocess.call(['rxvt', '-e', 'mpirun -np 1 hemelb -in /Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/config.xml'])

    # // I get this error 
#     xterm: Xt error: Can't open display: 
#    xterm: DISPLAY is not set



    # Should i just send this to the server via python
    # https://stackoverflow.com/questions/12354047/x11-forwarding-with-paramiko
    print "----- finished? -----\n\n"