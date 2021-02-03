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

hemelb_setup_exe = 'env PYTHONPATH=/home/vascrem/hemelb-dev/Tools:/home/vascrem/hemelb-dev/Tools/setuptool:/home/vascrem/hemelb-dev/Tools/hemeTools/converters:$PYTHONPATH /home/vascrem/hemelb-dev/Tools/setuptool/scripts/hemelb-setup-nogui'

def run_hemelb_setup():
    print "HemeLB setting up"
    heme_profile_filename = working_directory + 'config.pr2' 
    command = hemelb_setup_exe + ' ' + heme_profile_filename 
    subprocess.call(command, shell=True)
    print "HemeLB set up is DONE"



def write_pr2(outputDirectory, SimulationDuration, MinRadii):

    f = open(outputDirectory+"config.pr2", "w")

    mRadius = MinRadii

    V = 4 #Kinematic viscosity -- 4 mm^2/s  V = eta/rho
    deltaX = 2*mRadius/15
    deltaT = 0.1 * deltaX * deltaX/V

    InletPressure = 0.001
    OutletPressure = 0
    Duration = SimulationDuration*deltaT

    f.write("DurationSeconds: "+ str(Duration) +"\n")
    f.write("Iolets:\n"+ \
    "- Centre: {x: 0.005000, y: 0.000000, z: 0.000000}\n"+ \
    "  Name: Inlet1\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.005372\n"+ \
    "  Type: Inlet\n"+ \
    "- Centre: {x: 0.005000, y: -0.014000, z: 0.000000}\n"+ \
    "  Name: Inlet2\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.005372\n"+ \
    "  Type: Inlet\n"+ \
    "- Centre: {x: 0.072000, y: 0.000000, z: 0.000000}\n"+ \
    "  Name: Outlet1\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.005372\n"+ \
    "  Type: Outlet\n"+ \
    "- Centre: {x: 0.072000, y: -0.014000, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.005372\n"+ \
    "  Type: Outlet\n"+ \
    "OutputGeometryFile: config.gmy\n"+ \
    "OutputXmlFile: config.xml\n"+ \
    "SeedPoint: {x: 0.049260, y: -0.010895, z: 0.001275}\n"+ \
    "StlFile: config.stl\n"+ \
    "StlFileUnitId: 1\n"+ \
    "TimeStepSeconds: " +str(deltaT)+ "\n"+ \
    "VoxelSize: " +str(deltaX))

if __name__=="__main__":


    # chmod 700 RunFlowvtuBash
    # Currently this code does not generate pr2 or xml files :S 
    # subprocess.call("chmod 700 RunHemeLBCollapse", shell=True)
    Collapse = ['10']
    Parallel = 2
    SleepyTime = 60
    AvaliablePaths = range(Parallel)
    print AvaliablePaths
    for i in Collapse:

        TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/TestFlowThroughNonSymetricCollapse/'
        mHemeLBDirectory = TerminalOutputFolder+i+'/'
        command = 'mv ' + mHemeLBDirectory + 'HemeLBFluid/config.stl ' + mHemeLBDirectory + 'config.stl'
        subprocess.call(command, shell=True)

        command = 'mv ' + mHemeLBDirectory + 'HemeLBFluid/Chaste.vtu ' + mHemeLBDirectory + 'Chaste.vtu'
        subprocess.call(command, shell=True)

        
        write_pr2(TerminalOutputFolder, 40000, 0.00150*float(i)/10.0)
        run_hemelb_setup(mHemeLBDirectory )
        


    # Collapse = ['3','2']
    # Parallel = 2
    # SleepyTime = 60
    # AvaliablePaths = range(Parallel)
    # print AvaliablePaths
    # for i in Collapse:
    #     Core = AvaliablePaths[0]
    #     mHemeLBDirectory = TerminalOutputFolder+i+'/HemeLBFluid/'
    #     # Generate the flow vtus
    #     GenerateFlowVtus = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results2/Extracted/surface-pressure.xtr " + mHemeLBDirectory + "Results2/Extracted/wholegeometry-velocity.xtr"
    #     # Generate waitFile
    #     WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
    #     subprocess.Popen(['./RunFlowvtuBash', GenerateFlowVtus, WaitFileGeneration ])
    #     AvaliablePaths.remove(Core) 
    #     # # Check if all positions are taken
    #     while len(AvaliablePaths) ==0:
    #         time.sleep(SleepyTime)
    #         # print "Awake and checking for spare cores" 
    #         for P in range(Parallel):
    #             OutputFile = TerminalOutputFolder+'WaitFile'+str(P)+'.txt'
    #             if path.exists(OutputFile):
    #                 AvaliablePaths.append(P)
    #                 os.remove(OutputFile)
    #         if len(AvaliablePaths) >0:
    #             print AvaliablePaths, "Have found a spare core or two :-) " 
    #             print time.time() - t0, "seconds time"
