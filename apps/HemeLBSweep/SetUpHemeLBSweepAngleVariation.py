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
import tempfile 
import shlex


def write_pr2(outputDirectory, SimulationDuration, MinRadii, Seed):

    f = open(outputDirectory+"config.pr2", "w")

    V =  4 #Kinematic viscosity -- 4 X 10^-6 m^2/s  V = eta/rho Here v needs to be in the same dims as the input file!!!
    deltaX = MinRadii
    deltaT = float(0.1 * deltaX * deltaX/V)
    print ('DeltaX: ' , deltaX , 'DeltaT: ', deltaT, ' eta: ', V )

    InletPressure = 100
    OutletPressure = 0
    Duration = SimulationDuration*deltaT

    f.write("DurationSeconds: "+ str(Duration) +"\n")
    f.write("Iolets:\n"+ \
    # ------Inlets ----
    "- Centre: {x: 0.100, y: 0.000000, z: 0.000000}\n"+ \
    "  Name: Inlet1\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: 0.100, y: -0.27, z: 0.000000}\n"+ \
    "  Name: Inlet2\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: 0.100, y: 0.27, z: 0.000000}\n"+ \
    "  Name: Inlet3\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Inlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: 1.45, y: -0.27, z: 0.000000}\n"+ \
    "  Name: Outlet1\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Outlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: 1.45, y: 0, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Outlet\n"+ \
    # ------
    "- Centre: {x: 1.45, y: 0.270, z: 0.000000}\n"+ \
    "  Name: Outlet3\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Outlet\n"+ \
    "OutputGeometryFile: config.gmy\n"+ \
    "OutputXmlFile: config.xml\n"+ \
    "SeedPoint: {x: "+str(Seed[0])+", y: "+str(Seed[1])+", z: "+str(Seed[2])+"}\n"+
    "StlFile: config.stl\n"+ \
    "StlFileUnitId: 1\n"+ \
    "TimeStepSeconds: " +str(deltaT)+ "\n"+ \
    "VoxelSize: " +str(deltaX))


if __name__=="__main__":
    t0 = time.time()

    # chmod 700 RunHemeLBSweepBash
    # subprocess.call("chmod 700 RunHemeLBSweepBash", shell=True)
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/AngleVariation_3X3Network/'
    MeshDirectory = "/data/vascrem/MeshCollection/IdealisedNetwork/AngleVariation_3X3Network/"
    if path.isdir(TerminalOutputFolder)==0:
        os.mkdir(TerminalOutputFolder)

    Collapse = ['0.6','0.59','0.61','0.55','0.56','0.57','0.58','0.62','0.63','0.64','0.65']
   
    # Angle == 3: Seed = [2.19658, 1.97595, -0.155648]
    # Angle == 6:  Seed = [2.05519, 1.98863, 0.165396]
    # Angle == 2.2: Seed = [2.11408, 2.11277, 0.199536]

    counter = -1
    BifucationAngles = [ 'PI_2.2/','PI_3/','PI_6/' ] 
    Seeds = [[2.11408, 2.11277, 0.199536] , [2.19658, 1.97595, -0.155648],[2.05519, 1.98863, 0.165396] ]
    for i in Collapse:
        counter = -1
        for angle in BifucationAngles:
            counter = counter+1
            if path.isdir(TerminalOutputFolder+angle)==0:
                    os.mkdir(TerminalOutputFolder+angle)
        
            mHemeLBDirectory = TerminalOutputFolder+angle+str(float(i)*10)+'/'
            if path.isdir(mHemeLBDirectory)==0:
                os.mkdir(mHemeLBDirectory)
            if path.isdir(mHemeLBDirectory+'Results/')==1:
                os.rename(mHemeLBDirectory+'Results/',mHemeLBDirectory+'Results_PriorRun/')

            MeshFile = MeshDirectory+angle+"Mesh_"+i+".stl"

            shutil.copyfile(MeshFile, mHemeLBDirectory + 'config.stl')

            dX = 0.08/41.0

            write_pr2(mHemeLBDirectory, 4001, dX,  Seeds[counter])

            hemelb_setup_exe = '/home/vascrem/hemelb-dev/Tools/setuptool/scripts/hemelb-setup-nogui'

            heme_profile_filename = mHemeLBDirectory + 'config.pr2' 
            command = hemelb_setup_exe + ' ' + heme_profile_filename 
            args = shlex.split(command)
            subprocess.Popen(args)
        


            
            