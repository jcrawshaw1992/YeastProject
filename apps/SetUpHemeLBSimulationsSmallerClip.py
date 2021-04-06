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
import shlex

    # ---------------------------------


def write_pr2(outputDirectory, SimulationDuration, MinRadii):

    f = open(outputDirectory+"config.pr2", "w")

    V =  4 #Kinematic viscosity -- 4 X 10^-6 m^2/s  V = eta/rho Here v needs to be in the same dims as the input file!!!
    deltaX = MinRadii
    deltaT = float(0.1 * deltaX * deltaX/V)
    print deltaT
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
    "SeedPoint: {x: 0.686062, y: 0.243262, z: 0.025562}\n"+
    "StlFile: config.stl\n"+ \
    "StlFileUnitId: 1\n"+ \
    "TimeStepSeconds: " +str(deltaT)+ "\n"+ \
    "VoxelSize: " +str(deltaX))


if __name__=="__main__":
    t0 = time.time()

    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/TestDiscretisation/'

    MeshDirectory = "/data/vascrem/MeshCollection/IdealisedNetwork/CollapseOf3By3Network/UpperBranch/"
    # MeshDirectory = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/UpperBranch/'
   

    if path.isdir(TerminalOutputFolder)==0:
        os.mkdir(TerminalOutputFolder)

    Collapse = ['2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20']
    
    for i in Collapse:
        mHemeLBDirectory = TerminalOutputFolder+i+'/'
        if path.isdir(mHemeLBDirectory)==0:
            os.mkdir(mHemeLBDirectory)

        MeshFile = MeshDirectory+"ClippedAndScaledMesh."+i+".stl"
        # MeshFile = '/data/vascrem/MeshCollection/IdealisedNetwork/CollapseOf3By3Network/UpperBranch/Scaledmesh.'+str(i)+'.stl'
        shutil.copyfile(MeshFile, mHemeLBDirectory + 'config.stl')
        dX = 0.08/41.0
        print dX
        write_pr2(mHemeLBDirectory, 4001, dX)

        hemelb_setup_exe = '/home/vascrem/hemelb-dev/Tools/setuptool/scripts/hemelb-setup-nogui'

        heme_profile_filename = mHemeLBDirectory + 'config.pr2' 
        command = hemelb_setup_exe + ' ' + heme_profile_filename 
        subprocess.call(command, shell=True)
   
 
         # # # Run HemeLB

        GmyUnstructuredGridReader ="python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml "
        subprocess.call(GmyUnstructuredGridReader, shell=True)
   
   
   






        # # Update the xml file -- Cant do this untill all are set up 
        # update_xml_file(int(501*0.9), mHemeLBDirectory)
        

   