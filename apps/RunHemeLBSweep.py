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
from xml.etree import ElementTree

def update_xml_file(period, working_directory):

  
    print "update_xml_file"
    # Load automatically generated XML file
    filename = working_directory + 'config.xml'
    
    tree = ElementTree.parse(filename)
    root = tree.getroot()

    # Add monitoring of incompressibility and convergence
    monitoring = ElementTree.SubElement(root, 'monitoring')
    ElementTree.SubElement(monitoring, 'incompressibility') 
    convergence = ElementTree.SubElement(monitoring, 'steady_flow_convergence', {'tolerance': '1e-3', 'terminate': 'false'})
    ElementTree.SubElement(convergence, 'criterion', {'type': 'velocity', 'value': '0.01', 'units': 'm/s'})
    
    # Add definition of properties to be extracted
    extr = ElementTree.SubElement(root, 'properties')

 
    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'surface-traction.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction') 

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'surface-pressure.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='pressure')

    wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'wholegeometry-velocity.xtr'})
    ElementTree.SubElement(wholegeometry, 'geometry', type='whole')
    ElementTree.SubElement(wholegeometry, 'field', type='velocity')

    # Save XML file to disk
    tree.write(filename)

    # ---------------------------------


hemelb_setup_exe = 'env PYTHONPATH=/home/vascrem/hemelb-dev/Tools:/home/vascrem/hemelb-dev/Tools/setuptool:/home/vascrem/hemelb-dev/Tools/hemeTools/converters:$PYTHONPATH /home/vascrem/hemelb-dev/Tools/setuptool/scripts/hemelb-setup-nogui'

def run_hemelb_setup(working_directory):
    print "HemeLB setting up"
    heme_profile_filename = working_directory + 'config.pr2' 
    command = hemelb_setup_exe + ' ' + heme_profile_filename 
    subprocess.call(command, shell=True)
    print "HemeLB set up is DONE"



def write_pr2(outputDirectory, SimulationDuration, MinRadii):

    f = open(outputDirectory+"config.pr2", "w")

    V = 4 #Kinematic viscosity -- 4 mm^2/s  V = eta/rho
    deltaX = MinRadii
    deltaT = float(0.1 * deltaX * deltaX/V)
    print ("DeltaX: " , deltaX , "DeltaT: ", deltaT)

    InletPressure = 0.001
    OutletPressure = 0
    Duration = SimulationDuration*deltaT

    f.write("DurationSeconds: "+ str(Duration) +"\n")
    f.write("Iolets:\n"+ \
    

    "- Centre: {x: 0.500, y: 0.000000, z: 0.000000}\n"+ \
    "  Name: Inlet1\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Inlet\n"+ \
    "- Centre: {x: 0.500, y: -1.4000, z: 0.000000}\n"+ \
    "  Name: Inlet2\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Inlet\n"+ \
    "- Centre: {x: 7.2000, y: 0.000000, z: 0.000000}\n"+ \
    "  Name: Outlet1\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Outlet\n"+ \
    "- Centre: {x: 7.2000, y: -1.4000, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Outlet\n"+ \
    "OutputGeometryFile: config.gmy\n"+ \
    "OutputXmlFile: config.xml\n"+ \
    "SeedPoint: {x: 0.8, y: 0.2, z: 0}\n"+ \
    "StlFile: config.stl\n"+ \
    "StlFileUnitId: 1\n"+ \
    "TimeStepSeconds: " +str(deltaT)+ "\n"+ \
    "VoxelSize: " +str(deltaX))


if __name__=="__main__":
    t0 = time.time()

    # chmod 700 RunHemeLBSweepBash
    # Currently this code does not generate pr2 or xml files :S 
    # subprocess.call("chmod 700 RunHemeLBCollapse", shell=True)
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymmetricCollapse/Between5And6/'
    # command = 'mkdir ' + TerminalOutputFolder
    # subprocess.call(command, shell=True)
    Collapse = ['6','6.1','6.2','6.3','6.4','6.5']
   
    Parallel = 5
    SleepyTime = 150
    AvaliablePaths = range(Parallel)
    print AvaliablePaths
    for i in Collapse:

        Core = AvaliablePaths[0]
        mHemeLBDirectory = TerminalOutputFolder+i+'/'
        MakeNewFolder = 'mkdir '+ mHemeLBDirectory
        subprocess.call(MakeNewFolder, shell=True)
       
        # MeshFile = '/home/vascrem/MeshCollection/IdealisedNetwork/FlowThroughNonSymmetricGrowth/mesh.'+i+'.stl'
        MeshFile = '/home/vascrem/MeshCollection/IdealisedNetwork/CollapseBetween5And6/mesh.'+i+'.stl'
        MoveMeshFile = 'mv ' + MeshFile +' ' + mHemeLBDirectory + 'config.stl'
        subprocess.call(MoveMeshFile, shell=True)
        dX = 0.2/15
        # dX = min(0.2*float(i)/10/3, 0.2/15)
        # if i == '0':
        #     dX = 0.2/15
        # dX = 0.2/40
        write_pr2(mHemeLBDirectory, 4001, dX)
        run_hemelb_setup(mHemeLBDirectory)

        # Update the xml file
        update_xml_file(int(4001*0.8), mHemeLBDirectory)

        # # Run HemeLB
        RunHemeLB = 'mpirun -np 6 hemelb -in ' + mHemeLBDirectory+ 'config.xml -out '+mHemeLBDirectory +'Results2/'
        TerminalOutput = mHemeLBDirectory+'HemeLBTerminalOutput.txt'
        # # Generate the new config.vtu
        GmyUnstructuredGridReader ="python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml "
        # Generate the flow vtus
        GenerateFlowVtus = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results2/Extracted/surface-pressure.xtr " + mHemeLBDirectory + "Results2/Extracted/wholegeometry-velocity.xtr "
        # Generate waitFile
        WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
        subprocess.Popen(['./RunHemeLBSweepBash', RunHemeLB, TerminalOutput, GmyUnstructuredGridReader, GenerateFlowVtus, WaitFileGeneration ])
    
        # python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/GmyUnstructuredGridReader.py  /data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymmetricCollapse/10/config.xml
        # python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py  /data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymmetricCollapse/10/config.vtu /data/vascrem/testoutput/HemeLBSweep/FlowThroughNonSymmetricCollapse/10/Results22/Extracted/surface-pressure.xtr 

        AvaliablePaths.remove(Core) 
        # Check if all positions are taken
        while len(AvaliablePaths) ==0:
            time.sleep(SleepyTime)
            # print "Awake and checking for spare cores" 
            print "Sleep Time"
            for P in range(Parallel):
                OutputFile = TerminalOutputFolder+'WaitFile'+str(P)+'.txt'
                if path.exists(OutputFile):
                    AvaliablePaths.append(P)
                    os.remove(OutputFile)
            if len(AvaliablePaths) >0:
                print AvaliablePaths, "Have found a spare core or two :-) " 
                print time.time() - t0, "seconds time"