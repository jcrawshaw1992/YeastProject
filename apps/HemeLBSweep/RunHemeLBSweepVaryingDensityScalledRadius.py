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

    # Edit the initial condition to my satisfaction

    file = open(filename).readlines()
    #Create temp file
    file_path = filename
    pattern = '<uniform units="mmHg" value="0.0" />'
    subst = '<uniform units="mmHg" value="50.0" />'

    fh, abs_path = tempfile.mkstemp()
    with os.fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    #Copy the file permissions from the old file to the new file
    shutil.copymode(file_path, abs_path)
    #Remove original file
    os.remove(file_path)
    #Move new file
    shutil.move(abs_path, file_path)



    # ---------------------------------


hemelb_setup_exe = 'env PYTHONPATH=/home/vascrem/hemelb-dev/Tools:/home/vascrem/hemelb-dev/Tools/setuptool:/home/vascrem/hemelb-dev/Tools/hemeTools/converters:$PYTHONPATH /home/vascrem/hemelb-dev/Tools/setuptool/scripts/hemelb-setup-nogui'

def run_hemelb_setup(working_directory):
    print "HemeLB setting up"
    heme_profile_filename = working_directory + 'config.pr2' 
    command = hemelb_setup_exe + ' ' + heme_profile_filename 
    subprocess.call(command, shell=True)
    print "HemeLB set up is DONE"

def GetTheDetailsOfTheMesh(MeshFile):
    # Using an existing stl file:
    your_mesh = mesh.Mesh.from_file(MeshFile)

    # Or creating a new mesh (make sure not to overwrite the `mesh` import by
    # naming it `mesh`):
    MaxPoint = 0
    MinPoint =100

    for i in range(0,len(your_mesh.points)):
        vector = your_mesh.points[i]
        for i in [0,3,6]:
            # print vector[i]
            if vector[i]<MinPoint:
                MinPoint = vector[i]
            elif vector[i]>MaxPoint:
                MaxPoint = vector[i]

    Length = MaxPoint - MinPoint


    # So now the boundaries are

    LeftBoundary =  MinPoint + Length/100
    RightBoundary =  MaxPoint - Length/100
    Seed = your_mesh.points[200][0:3]

    Result = [LeftBoundary,RightBoundary, Seed[0],Seed[1],Seed[2]]

    # print "LeftBoundary ", LeftBoundary
    # print "RightBoundary ", RightBoundary
    # print "seed" , Seed
    return Result        
        



def write_pr21(outputDirectory, SimulationDuration, MinRadii, Seed, Boundaries):

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
    "- Centre: {x: "+str(0.06)+", y: 0.2, z: 0.000000}\n"+ \
    "  Name: Inlet1\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.0, z: 0.000000}\n"+ \
    "  Name: Inlet2\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.2, z: 0.000000}\n"+ \
    "  Name: Inlet2\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.4, z: 0.000000}\n"+ \
    "  Name: Inlet3\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: 0.2, z: 0.000000}\n"+ \
    "  Name: Outlet1\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: 0.0, z: 0.000000}\n"+ \
    "  Name: Outlet1\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: -0.2, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    # ------
    "- Centre: {x: "+str(1.49)+", y: -0.4, z: 0.000000}\n"+ \
    "  Name: Outlet3\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    "OutputGeometryFile: config.gmy\n"+ \
    "OutputXmlFile: config.xml\n"+ \
    "SeedPoint: {x: "+str(Seed[0])+", y: "+str(Seed[1])+", z: "+str(Seed[2])+"}\n"+
    "StlFile: config.stl\n"+ \
    "StlFileUnitId: 1\n"+ \
    "TimeStepSeconds: " +str(deltaT)+ "\n"+ \
    "VoxelSize: " +str(deltaX))



def write_pr22(outputDirectory, SimulationDuration, MinRadii, Seed, Boundaries):

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
    "- Centre: {x: "+str(0.06)+", y: 0.18, z: 0.000000}\n"+ \
    "  Name: Inlet1\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.0, z: 0.000000}\n"+ \
    "  Name: Inlet2\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.18, z: 0.000000}\n"+ \
    "  Name: Inlet2\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.34, z: 0.000000}\n"+ \
    "  Name: Inlet3\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.51, z: 0.000000}\n"+ \
    "  Name: Inlet3\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: 0.18, z: 0.000000}\n"+ \
    "  Name: Outlet1\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: 0.0, z: 0.000000}\n"+ \
    "  Name: Outlet1\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: -0.18, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    # ------
    "- Centre: {x: "+str(1.49)+", y: -0.34, z: 0.000000}\n"+ \
    "  Name: Outlet3\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
     # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: -0.51, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    "OutputGeometryFile: config.gmy\n"+ \
    "OutputXmlFile: config.xml\n"+ \
    "SeedPoint: {x: "+str(Seed[0])+", y: "+str(Seed[1])+", z: "+str(Seed[2])+"}\n"+
    "StlFile: config.stl\n"+ \
    "StlFileUnitId: 1\n"+ \
    "TimeStepSeconds: " +str(deltaT)+ "\n"+ \
    "VoxelSize: " +str(deltaX))


def write_pr23(outputDirectory, SimulationDuration, MinRadii, Seed, Boundaries):

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
    "- Centre: {x: "+str(0.06)+", y: 0.15, z: 0.000000}\n"+ \
    "  Name: Inlet1\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.0, z: 0.000000}\n"+ \
    "  Name: Inlet2\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.15, z: 0.000000}\n"+ \
    "  Name: Inlet2\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.28, z: 0.000000}\n"+ \
    "  Name: Inlet3\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.43, z: 0.000000}\n"+ \
    "  Name: Inlet3\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.57, z: 0.000000}\n"+ \
    "  Name: Inlet3\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: 0.15, z: 0.000000}\n"+ \
    "  Name: Outlet1\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: 0.0, z: 0.000000}\n"+ \
    "  Name: Outlet1\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: -0.15, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    # ------
    "- Centre: {x: "+str(1.49)+", y: -0.28, z: 0.000000}\n"+ \
    "  Name: Outlet3\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
     # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: -0.43, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
         # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: -0.57, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    "OutputGeometryFile: config.gmy\n"+ \
    "OutputXmlFile: config.xml\n"+ \
    "SeedPoint: {x: "+str(Seed[0])+", y: "+str(Seed[1])+", z: "+str(Seed[2])+"}\n"+
    "StlFile: config.stl\n"+ \
    "StlFileUnitId: 1\n"+ \
    "TimeStepSeconds: " +str(deltaT)+ "\n"+ \
    "VoxelSize: " +str(deltaX))


def write_pr24(outputDirectory, SimulationDuration, MinRadii, Seed, Boundaries):

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
    "- Centre: {x: "+str(0.06)+", y: 0.12, z: 0.000000}\n"+ \
    "  Name: Inlet1\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.0, z: 0.000000}\n"+ \
    "  Name: Inlet2\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.12, z: 0.000000}\n"+ \
    "  Name: Inlet2\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.25, z: 0.000000}\n"+ \
    "  Name: Inlet3\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.38, z: 0.000000}\n"+ \
    "  Name: Inlet3\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.50, z: 0.000000}\n"+ \
    "  Name: Inlet3\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(0.06)+", y: -0.63, z: 0.000000}\n"+ \
    "  Name: Inlet3\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Inlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: 0.12, z: 0.000000}\n"+ \
    "  Name: Outlet1\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: 0.0, z: 0.000000}\n"+ \
    "  Name: Outlet1\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: -0.12, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    # ------
    "- Centre: {x: "+str(1.49)+", y: -0.25, z: 0.000000}\n"+ \
    "  Name: Outlet3\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
     # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: -0.37, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
     # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: -0.50, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
         # -----Outlets -----
    "- Centre: {x: "+str(1.49)+", y: -0.63, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.06\n"+ \
    "  Type: Outlet\n"+ \
    "OutputGeometryFile: config.gmy\n"+ \
    "OutputXmlFile: config.xml\n"+ \
    "SeedPoint: {x: "+str(Seed[0])+", y: "+str(Seed[1])+", z: "+str(Seed[2])+"}\n"+
    "StlFile: config.stl\n"+ \
    "StlFileUnitId: 1\n"+ \
    "TimeStepSeconds: " +str(deltaT)+ "\n"+ \
    "VoxelSize: " +str(deltaX))





def EditRunTime(working_directory):

    print "update_xml_file"
    # Load automatically generated XML file
    filename = working_directory + 'config.xml'

    # Edit the initial condition to my satisfaction

    file = open(filename).readlines()
    #Create temp file
    file_path = filename
   
    subst = '<steps units="lattice" value="20000" />'
    pattern = '<steps units="lattice" value="2001" />'
    fh, abs_path = tempfile.mkstemp()
    with os.fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst)) 

    pattern = '<steps units="lattice" value="4001" />'
    fh, abs_path = tempfile.mkstemp()
    with os.fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))

    pattern = '<steps units="lattice" value="2000" />'
    fh, abs_path = tempfile.mkstemp()
    with os.fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))





if __name__=="__main__":
    t0 = time.time()
    time.sleep(5*60*60)

    # chmod 700 RunHemeLBSweepBash
    # subprocess.call("chmod 700 RunHemeLBSweepBash", shell=True)
    
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/IncreasingDensityWithUnchangingRadius/'
    MeshDirectory = "/data/vascrem/MeshCollection/IdealisedNetwork/IncreasingDensityScallingRadius/"
    if path.isdir(TerminalOutputFolder)==0:
        os.mkdir(TerminalOutputFolder)

    Scalling = ['2']     
    # Collapse = ['0','0.1227','0.2248', '0.3178', '0.417', '0.5124', '0.6119', '0.708', '0.8059', '0.9032', '1.0']
    Collapse = ['0.5124', '0.6119', '0.708', '0.8059', '0.9032', '1.0']

    dX1 = 0.08*0.75/41
    dX2 = 0.08*0.7/41

    dX3 = 0.08*0.51/41
    dX4 = 0.08*0.45/41

    Parallel = 3
    SleepyTime = 200
    AvaliablePaths = range(Parallel)
    for Level in Scalling:
        LevelsDirectory = MeshDirectory + Level
        if path.isdir(TerminalOutputFolder+'Density_'+Level)==0:
            os.mkdir(TerminalOutputFolder+'Density_'+Level)
        for i in Collapse:
            Core = AvaliablePaths[0] 

            if Level == '1' and Collapse == '0':
                continue
            
            mHemeLBDirectory = TerminalOutputFolder+'Density_'+Level+'/'+i+'/'
            if path.isdir(mHemeLBDirectory)==0:
                os.mkdir(mHemeLBDirectory)
            print mHemeLBDirectory
            if path.isdir(mHemeLBDirectory+'Results/')==1:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                os.rename(mHemeLBDirectory+'Results/',mHemeLBDirectory+'Results_'+str(current_time)+'/')
            ScaledMeshstl = LevelsDirectory+"/ScalledMesh"+str(i)+".stl"
            MeshFile = mHemeLBDirectory + 'config.stl' #MeshDirectory+""+Level+"/ScalledMesh"+i+".stl"
              #shutil.copyfile(ScaledMeshstl, mHemeLBDirectory + 'config.stl')
            MeshData = GetTheDetailsOfTheMesh(MeshFile)
            Boundaries =  MeshData[0:2]
            Seed = MeshData[2:5]
            
            if Level =='1':
                write_pr21(mHemeLBDirectory, 20000, dX1, Seed, Boundaries)
            elif Level =='2':
                write_pr22(mHemeLBDirectory, 20000, dX2, Seed, Boundaries)
            elif Level =='3':
                write_pr23(mHemeLBDirectory, 20000, dX3, Seed, Boundaries)
            elif Level =='4':
                write_pr24(mHemeLBDirectory, 20000, dX4, Seed, Boundaries)


            run_hemelb_setup(mHemeLBDirectory)

            # Update the xml file
            update_xml_file(int(20000*1), mHemeLBDirectory)
            # EditRunTime(mHemeLBDirectory)

            # # Run HemeLB
            RunHemeLB = 'mpirun -np 10 hemelb -in ' + mHemeLBDirectory+ 'config.xml -out '+mHemeLBDirectory +'Results/'
            TerminalOutput = mHemeLBDirectory+'HemeLBTerminalOutput.txt'
            # # Generate the new config.vtu
            GmyUnstructuredGridReader ="python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml "
            # Generate the flow vtus
            GenerateFlowVtus = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results/Extracted/surface-traction.xtr " + mHemeLBDirectory + "Results/Extracted/wholegeometry-velocity.xtr "
            # Generate waitFile
            WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
            print mHemeLBDirectory
            subprocess.Popen(['./RunHemeLBSweepBash', RunHemeLB, TerminalOutput, GmyUnstructuredGridReader, GenerateFlowVtus, WaitFileGeneration ])
        
        
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