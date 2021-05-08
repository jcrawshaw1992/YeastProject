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


def run_hemelb_setupParrallel( working_directory):
    command = 'python RunSetup.py  -working_directory ' + mHemeLBDirectory
    subprocess.Popen(command, shell=True)

    
def run_hemelb_setup(working_directory):
    print "HemeLB setting up"
    heme_profile_filename = working_directory + 'config.pr2' 
    command = hemelb_setup_exe + ' ' + heme_profile_filename 
    subprocess.call(command, shell=True)
    
    print "HemeLB set up is DONE"



def EditTimeStep(working_directory):

    # Load automatically generated XML file
    filename = working_directory + 'config.xml'

    # Edit the initial condition to my satisfaction

    file = open(filename).readlines()
    #Create temp file
    file_path = filename
   
    pattern = '<step_length units="s" value="9.51814396193e-09" />'
    subst = '<step_length units="s" value="9.51814396193e-08" />'

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
        



def write_pr2(outputDirectory, SimulationDuration, MinRadii, Seed, Boundaries, Pressure):

    f = open(outputDirectory+"config.pr2", "w")

    V =  4 #Kinematic viscosity -- 4 X 10^-6 m^2/s  V = eta/rho Here v needs to be in the same dims as the input file!!!
    deltaX = MinRadii
    deltaT = float(0.1 * deltaX * deltaX/V)
    print ('DeltaX: ' , deltaX , 'DeltaT: ', deltaT, ' eta: ', V )

    InletPressure = Pressure;#100
    OutletPressure = 0
    Duration = SimulationDuration*deltaT

    f.write("DurationSeconds: "+ str(Duration) +"\n")
    f.write("Iolets:\n"+ \
    # ------Inlets ----
    "- Centre: {x: "+str(Boundaries[0])+", y: 0.000000, z: 0.000000}\n"+ \
    "  Name: Inlet1\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(Boundaries[0])+", y: -0.27, z: 0.000000}\n"+ \
    "  Name: Inlet2\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Inlet\n"+ \
    # ------
    "- Centre: {x: "+str(Boundaries[0])+", y: 0.27, z: 0.000000}\n"+ \
    "  Name: Inlet3\n"+ \
    "  Normal: {x: 1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(InletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Inlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(Boundaries[1])+", y: -0.27, z: 0.000000}\n"+ \
    "  Name: Outlet1\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Outlet\n"+ \
    # -----Outlets -----
    "- Centre: {x: "+str(Boundaries[1])+", y: 0, z: 0.000000}\n"+ \
    "  Name: Outlet2\n"+ \
    "  Normal: {x: -1.000000, y: 0.000000, z: 0.000000}\n"+ \
    "  Pressure: {x: "+str(OutletPressure)+", y: 0.0, z: 0.0}\n"+ \
    "  Radius: 0.5\n"+ \
    "  Type: Outlet\n"+ \
    # ------
    "- Centre: {x: "+str(Boundaries[1])+", y: 0.270, z: 0.000000}\n"+ \
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



def EditPressure(working_directory):

    print "update_xml_file"
    # Load automatically generated XML file
    filename = working_directory + 'config.xml'

    # Edit the initial condition to my satisfaction

    file = open(filename).readlines()
    #Create temp file
    file_path = filename
   
    pattern = '<mean units="mmHg" value="200" />'
    subst = '<mean units="mmHg" value="100" />'

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


def EditRunTime(working_directory):

    print "update_xml_file"
    # Load automatically generated XML file
    filename = working_directory + 'config.xml'

    # Edit the initial condition to my satisfaction

    file = open(filename).readlines()
    #Create temp file
    file_path = filename
   
    pattern = '<steps units="lattice" value="3000" />'
    subst = '<steps units="lattice" value="10000" />'

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




if __name__=="__main__":
    t0 = time.time()
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/PressureVariation/'

    Collapse = ['3', '4','5', '6', '7', '8', '9', '10', '0', '1', '2']
    PressureSet = [50,100,200]

    Parallel =1
    SleepyTime = 200
    AvaliablePaths = range(Parallel)
    Setup = 0 
    if Setup == 1:
        for Pressure in PressureSet:
            if path.isdir(TerminalOutputFolder+'Pressure_'+str(Pressure))==0:
                os.mkdir(TerminalOutputFolder+'Pressure_'+str(Pressure))          
            for i in Collapse:
                Core = AvaliablePaths[0] 
                mHemeLBDirectory = TerminalOutputFolder+'Pressure_'+str(Pressure)+'/'+i+'/'
                if path.isdir(mHemeLBDirectory)==0:
                        os.mkdir(mHemeLBDirectory)
                print mHemeLBDirectory
                if path.isdir(mHemeLBDirectory+'Results/')==1:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    os.rename(mHemeLBDirectory+'Results/',mHemeLBDirectory+'Results_'+str(current_time)+'/')

                MeshFile = mHemeLBDirectory + 'config.stl'
                MeshData = GetTheDetailsOfTheMesh(MeshFile)
                Boundaries =  MeshData[0:2]
                Seed = MeshData[2:5]
        
                dX = 0.08/41.0
                write_pr2(mHemeLBDirectory, 3000, dX, Seed, Boundaries, Pressure)
                print "Pre"
                run_hemelb_setupParrallel(TerminalOutputFolder, mHemeLBDirectory, Core)
        
                WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
                subprocess.Popen('python RunSetup.py  -working_directory ' + mHemeLBDirectory, shell=True)

                AvaliablePaths.remove(Core) 
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
            
    Run = 0 
    if Run ==1:
        Parallel = 5
        SleepyTime = 100
        AvaliablePaths = range(Parallel-1)
        for Pressure in PressureSet: 
            for i in Collapse:
                Core = AvaliablePaths[0] 
                mHemeLBDirectory = TerminalOutputFolder+'Pressure_'+str(Pressure)+'/'+i+'/'
                update_xml_file(int(10000*0.5), mHemeLBDirectory)
                EditRunTime(mHemeLBDirectory)
                print mHemeLBDirectory
                RunHemeLB = 'mpirun -np 6 hemelb -in ' + mHemeLBDirectory+ 'config.xml -out '+mHemeLBDirectory +'Results/'
                TerminalOutput = mHemeLBDirectory+'HemeLBTerminalOutput.txt'

                WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
                subprocess.Popen(['./RunHemeLBSweepBash', RunHemeLB, TerminalOutput, " ",  " ", WaitFileGeneration ])
        
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


        OutputFile = TerminalOutputFolder+'WaitFile'+str(Parallel-1)+'.txt'
        if (os.path.isfile(OutputFile)==0):
            print 'Not exisit'
            time.sleep(SleepyTime)
            while (~path.exists(OutputFile)):
                time.sleep(SleepyTime)

        for P in range(Parallel):
            OutputFile = TerminalOutputFolder+'WaitFile'+str(P)+'.txt'
            if path.exists(OutputFile):
                os.remove(OutputFile)



    WriteOutVtus = 1
    if WriteOutVtus ==1:

        Parallel =31
        AvaliablePaths = range(Parallel)
        for Pressure in PressureSet: 
            for i in Collapse:
                Core = AvaliablePaths[0] 
                mHemeLBDirectory = TerminalOutputFolder+'Pressure_'+str(Pressure)+'/'+i+'/'
                TerminalOutput = mHemeLBDirectory+'HemeLBTerminalOutput.txt'
                GmyUnstructuredGridReader ="python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml "
                GenerateFlowVtus = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results/Extracted/surface-traction.xtr "  + mHemeLBDirectory + "Results/Extracted/surface-pressure.xtr " + mHemeLBDirectory + "Results/Extracted/wholegeometry-velocity.xtr "
                print mHemeLBDirectory
                WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
                subprocess.Popen(['./RunHemeLBSweepBash', " ", TerminalOutput, GmyUnstructuredGridReader, GenerateFlowVtus, WaitFileGeneration])


#  3319  sudo du -h | sort -rh | head -5
#
