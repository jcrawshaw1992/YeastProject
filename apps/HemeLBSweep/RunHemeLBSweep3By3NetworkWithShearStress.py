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
from datetime import datetime

def EditTimeStep(working_directory):

    print "update_xml_file"
    # Load automatically generated XML file
    filename = working_directory + 'config.xml'

    # Edit the initial condition to my satisfaction

    file = open(filename).readlines()
    #Create temp file
    file_path = filename
   
    pattern = '<step_length units="s" value="9.51814396193e-08" />'
    subst = '<step_length units="s" value="9.51814396193e-09" />'

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

 
    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'surface-shearstress.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='shearstress') # vonmisesstress

    # surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'surface-vonmisesstress.xtr'})
    # ElementTree.SubElement(surface, 'geometry', type='surface')
    # ElementTree.SubElement(surface, 'field', type='vonmisesstress')

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



def write_pr2(outputDirectory, SimulationDuration, MinRadii):

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
    "SeedPoint: {x: 0.267102, y: 0.118712, z: -0.0160343}\n"+
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
   
    pattern = '<steps units="lattice" value="8000" />'
    subst = '<steps units="lattice" value="7000" />'

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

   

    #Copy the file permissions from the old file to the new file
    shutil.copymode(file_path, abs_path)
    #Remove original file
    os.remove(file_path)
    #Move new file
    shutil.move(abs_path, file_path)


if __name__=="__main__":
    t0 = time.time()

    # chmod 700 RunHemeLBSweepBash
    # subprocess.call("chmod 700 RunHemeLBCollapse", shell=True)
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3CollapseWithShearStress/'
    Collapse =['0.5924' ]#, '0.5941', '0.5951' , '0.5961','0.5971' , '0.6098' , '0.6117', '0.6119', '0.6127', '0.6137', '0.6146']
 
    UpperBranchFolder = TerminalOutputFolder + 'UpperBranchFolder/AroundUpperBranchAround6/'
    if path.isdir(UpperBranchFolder)==0:
        os.mkdir(UpperBranchFolder)



    Setup = 0
    if Setup ==1:
        Parallel = 32
        SleepyTime = 10
        AvaliablePaths = range(Parallel)
        print AvaliablePaths
        for i in Collapse:
            Core = AvaliablePaths[0] 
        
            mHemeLBDirectory = UpperBranchFolder+i+'/'
           
            print mHemeLBDirectory
            if path.isdir(mHemeLBDirectory+'Results/')==1:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                os.rename(mHemeLBDirectory+'Results/',mHemeLBDirectory+'Results_'+str(current_time)+'/')

            dX = 0.08/41.0
            write_pr2(mHemeLBDirectory, 10000, dX)

            WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
            Command = 'python RunSetup.py  -working_directory ' + mHemeLBDirectory + ' -WaitFileGeneration ' + WaitFileGeneration
            subprocess.Popen(Command, shell=True)

    # Collapse =['0.5815', '0.5824', '0.5941', '0.5951' , '0.5961','0.5971' , '0.6098' , '0.6117', '0.6119', '0.6127', '0.6137', '0.6146','0.5924', '0.5932' ]
   
    # Parallel = 3
    # SleepyTime = 100
    # AvaliablePaths = range(Parallel)
    # print AvaliablePaths
    # for i in Collapse:
    #     Core = AvaliablePaths[0] 
    #     mHemeLBDirectory = UpperBranchFolder+i+'/'
    #     if path.isdir(mHemeLBDirectory)==0:
    #         os.mkdir(mHemeLBDirectory)
    #     if path.isdir(mHemeLBDirectory+'Results/')==1:
    #         now = datetime.now()
    #         current_time = now.strftime("%H:%M:%S")
    #         os.rename(mHemeLBDirectory+'Results/',mHemeLBDirectory+'Results_'+str(current_time)+'/')
              
    #     # dX = 0.08/41.0
    
    #     # write_pr2(mHemeLBDirectory, 4001, dX)
    #     # run_hemelb_setup(mHemeLBDirectory)

    #     # # Update the xml file
    #     update_xml_file(int( 10000*0.95), mHemeLBDirectory)
    #     EditRunTime(mHemeLBDirectory)

    #     print "About to run HemeLB"
    #     # # Run HemeLB
    #     RunHemeLB = 'mpirun -np 10 hemelb -in ' + mHemeLBDirectory+ 'config.xml -out '+mHemeLBDirectory +'Results/'
 
    #     TerminalOutput = mHemeLBDirectory+'HemeLBTerminalOutput.txt'
    #     # # Generate the new config.vtu
    #     GmyUnstructuredGridReader ="python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml "
    #     # Generate the flow vtus
    #     GenerateFlowVtus = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results/Extracted/surface-traction.xtr " + mHemeLBDirectory + "Results/Extracted/wholegeometry-velocity.xtr "
  
  
    #     WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
    #     subprocess.Popen(['./RunHemeLBSweepBash', RunHemeLB, TerminalOutput, GmyUnstructuredGridReader, GenerateFlowVtus, WaitFileGeneration ])
    
    
    #     AvaliablePaths.remove(Core) 
    #     # Check if all positions are taken
    #     while len(AvaliablePaths) ==0:
    #         time.sleep(SleepyTime)
    #         # print "Awake and checking for spare cores" 
    #         print "Sleep Time"
    #         for P in range(Parallel):
    #             OutputFile = TerminalOutputFolder+'WaitFile'+str(P)+'.txt'
    #             if path.exists(OutputFile):
    #                 AvaliablePaths.append(P)
    #                 os.remove(OutputFile)
    #         if len(AvaliablePaths) >0:
    #             print AvaliablePaths, "Have found a spare core or two :-) " 
    #             print len(AvaliablePaths)
    #             print time.time() - t0, "seconds time"




# -----------------------------------------------------------------------------
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3CollapseWithShearStress/'
    
    UpperBranchFolder = TerminalOutputFolder + 'UpperBranchFolder/'
    Collapse = ['1']#'4','2','3','5','6','7','8','9','10','0','5.9','6','6.1','5.5','5.6','5.7','5.8','6.2','6.3','6.4','6.5']
    
    Parallel = 1
    SleepyTime = 100
    AvaliablePaths = range(Parallel)
    print AvaliablePaths
    for i in Collapse:
        Core = AvaliablePaths[0] 
        mHemeLBDirectory = UpperBranchFolder+i+'/'
        if path.isdir(mHemeLBDirectory)==0:
            os.mkdir(mHemeLBDirectory)
        if path.isdir(mHemeLBDirectory+'Results/')==1:
            # shutil.rmtree(mHemeLBDirectory+'Results/', ignore_errors=False, onerror=None)
            # now = datetime.now()
            # current_time = now.strftime("%H:%M:%S")
            # os.rename(mHemeLBDirectory+'Results/',mHemeLBDirectory+'Results_'+current_time+'/')

        # EditRunTime(mHemeLBDirectory)

        # update_xml_file(7000*0.99, mHemeLBDirectory)

        print "About to run HemeLB"
        # # Run HemeLB
        # Copiedxml = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3CollapseWithShearStress/UpperBranchFolder/1/config.xml'
        # NewXml  =mHemeLBDirectory+ 'config.xml'
        # mv =  'cp ' +Copiedxml +' '+NewXml
        # subprocess.call(mv, shell=True)

        RunHemeLB = 'mpirun -np 31 hemelb -in ' + mHemeLBDirectory+ 'config.xml -out '+mHemeLBDirectory +'Results/'
        
        TerminalOutput = mHemeLBDirectory+'HemeLBTerminalOutput.txt'
        # # Generate the new config.vtu
        GmyUnstructuredGridReader =" " #python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml "
        # Generate the flow vtus                                                                                                                                                                           surface-shearstress.xtr                                                     
        GenerateFlowVtus = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results/Extracted/surface-shearstress.xtr " + mHemeLBDirectory + "Results/Extracted/wholegeometry-velocity.xtr "
  
  
        WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
        subprocess.Popen(['./RunHemeLBSweepBash', RunHemeLB, TerminalOutput, GmyUnstructuredGridReader, GenerateFlowVtus, WaitFileGeneration ])
    
    
        # AvaliablePaths.remove(Core) 
        # # Check if all positions are taken
        # while len(AvaliablePaths) ==0:
        #     time.sleep(SleepyTime)
        #     # print "Awake and checking for spare cores" 
        #     print "Sleep Time"
        #     for P in range(Parallel):
        #         OutputFile = TerminalOutputFolder+'WaitFile'+str(P)+'.txt'
        #         if path.exists(OutputFile):
        #             AvaliablePaths.append(P)
        #             os.remove(OutputFile)
        #     if len(AvaliablePaths) >0:
        #         print AvaliablePaths, "Have found a spare core or two :-) " 
        #         print len(AvaliablePaths)
        #         print time.time() - t0, "seconds time"

# -----------------------------------------------------------------------------

    # MiddleBranchFolder = TerminalOutputFolder + 'MiddleBranchFolder/'
    # Collapse = ['1','2','3','4','5','6','7','8','9','10','0','5.9','6','6.1','5.5','5.6','5.7','5.8','6.2','6.3','6.4','6.5']
    
    # Parallel = 10
    # SleepyTime = 100
    # AvaliablePaths = range(Parallel)
    # print AvaliablePaths
    # for i in Collapse:
    #     Core = AvaliablePaths[0] 
    #     mHemeLBDirectory = MiddleBranchFolder+i+'/'
    #     if path.isdir(mHemeLBDirectory)==0:
    #         os.mkdir(mHemeLBDirectory)
    #     if path.isdir(mHemeLBDirectory+'Results/')==1:
    #         # now = datetime.now()
    #         # current_time = now.strftime("%H:%M:%S")
    #         os.rename(mHemeLBDirectory+'Results/',mHemeLBDirectory+'Results_WallShearStress/')


    #     Copiedxml = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3CollapseWithShearStress/UpperBranchFolder/1/config.xml'
    #     NewXml  =mHemeLBDirectory+ 'config.xml'
    #     mv =  'cp ' +Copiedxml +' '+NewXml
    #     subprocess.call(mv, shell=True)


    #     print "About to run HemeLB"
       
    #     RunHemeLB = 'mpirun -np 3 hemelb -in ' + mHemeLBDirectory+ 'config.xml -out '+mHemeLBDirectory +'Results/'
        
 
    #     TerminalOutput = mHemeLBDirectory+'HemeLBTerminalOutput.txt'
    #     # # Generate the new config.vtu
    #     GmyUnstructuredGridReader =" " #python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml "
    #     # Generate the flow vtus
    #     GenerateFlowVtus = "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "Results/Extracted/surface-shearstress.xtr " + mHemeLBDirectory + "Results/Extracted/wholegeometry-velocity.xtr "
    
    #     WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
    #     subprocess.Popen(['./RunHemeLBSweepBash', RunHemeLB, TerminalOutput, GmyUnstructuredGridReader, GenerateFlowVtus, WaitFileGeneration ])
    
    
    #     AvaliablePaths.remove(Core) 
    #     # Check if all positions are taken
    #     while len(AvaliablePaths) ==0:
    #         time.sleep(SleepyTime)
    #         # print "Awake and checking for spare cores" 
    #         print "Sleep Time"
    #         for P in range(Parallel):
    #             OutputFile = TerminalOutputFolder+'WaitFile'+str(P)+'.txt'
    #             if path.exists(OutputFile):
    #                 AvaliablePaths.append(P)
    #                 os.remove(OutputFile)
    #         if len(AvaliablePaths) >0:
    #             print AvaliablePaths, "Have found a spare core or two :-) " 
    #             print len(AvaliablePaths)
    #             print time.time() - t0, "seconds time"