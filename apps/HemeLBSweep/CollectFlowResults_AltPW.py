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



if __name__=="__main__":

    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/IncreasingAlternativePathways4/'
    # Scalling = ['PI_2.2',]   
    Scalling = ['4']   
    Collapse = [ '0', '0.3178', '0.417','0.5124', '0.6119', '0.708', '0.8059', '0.9032' , '0.1227', '0.2248'  , '1.0' ]
    
    NewNumbering = [ '0','1','2', '3','4','5', '6', '7','8','9','10']

    # Collapse = [ '0','1','9','10']
    HeadNewPath = TerminalOutputFolder+'CollectedResults/'
    print HeadNewPath
            
    for j in Scalling:
        NewPath = HeadNewPath+"Levels_" +j
        if path.isdir(NewPath)==0:
            os.mkdir(NewPath)
        counter =-1
        for i in Collapse:
            counter = counter+1
            print counter
            Pressure = TerminalOutputFolder+"Levels_"+j+'/'+i+'/Results/Extracted/surface-pressure_9500.vtu'
            Velocity = TerminalOutputFolder+"Levels_"+j+'/'+ i+'/Results/Extracted/wholegeometry-velocity_9500.vtu'
            # # Traction = TerminalOutputFolder+j+'/'+ i+'/Results/Extracted/surface-traction_25000.vtu'
            # Mesh = TerminalOutputFolder+"Level_"+j+'/'+ i+'/config.stl'

            NewPressure = NewPath+'/surface-pressure_'+NewNumbering[counter]+'.vtu'
            NewVelocity = NewPath+'/wholegeometry-velocity_'+NewNumbering[counter]+'.vtu'
            # NewTraction = NewPath+'/surface-traction_'+i+'.vtu'
            # NewMesh = NewPath+'/mesh_'+i+'.vtu'

            mv =  'cp ' +Pressure +' '+NewPressure
            subprocess.call(mv, shell=True)

            mv =  'cp ' +Velocity +' '+NewVelocity
            subprocess.call(mv, shell=True)

            # mv =  'cp ' +Mesh +' '+NewMesh
            # subprocess.call(mv, shell=True)

            # mv =  'cp ' +Traction +' '+NewTraction
            # subprocess.call(mv, shell=True)






    # TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/UpperBranchFolder/'
    
    # # Collapse = [ '0', '0.1227', '0.2248', '0.3178', '0.417', '0.5124', '0.6119', '0.708', '0.8059', '0.9032', '1.0']
    
    
    # Collapse = ['1','2','3','4','5','6','7','8','9','10','0','5.9','6','6.1','5.5','5.6','5.7','5.8','6.2','6.3','6.4','6.5']
    
    #     # Collapse =['0.5815', '0.5824', '0.5941', '0.5951' , '0.5961','0.5971' , '0.6098' , '0.6117', '0.6119', '0.6127', '0.6137', '0.6146','0.5924', '0.5932' ]
   
    # # NewNumbering = ['1','2','3','4','5','6','7','8','9','10','0','5.9','6','6.1','5.5','5.6','5.7','5.8','6.2','6.3','6.4','6.5']
    # Levels = ['7']
    # NewPath = TerminalOutputFolder + 'CollectedResults/'
    # if path.isdir(NewPath)==0:
    #     os.mkdir(NewPath)
    # for j in Levels:
    #     # NewPath = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/AngleVariation_3X3Network/CollectedResults/Alpha_'+j
    #     NewPath = TerminalOutputFolder + 'CollectedResults/Levels_'+j
    #     if path.isdir(NewPath)==0:
    #         os.mkdir(NewPath)
    #     counter =-1
    #     for i in Collapse:
    #         counter = counter+1
    #         print counter
    #         # Pressure = TerminalOutputFolder +'Levels_' +j+'/'+i+'/Results/Extracted/surface-pressure_3600.vtu'
    #         # Velocity = TerminalOutputFolder +'Levels_' +j+'/'+ i+'/Results/Extracted/wholegeometry-velocity_3600.vtu'
    #         Pressure = TerminalOutputFolder +'Levels_' +j+'/'+i+'/Results/Extracted/surface-pressure_9000.vtu'
    #         Velocity = TerminalOutputFolder +'Levels_' +j+'/'+ i+'/Results/Extracted/wholegeometry-velocity_9000.vtu'
    #         MeshFile = TerminalOutputFolder +'Levels_' +j+'/'+ i+'/config.stl'

    #         NewPressure = NewPath+'/surface-pressure_'+NewNumbering[counter]+'.vtu'
    #         NewVelocity = NewPath+'/wholegeometry-velocity_'+NewNumbering[counter]+'.vtu'
    #         NewMesh = NewPath+'/mesh_'+NewNumbering[counter]+'.stl'

    #         mv =  'cp ' +Pressure +' '+NewPressure
    #         subprocess.call(mv, shell=True)
    #         mv =  'cp ' +Velocity +' '+NewVelocity
    #         subprocess.call(mv, shell=True)
    #         mv =  'cp ' +MeshFile +' '+NewMesh
    #         subprocess.call(mv, shell=True)

    









    # t0 = time.time()

    # # chmod 700 RunFlowvtuBash
    # # Currently this code does not generate pr2 or xml files :S 
    # # subprocess.call("chmod 700 RunHemeLBCollapse", shell=True)
    # TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/AngleVariation_3X3Network/'
    # TerminalOutputFolder = '/Volumes/Hardrive/Projects/FSI/VaryingPressure/FlowFiles/Pressure_'


    # # Scalling = ['1']#'0.2', '0.4', '0.6', '0.8', '1.2', '1.4', '1.6', '1.8', '2', '2.2', '2.4', '2.6', '2.8', '3']   
    # Pressure = ['50', '200']
    # Collapse = [ '0.0','10.0','20.0', '30.0','40.0','50.0', '60.0', '70.0','80.0','90.0']
    # # Scalling = ['PI_3','PI_5']

    # NewPath = '/Volumes/Hardrive/Projects/FSI/VaryingPressure/FlowFiles2/'
    # if path.isdir(NewPath)==0:
    #     os.mkdir(NewPath)
             

    # for j in Pressure:
    #     # NewPath = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/AngleVariation_3X3Network/CollectedResults/Alpha_'+j
    #     NewPath = '/Volumes/Hardrive/Projects/FSI/VaryingPressure/FlowFiles2/Pressure_'+j

    #     if path.isdir(NewPath)==0:
    #         os.mkdir(NewPath)
    #     for i in Collapse:
    #         # i = str(float(I)*10)
    #         Pressure = TerminalOutputFolder +j+'/'+i+'/Results/Extracted/surface-pressure_2700.vtu'
    #         Velocity = TerminalOutputFolder +j+'/'+ i+'/Results/Extracted/wholegeometry-velocity_2700.vtu'
    #         I = str(float(i)/100)
    #         NewPressure = NewPath+'/surface-pressure_'+I+'.vtu'
    #         NewVelocity = NewPath+'/wholegeometry-velocity_'+I+'.vtu'

    #         mv =  'cp ' +Pressure +' '+NewPressure
    #         subprocess.call(mv, shell=True)

    #         mv =  'cp ' +Velocity +' '+NewVelocity
    #         subprocess.call(mv, shell=True)

 

 