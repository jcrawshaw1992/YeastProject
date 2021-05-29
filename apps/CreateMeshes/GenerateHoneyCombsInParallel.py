#!/usr/bin/env python
#
import subprocess
import vtk
import os
import pdb
import string
import math as m
import vmtk
from vmtk import pypes
import array as A
# from vtk import vtkXMLUnstructuredGridWriter as Writer, VTK_VERSION
# from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
import numpy as np
import meshio
from vtk.util import numpy_support
import vtk
import clip
import math
import ConvertVTKtoSTL
from os import path
import time
 

def CreateNewFolder(Directory):
    if path.isdir(Directory)==0:
        os.mkdir(Directory)
        
def Waiting(AvaliablePaths,SleepyTime,Parallel, TerminalOutputFolder, t0): 
    while len(AvaliablePaths) ==0:
        time.sleep(SleepyTime)
        print "Sleep Time"
        for P in range(Parallel):
            OutputFile = TerminalOutputFolder+'WaitFile'+str(P)+'.txt'
            if path.exists(OutputFile):
                AvaliablePaths.append(P)
                os.remove(OutputFile)
        if len(AvaliablePaths) >0:
            print AvaliablePaths, "Have found a spare core or two :-) " 
            print time.time() - t0, "seconds time"

    
        

if __name__=="__main__":
    t0 = time.time()
    Directory = "/data/vascrem/MeshCollection/IdealisedNetwork/VaryingLengthAndAngleNew/"
    CreateNewFolder(Directory) 
    TerminalOutputFolder = Directory + "TerminalOutputFolder/"
    CreateNewFolder(TerminalOutputFolder) 

    Parallel = 25
    SleepyTime = 50
    AvaliablePaths = range(Parallel)
        
    FileLabels = ['PI_2_2', 'PI_3', 'PI_4', 'PI_5', 'PI_6'] 
    Alpha = [m.pi/2.2, m.pi/3, m.pi/4,m.pi/5, m.pi/6] 

    # FileLabels = [ 'PI_5'] 
    # Alpha = [m.pi/5] 
    Length = [ 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2]
    Collapse = [0.2248, 0.3178, 0.4170 , 0.5124, 0.6119, 0.7080, 0.8059, 0.9032]
    Collapse = [0.4170]
    counter = -1
    for A in Alpha:
        counter = counter+1
        AngleDirectory = Directory+FileLabels[counter]+'/'
        CreateNewFolder(AngleDirectory) 
        Lengthcounter = -1
        for L in Length:
            LengthDirectory = AngleDirectory+'HorizontalLength_'+str(L)+'/'
            CreateNewFolder(LengthDirectory) 
    
            for i in Collapse:
                Core = AvaliablePaths[0] 
                AvaliablePaths.remove(Core) 

                WaitFileGeneration = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
                Command = 'python CreateHoneycombMesh.py -Directory ' + LengthDirectory + ' -Collapse ' +str(i) + ' -Length '+str(L) + ' -Angle ' + str(A)  + ' -WaitFileGeneration ' + WaitFileGeneration 
                subprocess.Popen(Command, shell=True)

                Waiting(AvaliablePaths,SleepyTime,Parallel, TerminalOutputFolder, t0)