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

    TerminalOutputFolder = '/data/vascrem/MeshCollection/IdealisedNetwork/VaryingLengthAndAngle/'
    # Scalling = ['PI_2.2',]   
    Scalling = ['4']   
    Collapse = [ '0', '0.3178', '0.417','0.5124', '0.6119', '0.708', '0.8059', '0.9032' , '0.1227', '0.2248'  , '1.0' ]
    Angle = [ 'PI_2_2' , 'PI_3' , 'PI_4' , 'PI_5' , 'PI_6' ]
    Length = [ '0.2' , '0.4' , '0.6' , '0.8', '1' ,  '1.2' , '1.4' , '1.6' , '1.8',  '2' ]

    
    NewNumbering = [ '0','1','2', '3','4','5', '6', '7','8','9','10']

    # Collapse = [ '0','1','9','10']
    HeadNewPath = TerminalOutputFolder+'CollectedResults/'
    if path.isdir(HeadNewPath)==0:
        os.mkdir(HeadNewPath)

    for k in   Angle:
        AnglePath = HeadNewPath + k
        if path.isdir(AnglePath)==0:
            os.mkdir(AnglePath)
        for j in Length:
            NewPath = AnglePath+"/HorizontalLength_" +j
            if path.isdir(NewPath)==0:
                os.mkdir(NewPath)
            for i in Collapse:
                # counter = counter+1
                # print counter
                Mesh = TerminalOutputFolder + k + "/HorizontalLength_" + j + "/ScaledMesh." + i + ".stl"
                NewMesh = NewPath+'/mesh_'+i+'.stl'

                mv =  'cp ' +Mesh +' '+NewMesh
                subprocess.call(mv, shell=True)
