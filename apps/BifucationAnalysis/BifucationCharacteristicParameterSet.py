#!/usr/bin/env python
# vasc_res_sim.py
import subprocess
import vtk
import shutil
import os
from xml.etree import ElementTree
import glob
from argparse import ArgumentParser
import numpy as np
from vmtk import vmtkscripts
import pdb


from os import path
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

def vmtk_compute_stl_radii(Input, output_filename, counter):
    devnull = open(os.devnull, 'w')
    command = 'vmtk vmtknetworkextraction -ifile ' + Input + ' -ofile ' + output_filename
    subprocess.call(command, shell=True)


    centerlineReader = vmtkscripts.vmtkSurfaceReader()
    centerlineReader.InputFileName = output_filename
    centerlineReader.Execute()

    clNumpyAdaptor = vmtkscripts.vmtkCenterlinesToNumpy()
    clNumpyAdaptor.Centerlines = centerlineReader.Surface
    clNumpyAdaptor.Execute()
    
    numpyCenterlines = clNumpyAdaptor.ArrayDict

    # pdb.set_trace()
    
    Radius = numpyCenterlines['PointData']['Radius']
    DataPoints = numpyCenterlines['Points']

    if np.size(Radius)<2:
        Radius = [float("NaN")]
        DataPoints = [float("NaN"), float("NaN"),float("NaN")]


    np.savetxt(OutputDirectory+"CenterPoints_"+counter+".txt", DataPoints, delimiter=' ')   # Save CenterPoints
    np.savetxt(OutputDirectory+"Radius_"+counter+".txt", Radius, delimiter=',')   # Save Radius


if __name__=="__main__":

    parser = ArgumentParser(description='AStuff')
    parser.add_argument('-Area', type=str, help='Need to supply a input folder')
    parser.add_argument('-Dil', type=str, help='Need to supply a input folder')
    parser.add_argument('-Shear', type=str, help='Need to supply a input folder')
    parser.add_argument('-Bending', type=str, help='Need to supply a input folder')
    parser.add_argument('-input', type=str, help='Need to supply a output folder')
    parser.add_argument('-output', type=str, help='Need to supply a output folder')
    args = parser.parse_args()

    i = args.Area
    j = args.Dil
    k =  args.Shear
    l =  args.Bending

    Folder = args.input
    NewDirectory = args.output

    Directory = Folder +"Area_"+ i+"_Dil_"+ j+"_Shear_"+ k+"_Bend_"+ l+"/results_from_time_0/"
    OutputDirectory = NewDirectory + "Area_"+ i+"_Dil_"+ j+"_Shear_"+ k+"_Bend_"+ l+"/"  
    MeshOutputDirectory =   Folder +"Meshes_Area_"+ i+"_Dil_"+ j+"_Shear_"+ k+"_Bend_"+ l+"/"

    if path.isdir(OutputDirectory)==0:
        os.mkdir(OutputDirectory)

    if path.isdir(MeshOutputDirectory)==0:
        os.mkdir(MeshOutputDirectory)

    Number = np.array([ ])

    for file in os.listdir(Directory):
        if file.endswith(".vtu") & file.startswith("mesh"):
            Numb =file[5:-4]

            ConvertToSt = "meshio-convert " + Directory +  file + " " + MeshOutputDirectory + "mesh_"+ Numb +".stl"
            CenterlinesFile = OutputDirectory + "Centerlines_"+ Numb +".vtp"
            subprocess.call(ConvertToSt, shell=True)
            vmtk_compute_stl_radii(MeshOutputDirectory + "mesh_"+Numb+".stl",CenterlinesFile, Numb)

            Number = np.append(Number, [int(Numb)])
    Number = np.sort(Number, axis=0) 
    np.savetxt(OutputDirectory+"WriteOutTimes.txt", Number, delimiter=' ')   # Save Radius

            




    
    

                


            # 1) Load the folder 
            # 2) Loop over all exisiting vtu files and convert to stl 
            # 3) Turn all stls into centerlines files 
            # 4) Try find the radi in the centerlines files 
            #   a) get the max radius
        #   b) get the radius at the end of each branch
        #   c) get gradrent of radius over first branch
        #  

        #
        # Step 5: Convert Chaste output (vtu file) into the input of the next iteration (stl file required by hemelb setup tool)
        # #
        # vtu2stl()

        # if args.compute_radii:
        #     vmtk_compute_stl_radii(iter+1)

        
            
