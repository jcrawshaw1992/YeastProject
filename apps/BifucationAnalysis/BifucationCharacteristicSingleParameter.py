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
    parser.add_argument('-Shear', type=str, help='Need to supply a input folder')
    parser.add_argument('-input', type=str, help='Need to supply a output folder')
    parser.add_argument('-output', type=str, help='Need to supply a output folder')
    args = parser.parse_args()

    k =  args.Shear
    print k
    Folder = args.input
    NewDirectory = args.output


    AreaParameter = [6]
    DilationParameter =  [6.5]
    # DeformationParamter = [5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    BendingParamter = [7,8,9,10,11,12,13,14,15,16]


    for i in AreaParameter:
        for j in DilationParameter:
            for l in BendingParamter:
                Directory = Folder +"Area_"+ str(i)+"_Dil_"+ str(j)+"_Shear_"+ k+"_Bend_"+ str(l)+"/results_from_time_0/"
                OutputDirectory = NewDirectory + "Area_"+ str(i)+"_Dil_"+ str(j)+"_Shear_"+ k+"_Bend_"+ str(l)+"/"            
            
                if path.isdir(OutputDirectory)==0:
                    os.mkdir(OutputDirectory)
        
                Number = np.array([ ])

                for file in os.listdir(Directory):
                    if file.endswith(".vtu") & file.startswith("mesh"):
                        Numb =file[5:-4]

                        ConvertToSt = "meshio-convert " + Directory +  file + " " + OutputDirectory + "mesh_"+ Numb +".stl"
                        CenterlinesFile = OutputDirectory + "Centerlines_"+ Numb +".vtp"
                        subprocess.call(ConvertToSt, shell=True)
                        
                        vmtk_compute_stl_radii(OutputDirectory + "mesh_"+Numb+".stl",CenterlinesFile, Numb)

                        # Numb = np.array([Numb])
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

        
            
