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



# Start with cylinder 

working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/HetrogeneousCylinder/CenterCollapse' 
data_path = working_directory  + '/SetUpData/'
chaste_setup_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/TestSetupMeshRunner'
# chaste_run_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/TestRunFlowInPipeRunner'


print "Setting Up Mesh"

   
if __name__=="__main__":


# Create the directory 
    if not os.path.exists(data_path):
        os.makedirs(data_path)

# **** Set up what you need for the collected chaste mesh folder
  
    
 
    os.environ['CHASTE_TEST_OUTPUT'] = working_directory
    
    #
    # Step 1: Run preliminary Chaste setup
    #
    vtu_filename = working_directory +'/config.vtu' # Use proper path concatentation


    print '********** ------- Setting up Chaste simulation ------- **********'
 
    chaste_setup_call = chaste_setup_exe #+ ' -mesh ' + mesh_filename + ' -xml ' + xml_filename + ' -div_threshold ' +  str(args.div_threshold) + ' -mesh_scale ' +  str(args.mesh_scale) + ' -wss_threshold ' + str(0.025)
    subprocess.call(chaste_setup_call, shell=True)

    print '********** ------- Chaste is set up ------- **********'
    shutil.copy( vtu_filename , data_path+'config.vtu'  )
  
    
    # Step 1.a: Prepare mesh for HemeLB setup tool and compute radii along the axes of the domain (if requested)
    
    

    print "**************   Convert vtu to stl   **************  "
    # Read the VTU file from disk
    vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    vtu_reader.SetFileName(data_path+ 'config.vtu')
    # cyl_581_nodes

    # vtkAppendFilter appends one or more datasets together into a single unstructured grid
    extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

    # Write out the data in unstructured grid format
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    stl_writer.SetFileName(data_path + 'config.stl')
    print data_path + 'config.stl' 
    stl_writer.Write()

    # shutil.copy( xml_filename , working_directory +'/config.xml'  )
  


    print "**************   vtu and stl created   **************  "

    #  shutil.rmtree(path)





        
    print '\n ********* ------ Completed ------ ********* \n'

# 














