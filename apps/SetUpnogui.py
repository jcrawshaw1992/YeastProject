import subprocess
import vtk
import shutil
import os
from xml.etree import ElementTree
import glob
from argparse import ArgumentParser
import numpy as np
import time
import matplotlib.pyplot as plt
import csv
import pdb
import string
import math

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")


t0 = time.time()
hemelb_setup_exe = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools:/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'
hemelb_setup_exe1 = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'



# voxel_size =  0.25e-4 #0.1953e-06
# subfolder = [ 'Growing_1_5','Growing_2', 'Growing_2_5']#, 'Growing_3', 'Growing_3_5', 'Growing_3', 'Growing_3_5', 'Growing_4', 'Growing_4_5', 'Growing_5']

# for i in subfolder :
#     working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/TestingShapes1/' + i +'/' #TestCylinderRadi_voxal/NotScalledUp/Radi_e-3/Voxal1.4e-8/' 

#     print('here')
#     # _______________________
#     #  COnvert vtu to the stl i need 
#     #  _______________________



#     # print "  Convert vtu to stl    "
#     #     # Read the VTU file from disk
#     # vtu_reader = vtk.vtkXMLUnstructuredGridReader()
#     # vtu_reader.SetFileName(working_directory + 'config.vtu')

#     # extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
#     # extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

#     # # Write out the data in unstructured grid format
#     # stl_writer = vtk.vtkSTLWriter()
#     # stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
#     # stl_writer.SetFileName(working_directory + 'config.stl')
#     # stl_writer.Write()


#     # _______________________
#     # Set up the inital config.xml file 
#     # _______________________


#     print "HemeLB setting up"
#     # command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/GmyUnstructuredGridReader.py ~/Documents/ChasteWorkingDirectory/config.xml'
#     # subprocess.call(command, shell=True)

#     # TODO: It seems that when you use the --stl option below, the voxel size specified in .pr2 is ignored and the setup tool guesses one for you



#     # VoxelSize = 0.48330092430114746

#     heme_profile_filename = working_directory + 'config.pr2' 
#     ConfigDirectory = working_directory + 'config.stl'
#     # I think the problem is in this command 
#     # command = hemelb_setup_exe1 + ' ' + heme_profile_filename + ' --voxel ' + str(voxel_size)  + ' --geometry ' + working_directory + '/config.gmy' + ' --xml ' + working_directory + '/config.xml'

#     command = hemelb_setup_exe1 + ' ' + heme_profile_filename + ' --voxel ' + str(voxel_size) # + ' --geometry ' + working_directory + '/config.gmy' + ' --xml ' + working_directory + '/config.xml'

#     subprocess.call(command, shell=True)
#     print "HemeLB set up is DONE"

#     # pause()

#     print "update_xml_file"


#     # pause()

#     # _______________________
#     # Now update the xml file 
#     # _______________________


#     # Load automatically generated XML file
#     filename = working_directory + 'config.xml'

#     tree = ElementTree.parse(filename)
#     root = tree.getroot()

#     # Add monitoring of incompressibility and convergence
#     monitoring = ElementTree.SubElement(root, 'monitoring')
#     ElementTree.SubElement(monitoring, 'incompressibility') 
#     convergence = ElementTree.SubElement(monitoring, 'steady_flow_convergence', {'tolerance': '1e-3', 'terminate': 'false'})
#     ElementTree.SubElement(convergence, 'criterion', {'type': 'velocity', 'value': '0.01', 'units': 'm/s'})

#     # Add definition of properties to be extracted
#     extr = ElementTree.SubElement(root, 'properties')

#     surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-tangentialprojectiontraction.xtr'})
#     ElementTree.SubElement(surface, 'geometry', type='surface')
#     ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')

#     surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-traction.xtr'})
#     ElementTree.SubElement(surface, 'geometry', type='surface')
#     ElementTree.SubElement(surface, 'field', type='traction')

#     surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-tractions.xtr'})
#     ElementTree.SubElement(surface, 'geometry', type='surface')
#     ElementTree.SubElement(surface, 'field', type='traction')
#     ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')    

#     surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-pressure.xtr'})
#     ElementTree.SubElement(surface, 'geometry', type='surface')
#     ElementTree.SubElement(surface, 'field', type='pressure')

#     wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'wholegeometry-velocity.xtr'})
#     ElementTree.SubElement(wholegeometry, 'geometry', type='whole')
#     ElementTree.SubElement(wholegeometry, 'field', type='velocity')

#     # Save XML file to disk
#     tree.write(filename)

#     print('here')

#     print "done"


#     # _______________________
#     # Run the HemeLB simulation  
#     # _______________________

#     print "**************   hemelb step  **************   "
#     shutil.rmtree(working_directory + 'results/', ignore_errors=True)
#     xml_config_file = working_directory + 'config.xml'
#     command = 'mpirun -np 1 hemelb -in ' + xml_config_file
#     subprocess.call(command, shell=True)
#     # mpirun -np 1 hemelb -in "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/config.xml"


#     # _______________________
#     # Generate the flow files  
#     # _______________________


#     print "  generate_flow_vtus   "
#         # Make the vtu 
#     command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/GmyUnstructuredGridReader.py' + ' ' + working_directory + 'config.xml'
#     subprocess.call(command, shell=True)

#     # this need to be generalised 
#     command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py' + ' ' +  working_directory +'config.vtu' +' ' +  working_directory+ 'results/Extracted/surface-pressure.xtr' +' ' +  working_directory+ 'results/Extracted/wholegeometry-velocity.xtr'+' ' +  working_directory+ 'results/Extracted/surface-traction.xtr'+' ' +  working_directory+ 'results/Extracted/surface-tangentialprojectiontraction.xtr'
#     subprocess.call(command, shell=True)

#         # output_folders = glob.glob(working_directory + '/results_from_time_*')

#         # for folder in output_folders:
#         #     folder_timestep = folder.split('results_from_time_')[-1]
#         #     if float(folder_timestep) == timestep:
#         #         shutil.copyfile(working_directory + '/'+ 'results/Extracted/surface-pressure_2000.vtu', folder + '/surface-pressure.vtu')
#         #         shutil.copyfile(working_directory + '/'+ 'results/Extracted/surface-traction_2000.vtu', folder + '/surface-traction.vtu')
#         #         shutil.copyfile(working_directory + '/'+ 'results/Extracted/wholegeometry-velocity_2000.vtu', folder + '/wholegeometry-velocity.vtu')
#         #         shutil.copyfile(working_directory + '/'+ 'results/Extracted/surface-tangentialprojectiontraction_2000.vtu', folder + '/surface-tangentialprojectiontraction.vtu')
    
#     print "generate_flow_vtus DONE    "
#     t1 = time.time()

#     total = t1-t0
#     print "Time taken for simulaiton + " + str(total)





voxel_size =  0.390625e-6 #0.1953e-06
subfolder = [ 'Growing_1_5','Growing_2', 'Growing_2_5']#, 'Growing_3', 'Growing_3_5', 'Growing_3', 'Growing_3_5', 'Growing_4', 'Growing_4_5', 'Growing_5']

for i in subfolder :
    working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/TestingShapes3/' + i +'/' #TestCylinderRadi_voxal/NotScalledUp/Radi_e-3/Voxal1.4e-8/' 

    print('here')
    # _______________________
    #  COnvert vtu to the stl i need 
    #  _______________________



    # print "  Convert vtu to stl    "
    #     # Read the VTU file from disk
    # vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    # vtu_reader.SetFileName(working_directory + 'config.vtu')

    # extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    # extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

    # # Write out the data in unstructured grid format
    # stl_writer = vtk.vtkSTLWriter()
    # stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    # stl_writer.SetFileName(working_directory + 'config.stl')
    # stl_writer.Write()


    # _______________________
    # Set up the inital config.xml file 
    # _______________________


    print "HemeLB setting up"
    # command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/GmyUnstructuredGridReader.py ~/Documents/ChasteWorkingDirectory/config.xml'
    # subprocess.call(command, shell=True)

    # TODO: It seems that when you use the --stl option below, the voxel size specified in .pr2 is ignored and the setup tool guesses one for you



    # VoxelSize = 0.48330092430114746

    heme_profile_filename = working_directory + 'config.pr2' 
    ConfigDirectory = working_directory + 'config.stl'
    # I think the problem is in this command 
    # command = hemelb_setup_exe1 + ' ' + heme_profile_filename + ' --voxel ' + str(voxel_size)  + ' --geometry ' + working_directory + '/config.gmy' + ' --xml ' + working_directory + '/config.xml'

    command = hemelb_setup_exe1 + ' ' + heme_profile_filename + ' --voxel ' + str(voxel_size) # + ' --geometry ' + working_directory + '/config.gmy' + ' --xml ' + working_directory + '/config.xml'

    subprocess.call(command, shell=True)
    print "HemeLB set up is DONE"

    # pause()

    print "update_xml_file"


    # pause()

    # _______________________
    # Now update the xml file 
    # _______________________


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

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-tangentialprojectiontraction.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-traction.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-tractions.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')
    ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')    

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-pressure.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='pressure')

    wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'wholegeometry-velocity.xtr'})
    ElementTree.SubElement(wholegeometry, 'geometry', type='whole')
    ElementTree.SubElement(wholegeometry, 'field', type='velocity')

    # Save XML file to disk
    tree.write(filename)

    print('here')

    print "done"


    # _______________________
    # Run the HemeLB simulation  
    # _______________________

    print "**************   hemelb step  **************   "
    shutil.rmtree(working_directory + 'results/', ignore_errors=True)
    xml_config_file = working_directory + 'config.xml'
    command = 'mpirun -np 1 hemelb -in ' + xml_config_file
    subprocess.call(command, shell=True)
    # mpirun -np 1 hemelb -in "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/config.xml"


    # _______________________
    # Generate the flow files  
    # _______________________


    print "  generate_flow_vtus   "
        # Make the vtu 
    command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/GmyUnstructuredGridReader.py' + ' ' + working_directory + 'config.xml'
    subprocess.call(command, shell=True)

    # this need to be generalised 
    command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py' + ' ' +  working_directory +'config.vtu' +' ' +  working_directory+ 'results/Extracted/surface-pressure.xtr' +' ' +  working_directory+ 'results/Extracted/wholegeometry-velocity.xtr'+' ' +  working_directory+ 'results/Extracted/surface-traction.xtr'+' ' +  working_directory+ 'results/Extracted/surface-tangentialprojectiontraction.xtr'
    subprocess.call(command, shell=True)

        # output_folders = glob.glob(working_directory + '/results_from_time_*')

        # for folder in output_folders:
        #     folder_timestep = folder.split('results_from_time_')[-1]
        #     if float(folder_timestep) == timestep:
        #         shutil.copyfile(working_directory + '/'+ 'results/Extracted/surface-pressure_2000.vtu', folder + '/surface-pressure.vtu')
        #         shutil.copyfile(working_directory + '/'+ 'results/Extracted/surface-traction_2000.vtu', folder + '/surface-traction.vtu')
        #         shutil.copyfile(working_directory + '/'+ 'results/Extracted/wholegeometry-velocity_2000.vtu', folder + '/wholegeometry-velocity.vtu')
        #         shutil.copyfile(working_directory + '/'+ 'results/Extracted/surface-tangentialprojectiontraction_2000.vtu', folder + '/surface-tangentialprojectiontraction.vtu')
    
    print "generate_flow_vtus DONE    "
    t1 = time.time()

    total = t1-t0
    print "Time taken for simulaiton + " + str(total)






voxel_size =  0.1953125e-6 #0.1953e-06
subfolder = [ 'Growing_1_5','Growing_2', 'Growing_2_5']#, 'Growing_3', 'Growing_3_5', 'Growing_3', 'Growing_3_5', 'Growing_4', 'Growing_4_5', 'Growing_5']

for i in subfolder :
    working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/TestingShapes4/' + i +'/' #TestCylinderRadi_voxal/NotScalledUp/Radi_e-3/Voxal1.4e-8/' 

    print('here')
    # _______________________
    #  COnvert vtu to the stl i need 
    #  _______________________



    # print "  Convert vtu to stl    "
    #     # Read the VTU file from disk
    # vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    # vtu_reader.SetFileName(working_directory + 'config.vtu')

    # extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    # extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

    # # Write out the data in unstructured grid format
    # stl_writer = vtk.vtkSTLWriter()
    # stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    # stl_writer.SetFileName(working_directory + 'config.stl')
    # stl_writer.Write()


    # _______________________
    # Set up the inital config.xml file 
    # _______________________


    print "HemeLB setting up"
    # command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/GmyUnstructuredGridReader.py ~/Documents/ChasteWorkingDirectory/config.xml'
    # subprocess.call(command, shell=True)

    # TODO: It seems that when you use the --stl option below, the voxel size specified in .pr2 is ignored and the setup tool guesses one for you



    # VoxelSize = 0.48330092430114746

    heme_profile_filename = working_directory + 'config.pr2' 
    ConfigDirectory = working_directory + 'config.stl'
    # I think the problem is in this command 
    # command = hemelb_setup_exe1 + ' ' + heme_profile_filename + ' --voxel ' + str(voxel_size)  + ' --geometry ' + working_directory + '/config.gmy' + ' --xml ' + working_directory + '/config.xml'

    command = hemelb_setup_exe1 + ' ' + heme_profile_filename + ' --voxel ' + str(voxel_size) # + ' --geometry ' + working_directory + '/config.gmy' + ' --xml ' + working_directory + '/config.xml'

    subprocess.call(command, shell=True)
    print "HemeLB set up is DONE"

    # pause()

    print "update_xml_file"


    # pause()

    # _______________________
    # Now update the xml file 
    # _______________________


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

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-tangentialprojectiontraction.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-traction.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-tractions.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')
    ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')    

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-pressure.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='pressure')

    wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'wholegeometry-velocity.xtr'})
    ElementTree.SubElement(wholegeometry, 'geometry', type='whole')
    ElementTree.SubElement(wholegeometry, 'field', type='velocity')

    # Save XML file to disk
    tree.write(filename)

    print('here')

    print "done"


    # _______________________
    # Run the HemeLB simulation  
    # _______________________

    print "**************   hemelb step  **************   "
    shutil.rmtree(working_directory + 'results/', ignore_errors=True)
    xml_config_file = working_directory + 'config.xml'
    command = 'mpirun -np 1 hemelb -in ' + xml_config_file
    subprocess.call(command, shell=True)
    # mpirun -np 1 hemelb -in "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/config.xml"


    # _______________________
    # Generate the flow files  
    # _______________________


    print "  generate_flow_vtus   "
        # Make the vtu 
    command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/GmyUnstructuredGridReader.py' + ' ' + working_directory + 'config.xml'
    subprocess.call(command, shell=True)

    # this need to be generalised 
    command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py' + ' ' +  working_directory +'config.vtu' +' ' +  working_directory+ 'results/Extracted/surface-pressure.xtr' +' ' +  working_directory+ 'results/Extracted/wholegeometry-velocity.xtr'+' ' +  working_directory+ 'results/Extracted/surface-traction.xtr'+' ' +  working_directory+ 'results/Extracted/surface-tangentialprojectiontraction.xtr'
    subprocess.call(command, shell=True)

        # output_folders = glob.glob(working_directory + '/results_from_time_*')

        # for folder in output_folders:
        #     folder_timestep = folder.split('results_from_time_')[-1]
        #     if float(folder_timestep) == timestep:
        #         shutil.copyfile(working_directory + '/'+ 'results/Extracted/surface-pressure_2000.vtu', folder + '/surface-pressure.vtu')
        #         shutil.copyfile(working_directory + '/'+ 'results/Extracted/surface-traction_2000.vtu', folder + '/surface-traction.vtu')
        #         shutil.copyfile(working_directory + '/'+ 'results/Extracted/wholegeometry-velocity_2000.vtu', folder + '/wholegeometry-velocity.vtu')
        #         shutil.copyfile(working_directory + '/'+ 'results/Extracted/surface-tangentialprojectiontraction_2000.vtu', folder + '/surface-tangentialprojectiontraction.vtu')
    
    print "generate_flow_vtus DONE    "
    t1 = time.time()

    total = t1-t0
    print "Time taken for simulaiton + " + str(total)
