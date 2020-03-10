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

working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/Cylinder2' # dont know what this is.... testoutput????ChasteWorkingDirectory

data_path ='/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/test/data/Clyinder/'
# Have changed the following to the debug folder 
chaste_setup_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/debug/TestSetupFlowInPipeRunner'
chaste_run_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/debug/TestRunFlowInPipeRunner'
hemelb_setup_exe = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools:/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'



# data_path = '/Users/mobernabeu/Documents/workspace/Chaste/projects/VascularRemodelling/test/data/'
# working_directory = '/Users/mobernabeu/Documents/workspace/ChasteWorkingDirectory'

# chaste_setup_exe = '/Users/mobernabeu/Documents/workspace/Chaste/projects/VascularRemodelling/build/optimised/TestSetupFlowInPipeRunner'
# chaste_run_exe = '/Users/mobernabeu/Documents/workspace/Chaste/projects/VascularRemodelling/build/optimised/TestRunFlowInPipeRunner'
# hemelb_setup_exe = 'env PYTHONPATH=hemelb/Tools/setuptool/:$PYTHONPATH hemelb-dev/Tools/setuptool/scripts/hemelb-setup-nogui'

radii_over_time = []

def vmtk_compute_stl_radii(iter_number):
    output_filename = working_directory + 'centerlines' +  str(iter_number) + '.vtp'
    command = 'vmtk vmtknetworkextraction -ifile ' + working_directory + 'config.stl' + ' -ofile ' + output_filename
    subprocess.call(command, shell=True)

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(output_filename)
    reader.Update()
    point_data = reader.GetOutput().GetPointData()
    assert point_data.GetArrayName(0) == 'Radius', "VTP file doesn't contain array Radius"
    radii = point_data.GetArray(0)
    radii_over_time.append(np.array([radii.GetValue(i) for i in range(radii.GetSize())]))

# def generate_flow_vtus(timestep):
#     print 'Turning hemelb results into vtus'
#    # command = 'python -m hemeTools.converters.ExtractedPropertyUnstructuredGridReader ' + working_directory + 'config.xml ' + working_directory + 'results/Extracted/surface-pressure.xtr ' + working_directory + 'results/Extracted/wholegeometry-velocity.xtr ' + working_directory + 'results/Extracted/surface-traction.xtr ' + working_directory + 'results/Extracted/surface-tangentialprojectiontraction.xtr'
#     command = python hemeTools/converters/ExtractedPropertyUnstructuredGridReader Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/test/data/Untitled.gmy  Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/test/data/results/Extracted/surface-pressure.xtr  Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/test/data/results/Extracted/wholegeometry-velocity.xtr Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/test/data/results/Extracted/surface-traction.xtr  Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/test/data/results/Extracted/surface-tangentialprojectiontraction.xtrpython -m hemeTools.converters.ExtractedPropertyUnstructuredGridReader Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/test/data/Untitled.gmy  Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/test/data/results/Extracted/surface-pressure.xtr  Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/test/data/results/Extracted/wholegeometry-velocity.xtr Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/test/data/results/Extracted/surface-traction.xtr  Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/test/data/results/Extracted/surface-tangentialprojectiontraction.xtr
#     Working direct = Users.jcrawshaw/Documents.Chaste/projects/VascularRemodelling/test/data/
#     subprocess.call(command, shell=True)

#     output_folders = glob.glob(working_directory + '/results_from_time_*')

#     for folder in output_folders:
#         folder_timestep = folder.split('results_from_time_')[-1]
#         if float(folder_timestep) == timestep:
#             shutil.copyfile(working_directory + 'results/Extracted/surface-pressure_2000.vtu', folder + '/surface-pressure.vtu')
#             shutil.copyfile(working_directory + 'results/Extracted/surface-traction_2000.vtu', folder + '/surface-traction.vtu')
#             shutil.copyfile(working_directory + 'results/Extracted/wholegeometry-velocity_2000.vtu', folder + '/wholegeometry-velocity.vtu')
#             shutil.copyfile(working_directory + 'results/Extracted/surface-tangentialprojectiontraction_2000.vtu', folder + '/surface-tangentialprojectiontraction.vtu')

def vtu2stl():
    # Read the VTU file from disk
    vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    vtu_reader.SetFileName(working_directory + 'config.vtu')

    # vtkAppendFilter appends one or more datasets together into a single unstructured grid
    extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

    # Write out the data in unstructured grid format
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    stl_writer.SetFileName(working_directory + 'config.stl')
    stl_writer.Write()

def run_hemelb_setup():
    # TODO: It seems that when you use the --stl option below, the voxel size specified in .pr2 is ignored and the setup tool guesses one for you
    if  input_folder == 'embryo_plexus/': 
        voxel_size = 1.222e-6 
    else: 
        voxel_size = 0.15e-3

    heme_profile_filename = data_path + input_folder + 'config.pr2' 


def update_xml_file(iter_num, num_iters):
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

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(1000), 'file': 'surface-tangentialprojectiontraction.xtr'})
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


def run_hemelb():
    shutil.rmtree(working_directory + '/results/', ignore_errors=True)
    xml_config_file = working_directory + 'config.xml'
    command = 'mpirun -np 4 hemelb -in ' + xml_config_file
    subprocess.call(command, shell=True)


if __name__=="__main__":

    # These simulations need updating
    old_simulations = ['straight_vessel', 'cylinder_validation', 'bifurcation_cut_changing_pressure', 'bifurcation_weak_region', 'bifurcation_shrink', 'capillary_network']

    existing_simulations = ['bifurcation_cut', 'embryo_plexus']
        
    # Define arguments to be parsed
    parser = ArgumentParser(description='Run a vascular remodelling simulation')
    parser.add_argument('-o', dest='overwrite_results', action='store_true', help='Allow overwriting of existing results folder.')
    parser.add_argument('--hemelb_vtu', dest='flow_to_vtu', action='store_true', help='Create VTU files with flow results.')
    parser.add_argument('--distributed', dest='run_distributed_with_muscle', action='store_true', help='Enables distributed execution via MAPPER toolkit.')
    # parser.add_argument('simulation_type', choices=existing_simulations, help='Type of simulation to be run.')
    parser.add_argument('num_iterations', nargs='?', default=5, type=int, help='Number of Hemelb/Chaste iterations to be run (optional, default is 5).')
    parser.add_argument('--compute_radii', dest='compute_radii', action='store_true', help='Use VMTK to compute axis radii.')
    parser.add_argument('--output_postfix', dest='output_postfix', default='', help='This string will be added to ChasteWorkingDirectory to get the output folder.')
    parser.add_argument('--div_threshold', dest='div_threshold', default=1e10, help='This specifies the length that edges will divide at. (Defaults to no division, but 6e-4 is good for mm meshes')
    parser.add_argument('--mesh_scale', dest='mesh_scale', default=1e-3, help='This specifies what to scale the mesh by so that all distances are in meters (defaults to mm).')
    
    # Parse arguments (this will create args.flow_to_vtu etc. variables)
    args = parser.parse_args()

    # Define input folder from simulation type chosen
    # input_folder = args.simulation_type + '/'
    
    # Sort out working directory (create, overwrite if allowed, etc.)
    if len(args.output_postfix) == 0:
        working_directory = working_directory + '/'
    # else:
    #     working_directory = working_directory + '-' + args.output_postfix + '/'

    # if os.path.isdir(working_directory):
    #     if args.overwrite_results:
    #         shutil.rmtree(working_directory)
    #         os.mkdir(working_directory)
    #     else:
    #         raise Exception('Results folder {} exists. Enable results overwriting with -o. '.format(working_directory))
    # else:
    #     os.mkdir(working_directory)

    os.environ['CHASTE_TEST_OUTPUT'] = working_directory

    print 'Working directory = ' + working_directory


    #
    # Step 1: Run preliminary Chaste setup
    #
    # mesh_filename = data_path + input_folder + 'config.vtu' # Use proper path concatentation
    # xml_filename = data_path + input_folder + 'config.xml' # Use proper path concatentation
    mesh_filename = data_path + 'config.vtu' # Use proper path concatentation
    xml_filename = data_path + 'config.xml' # Use proper path concatentation
    print mesh_filename
    # print 'Division threshold ', str(args.div_threshold)

    # chaste_setup_call = chaste_setup_exe + ' -mesh ' + mesh_filename + ' -xml ' + xml_filename + ' -div_threshold ' +  str(args.div_threshold) + ' -mesh_scale ' +  str(args.mesh_scale) + ' -wss_threshold ' + str(0.025)
    # # if args.simulation_type == 'bifurcation_weak_region':
    # #     chaste_setup_call += ' -label_weak_region'
    # # if args.simulation_type == 'bifurcation_shrink':
    # #     chaste_setup_call += ' -shrink_below_stress_threshold'
    # subprocess.call(chaste_setup_call, shell=True)

    chaste_setup_call = chaste_setup_exe + ' -mesh ' + mesh_filename + ' -xml ' + xml_filename + ' -div_threshold ' +  str(args.div_threshold) + ' -mesh_scale ' +  str(args.mesh_scale) + ' -wss_threshold ' + str(0.025)
    subprocess.call(chaste_setup_call, shell=True)


    #
    # Step 1.a: Prepare mesh for HemeLB setup tool and compute radii along the axes of the domain (if requested)
    #
    vtu2stl()
    if args.compute_radii:
        vmtk_compute_stl_radii(0)


    start_time = 0
    duration = 1
    for iter in range(args.num_iterations):
        print 'Iteration ' + str(iter)
        #
        # Step 2: Run HemeLB setup
        #
        run_hemelb_setup()

        #
        # Step 3: Use HemeLB to calculate tractions
        #
        update_xml_file(iter, args.num_iterations)
        if args.run_distributed_with_muscle:
            # Run the simulation described by /home/mobernabeu/workspace/ChasteWorkingDirectory/config.xml and /home/mobernabeu/workspace/ChasteWorkingDirectory/config.gmy
            raise NotImplementedError
            # Ensure that /home/mobernabeu/workspace/results/Extracted/surface-traction.xtr' exists for Chaste to carry on.
        else:
            run_hemelb()
        
        #
        # Step 4: Run Chaste with the tractions previously computed
        #
        traction_filename = working_directory + 'results/Extracted/surface-tractions.xtr'
        command = chaste_run_exe + ' -start_time ' + str(start_time) + ' -duration ' + str(duration) + ' -traction_file ' +  traction_filename + ' -mesh_scale ' +  str(args.mesh_scale)
        subprocess.call(command, shell=True)
        start_time += duration

        # This has to be done after running Chaste, otherwise it will overwrite the folder
        # if args.flow_to_vtu:
        #     generate_flow_vtus(iter*duration)
                
        #
        # Step 5: Convert Chaste output (vtu file) into the input of the next iteration (stl file required by hemelb setup tool)
        #
        vtu2stl()

        #
        # Step 6: Compute radii along the axes of the domain (if requested)
        if args.compute_radii:
            vmtk_compute_stl_radii(iter+1)

    if args.compute_radii:
        import matplotlib.pyplot as plt
        import csv
        csvfile = open(working_directory + 'radii_over_time.csv', 'wb')
        writer = csv.writer(csvfile)
        for iter_num, radii in enumerate(radii_over_time):
            time = iter_num*duration # The end of the n-th iteration
            plt.plot(radii, label='t={}'.format(time))
            line = [time]
            line.extend(radii)
            writer.writerow(line)
        plt.legend()
        plt.savefig(working_directory + 'radii_over_time.png')
        csvfile.close()
        
            
