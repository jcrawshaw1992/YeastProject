#!/usr/bin/env python
# vasc_res_sim.py
import subprocess
# import vtk
import shutil
import os
from xml.etree import ElementTree
from argparse import ArgumentParser
import string
import math


if __name__=="__main__":

    # Define arguments to be parsed
    parser = ArgumentParser(description='Update the necessary xml file for running HemeLB')
    parser.add_argument('-period', dest='period', default='1000', help='Period to print out results, which will be turned into vtus --- we really only need the last output.')
    parser.add_argument('-directory', dest='directory', help='Need to provide a destination to find the xml file to edit ')
    parser.add_argument('-InitalConditions', dest='IC', default=0, help='Need to provide a destination to find the xml file to edit ')
    parser.add_argument('-ConvergenceTermination', dest='ConvergenceTermination', type=str, default='false', help='To terminate when the simulation reaches a steady state')
    parser.add_argument('-AverageVelocity', dest='AverageVelocity', type=float, default='false', help='To terminate when the simulation reaches a steady state')
## Need to play with this Jess

    

    # Only the final time output is needed, so this should be the period --- this is because it is time costly
    # to write out files, aaaand until the last time point, I cant be sure 
    # that the simulation is at equalibirum --- prior to equi I just have rubbish
   
    args = parser.parse_args()
    period = args.period
    filename = args.directory + 'config.xml'
    Terminate= args.ConvergenceTermination # Set to true or false :) 
    IC =  args.IC
    AverageVelocity= args.AverageVelocity
    tolerance= AverageVelocity*1e-5

    tree = ElementTree.parse(filename)
    root = tree.getroot()

    # Add monitoring of incompressibility and convergence
    monitoring = ElementTree.SubElement(root, 'monitoring')
    ElementTree.SubElement(monitoring, 'incompressibility') 
    # convergence = ElementTree.SubElement(monitoring, 'steady_flow_convergence', {'tolerance': '1e-3', 'terminate': 'false'})
    convergence = ElementTree.SubElement(monitoring, 'steady_flow_convergence', {'tolerance': str(tolerance), 'terminate': Terminate})
    # ElementTree.SubElement(convergence, 'criterion', {'type': 'velocity', 'value': '0.01', 'units': 'm/s'})
    ElementTree.SubElement(convergence, 'criterion', {'type': 'velocity', 'value': str(AverageVelocity), 'units': 'm/s'})
    
    # Add definition of properties to be extracted
    extr = ElementTree.SubElement(root, 'properties')

    # surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-tangentialprojectiontraction.xtr'})
    # ElementTree.SubElement(surface, 'geometry', type='surface')
    # ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'surface-traction.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'surface-tractions.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')
    ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')    

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'surface-pressure.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='pressure')

    wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'wholegeometry-velocity.xtr'})
    ElementTree.SubElement(wholegeometry, 'geometry', type='whole')
    ElementTree.SubElement(wholegeometry, 'field', type='velocity')

    # Save XML file to disk
    tree.write(filename)
    print "Finished updateing xml file "


    xml_file = open(filename, "r")
    list_of_lines = xml_file.readlines()
    NumberOfLines = len(list_of_lines)
    print len(list_of_lines)
    print list_of_lines[NumberOfLines-4][7:15]

    list_of_lines[NumberOfLines-4] ='\t\t\t <uniform units="mmHg" value="' +str(IC)+'" />\n'
    xml_file = open(filename, "w")
    xml_file.writelines(list_of_lines)
    xml_file.close()





