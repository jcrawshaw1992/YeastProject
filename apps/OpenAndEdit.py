# Formate .pr2 files -- this code will put a new line between each line, 
# elsewhise when updating the file, some lines will be thrown away 
import subprocess
import vtk
import shutil
import os
from xml.etree import ElementTree
import glob
from argparse import ArgumentParser
import numpy as np
import time
import csv
import pdb
import string
import math

from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove


# def Format_pr2_files(working_directory):
working_directory = '/Users/jcrawshaw/Downloads/'
filename = working_directory + 'config.xml'
file = open(filename).readlines()
#Create temp file
file_path = filename
pattern = '<uniform units="mmHg" value="0.0" />'
subst = '<uniform units="mmHg" value="100.0" />'

fh, abs_path = mkstemp()
with os.fdopen(fh,'w') as new_file:
    with open(file_path) as old_file:
        for line in old_file:
            new_file.write(line.replace(pattern, subst))
#Copy the file permissions from the old file to the new file
copymode(file_path, abs_path)
#Remove original file
os.remove(file_path)
#Move new file
move(abs_path, file_path)