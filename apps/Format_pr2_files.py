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
import matplotlib.pyplot as plt
import csv
import pdb
import string
import math


# def Format_pr2_files(working_directory):

filename = working_directory + 'config.pr2'

file = open(filename).readlines()

for line in file:
    s = open(filename).read()
    s = s.replace(line, line +'n')
    f = open(filename, 'w')
    f.write(s)
    f.close()



print "Formated .pr2 file"

    # ---------------------------------
    # ---------------------------------
    # ---------------------------------



