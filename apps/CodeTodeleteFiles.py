# File to delete on mass

import subprocess
import vtk
import shutil
import os


directory = '/Users/jcrawshaw/Documents/Projects/2020SupervisorMeetings/April/' 

os.walk(directory)

# for x in os.walk(directory):
#     print x

for file in os.listdir(directory):
    if file.endswith(".log"):
        os.remove(directory + file)
    if file.endswith(".aux"):
        os.remove(directory + file)
    if file.endswith(".fls"):
        os.remove(directory + file)
    if file.endswith(".synctex.gz"):
        os.remove(directory + file)
    if file.endswith(".fdb_latexmk"):
        os.remove(directory + file)
    if file.endswith(".out"):
        os.remove(directory + file)
    if file.endswith(".toc"):
        os.remove(directory + file)
    if file.endswith(".blg"):
        os.remove(directory + file)
    if file.endswith(".bbl"):
        os.remove(directory + file)
    if file.endswith(".bcf"):
        os.remove(directory + file)
    if file.endswith(".xml"):
        os.remove(directory + file)
    if file.endswith(".tdo"):
        os.remove(directory + file)