# File to delete on mass

import subprocess
import vtk
import shutil
import os


def DeleteBuildtextFiles(Directory):
    os.walk(directory)

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


if __name__=="__main__":
  
    directory = '/Users/jcrawshaw/Documents/Projects/Writing/JobApplications/Oxford_149153/' 
    DeleteBuildtextFiles(directory)















    