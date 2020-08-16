import math
import subprocess
# import vtk
import shutil
import os
from argparse import ArgumentParser
import string
import math


if __name__=="__main__":

    # Define arguments to be parsed
    parser = ArgumentParser(description='Inputs for log calculation')
    parser.add_argument('-base', dest='Base', type =float, help='log base')
    parser.add_argument('-number', dest='Number', type =float, help='Number getting log of')
   
    args = parser.parse_args()
    base = args.Base
    number = args.Number
    exponent = math.log(number, base)  # = 3
    print exponent
        


