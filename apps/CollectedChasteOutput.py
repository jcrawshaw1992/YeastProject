import os
import shutil
from argparse import ArgumentParser



Iterations =19
parser = ArgumentParser(description='Run a vascular remodelling simulation')
parser.add_argument('--num_iterations', default=Iterations, type=int, help='Number of Hemelb/Chaste iterations to be run (optional, default is 5).')
args = parser.parse_args()
working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory' 


# newpath = r'/Users/jcrawshaw/Documents/ChasteWorkingDirectory/NewCylinder/ChasteMeshes/' 
newpath =  '/Users/jcrawshaw/Documents/testoutput/PottsSimlation/' 
if not os.path.exists(newpath):
    os.makedirs(newpath)

oldpath = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/results_from_time_'

NumberOfFiles = 11
list = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5,5,5.5,6]
print list
i=0

Directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/TestCollection/'
for j in list:
    directory = Directory + 'results_from_time_'+ str(j)+ '/'
    print "Jess is great"
    for file in os.listdir(directory):
        if file.endswith(".vtu"):
            if file.startswith("mesh"):
             shutil.move(directory + file , newpath+'mesh_'+ str(i) +'.vtu'  )
            # if file.startswith("results"):
            #  shutil.move(directory + file , newpath+'results_'+ str(i) +'.vtu'  )
            i+=1

print "done"
            