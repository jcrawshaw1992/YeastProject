import os.path
import time


file_path = 
while not os.path.exists(file_path):
    time.sleep(1)

if os.path.isfile(file_path):
    # read file
else:
    raise ValueError("%s isn't a file!" % file_path)