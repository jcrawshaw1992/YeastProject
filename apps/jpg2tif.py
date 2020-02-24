import subprocess
import vtk
import shutil
import os
from xml.etree import ElementTree
import glob
import numpy as np
import time
import pdb
import string
import math
import sys

from PIL import Image
im = Image.open('/Users/jcrawshaw/docker-polnet-master/PlexusWithLongInlets/PlexusWithInlets.jpg')
im.save('/Users/jcrawshaw/docker-polnet-master/PlexusWithLongInlets/test.tiff')  # or 'test.tif'
