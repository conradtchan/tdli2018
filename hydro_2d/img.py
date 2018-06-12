"""
TDLI Summer School 2018. Plotting routine for 1D hydro output.
Output should be located in a directory ./output/, with filename format 'out0001.dat', 'out0002.dat', etc.
"""

import numpy as np  # NumPy library for array operations
import matplotlib.pyplot as plt  # Matplotlib library for plotting
import glob # Glob module for filename pattern matching
import sys  # Sys module for accessing methods/constants/functions from the python interpreter
import os  # Os module for miscellaneous operating system interfaces
from PIL import Image

file_template = 'output/out*.dat'   # Set the template to search for files
files = glob.glob(file_template)    # Collect all of the files which match the template
files = sorted(files)               # Sort the files so that we can plot them in the correct sequence

if not files:   # Print an error message if no files are found (files array is empty), and exit the script
    print()
    print('Output data not found, searching for files output/*.dat')
    print('Ensure your output is in this location with compatible filename, or edit the python script')
    print()
    sys.exit()

x = np.loadtxt('output/x.dat')
y = np.loadtxt('output/y.dat')

n = 0
for output in files:  # Loop over each of the output files
    print('plotting', n)
    data = np.loadtxt(output)  # Load the file into a numpy array
    data = data.reshape(len(y), len(x), 7)

    plt.clf()  # Clear the axis so that we can plot a new line
    img_data = data[:,:,1:4]
    img_data = (img_data + 0.5)*256
    img_data = img_data.astype('uint8')
    img = Image.fromarray(img_data)
    img.save(f'plot/alex_{n:0>5d}.jpg')
    n += 1

input()  # This requires some input before the code will complete, i.e., the animation will freeze on
         # the last frame until you hit enter at the command line.
