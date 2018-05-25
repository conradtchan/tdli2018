"""
TDLI Summer School 2018. Plotting routine for 1D hydro output.
Output should be located in a directory ./output/, with filename format 'out0001.dat', 'out0002.dat', etc.
"""

import numpy as np  # NumPy library for array operations
import matplotlib.pyplot as plt  # Matplotlib library for plotting
import glob # Glob module for filename pattern matching
import sys  # Sys module for accessing methods/constants/functions from the python interpreter
import os  # Os module for miscellaneous operating system interfaces

var_dict = {    # Dictionary of variables : column number in output files. Adjust this to
            'rho' : 0,  # match your output format.
            'v1' : 1,
            'v2' : 2,
            'v3' : 3,
            'eps' : 4,
            'pre' : 5,
            'cs' : 6
            }

def help():
    print('Pass the variable you wish to plot as an argument (u, rho, pre, etc..)')
    print('e.g.; $./plot.py rho')
    print()
    print('Variable must be one of {}, otherwise, edit the python script'.format(list(var_dict.keys())))
    print()

if len(sys.argv) == 2:
    var = str(sys.argv[1])  # The variable to plot is a command line argument, e.g.; "$./plot.py rho"
else:
    help()
    sys.exit()

if var not in var_dict:   # Print an error message if the variable passed at the command line is
    print()               # not found in the dictionary of output variables.
    print('Variable {} not found.')
    help()
    sys.exit()

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
    plt.imshow(data[:,:,var_dict[var]])
    plt.colorbar()
    # plt.show(block=False)  # Show the plot. block=False allows the code to continue, if block=True the code will stop here
    # plt.pause(0.01)  # Pause the animation for a short time between frames
    plt.savefig(f'plot/{var}_{n:0>5d}.png')
    n += 1

input()  # This requires some input before the code will complete, i.e., the animation will freeze on
         # the last frame until you hit enter at the command line.
