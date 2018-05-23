"""
TDLI Summer School 2018. Plotting routine for 1D hydro output.
Output should be located in a directory ./output/, with filename format 'out0001.dat', 'out0002.dat', etc.
"""

import numpy as np  # NumPy library for array operations
import matplotlib.pyplot as plt  # Matplotlib library for plotting
import glob # Glob module for filename pattern matching
import sys  # Sys module for accessing methods/constants/functions from the python interpreter
import os  # Os module for miscellaneous operating system interfaces

var_dict = {'x' : 0,    # Dictionary of variables : column number in output files. Adjust this to
            'rho' : 1,  # match your output format.
            'v1' : 2,
            'v2' : 3,
            'v3' : 4,
            'eps' : 5,
            'pre' : 6,
            'cs' : 7
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

for output in files:  # Loop over each of the output files
    data = np.loadtxt(output)  # Load the file into a numpy array

    plt.cla()  # Clear the axis so that we can plot a new line
    plt.xlabel ('x')  # Set the x-axis label
    plt.ylabel (var)  # Set the y-axis label
    plt.ylim((-1.,1.1))  # Set the axis limits for the y-axis
    plt.plot(data[:,0],data[:,var_dict[var]])  # Plot the data
    plt.show(block=False)  # Show the plot. block=False allows the code to continue, if block=True the code will stop here
    plt.pause(0.1)  # Pause the animation for a short time between frames

input()  # This requires some input before the code will complete, i.e., the animation will freeze on
         # the last frame until you hit enter at the command line.
