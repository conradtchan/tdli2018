#!/usr/bin/python3.5

'''
TDLI Summer School 2018
Lab Exercise 3a - Solving a simple differential equation with the forward Euler method
'''

import numpy as np                  # NumPy library for arrays
import matplotlib.pyplot as plt     # Matplotlib library for plotting

tmax = 10.0                         # Time to perform integration to
dt = 0.01                           # Timestep
time = np.arange(0.0, tmax, dt)     # Initialise the time array; from t = 0 -> t = tmax, in steps of size dt

# Set the constants and initial conditions in the equation to solve; du/dt = -lambda * (u - sin(t))
lambd = 100.0                       # Set lambda as a constant (= 100)
u = np.zeros(len(time))             # Initialise the array; u = 0 for each timestep
u[0] = 1.0                          # Intial condition, u(0) = 1

for i in range(len(time)-1):       # Loop over i in 0,1,2,..,(total number of timesteps - 1)
    # We solve for u(i+1);
    # du/dt = -lambda * (u - sin(t))
    # So,
    # u(i+1) = u(i) + du/dt * dt

    dudt = -lambd * (u[i] - np.sin(time[i]))    # We use the numpy sine function (np.sin)  
    u[i+1] = u[i]  + dt * dudt

    plt.cla()                       # Clear the axis to prevent repeatedly plotting on top of old lines
    plt.xlabel('t')                # x-label for the plot
    plt.ylabel('u')                # y-label for the plot
    plt.plot(time,u,label='numerical solution')   # Plot our numerical solution for the current timestep
    plt.plot(time,np.sin(time),ls='dashed',label=r'$\sin t$') # Plot sin(t) for the current timestep as reference
                                                              # with a dashed linestyle (ls='dashed')
    plt.title(f'time={time[i]:.2f}')# Update the title as the current time to two decimal places    
    plt.legend()                    # Display the legend, where the labels are set in the two lines above
    plt.show(block = False)         # Tell matplotlib to show the plot (block=False allows the code to continue)
    plt.pause(0.02)                 # Pause for a short time for animation
