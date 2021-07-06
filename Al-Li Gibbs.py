#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:30:26 2021

@author: seblohier
"""
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import background
import csv
import os
import pandas as pd
import matplotlib.pyplot as plt
from pycalphad import Database, calculate, variables as v
from pycalphad.plot.utils import phase_legend
import numpy as np

# Load database and choose the phases that will be plotted
db_nbre = Database('Al-Li.tdb')
my_phases_nbre = ['LIQUID','BCC_A2','FCC_A1','BCC_DIS','BCC_B2','AL2LI3','AL4LI9']

# Get the colors that map phase names to colors in the legend
legend_handles, color_dict = phase_legend(my_phases_nbre)

fig = plt.figure(figsize=(9,6))
ax = fig.gca()
DataValues=[]
temper=600
df = pd.DataFrame(list())
phase_potentials={}
hullval=[]
# Loop over phases, calculate the Gibbs energy, and scatter plot GM vs. X(RE)
for phase_name in my_phases_nbre:
    result = calculate(db_nbre, ['AL', 'LI','VA'], phase_name, P=101325, T=temper, output='GM')
    print(phase_name)

    temp=result.GM.values
    tempx=result.X.sel(component='LI').values

   
    x = tempx.squeeze() 
    y = temp.squeeze()
    #Here the (x,y) pairs of points are but into a more simple and callable form
    #We use a try except case here because one of the phases has only one data point. Meaning when
    #try to ge the length it gives a type eror
    try:
        #This portion of the code needs to be changed but all it does is get the hull of the gibbs surface
        if phase_name=='BCC_B2':
            minihull=[]
        for i in range(len(x)):
            app=[x[i],y[i]]
            if phase_name == 'BCC_B2':
                minihull.append(app)
            hullval.append(app)
        if phase_name=='BCC_B2':
            minihull=np.array(minihull)
            liquidx=[]
            liquidy=[]
            lhull=ConvexHull(minihull)
            for simplex in lhull.simplices:
                liquidx.extend(minihull[simplex,0].tolist())
                
                liquidy.extend(minihull[simplex,1].tolist())
            liquidx.extend(minihull[lhull.vertices,0].tolist())
            liquidy.extend(minihull[lhull.vertices,1].tolist())

            x=liquidx
            y=liquidy
        #From this point all the code makes a polynomial estimate of the specific phase curve 
        #then find the derivative
        estimate=np.polyfit(x,y,8)
        derivative=background.find_derivative(estimate)
        alphapot=[]
        betapot=[]
        for i in range(len(x)):
            #This derivative is used to find the potential of the phase at each composition
            slope=derivative(x[i])
            potential=background.find_potentials(slope, x[i], y[i], True)
            alphapot.append(potential[0])
            betapot.append(potential[1])
        print(len(x))
        phase_potentials[phase_name]={'alpha':alphapot,'beta':betapot,'x':x}
    except TypeError:
        print("We good")
    #The points are the plotted
    ax.scatter(x,y, marker='.', s=5, color=color_dict[phase_name])

    #This commented portion takes the values attained and turns them into a csv file labeling what it is
    # background.csvform("Al", "Li", phase_name, temper)
          

# Format the plot
ax.set_xlabel('X(NI)')
ax.set_ylabel('GM')
ax.set_xlim((0, 1))
ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1, 0.6))

pplot=plt.figure(2)
ax2=pplot.gca()
phase='LIQUID'
ax2.scatter(phase_potentials[phase]['x'],phase_potentials[phase]['alpha'],color='red')
ax2.scatter(phase_potentials[phase]['x'],phase_potentials[phase]['beta'],color='blue')

plt.show()
