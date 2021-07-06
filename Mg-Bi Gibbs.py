#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:30:26 2021

@author: seblohier
"""
import background
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import os
import csv
import matplotlib.pyplot as plt
from pycalphad import Database, calculate, variables as v
from pycalphad.plot.utils import phase_legend
import numpy as np

# Load database and choose the phases that will be plotted
db_nbre = Database('Mg-Bi.tdb')
my_phases_nbre = ['LIQUID', 'BI2MG3_H', 'BI2MG3_L']

# Get the colors that map phase names to colors in the legend
legend_handles, color_dict = phase_legend(my_phases_nbre)
hullval=[]
fig = plt.figure(num=1,figsize=(9,6))
ax = fig.gca()
temper=1000
# Loop over phases, calculate the Gibbs energy, and scatter plot GM vs. X(RE)
phase_potentials={}
for phase_name in my_phases_nbre:
    result = calculate(db_nbre, ['BI','MG','VA'], phase_name, P=101325, T=temper, output='GM')
    
       
        #ax.scatter(result.X.sel(component='MG'), result.GM, marker='.', s=5, color=color_dict[phase_name])
    
    print(phase_name)
    temp=result.GM.values
    tempx=result.X.sel(component='MG').values
    
    x = tempx.squeeze()     
    y = temp.squeeze()

    #Here the (x,y) pairs of points are but into a more simple and callable form
     #This portion of the code needs to be changed but all it does is get the hull of the gibbs surface
    if phase_name=='LIQUID':
        minihull=[]
    for i in range(len(x)):
        app=[x[i],y[i]]
        if phase_name == 'LIQUID':
            minihull.append(app)
        hullval.append(app)
    if phase_name=='LIQUID':
        minihull=np.array(minihull)
        liquidx=[]
        liquidy=[]
        lhull=ConvexHull(minihull)
        for simplex in lhull.simplices:
            liquidx.extend(minihull[simplex,0].tolist())
            liquidy.extend(minihull[simplex,1].tolist())
        x=liquidx
        y=liquidy
    #From this point all the code makes a polynomial estimate of the specific phase curve 
    #then find the derivative
    estimate=np.polyfit(x,y,8)
    derivative=background.find_derivative(estimate)
    alphapot=[]
    betapot=[]
    for i in range(len(x)):
        slope=derivative(x[i])
        potential=background.find_potentials(slope, x[i], y[i], True)
        alphapot.append(potential[0])
        betapot.append(potential[1])
    print(len(x))
    phase_potentials[phase_name]={'alpha':alphapot,'beta':betapot,'x':x}
  

    ax.scatter(x,y, marker='.', s=5, color=color_dict[phase_name])

        
    #This commented portion takes the values attained and turns them into a csv file labeling what it is
    # background.csvform("Mg", "Bi", phase_name, temper)
#This finds the convex hull of the plot
hullval=np.array(hullval)
hull = ConvexHull(hullval)

# plt.plot(hullval[:][0], hullval[:][1], 'o')
#for simplex in hull.simplices:

    #plt.plot(hullval[simplex,0], hullval[simplex,1], 'k-')
#plt.plot(hullval[hull.vertices,0], hullval[hull.vertices,1], 'r--', lw=2)
#plt.plot(hullval[hull.vertices[0],0], hullval[hull.vertices[0],1], 'ro')

# plt.show()
# Format the plot
ax.set_xlabel('X(Mg)')
ax.set_ylabel('GM')
ax.set_xlim((0, 1))
ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1, 0.6))

pplot=plt.figure(2)
ax2=pplot.gca()
phase='LIQUID'
ax2.set_xlabel('Mole Percent of Mg',fontsize=16)
ax2.set_ylabel('Chemical Potential')
ax2.set_title('Chemical Potential of the '+phase+' phase at '+str(temper)+" K")
ax2.scatter(phase_potentials[phase]['x'],phase_potentials[phase]['alpha'],color='red',label='Pure Mg potential')
ax2.scatter(phase_potentials[phase]['x'],phase_potentials[phase]['beta'],color='blue',label='Pure Bi potential')
ax2.legend(loc='best', scatterpoints=100)

plt.show()
