import sympy as sym
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import os
import csv
import matplotlib.pyplot as plt
from pycalphad import Database, calculate, variables as v
from pycalphad.plot.utils import phase_legend
import numpy as np

def line(m,x,y,comp):
    val=m*(comp-x)+y
    return val
def find_derivative(curve,second=False):
    order=len(curve)-1
    x=sym.Symbol('x')
    function_form=curve[0]*(x**order)
    for i in range(1,len(curve)):
        order-=1
        function_form+=curve[i]*(x**order)
    der=sym.diff(function_form,x)
    if second:
        
        der=sym.diff(function_form,x,2)
        print(der)

    return sym.utilities.lambdify(x,expr=der)
def find_potentials(m,x,y,test,comp=0):
    if test:
        zero_pot=m*(0-x)+y
        one_pot=m*(1-x)+y
        return (zero_pot,one_pot)
    else:
        val=m*(comp-x)+y
        return val
