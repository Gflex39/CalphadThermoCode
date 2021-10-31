import sympy as sym
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import os
import csv
import matplotlib.pyplot as plt
from pycalphad import Database, calculate, variables as v
from pycalphad.plot.utils import phase_legend
import numpy as np
import mpmath as mp
#This file is just a dump for the functions used for all of the Gibbs curve files.
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
def csvform(e1,e2,phase_name,temper,data,keys):
    print("writen")
    name=e1+"-"+e2+" |"+phase_name+"| T="+str(temper)+".csv"
    csv_columns = keys
    # print(data)
    cur=os.path.dirname(__file__)
    os.chdir(cur+"/Data_Files_"+e1+"_"+e2)
    with open(name, "w+") as f:
        writer = csv.writer(f,delimiter = "\t")
        writer.writerow(csv_columns)
        # print(list(zip(*[data[key] for key in csv_columns]))[0])
        writer.writerows(zip(*[data[key] for key in csv_columns]))
    os.chdir(cur)
        # row_writer=csv.writer(f,delimiter=',')
        # try:
        #     for i in range(len(data)):
        #         row_writer.writerow([str(x[i]),str(y[i])])
        # except TypeError:
        #         row_writer.writerow([str(x),str(y)])
def checkdata(e1,e2,phase_name,temper):
    data={}
    transfer={}
    name=e1+"-"+e2+" |"+phase_name+"| T="+str(temper)+".csv"
    cur=os.path.dirname(__file__)
    os.chdir(cur+"/Data_Files_"+e1+"_"+e2)
    # if True:
    try:
        with open(name,"r") as f:
            reader = csv.reader(f, delimiter='\t')
            first=True
            
            for row in reader:
                if first:
                    first=False
                    count=0
                    keys=row
                    # print(keys)
                    for j in keys:
                        data[j]=[]
                        transfer[count]=j
                        count+=1
                    
                else:
                    for j in range(len(row)):
                        data[transfer[j]].append(float(row[j]))
        os.chdir(cur)
    except:
        os.chdir(cur)
        assert False
    
    return data

            
            


    os.chdir(cur)
    
def ideal_entropy(x):
    return mp.fmul(-8.314,(mp.fmul(x,mp.log(x))+mp.fmul((1-x),mp.log((1-x)))))
def ideal_gibbs(x,t):
    return mp.fmul(mp.fmul(8.314,(mp.fmul(x,mp.log(x))+mp.fmul((1-x),mp.log((1-x))))),t)


# c1, c2,c3,c4,c5,c6 = parameters('c1, c2','c3','c4','c5','c6')
# print(Poly( {(2, 0): c1, (0, 2): c3, (1, 1): c4,(1, 0): c2,(0, 1): c5,(3, 0): c6}, x ,y).as_expr())
# Generate example data