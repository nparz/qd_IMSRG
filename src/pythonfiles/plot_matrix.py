from numpy import *
from pylab import *
from matplotlib import rc, rcParams
import matplotlib.units as units
import matplotlib.ticker as ticker
import sys
import os
import time

# set tickers etc
#rc('text',usetex=True)
#rc('font',**{'family':'serif','serif':['Binding energies']})
#title(r'{\bf Neutron Seperation energies}', fontsize=20)
title(r'Decoupled from Ground State',fontsize=20)     
# read in data from file

fl = open('matrix','r')

ix=[]
jy=[]
zx=[]

i=0

for line in fl:
    i=i+1
    a=line.strip().split()
    g=len(a) 
    for j in range(1,g+1):
        ix.append(i)
        jy.append(j)
        zx.append(float(a[j-1])) 
    
jx = jy[::-1]


axis([0.5,g+.5,0.5,g+.5])
xlabel(r'',fontsize=20)
ylabel(r'',fontsize=20)

sze= g*g/.4266

scatter(ix, jx , c=zx, s=sze, vmin=0.0, vmax=0.5,marker='s',edgecolors='none')

colorbar()


# Save the figure in a separate file
savefig('matrixview.pdf', format='pdf')

# Draw the plot to screen
show()
