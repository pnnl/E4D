#!/bin/python

# This python script builds the conductivity file for to be used in this example

import numpy as np


# Load the mesh elements
fl = 'sp'
ele = np.loadtxt(fl+'1.ele',skiprows=1,dtype='i4')
nele=len(ele[:,0]) # number of elements in mesh

#load the nodes and translate
nods = np.genfromtxt(fl+'1.node',skip_header=1,skip_footer=1)
trn = np.loadtxt(fl+'trn')
nods[:,1] = nods[:,1] + trn[0]
nods[:,2] = nods[:,2] + trn[1]
nods[:,3] = nods[:,3] + trn[2]
nnods = len(nods[:,0])

#calculate element centroids
mids = np.zeros((nele,3))
for i in range(nele):
   for j in range(3):
      mids[i,j]=0.25*sum(nods[ele[i,1:5]-1,j+1])



sig=0.001 + np.zeros((nele,1)) # allocate a vector of conductivities

#zt is the top of the dipping plane at this centroids x-value
zt=.25*mids[:,0] - 1.75

#zb is the bottom of the dipping plane at this centroids x-value
zb=.25*mids[:,0] - 2.25

#Check each element to determine if the centroid is within the 
#dipping plane, and if so set the conductivity to 0.01 S/m.
for i in range(nele):
  #if the z-value of this centroid is between the top and bottom
  #of the dipping plane set the value to 0.01 S/m.
  if (mids[i,2]<=zt[i] and mids[i,2]>=zb[i]):
     sig[i]=0.01;
  


np.savetxt('mbsl_truesig.sig',sig,newline='\n',fmt='%-10.5e',header=str(nele)+' 1',comments='')



