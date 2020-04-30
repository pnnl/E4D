#!/bin/python

# This python script builds the conductivity file for to be used in this example

import numpy as np


# Load the mesh elements
fl = 'sp'
ele = np.loadtxt(fl+'.1.ele',skiprows=1,dtype='i4')
nele=len(ele[:,0]) # number of elements in mesh

#load the nodes and translate
nods = np.genfromtxt(fl+'.1.node',skip_header=1,skip_footer=1)
trn = np.loadtxt(fl+'.trn')
nods[:,1] = nods[:,1] + trn[0]
nods[:,2] = nods[:,2] + trn[1]
nods[:,3] = nods[:,3] + trn[2]
nnods = len(nods[:,0])

#calculate element centroids
# The centroid is the average position of the 4 nodes making the element.
mids = np.zeros((nele,3))
for i in range(nele):
   for j in range(3):
      mids[i,j]=0.25*sum(nods[ele[i,1:5]-1,j+1])


# baseline conductivity
sig=0.001 + np.zeros((nele,1)) # allocate a vector of conductivities
np.savetxt('baseline.sig',sig,newline='\n',fmt='%-10.5e',header=str(nele)+' 1',comments='')

xmid=0
ymid=0
R=2
nfile=0
for zmid in np.arange(0,-10.5,-0.5):
    sig=0.001 + np.zeros((nele,1)) # initialize a vector of conductivities
    r=np.sqrt((mids[:,0]-xmid)**2 + (mids[:,1]-ymid)**2 + (mids[:,2]-zmid)**2)
    r=np.reshape(r, (len(mids), 1))
    #find the elements where r<=R
    inds = np.argwhere(r<=R)
    sig[inds]=0.1

    #find elements where 1.0<=r<= R
    inds=np.all((r>=1, r<=R), axis=0)
    sig[inds] = 0.1*(r[inds]**(-4))

    fout='sig_'+str(abs(zmid))+'.sig'
    np.savetxt(fout,sig,newline='\n',fmt='%-10.5e',header=str(nele)+' 1',comments='')
    nfile+=1
    
#%% Build visualization file
import subprocess

cmd='px -f sp baseline.sig models -1'
subprocess.run(cmd, shell=True,stdout=subprocess.PIPE)
    
inc=0
for i in range(nfile):
    fn='sig_'+str(inc)+'.sig'
    inc+=0.5
    print (fn)

    cmd='px -af sp '+fn+' models ' + str(inc)
    subprocess.run(cmd, shell=True,stdout=subprocess.PIPE)