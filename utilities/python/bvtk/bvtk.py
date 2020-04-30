#!/usr/bin/python
import sys
import io
import os
import time
import datetime
import glob
import numpy as np

#Function definitions--------------------------------------------------------------
def check_args(strv):
    #this function makes sure there are enough
    #arguments on the command line
    passed = True
    if (len(sys.argv) <= 4):
        print ("\n-Not enough arguments on the command line")
        passed = False
        
    if (len(sys.argv) > 1):
        if (sys.argv[1] != '-f' and  sys.argv[1] !='-af'):
            print ("\n-The specified flag is "+sys.argv[1])
            passed = False

    if (len(sys.argv) > 2):
        mfcheck = True
        if(not os.path.isfile(sys.argv[2]+".1.node")):
            mfcheck = False
        if(not os.path.isfile(sys.argv[2]+".1.ele")):
            mfcheck = False
        if(not os.path.isfile(sys.argv[2]+".trn")):
            mfcheck = False
        if(not mfcheck):
            print ("\n-Specified mesh file prefix is: "+sys.argv[2])
            print (" Cannot find one of the following mesh files ...")
            print (" "+sys.argv[2]+".1.node")
            print (" "+sys.argv[2]+".1.ele")
            print (" "+sys.argv[2]+".trn")
            passed = False

    if(len(sys.argv) > 3):
    	if glob.glob(sys.argv[3]+'*')==[]:
    	    print ("\n-Cannot find any parameter files starting with "+sys.argv[3])
    	    passed = False
 
    if(len(sys.argv) > 4 and sys.argv[1]=='-af' and passed):
        if(not os.path.isfile(sys.argv[4])):
            print ("\n-Append options -f specified, but I")
            print (" cannot find the specified vtk file: "+sys.argv[4])
            passed = False
    
    if(len(sys.argv) == 5 and passed):
        print ("\n-No time stamp specified, assuming starting time 0 and recording all zones")
    elif(len(sys.argv) > 5 and passed):
    	print ("\n-Starting Time: "+sys.argv[5])
    
    if(len(sys.argv) > 6 and passed):
        print ("Recording only the following zones ...")
    	for i in range(6,len(sys.argv)):
    	    print ("  Zone: "+sys.argv[i])
    	      
    if (not passed):
        print ("\nCalling sequence: bvtk -flag mesh_pre param_file vtk_out tstamp")
        print ("where ...")
        print ("  flag = f or af (f for new vtk file, af for append to existing)")
        print ("  mesh_pre = mesh files prefix")
        print ("  param_file = mesh or node centered e4d file name")
        print ("  vtk_out = name of output file (must exist for flag = af)")
        print ("  tstamp = parameter file time stamp")
        print ("\n")
        
    return passed
#end of check args---------------------------------------------------------

def load_nodes(mesh_pre):
    f1 = open(mesh_pre+'.1.node','r')
    str = f1.readline().split()
    nnods = int(str[0])
    f1.close
    trn = np.genfromtxt(mesh_pre+'.trn')
    nods = np.genfromtxt(mesh_pre+'.1.node',skip_header=1,skip_footer=1,usecols=(1,2,3,5))
    #nod_flags = np.genfromtxt(mesh_pre+'.1.node',skip_header=1,skip_footer=1,usecols=(5))
    nods[:,0]=nods[:,0]+trn[0]
    nods[:,1]=nods[:,1]+trn[1]
    nods[:,2]=nods[:,2]+trn[2] 
    return nnods,nods
#end of load_nodes---------------------------------------------------------

def load_elements(mesh_pre):
    f1 = open(mesh_pre+'.1.ele','r')
    str = f1.readline().split()
    nele = int(str[0])
    f1.close()
    ele=np.genfromtxt(mesh_pre+'.1.ele',skip_header=1,skip_footer=0,usecols=(1,2,3,4,5))
    return nele,ele
#end of load_elements--------------------------------------------------------- 

def load_parameters(fname):
    f1 = open(fname,'r')
    str = f1.readline().split()
    nparm = int(str[0])
    ncol = int(str[1])
    f1.close()

    ucols = (0)
    if (ncol==2):
        ucols=(0,1)

    parm=np.genfromtxt(fname,skip_header=1,usecols=ucols)
    return nparm,parm
#end of load_parameters--------------------------------------------------------- 
 
def build_new_vtk(nnods,nele,nods,ele,fnout):
    
    #determine which zones to print
    if(len(sys.argv)>6):
    	nzp = len(sys.argv)-6
    	cnt=0
    	zp = np.zeros(nzp)
    	for zn in range(6,len(sys.argv)):
    	   zp[cnt] =  sys.argv[zn]
    	   cnt=cnt+1
    else:
    	nzp = int(max(ele[:,4]))
    	zp = np.zeros(nzp)
    	for i in range(nzp):
    	   zp[i] = i+1
    	   print zp[i]
    
    #loop over the parameter files
    cnt = 0
    
    f1=open('.visit','w')
    f1.write('!NBLOCKS %-10.0f\n' % nzp)
    for pfile in glob.glob(sys.argv[3]+'*'):
    	nparms,parms = load_parameters(pfile)
    	
    	if(nparms == nele or nparms == nnods):
    	   #build the filenames for each saved zone of this parameter file
    	   filenames = [sys.argv[4]+'_zn'+str(int(zind))+'_t'+str(cnt)+'.vtk' for zind in zp]
    	   for fn in filenames:
    	   	f1.write(fn+'\n')
    	   #open the files
           file_descriptors = [open(filename, 'w') for filename in filenames]
           zcnt = 0
           for fn in file_descriptors:
           	fn.write('# vtk DataFile Version 3.0\n')
	   	fn.write('This is an E4D vis file constructed using bvtk.py\n')
	   	fn.write('ASCII\n')
	   	fn.write('DATASET UNSTRUCTURED_GRID\n')
	   	fn.write('FIELD FieldData 1\n')
		fn.write('TIME 1 1 double\n')
	        fn.write('%f\n' % float(sys.argv[5]))
	   	fn.write('POINTS %-10.0f float\n' % nnods)
	   	for i in range(nnods):
	   	    fn.write('%f %f %f\n' % (nods[i,0],nods[i,1],nods[i,2]))
                fn.write('\n')
           	
           	zinds = np.where(ele[:,4]==zp[zcnt])
           	nelez = len(zinds[0])
           	nelez5 = 5*int(nelez)
           	fn.write('CELLS %-10.0f %-10.0f \n' % (nelez,nelez5) )
       	
       		ecnt = 0
           	for i in zinds[0]:
        	    fn.write('4 %-10.0f %-10.0f %-10.0f %-10.0f\n' %(ele[i,0]-1,ele[i,1]-1,ele[i,2]-1,ele[i,3]-1))
           	
           	fn.write('\nCELL_TYPES %-10.0f\n' % nelez)
           	fn.write('%-3.0f' %10*nelez)
           	
           	fn.write('\nCELL_DATA %-10.0f\n' %nelez)
           	fn.write('SCALARS CellData float 1\n')
           	fn.write('LOOKUP_TABLE default\n')
           	for i in zinds[0]:
           	    fn.write('%f\n'%parms[i])
           	fn.close()
           	zcnt = zcnt+1
        cnt=cnt+1
    f1.close()
    quit()
    
  
lpass = True
check_args(lpass)
if(not lpass):
   quit()
nnods,nods=load_nodes(sys.argv[2])
nele,ele=load_elements(sys.argv[2])
build_new_vtk(nnods,nele,nods,ele,sys.argv[4])