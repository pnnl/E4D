#!/bin/python


# Imports mesh files and conductivity output files and outputs two files 
# for viewing in ViSit or ParaView
# 1) h5 file
# 2) xdmf file (this file references h5 file)


# input variables
# - mesh prefix 
# - filename of text file containing time lapse files OR a single conductivity file
# - output prefix for xmf and h5 file
import numpy as np
import h5py
import re
import sys

from input_validation import check_command_input
from input_validation import check_file_input
from input_validation import check_mode

from h5_outputs import h5_init

from utils import file_num
from utils import file_text_insert

from xdmf_outputs import xmlStart
from xdmf_outputs import write_xdmf_element
from xdmf_outputs import write_xdmf_node


# *****INPUT VALIDATION (input_validation.py)*****
# check command line input
[mode, grid, visF, outF, timeStamp, varName]= check_command_input(sys.argv)

# check input - make sure files exist
check_file_input(grid, visF)

# check if creating a new file of adding to existing file
newF=check_mode(mode, outF)

# *****

# *****MESH IMPORT AND VALIDATION IN h5 FILE ******
# read the mesh translation file
trn = np.genfromtxt(grid+'.trn',dtype='float')

# load the node coordinates and translate
# get the first line of the node file to check number of nodes loaded in matrix
f=open(grid+'.1.node')
first=f.readline()
first=first.split()
nnods_header=int(first[0])

nods = np.genfromtxt(grid+'.1.node',skip_header=1,skip_footer=1,usecols=(1,2,3),dtype='float')
nnods = np.size(nods[:,0])

# this is necessary if node file was generated through a .vol file refinement in tetgen
if nnods != nnods_header:
    nods = np.genfromtxt(grid+'.1.node',skip_header=1,usecols=(1,2,3),dtype='float')
    nnods = np.size(nods[:,0])

nods[:,0]=nods[:,0]+trn[0]
nods[:,1]=nods[:,1]+trn[1]
nods[:,2]=nods[:,2]+trn[2]

print("Reading "+grid+ '.1.node')

#read and load the elements
ele = np.genfromtxt(grid+'.1.ele',skip_header=1,usecols=(1,2,3,4),dtype='int')
ele = ele - 1
nele = np.size(ele[:,0])
zone = np.genfromtxt(grid+'.1.ele',skip_header=1,usecols=5,dtype='i')
minzone = np.min(zone)
maxzone = np.max(zone)
print("Reading "+grid+ '.1.ele')

print("THE NUMBER OF NODES IS:" + str(nnods))
print("THE NUMBER OF ELEMENTS IS:" + str(nele))
print("THE NUMBER OF ELEMENT ZONES IS:" + str((maxzone-minzone)+1))

# initialize h5 or validate existing h5 file (h5_outputs.py)
h5_init(newF, outF, nele, ele, nnods, nods,grid)
# ***** 

# List or single file?  Potential or Conductivity file? (utils.py)
[eleF, nodeF, file_list, timeStamp]=file_num(visF, nele, nnods, timeStamp)

flines=list()

# Create and/or update xml and complete h5 files in tandem  
# associated xmf file
tmp=[]
if newF: # create new file
    flines=xmlStart(flines)
    flines.append('<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n')
else: # update existing h5 and xml files
    f2=open(outF+'.xmf', "r")
    
    # new xml lines belong after last time step
    line_cnt=0
    time_rows=[]
    print ("Looping through previous time steps:")
    for line in f2:
        if 'TimeType="Single"' in line:
               lineList=line.split()     
               val=re.findall(r'\"(.+?)\"', lineList[2])
               if tmp!=val[0]:
                   print ('Time step: ' + val[0])
               time_rows.append(line_cnt)
               
               for i in np.arange(len(timeStamp)):
                   if str(val[0])==timeStamp[i]:
                       print('The time stamp definition has already been defined for:' + val[0])
                       print ('Aborting....')
                       sys.exit()
               tmp=val[0]
        line_cnt+=1
    f2.close()

    f2=open(outF+'.xmf', "r") 
    contents=f2.readlines()
    f2.close()
    size=len(contents)
    insert_row=size-3 # insert before #finish up xmf write, 3 lines before end of file

    # This file must have the same variable name definition as the previous file
    listToStr = ' '.join([str(elem) for elem in contents]) 
    string = ''.join(listToStr.split())
    
    # first variable name defintion in xmf file
    attr_str = re.search('</Geometry><AttributeName="(.*?)"', string)

    
    # second variable name definition, optional
    attr_str2 = re.search('</Attribute><AttributeName="(.*?)"', string)

    if attr_str2:
        attr=attr_str.group(1)+ ',' + attr_str2.group(1)
    else:
        attr=attr_str.group(1)
        

    if varName!='':  # see if there was a user command line entry      
        if attr!=varName:
            print ('The variable name definition of this file must be ' + attr + ' NOT ' + varName)
            print ('Setting variable name defnition to: ' + attr)
            varName=attr
    else:
        varName=attr # set variable name to the same as previous xmf file

        
# append h5 file
f5 = h5py.File(outF+'.h5','a')
for i in range(len(file_list)):
    var = np.genfromtxt(file_list[i], skip_header=1,dtype='float')  
    if nodeF:    
        # h5 file
        if varName=='':
           vn1='Real_potential'
           vn2='Imag_potential'
        else:
          try:
             tmp=varName.split(',') 
             vn1=tmp[0]
             vn2=tmp[1]
          except:
             vn1=varName
             vn2=''
        
        
        if var.ndim==1:
            if len(var)!=nnods:
                print("Error in node file.  The number of node values does not match the number of nodes in the mesh.")
                sys.exit()
            
            print ('Recording node file ' +file_list[i] + ' at time stamp:'+timeStamp[i])
            f5.create_dataset('/Mesh/Time '+ timeStamp[i]+' Real',data=var)
        elif var.ndim==2:    # contains complex potentials
            if len(var)!=nnods:
                print("Error in complex node file.  The number of complex node values does not match the number of nodes in the mesh.")
                sys.exit()

            print ('Recording complex node file ' + file_list[i] + ' at time stamp:'+timeStamp[i])
            f5.create_dataset('/Mesh/Time '+ timeStamp[i]+' Real',data=var[:,0])
            f5.create_dataset('/Mesh/Time '+timeStamp[i]+' Imag',data=var[:,1])
            

        flines.append('  <Time TimeType="Single" Value="'+timeStamp[i]+'"/>\n')    
        # associated xdmf file
        flines=write_xdmf_node(flines,nele, nnods, outF, var.ndim, timeStamp[i],vn1,vn2)   

    elif eleF: # element file
        # h5 file
        if len(var)!=nele:
            print("Error in element file.  The number of element values does not match the number of elements in the mesh.")
            sys.exit()
        zone_size = np.zeros(maxzone-minzone+1,dtype=int)
        
        if varName=='':
           vn1='Real_conductivity'
           vn2='Imag_conductivity'
        else:
          try:
             tmp=varName.split(',') 
             vn1=tmp[0]
             vn2=tmp[1]
          except:
             vn1=varName
             vn2=''
                 
        cnt = 0
        if var.ndim==1:
            print ('Recording element file ' + file_list[i] + ' at time stamp:'+ timeStamp[i])
            f5.create_dataset('/Mesh/Time '+timeStamp[i]+' AllZone', data=ele)
            f5.create_dataset('/Mesh/Time '+timeStamp[i]+' AllReal', data=var)

            for zn in range(minzone,maxzone+1):
                inds=np.where(zone==zn)
                zone_size[cnt] = np.size(inds)
                if(zone_size[cnt]>0):
                    f5.create_dataset('/Mesh/Time '+timeStamp[i]+' Zone '+str(zn),data=ele[inds])
                    f5.create_dataset('/Mesh/Time '+timeStamp[i]+' Real '+str(zn),data=var[inds])
                    
                    
                    cnt = cnt+1
        elif var.ndim==2: # contains complex conductivities
            if len(var)!=nele:
                print("Error in complex element file.  The number of complex element values does not match the number of elements in the mesh.")
                sys.exit()
            print ('Recording complex element file ' + file_list[i] + ' at time stamp:'+timeStamp[i])
            f5.create_dataset('/Mesh/Time '+timeStamp[i]+' AllZone', data=ele)
            f5.create_dataset('/Mesh/Time '+timeStamp[i]+' AllReal', data=var[:,0])
            f5.create_dataset('/Mesh/Time '+timeStamp[i]+' AllImag', data=var[:,1])
            for zn in range(minzone,maxzone+1):
                inds=np.where(zone==zn)
                zone_size[cnt] = np.size(inds)
                if(zone_size[cnt]>0):
                    f5.create_dataset('/Mesh/Time '+timeStamp[i]+' Zone '+str(zn),data=ele[inds])
                    f5.create_dataset('/Mesh/Time '+timeStamp[i]+' Real '+str(zn),data=var[inds,0])
                    f5.create_dataset('/Mesh/Time '+timeStamp[i]+' Imag '+str(zn),data=var[inds,1])
                    cnt = cnt+1
                        
        
        # xmf file
        if maxzone!=minzone:
            flines.append('<Grid Name = "mesh" GridType="Collection" CollectionType="Spatial">\n\n') # main grid
        flines.append('  <Time TimeType="Single" Value="'+timeStamp[i]+'"/>\n')
        cnt = 0
        for zn in range(minzone,maxzone+1):
            if(zone_size[cnt]>0):
                lines=list()
                flines=write_xdmf_element (flines, zn, zone_size[cnt], nnods, outF, var.ndim, timeStamp[i], vn1, vn2)
            cnt = cnt+1
            
        if minzone-maxzone!=0:
            flines.append('</Grid>\n')
            
    else:
        print("The number of elements in "+ file_list[i]+ " does not match the mesh (.ele) file")
        print("Skipping file.....")
        continue

# done with h5
f5.close()
print("Done writing "+outF+ '.h5 file')


#finish up xmf write
if newF:
    flines.append('</Grid>\n')
    flines.append('</Domain>\n')
    flines.append('</Xdmf>\n')
    
    f1=open(outF+".xmf","w")
    file_text_insert(f1, flines)
    f1.close()
    
    print("Finished build of "+outF+ '.xmf file')
    
else: # insert new lines
    for item in flines:
        contents.insert(insert_row, item)
        insert_row+=1    

    f1=open(outF+'.xmf', "w")
    contents = "".join(contents)
    f1.write(contents)
    f1.close()

    print("Finished rebuild of "+outF+ '.xmf file')
    
    
    
   