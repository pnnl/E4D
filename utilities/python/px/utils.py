
from itertools import islice
import sys
import os
import re
import numpy as np

def file_num(visF, nele, nnods, time_tmp):

    with open(visF) as inputfile:
        firstLine = list(islice(inputfile, 1))
        listFirstLine=re.findall('\d+', str(firstLine))
        num=int(listFirstLine[0])

    sigF=False
    potF=False
    file_list=[]
    file_list_tmp=[]
    timeStamp=[]
    if num==nele or num==nnods: # single file
       if num==nele:  # conductivity/complex conductivity file
           sigF=True 
       if num==nnods: # conductivity/complex potential file
           potF=True
       if time_tmp is None:
           timeStamp.append('0') 
       else:
           timeStamp.append(time_tmp)
       file_list.append(visF)
       
    else:
        # check input file list
        cnt_skip=0
        with open(visF) as inputfile:
            line=inputfile.readline()
            print (line)
            cnt=0
            while line:
                if cnt>0:
                    lineList=line.split()
                    if lineList:
                        file_list_tmp.append(lineList[0])
                        timeStamp.append(lineList[1])                
    
                        #keep a count of missing files
                        if not os.path.isfile(lineList[0]):
                            cnt_skip+=1
                            
                line=inputfile.readline()
                cnt+=1
        
        if cnt_skip==len(file_list_tmp): # make sure at least one file exists
            print("All files listed in "+visF+ " do not exist.")
            sys.exit(1)
            
        
        if num!=len(file_list_tmp):
            print("The files listed in "+visF+ " does not match the number reported.")
            print("Aborting.....")
            sys.exit(1)
    
        for i in range(len(file_list_tmp)):
            if not os.path.isfile(file_list_tmp[i]):
               print ("Cannot find the input file: "+file_list_tmp[i])
               print ("Skipping over file: "+file_list_tmp[i])
               continue
    
            var = np.genfromtxt(file_list_tmp[i], skip_header=1,dtype='float')  
       
            if len(var)==nnods: # potential file     
                potF=True
                file_list.append(file_list_tmp[i])
            elif len(var)==nele: # conductivity file
                sigF=True
                file_list.append(file_list_tmp[i])
            else:
                print ("The file: " + file_list_tmp[i] + ' does not have the correct number of values')
                print ("Skipping over file: "+file_list_tmp[i])
                            
    return (sigF, potF, file_list, timeStamp )


def file_text_insert(f1, flines):
    # add contents to file
    for item in flines:
        f1.write(item)

   
    
