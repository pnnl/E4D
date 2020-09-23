# validation of inputs for px visualizaton module


import os
import sys
import re


def check_command_input(val):
    varName=''
    if len(val)<6:
        print ("Not enough input variables entered\n")
        print ("px mode[-af or f] meshprefix visfile(s) output[.h5 .xml] timestamp variable_name (optional)\n")
        sys.exit(1)

    mode, grid, visF, outF, timeStamp =([] for i in range(5))
        
    mode=val[1]
    if not (re.search(mode, '-f')) and not (re.search(mode, '-af')):
        print ("mode needs to be entered as -f or -af")
        print ("px mode[-af or f] meshprefix visfile(s) output[.h5 .xml] timestamp(optional) variable_name (optional)\n")
        sys.exit()
    grid = val[2]
    visF = val[3]
    outF = val[4]
      
    if len(val)>=6: # user entered timestamp
        timeStamp=str(val[5])
    else:
        timeStamp='0'

    if len(val)>=7: # user entered variable name
        varName=str(val[6])
        if len(val)==8: # user entered two variable names
                varName=varName+ ',' + (str(val[7]))
        print ('Variable name entered:' + varName)
    else:
        varName=''         
        
    return mode, grid, visF, outF, timeStamp, varName
    
def check_file_input(grid, visF):
    if not os.path.isfile(grid+".1.node"):
       print ("Cannot find the node mesh file: ",grid, ".1.node")
       sys.exit(1)

    if not os.path.isfile(grid+".1.ele"):
       print ("Cannot find the element mesh file: ",grid, ".1.ele")
       sys.exit(1)

    if not os.path.isfile(grid+".trn"):
       print ("Cannot find the mesh translation file: ",grid, ".trn")
       sys.exit(1)
    
    if not os.path.isfile(visF):
       print ("Cannot find the input file: ",visF)
       sys.exit(1)

def check_mode(mode, outF):
    newF=False
    if re.search(mode,'-f'):
        print ('Creating new xmf and h5 file.....')
        newF=True
    elif re.search(mode,'-af'):
        if not os.path.isfile(outF +'.h5'):
            print (outF+'.xmf/.h5 does not exist.... creating new file')
            newF=True
        else:
            print ('Updating '+ outF+'.xmf and '+ outF+'.h5.....')
    return newF            


