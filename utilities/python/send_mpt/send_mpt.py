import sys
import io
import os
import time
import datetime

iostat = 0
# see if an options filename was specified
if len(sys.argv) < 2:
    print 'On the command line, please specify an input file name that'
    print 'has the following format ...\n'
    print 'line 1: The prefix of the MPT data file'
    print 'line 2: The name of the remote computer'
    print 'line 3: The absolute path name of the directy where the file'
    print '        should be copied onto the remote computer'
    iostat = -1
else:
    #get the input file name
    inp_file = sys.argv[1]
    
    #see if the input file exists
    if(os.path.isfile(inp_file)):
       print '\nReading options from the input file:',inp_file
    else:
       print '\nCannot find the inpute file:',inp_file
       iostat = -1
       
if(iostat != 0):
    exit(-1)

#open and read the input file
f = open(inp_file,'r')
mpt_pre = f.readline()
dest_comp = f.readline()
dest_dir = f.readline()
f.close

#strip the newline characters
n=len(mpt_pre)
mpt_pre = mpt_pre[:n-1]
npre = len(mpt_pre)

n=len(dest_comp)
dest_comp = dest_comp[:n-1]
ncp = len(dest_comp)

n=len(dest_dir)
dest_dir = dest_dir[:n-1]
ndr = len(dest_dir)

print '\nMPT prefix: ',mpt_pre
print 'Destination computer: ',dest_comp
print 'Destination directory: ',dest_dir,'\n'

#make a directory for the data files that have been trasnferred
exp_dir = 'Exported_Data'
if not os.path.exists(exp_dir):
     os.mkdir(exp_dir)

#loop to check for new MPT data files
while True:
   print (time.strftime("%I:%M:%S")),'checking for new data files'
   for fil in os.listdir('.'):
     if os.path.isfile(fil):
     	pre = fil[:npre-1]
        if (pre == mpt_pre[:npre-1]) & fil.endswith('.Data'):
           #make sure the file is done
           done = False
           if os.access(fil,os.R_OK):
              with open(fil,'r') as f:
                 for line in f:
                     if 'Complete' in line:
          	        print '\nsending',fil,'to',dest_comp
                        done = True
                        
                        #send this file to the remote computer
                        str1 = 'scp '+fil+' '+dest_comp[:ncp-1]+':'+dest_dir[:ndr-1]
                        #print(str1)
                        os.system(str1)
                      
                        #move this file to exp_dir
                        os.system('mv '+fil+' '+exp_dir)
                 if not done:
                     print fil,'is incomplete'
           else:
               print 'No read access for ',fil
   time.sleep(5)   
      
