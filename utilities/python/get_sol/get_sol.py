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
    print 'line 1: The prefix of the sigma files'
    print 'line 2: The name of the remote computer'
    print 'line 3: The absolute path name of the directy where the sigma files'
    print '        should be pulled from remote computer'
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

print '\nFile prefix: ',mpt_pre
print 'Origin computer: ',dest_comp
print 'Origing directory: ',dest_dir,'\n'

#make a directory for the data files that have been trasnferred
imp_dir = 'Imported'
if not os.path.exists(imp_dir):
     os.mkdir(imp_dir)

#reused strings
str1 = 'ssh '+dest_comp[:ncp-1]
str2 = ' ls -1 '+dest_dir[:ndr-1]
str3 = 'scp '+dest_comp[:ncp-1]+':'+dest_dir[:ndr-1]
etm=0;
dt=4;
#loop to check for new MPT data files
while True:
   print (time.strftime("%I:%M:%S")),'checking for sigma files'
   
   #use ssh to query the desination directory and make a list of the files in the destination directory
   os.system(str1+str2+' > origin_files.txt') 
   ofl=open('origin_files.txt','r')
   olist = ofl.readlines()
   ofl.close()
   
   #make a list of the files we have
   os.system('/bin/ls -1 Imported > imported.txt')
   ifl = open('imported.txt','r')
   ilist = ifl.readlines()
   ifl.close
   
   #see if there are any new files
   for ofl in olist:
   	if ofl not in ilist:
   	   print('\nDownloading '+ofl.rstrip())
   	   #compress the file on the remote system
   	   #os.system(str1+' tar -cvzf '+dest_dir.rstrip()+ofl.rstrip()+'.tgz '+dest_dir.rstrip()+'/'+ofl.rstrip())
   	   
   	   #download the compressed file
  	   #os.system(str3+'/'+ofl.rstrip()+'.tgz .')
  	   
  	   #remove the compressed file from the remote system
  	   #os.system(str1+' rm '+dest_dir.rstrip()+'/'+ofl.rstrip()+'.tgz ')
  	   
  	   #decompress the file on the locale system
  	   #print('tar -xvzf '+ofl.rstrip()+'.tgz')
           #os.system('tar -xvzf '+ofl.rstrip()+'.tgz')
  	   
  	   #download the sigma file
  	   os.system(str3+'/'+ofl.rstrip()+' .')
  	   
  	   #see if the exodus file exists and build it if it doesn't
  	   exofil = mpt_pre.rstrip()+'.exo'
  	   if not os.path.isfile(exofil):
  	   	print('\nBuilding exodus file '+exofil)
  	   	os.system('bx -f tank2.1 '+ofl.rstrip()+' tmp.exo '+str(etm))
  	   	#os.system('mv tmp.exo '+exofil)
  	   	os.system('mv tmp.exo tmp2.exo')
  	   else:
  	        os.system('cp '+exofil+' tmp.exo')
  	        etm=etm+dt;
  	        os.system('bx -af tank2.1 '+ofl.rstrip()+' tmp.exo '+str(etm))
  	        #os.system('mv tmp.exo '+exofil)
  	        os.system('mv tmp.exo tmp2.exo')
  	   #mv the file to Imported
  	   os.system('mv '+ofl.rstrip()+' Imported');
  	   
  #sleep
   time.sleep(5)
  	   
  
   #for fil in os.listdir('.'):
   #  if os.path.isfile(fil):
   #  	pre = fil[:npre-1]
   #     if (pre == mpt_pre[:npre-1]) & fil.endswith('.Data'):
   #        #make sure the file is done
   #        done = False
   #        if os.access(fil,os.R_OK):
   #           with open(fil,'r') as f:
   #              for line in f:
   #                  if 'Complete' in line:
   #       	        print '\nsending',fil,'to',dest_comp
   #                     done = True
   #                     
   #                     #send this file to the remote computer
   #                     str1 = 'scp '+fil+' '+dest_comp[:ncp-1]+':'+dest_dir[:ndr-1]
   #                     #print(str1)
   #                     os.system(str1)
   #                   
   #                     #move this file to exp_dir
   #                     os.system('mv '+fil+' '+exp_dir)
   #              if not done:
   #                  print fil,'is incomplete'
   
   #else:
   #            print 'No read access for ',fil
   #time.sleep(5)   
      
