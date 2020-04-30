#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 11:36:30 2019

@author: robi526
"""

# view the files in visit
import os
import subprocess



#cmd='visit' 
#result = subprocess.run(cmd, shell=True, check=True)

#cmd='/home/robi526/codes/visit2_13_0.linux-x86_64/bin/visit -cli -s script.py' 


cmd='/home/robi526/codes/visit2_13_0.linux-x86_64/bin/visit' 
subprocess.run(['/home/robi526/codes/visit2_13_0.linux-x86_64/bin/visit' , 'cli', '-s', 'script.py'])



#result = subprocess.run(cmd, check=True, shell=True)