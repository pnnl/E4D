# xdmf outputs for px visualizaton module

import h5py
import sys
import numpy as np



def h5_init(newF, outF, nele, ele, nnods, nods, grid):
    if newF:
        #open an hdf5 file
        f5 = h5py.File(outF+'.h5','w')
        print("Creating "+outF+ '.h5 file')
        
        #create the mesh group and items within this group
        mesh_group = f5.create_group('Mesh')
        mesh_group.create_dataset("Nodes",data=nods)        
        mesh_group.create_dataset('Elements',data=ele)                    
        
        f5.close()
    
    else:
        f5 = h5py.File(outF+'.h5','r+')    
        grp=f5.get('Mesh')
        item=grp.get('Elements')
        file_ele=np.array(item).shape
        
        if file_ele[0]!=nele:
            print ("The number of elements in " + outF+'.h5 does not match '+ grid + '.1.ele')
            print ("Aborting....")
            sys.exit()
            
        item2=grp.get('Nodes')
        file_nods=np.array(item2).shape
    
        if file_nods[0]!=nnods:
            print ("The number of nodes in " + outF+'.h5 does not match '+grid+'.1.node')
            print ("Aborting....")
            sys.exit()
        
        print ("All okay with existing file: "+ outF+'.h5.....')

        f5.close()
