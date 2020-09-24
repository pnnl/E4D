# xdmf outputs for px visualizaton module


def xmlStart(lines):
    lines.append( '<?xml version="1.0" ?>\n' )
    lines.append('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
    lines.append('<Xdmf Version="2.0">\n')
    lines.append('<Domain>\n') 
    
    return lines



def write_xdmf_element (lines, zn, zn_size, nnods, outFilepre, dim, timeStamp,vn1,vn2):
    lines.append('  <Grid Name = "'+str(zn)+'">\n')
    lines.append('    <Topology TopologyType="Tetrahedron" NumberOfElements="'+str(zn_size)+'">\n')
    lines.append('      <DataItem Dimensions="'+str(zn_size)+' 4" Format="HDF">\n')    
    lines.append('      '+outFilepre+'.h5:/Mesh/Time '+ timeStamp+ ' Zone '+str(zn)+'\n')        
    lines.append('      </DataItem>\n')
    lines.append('    </Topology>\n\n')
    
    #record the nodes
    lines.append('  <Geometry GeometryType="XYZ">\n')
    lines.append('    <DataItem Dimensions= "'+str(nnods)+' 3" Format="HDF">\n')
    lines.append('    '+outFilepre+'.h5:/Mesh/Nodes\n')
    lines.append('    </DataItem>\n')
    lines.append('  </Geometry>\n\n')
    
    lines.append('    <Attribute Name="'+vn1+'" Center="Cell">\n')
    lines.append('      <DataItem Dimensions="'+str(zn_size)+'" Format="HDF">\n')    
    lines.append('      '+outFilepre+'.h5:/Mesh/Time '+timeStamp+' Real '+str(zn)+'\n')
        
    lines.append('      </DataItem>\n')
    lines.append('    </Attribute>\n\n')           
    
    if dim==2:
        lines.append('    <Attribute Name="'+vn2+'" Center="Cell">\n')
        lines.append('      <DataItem Dimensions="'+str(zn_size)+'" Format="HDF">\n')        
        lines.append('      '+outFilepre+'.h5:/Mesh/Time '+timeStamp+' Imag '+str(zn)+'\n')
        lines.append('      </DataItem>\n')
        lines.append('    </Attribute>\n\n')
    
    lines.append('</Grid>\n')
    return lines
                
                
def write_xdmf_node(lines, nele, nnods, outFilepre, dim, timeStamp, vn1, vn2):    
    lines.append('<Grid Name = "mesh" GridType="Collection" CollectionType="Spatial">\n\n') # main grid
    lines.append('  <Time TimeType="Single" Value="'+timeStamp+'"/>\n')
    lines.append('  <Grid Name = "'+str(timeStamp)+'">\n')
    lines.append('    <Topology TopologyType="Tetrahedron" NumberOfElements="'+str(nele)+'">\n')
    lines.append('      <DataItem Dimensions="'+str(nele)+' 4" Format="HDF">\n')
    lines.append('      '+outFilepre+'.h5:/Mesh/Elements\n')
    lines.append('      </DataItem>\n')
    lines.append('    </Topology>\n\n')
    
    #record the nodes
    lines.append('  <Geometry GeometryType="XYZ">\n')
    lines.append('    <DataItem Dimensions= "'+str(nnods)+' 3" Format="HDF">\n')
    lines.append('    '+outFilepre+'.h5:/Mesh/Nodes\n')
    lines.append('    </DataItem>\n')
    lines.append('  </Geometry>\n\n')
    
    lines.append('    <Attribute Name="'+vn1+'" Center="Node">\n')
    lines.append('      <DataItem Dimensions="'+str(nnods)+'" Format="HDF">\n')
    lines.append('      '+outFilepre+'.h5:/Mesh/'+ 'Time '+str(timeStamp)+' Real\n')
    lines.append('      </DataItem>\n')
    lines.append('    </Attribute>\n\n')           
    
    if dim==2:
        lines.append('    <Attribute Name="'+vn2+'" Center="Node">\n')
        lines.append('      <DataItem Dimensions="'+str(nnods)+'" Format="HDF">\n')
        lines.append('      '+outFilepre+'.h5:/Mesh/Time ' + str(timeStamp)+' Imag\n')
        lines.append('      </DataItem>\n')
        lines.append('    </Attribute>\n\n')
    
    lines.append('</Grid>\n')
    lines.append('</Grid>\n')
    return lines
    