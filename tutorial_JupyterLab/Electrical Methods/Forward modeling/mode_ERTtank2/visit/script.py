#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 13:22:10 2019

This script visualizes the potential file embedded in the .xmf and .h5 file 
in this directory


Open visit from the command line and in the menu select
Controls -> Launch CLI

At the Visit command interface, type
Source("script.py")


@author: robi526
"""
# Name of file to open
fn = "pot1_tank.xmf"

 
OpenDatabase(fn)
AddPlot("Pseudocolor", "Real_potential")
iso_atts = IsosurfaceAttributes()
iso_atts.contourMethod = iso_atts.Value
AddOperator("Isosurface")

DrawPlots()
v = GetView3D()
print ("The view is: ", v)
v.viewNormal = (0, -1, 0)
v.viewUp = (0, 0, 1)
SetView3D(v)

aatts = AnnotationAttributes()
#aatts.axes3D.xAxis.title.title='X'
aatts.axes3D.xAxis.title.userTitle=1
aatts.axes3D.xAxis.title.title='X (m)'

aatts.axes3D.yAxis.title.userTitle=1
aatts.axes3D.yAxis.title.title='Y (m)'

aatts.axes3D.zAxis.title.userTitle=1
aatts.axes3D.zAxis.title.title='Z (m)'

SetAnnotationAttributes(aatts)

