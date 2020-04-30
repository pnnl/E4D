#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 13:22:10 2019

This script visualizes the mesh file embedded in the .xmf and .h5 file 
in this directory

Open visit from the command line and in the menu select
Controls -> Launch CLI

At the Visit command interface, type
Source("script.py")


@author: robi526
"""
# Name of file to open
fn = "two_blocks.xmf"
 
OpenDatabase(fn)
AddPlot("Pseudocolor", "Real_conductivity")
DrawPlots()
d = GetDomains()
TurnDomainsOff(d[3])

AddPlot("Mesh", "mesh")
DrawPlots()
d = GetDomains()
TurnDomainsOff(d[3])

# add the electrodes
OpenDatabase("electrodes.txt", 0)
ScatterAtts = ScatterAttributes()
ScatterAtts.var1 = "var01"
ScatterAtts.var1Role = ScatterAtts.Coordinate0  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var1MinFlag = 0
ScatterAtts.var1MaxFlag = 0
ScatterAtts.var1Min = 0
ScatterAtts.var1Max = 1
ScatterAtts.var1Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var1SkewFactor = 1
ScatterAtts.var2Role = ScatterAtts.Coordinate1  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var2 = "var02"
ScatterAtts.var2MinFlag = 0
ScatterAtts.var2MaxFlag = 0
ScatterAtts.var2Min = 0
ScatterAtts.var2Max = 1
ScatterAtts.var2Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var2SkewFactor = 1
ScatterAtts.var3Role = ScatterAtts.Coordinate2  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var3 = "var03"
ScatterAtts.var3MinFlag = 0
ScatterAtts.var3MaxFlag = 0
ScatterAtts.var3Min = 0
ScatterAtts.var3Max = 1
ScatterAtts.var3Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var3SkewFactor = 1
ScatterAtts.var4Role = ScatterAtts.None  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var4 = "default"
ScatterAtts.var4MinFlag = 0
ScatterAtts.var4MaxFlag = 0
ScatterAtts.var4Min = 0
ScatterAtts.var4Max = 1
ScatterAtts.var4Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var4SkewFactor = 1
ScatterAtts.pointSize = 0.05
ScatterAtts.pointSizePixels = 1
ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.scaleCube = 1
ScatterAtts.colorType = ScatterAtts.ColorByForegroundColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (255, 0, 0, 255)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 1
SetDefaultPlotOptions(ScatterAtts)
AddPlot("Scatter", "var01", 1, 1)
DrawPlots()

SetActivePlots(0)
ScatterAtts.pointSize = 0.05
ScatterAtts.pointSizePixels = 10
ScatterAtts.pointType = ScatterAtts.Sphere  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.scaleCube = 0
ScatterAtts.colorType = ScatterAtts.ColorBySingleColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (255, 255, 255, 255)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 1
SetPlotOptions(ScatterAtts)

# set the view
v = GetView3D()
print ("The view is: ", v)
v.viewNormal = (0, 0, 1)
v.viewUp = (0, 1, 0)
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

