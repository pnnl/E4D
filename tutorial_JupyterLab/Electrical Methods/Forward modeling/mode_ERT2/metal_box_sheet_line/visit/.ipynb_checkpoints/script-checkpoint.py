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
SetPlotOptions(ScatterAtts)
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



# add the box control points
OpenDatabase("box_control_points.txt", 0)
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
ScatterAtts.pointSizePixels = 10
ScatterAtts.pointType = ScatterAtts.Sphere  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.scaleCube = 0
ScatterAtts.colorType = ScatterAtts.ColorByForegroundColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (0, 0, 0, 255)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 1
SetDefaultPlotOptions(ScatterAtts)
AddPlot("Scatter", "var01", 1, 0)
DrawPlots()
SetActivePlots(4)
SetActivePlots(4)

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
ScatterAtts.pointSizePixels = 10
ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.scaleCube = 0
ScatterAtts.colorType = ScatterAtts.ColorBySingleColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (0, 0, 0, 0)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 0
SetPlotOptions(ScatterAtts)

# sheet control points
OpenDatabase("sheet_control_points.txt", 0)
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
ScatterAtts.pointSizePixels = 10
ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.scaleCube = 0
ScatterAtts.colorType = ScatterAtts.ColorByForegroundColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (0, 0, 255, 255)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 0
SetDefaultPlotOptions(ScatterAtts)
AddPlot("Scatter", "var01", 1, 0)
DrawPlots()
SetActivePlots(5)
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
ScatterAtts.pointSizePixels = 10
ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.scaleCube = 0
ScatterAtts.colorType = ScatterAtts.ColorBySingleColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (0, 0, 0, 0)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 0
SetPlotOptions(ScatterAtts)

# line control points
OpenDatabase("line_control_points.txt", 0)
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
ScatterAtts.pointSizePixels = 10
ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.scaleCube = 0
ScatterAtts.colorType = ScatterAtts.ColorByForegroundColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (0, 0, 255, 255)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 0
SetDefaultPlotOptions(ScatterAtts)
AddPlot("Scatter", "var01", 1, 0)
DrawPlots()
SetActivePlots(5)
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
ScatterAtts.pointSizePixels = 10
ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.scaleCube = 0
ScatterAtts.colorType = ScatterAtts.ColorBySingleColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (0, 0, 0, 0)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 0
SetPlotOptions(ScatterAtts)


# Name of file to open
fn = "pot2071.xmf"
OpenDatabase(fn)
AddPlot("Pseudocolor", "Real_potential")

PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 1
PseudocolorAtts.min = -50
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.max = 50
PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "hot_desaturated"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
PseudocolorAtts.opacityVariable = ""
PseudocolorAtts.opacity = 1
PseudocolorAtts.opacityVarMin = 0
PseudocolorAtts.opacityVarMax = 1
PseudocolorAtts.opacityVarMinFlag = 0
PseudocolorAtts.opacityVarMaxFlag = 0
PseudocolorAtts.pointSize = 0.05
PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
PseudocolorAtts.pointSizeVarEnabled = 0
PseudocolorAtts.pointSizeVar = "default"
PseudocolorAtts.pointSizePixels = 2
PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
PseudocolorAtts.lineWidth = 0
PseudocolorAtts.tubeResolution = 10
PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.tubeRadiusAbsolute = 0.125
PseudocolorAtts.tubeRadiusBBox = 0.005
PseudocolorAtts.tubeRadiusVarEnabled = 0
PseudocolorAtts.tubeRadiusVar = ""
PseudocolorAtts.tubeRadiusVarRatio = 10
PseudocolorAtts.tailStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.headStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.endPointRadiusAbsolute = 0.125
PseudocolorAtts.endPointRadiusBBox = 0.05
PseudocolorAtts.endPointResolution = 10
PseudocolorAtts.endPointRatio = 5
PseudocolorAtts.endPointRadiusVarEnabled = 0
PseudocolorAtts.endPointRadiusVar = ""
PseudocolorAtts.endPointRadiusVarRatio = 10
PseudocolorAtts.renderSurfaces = 1
PseudocolorAtts.renderWireframe = 0
PseudocolorAtts.renderPoints = 0
PseudocolorAtts.smoothingLevel = 0
PseudocolorAtts.legendFlag = 1
PseudocolorAtts.lightingFlag = 1
PseudocolorAtts.wireframeColor = (0, 0, 0, 0)
PseudocolorAtts.pointColor = (0, 0, 0, 0)
SetPlotOptions(PseudocolorAtts)


iso_atts = IsosurfaceAttributes()
iso_atts.contourMethod = iso_atts.Value
AddOperator("Isosurface")
IsosurfaceAtts = IsosurfaceAttributes()
IsosurfaceAtts.contourNLevels = 10
IsosurfaceAtts.contourValue = ()
IsosurfaceAtts.contourPercent = ()
IsosurfaceAtts.contourMethod = IsosurfaceAtts.Level  # Level, Value, Percent
IsosurfaceAtts.minFlag = 0
IsosurfaceAtts.min = 0
IsosurfaceAtts.maxFlag = 0
IsosurfaceAtts.max = 1
IsosurfaceAtts.scaling = IsosurfaceAtts.Linear  # Linear, Log
IsosurfaceAtts.variable = "default"
SetOperatorOptions(IsosurfaceAtts, 1)

AddOperator("Box", 0)
DrawPlots()
BoxAtts = BoxAttributes()
BoxAtts.amount = BoxAtts.All  # Some, All
BoxAtts.minx = 0
BoxAtts.maxx = 1
BoxAtts.miny = 0
BoxAtts.maxy = 1
BoxAtts.minz = 0
BoxAtts.maxz = 1
BoxAtts.inverse = 0
SetOperatorOptions(BoxAtts, 1)
BoxAtts = BoxAttributes()
BoxAtts.amount = BoxAtts.All  # Some, All
BoxAtts.minx = -8
BoxAtts.maxx = 8
BoxAtts.miny = -6
BoxAtts.maxy = 6
BoxAtts.minz = -10
BoxAtts.maxz = 0
BoxAtts.inverse = 0
SetOperatorOptions(BoxAtts, 1)


# change view
DrawPlots()
v = GetView3D()
print ("The view is: ", v)
v.viewNormal = (0, 0, 1)
v.viewUp = (0, 1, 0)
SetView3D(v)

# change axis labels
aatts = AnnotationAttributes()
aatts.axes3D.xAxis.title.userTitle=1
aatts.axes3D.xAxis.title.title='X (m)'

aatts.axes3D.yAxis.title.userTitle=1
aatts.axes3D.yAxis.title.title='Y (m)'

aatts.axes3D.zAxis.title.userTitle=1
aatts.axes3D.zAxis.title.title='Z (m)'

SetAnnotationAttributes(aatts)




