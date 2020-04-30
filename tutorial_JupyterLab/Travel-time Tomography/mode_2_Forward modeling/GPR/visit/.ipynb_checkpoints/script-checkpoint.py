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
# add the source and receiver locations


OpenDatabase("sources.txt", 0)

AddPlot("Scatter", "var01", 1, 1)
DrawPlots()

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
ScatterAtts.colorType = ScatterAtts.ColorBySingleColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (255, 255, 0, 255)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 0
SetPlotOptions(ScatterAtts)

OpenDatabase("receivers.txt", 0)
AddPlot("Scatter", "var01", 1, 1)
DrawPlots()
SetActivePlots(1)
SetActivePlots(1)
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
ScatterAtts.colorType = ScatterAtts.ColorBySingleColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (0, 255, 0, 255)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 0
SetPlotOptions(ScatterAtts)


OpenDatabase("tt100.xmf", 0)
AddPlot("Pseudocolor", "Real_potential", 1, 0)
DrawPlots()
SetActivePlots(2)
SetActivePlots(2)
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 1
PseudocolorAtts.min = 0
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.max = 200
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
PseudocolorAtts.lightingFlag = 0
PseudocolorAtts.wireframeColor = (0, 0, 0, 0)
PseudocolorAtts.pointColor = (0, 0, 0, 0)
SetPlotOptions(PseudocolorAtts)

SetActivePlots((0, 2))
SetActivePlots((0, 2))
SetActivePlots(2)
AddOperator("ThreeSlice", 0)
DrawPlots()
ThreeSliceAtts = ThreeSliceAttributes()
ThreeSliceAtts.x = 0
ThreeSliceAtts.y = 1.25
ThreeSliceAtts.z = -40
ThreeSliceAtts.interactive = 1
SetOperatorOptions(ThreeSliceAtts, 0)

# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (-0.6887037510620115, -0.6442120524288911, 0.3326829944233257)

View3DAtts.focus = (-1.1708, 0, -33.0645)
View3DAtts.viewUp = (0.2474604520699948, 0.2224410054681357, 0.9430181990543125)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 36.2656
View3DAtts.nearPlane = -72.5311
View3DAtts.farPlane = 72.5311
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 0.826446
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (-1.1708, 0, -33.0645)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state





