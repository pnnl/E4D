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
OpenDatabase("sources.txt", 0)
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
ScatterAtts.legendFlag = 1
SetPlotOptions(ScatterAtts)




OpenDatabase("receivers.txt", 0)
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
ScatterAtts.singleColor = (255, 255, 0, 255)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 1
SetDefaultPlotOptions(ScatterAtts)
AddPlot("Scatter", "var01", 1, 1)
DrawPlots()
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
ScatterAtts.legendFlag = 1
SetPlotOptions(ScatterAtts)


OpenDatabase("bcd.xmf", 0)
AddPlot("Pseudocolor", "Real_conductivity", 1, 1)
DrawPlots()
SetActivePlots(2)
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 0
PseudocolorAtts.min = 0
PseudocolorAtts.maxFlag = 0
PseudocolorAtts.max = 1
PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "hot"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.opacityType = PseudocolorAtts.Constant  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
PseudocolorAtts.opacityVariable = ""
PseudocolorAtts.opacity = 0.258824
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


AddPlot("Mesh", "mesh", 1, 1)
DrawPlots()
SetActivePlots(3)
MeshAtts = MeshAttributes()
MeshAtts.legendFlag = 1
MeshAtts.lineStyle = MeshAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
MeshAtts.lineWidth = 0
MeshAtts.meshColor = (51, 51, 51, 255)
MeshAtts.meshColorSource = MeshAtts.MeshCustom  # Foreground, MeshCustom
MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom
MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off
MeshAtts.pointSize = 0.05
MeshAtts.opaqueColor = (255, 255, 255, 255)
MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
MeshAtts.pointSizeVarEnabled = 0
MeshAtts.pointSizeVar = "default"
MeshAtts.pointType = MeshAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
MeshAtts.showInternal = 0
MeshAtts.pointSizePixels = 2
MeshAtts.opacity = 0.247059
SetPlotOptions(MeshAtts)

AnnotationAtts = AnnotationAttributes()
AnnotationAtts.axes2D.visible = 1
AnnotationAtts.axes2D.autoSetTicks = 1
AnnotationAtts.axes2D.autoSetScaling = 1
AnnotationAtts.axes2D.lineWidth = 0
AnnotationAtts.axes2D.tickLocation = AnnotationAtts.axes2D.Outside  # Inside, Outside, Both
AnnotationAtts.axes2D.tickAxes = AnnotationAtts.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
AnnotationAtts.axes2D.xAxis.title.visible = 1
AnnotationAtts.axes2D.xAxis.title.font.font = AnnotationAtts.axes2D.xAxis.title.font.Courier  # Arial, Courier, Times
AnnotationAtts.axes2D.xAxis.title.font.scale = 1
AnnotationAtts.axes2D.xAxis.title.font.useForegroundColor = 1
AnnotationAtts.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
AnnotationAtts.axes2D.xAxis.title.font.bold = 1
AnnotationAtts.axes2D.xAxis.title.font.italic = 1
AnnotationAtts.axes2D.xAxis.title.userTitle = 0
AnnotationAtts.axes2D.xAxis.title.userUnits = 0
AnnotationAtts.axes2D.xAxis.title.title = "X-Axis"
AnnotationAtts.axes2D.xAxis.title.units = ""
AnnotationAtts.axes2D.xAxis.label.visible = 1
AnnotationAtts.axes2D.xAxis.label.font.font = AnnotationAtts.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
AnnotationAtts.axes2D.xAxis.label.font.scale = 1
AnnotationAtts.axes2D.xAxis.label.font.useForegroundColor = 1
AnnotationAtts.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
AnnotationAtts.axes2D.xAxis.label.font.bold = 1
AnnotationAtts.axes2D.xAxis.label.font.italic = 1
AnnotationAtts.axes2D.xAxis.label.scaling = 0
AnnotationAtts.axes2D.xAxis.tickMarks.visible = 1
AnnotationAtts.axes2D.xAxis.tickMarks.majorMinimum = 0
AnnotationAtts.axes2D.xAxis.tickMarks.majorMaximum = 1
AnnotationAtts.axes2D.xAxis.tickMarks.minorSpacing = 0.02
AnnotationAtts.axes2D.xAxis.tickMarks.majorSpacing = 0.2
AnnotationAtts.axes2D.xAxis.grid = 0
AnnotationAtts.axes2D.yAxis.title.visible = 1
AnnotationAtts.axes2D.yAxis.title.font.font = AnnotationAtts.axes2D.yAxis.title.font.Courier  # Arial, Courier, Times
AnnotationAtts.axes2D.yAxis.title.font.scale = 1
AnnotationAtts.axes2D.yAxis.title.font.useForegroundColor = 1
AnnotationAtts.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
AnnotationAtts.axes2D.yAxis.title.font.bold = 1
AnnotationAtts.axes2D.yAxis.title.font.italic = 1
AnnotationAtts.axes2D.yAxis.title.userTitle = 0
AnnotationAtts.axes2D.yAxis.title.userUnits = 0
AnnotationAtts.axes2D.yAxis.title.title = "Y-Axis"
AnnotationAtts.axes2D.yAxis.title.units = ""
AnnotationAtts.axes2D.yAxis.label.visible = 1
AnnotationAtts.axes2D.yAxis.label.font.font = AnnotationAtts.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
AnnotationAtts.axes2D.yAxis.label.font.scale = 1
AnnotationAtts.axes2D.yAxis.label.font.useForegroundColor = 1
AnnotationAtts.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
AnnotationAtts.axes2D.yAxis.label.font.bold = 1
AnnotationAtts.axes2D.yAxis.label.font.italic = 1
AnnotationAtts.axes2D.yAxis.label.scaling = 0
AnnotationAtts.axes2D.yAxis.tickMarks.visible = 1
AnnotationAtts.axes2D.yAxis.tickMarks.majorMinimum = 0
AnnotationAtts.axes2D.yAxis.tickMarks.majorMaximum = 1
AnnotationAtts.axes2D.yAxis.tickMarks.minorSpacing = 0.02
AnnotationAtts.axes2D.yAxis.tickMarks.majorSpacing = 0.2
AnnotationAtts.axes2D.yAxis.grid = 0
AnnotationAtts.axes3D.visible = 1
AnnotationAtts.axes3D.autoSetTicks = 1
AnnotationAtts.axes3D.autoSetScaling = 0
AnnotationAtts.axes3D.lineWidth = 0
AnnotationAtts.axes3D.tickLocation = AnnotationAtts.axes3D.Inside  # Inside, Outside, Both
AnnotationAtts.axes3D.axesType = AnnotationAtts.axes3D.ClosestTriad  # ClosestTriad, FurthestTriad, OutsideEdges, StaticTriad, StaticEdges
AnnotationAtts.axes3D.triadFlag = 1
AnnotationAtts.axes3D.bboxFlag = 0
AnnotationAtts.axes3D.xAxis.title.visible = 1
AnnotationAtts.axes3D.xAxis.title.font.font = AnnotationAtts.axes3D.xAxis.title.font.Arial  # Arial, Courier, Times
AnnotationAtts.axes3D.xAxis.title.font.scale = 1
AnnotationAtts.axes3D.xAxis.title.font.useForegroundColor = 1
AnnotationAtts.axes3D.xAxis.title.font.color = (0, 0, 0, 255)
AnnotationAtts.axes3D.xAxis.title.font.bold = 0
AnnotationAtts.axes3D.xAxis.title.font.italic = 0
AnnotationAtts.axes3D.xAxis.title.userTitle = 1
AnnotationAtts.axes3D.xAxis.title.userUnits = 0
AnnotationAtts.axes3D.xAxis.title.title = "X (m)"
AnnotationAtts.axes3D.xAxis.title.units = ""
AnnotationAtts.axes3D.xAxis.label.visible = 1
AnnotationAtts.axes3D.xAxis.label.font.font = AnnotationAtts.axes3D.xAxis.label.font.Arial  # Arial, Courier, Times
AnnotationAtts.axes3D.xAxis.label.font.scale = 2
AnnotationAtts.axes3D.xAxis.label.font.useForegroundColor = 1
AnnotationAtts.axes3D.xAxis.label.font.color = (0, 0, 0, 255)
AnnotationAtts.axes3D.xAxis.label.font.bold = 0
AnnotationAtts.axes3D.xAxis.label.font.italic = 0
AnnotationAtts.axes3D.xAxis.label.scaling = 0
AnnotationAtts.axes3D.xAxis.tickMarks.visible = 0
AnnotationAtts.axes3D.xAxis.tickMarks.majorMinimum = 0
AnnotationAtts.axes3D.xAxis.tickMarks.majorMaximum = 1
AnnotationAtts.axes3D.xAxis.tickMarks.minorSpacing = 0.02
AnnotationAtts.axes3D.xAxis.tickMarks.majorSpacing = 0.2
AnnotationAtts.axes3D.xAxis.grid = 0
AnnotationAtts.axes3D.yAxis.title.visible = 1
AnnotationAtts.axes3D.yAxis.title.font.font = AnnotationAtts.axes3D.yAxis.title.font.Arial  # Arial, Courier, Times
AnnotationAtts.axes3D.yAxis.title.font.scale = 1
AnnotationAtts.axes3D.yAxis.title.font.useForegroundColor = 1
AnnotationAtts.axes3D.yAxis.title.font.color = (0, 0, 0, 255)
AnnotationAtts.axes3D.yAxis.title.font.bold = 0
AnnotationAtts.axes3D.yAxis.title.font.italic = 0
AnnotationAtts.axes3D.yAxis.title.userTitle = 1
AnnotationAtts.axes3D.yAxis.title.userUnits = 0
AnnotationAtts.axes3D.yAxis.title.title = "Y (m)"
AnnotationAtts.axes3D.yAxis.title.units = ""
AnnotationAtts.axes3D.yAxis.label.visible = 1
AnnotationAtts.axes3D.yAxis.label.font.font = AnnotationAtts.axes3D.yAxis.label.font.Arial  # Arial, Courier, Times
AnnotationAtts.axes3D.yAxis.label.font.scale = 2
AnnotationAtts.axes3D.yAxis.label.font.useForegroundColor = 1
AnnotationAtts.axes3D.yAxis.label.font.color = (0, 0, 0, 255)
AnnotationAtts.axes3D.yAxis.label.font.bold = 0
AnnotationAtts.axes3D.yAxis.label.font.italic = 0
AnnotationAtts.axes3D.yAxis.label.scaling = 0
AnnotationAtts.axes3D.yAxis.tickMarks.visible = 0
AnnotationAtts.axes3D.yAxis.tickMarks.majorMinimum = 0
AnnotationAtts.axes3D.yAxis.tickMarks.majorMaximum = 1
AnnotationAtts.axes3D.yAxis.tickMarks.minorSpacing = 0.02
AnnotationAtts.axes3D.yAxis.tickMarks.majorSpacing = 0.2
AnnotationAtts.axes3D.yAxis.grid = 0
AnnotationAtts.axes3D.zAxis.title.visible = 1
AnnotationAtts.axes3D.zAxis.title.font.font = AnnotationAtts.axes3D.zAxis.title.font.Arial  # Arial, Courier, Times
AnnotationAtts.axes3D.zAxis.title.font.scale = 1
AnnotationAtts.axes3D.zAxis.title.font.useForegroundColor = 1
AnnotationAtts.axes3D.zAxis.title.font.color = (0, 0, 0, 255)
AnnotationAtts.axes3D.zAxis.title.font.bold = 0
AnnotationAtts.axes3D.zAxis.title.font.italic = 0
AnnotationAtts.axes3D.zAxis.title.userTitle = 1
AnnotationAtts.axes3D.zAxis.title.userUnits = 0
AnnotationAtts.axes3D.zAxis.title.title = "Z (m)"
AnnotationAtts.axes3D.zAxis.title.units = ""
AnnotationAtts.axes3D.zAxis.label.visible = 1
AnnotationAtts.axes3D.zAxis.label.font.font = AnnotationAtts.axes3D.zAxis.label.font.Arial  # Arial, Courier, Times
AnnotationAtts.axes3D.zAxis.label.font.scale =2
AnnotationAtts.axes3D.zAxis.label.font.useForegroundColor = 1
AnnotationAtts.axes3D.zAxis.label.font.color = (0, 0, 0, 255)
AnnotationAtts.axes3D.zAxis.label.font.bold = 0
AnnotationAtts.axes3D.zAxis.label.font.italic = 0
AnnotationAtts.axes3D.zAxis.label.scaling = 0
AnnotationAtts.axes3D.zAxis.tickMarks.visible = 1
AnnotationAtts.axes3D.zAxis.tickMarks.majorMinimum = 0
AnnotationAtts.axes3D.zAxis.tickMarks.majorMaximum = 1
AnnotationAtts.axes3D.zAxis.tickMarks.minorSpacing = 0.02
AnnotationAtts.axes3D.zAxis.tickMarks.majorSpacing = 0.2
AnnotationAtts.axes3D.zAxis.grid = 0
AnnotationAtts.axes3D.setBBoxLocation = 0
AnnotationAtts.axes3D.bboxLocation = (0, 1, 0, 1, 0, 1)
AnnotationAtts.userInfoFlag = 1
AnnotationAtts.userInfoFont.font = AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
AnnotationAtts.userInfoFont.scale = 1
AnnotationAtts.userInfoFont.useForegroundColor = 1
AnnotationAtts.userInfoFont.color = (0, 0, 0, 255)
AnnotationAtts.userInfoFont.bold = 0
AnnotationAtts.userInfoFont.italic = 0
AnnotationAtts.databaseInfoFlag = 0
AnnotationAtts.timeInfoFlag = 1
AnnotationAtts.databaseInfoFont.font = AnnotationAtts.databaseInfoFont.Arial  # Arial, Courier, Times
AnnotationAtts.databaseInfoFont.scale = 1
AnnotationAtts.databaseInfoFont.useForegroundColor = 1
AnnotationAtts.databaseInfoFont.color = (0, 0, 0, 255)
AnnotationAtts.databaseInfoFont.bold = 0
AnnotationAtts.databaseInfoFont.italic = 0
AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
AnnotationAtts.databaseInfoTimeScale = 1
AnnotationAtts.databaseInfoTimeOffset = 0
AnnotationAtts.legendInfoFlag = 1
AnnotationAtts.backgroundColor = (255, 255, 255, 255)
AnnotationAtts.foregroundColor = (0, 0, 0, 255)
AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
AnnotationAtts.gradientColor1 = (0, 0, 255, 255)
AnnotationAtts.gradientColor2 = (0, 0, 0, 255)
AnnotationAtts.backgroundMode = AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
AnnotationAtts.backgroundImage = ""
AnnotationAtts.imageRepeatX = 1
AnnotationAtts.imageRepeatY = 1
AnnotationAtts.axesArray.visible = 1
AnnotationAtts.axesArray.ticksVisible = 1
AnnotationAtts.axesArray.autoSetTicks = 1
AnnotationAtts.axesArray.autoSetScaling = 1
AnnotationAtts.axesArray.lineWidth = 0
AnnotationAtts.axesArray.axes.title.visible = 1
AnnotationAtts.axesArray.axes.title.font.font = AnnotationAtts.axesArray.axes.title.font.Arial  # Arial, Courier, Times
AnnotationAtts.axesArray.axes.title.font.scale = 1
AnnotationAtts.axesArray.axes.title.font.useForegroundColor = 1
AnnotationAtts.axesArray.axes.title.font.color = (0, 0, 0, 255)
AnnotationAtts.axesArray.axes.title.font.bold = 0
AnnotationAtts.axesArray.axes.title.font.italic = 0
AnnotationAtts.axesArray.axes.title.userTitle = 0
AnnotationAtts.axesArray.axes.title.userUnits = 0
AnnotationAtts.axesArray.axes.title.title = ""
AnnotationAtts.axesArray.axes.title.units = ""
AnnotationAtts.axesArray.axes.label.visible = 1
AnnotationAtts.axesArray.axes.label.font.font = AnnotationAtts.axesArray.axes.label.font.Arial  # Arial, Courier, Times
AnnotationAtts.axesArray.axes.label.font.scale = 1
AnnotationAtts.axesArray.axes.label.font.useForegroundColor = 1
AnnotationAtts.axesArray.axes.label.font.color = (0, 0, 0, 255)
AnnotationAtts.axesArray.axes.label.font.bold = 0
AnnotationAtts.axesArray.axes.label.font.italic = 0
AnnotationAtts.axesArray.axes.label.scaling = 0
AnnotationAtts.axesArray.axes.tickMarks.visible = 1
AnnotationAtts.axesArray.axes.tickMarks.majorMinimum = 0
AnnotationAtts.axesArray.axes.tickMarks.majorMaximum = 1
AnnotationAtts.axesArray.axes.tickMarks.minorSpacing = 0.02
AnnotationAtts.axesArray.axes.tickMarks.majorSpacing = 0.2
AnnotationAtts.axesArray.axes.grid = 0
SetAnnotationAttributes(AnnotationAtts)



# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0.17, -0.87, 0.47)
View3DAtts.focus = (0,0,-35)
View3DAtts.viewUp = (-0.07, 0.46, 0.88)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 40.1203
View3DAtts.nearPlane =-80.2407
View3DAtts.farPlane = 80.2407
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (1.15217, 1.14826, 1806.5)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state