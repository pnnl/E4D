# Visit 3.1.0 log file
ScriptVersion = "3.1.0"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
ShowAllWindows()
Source("script.py")
# MAINTENANCE ISSUE: SetSuppressMessagesRPC is not handled in Logging.C. Please contact a VisIt developer.
SaveSession("/home/robi526/.visit/crash_recovery.19924.session")
# MAINTENANCE ISSUE: SetSuppressMessagesRPC is not handled in Logging.C. Please contact a VisIt developer.
# MAINTENANCE ISSUE: SetSuppressMessagesRPC is not handled in Logging.C. Please contact a VisIt developer.
SaveSession("/home/robi526/.visit/crash_recovery.19924.session")
# MAINTENANCE ISSUE: SetSuppressMessagesRPC is not handled in Logging.C. Please contact a VisIt developer.
OpenDatabase("localhost:/home/robi526/codes/e4d_dev/e4d_dev/tutorial_JupyterLab/Electrical Methods/Time-lapse inverse modeling/mode_ERT4/sinking_plume/forward/models.xmf", 0)
# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (-0.4, -0.9, 0.4)
View3DAtts.focus = (0, 3.43019, -5)
View3DAtts.viewUp = (0.1, 0.37, 0.92)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 9.66962
View3DAtts.nearPlane = -21.8403
View3DAtts.farPlane = 21.8403
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 0.683013
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (0, 0.5, -5)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

AddPlot("Pseudocolor", "Real_conductivity", 1, 1)
DrawPlots()
SetActivePlots(2)
SetActivePlots(2)
SetActivePlots((1, 2))
SetActivePlots((1, 2))
SetActivePlots(1)
HideActivePlots()
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
SetTimeSliderState(6)
SetTimeSliderState(15)
SetActivePlots((1, 2))
SetActivePlots((1, 2))
SetActivePlots(2)
SetActivePlots((1, 2))
SetActivePlots((1, 2))
SetActivePlots(1)
SetActivePlots((1, 2))
SetActivePlots((1, 2))
SetActivePlots(2)
HideActivePlots()
HideActivePlots()
SetActivePlots((1, 2))
SetActivePlots(1)
HideActivePlots()
SetActivePlots((1, 2))
SetActivePlots(2)
HideActivePlots()
HideActivePlots()
HideActivePlots()
HideActivePlots()
silr = SILRestriction()
silr.TurnOnAll()
silr.TurnOffSet(2)
SetPlotSILRestriction(silr ,1)
SetActivePlots((1, 2))
SetActivePlots(1)
HideActivePlots()
SetActivePlots((1, 2))
SetActivePlots(2)
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Log  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 0
PseudocolorAtts.min = 0
PseudocolorAtts.useBelowMinColor = 0
PseudocolorAtts.belowMinColor = (0, 0, 0, 255)
PseudocolorAtts.maxFlag = 0
PseudocolorAtts.max = 1
PseudocolorAtts.useAboveMaxColor = 0
PseudocolorAtts.aboveMaxColor = (0, 0, 0, 255)
PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "Default"
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
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
SetTimeSliderState(0)
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
# MAINTENANCE ISSUE: SetSuppressMessagesRPC is not handled in Logging.C. Please contact a VisIt developer.
SaveSession("/home/robi526/.visit/crash_recovery.19924.session")
# MAINTENANCE ISSUE: SetSuppressMessagesRPC is not handled in Logging.C. Please contact a VisIt developer.
OpenDatabase("localhost:/home/robi526/codes/e4d_dev/e4d_dev/tutorial_JupyterLab/Electrical Methods/Time-lapse inverse modeling/mode_ERT4/sinking_plume/forward/models.xmf", 0)
SetTimeSliderState(5)
SetTimeSliderState(6)
SetTimeSliderState(7)
SetTimeSliderState(8)
SetTimeSliderState(9)
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
OpenDatabase("localhost:/home/robi526/codes/e4d_dev/e4d_dev/tutorial_JupyterLab/Electrical Methods/Time-lapse inverse modeling/mode_ERT4/sinking_plume/forward/models.xmf", 0)
OpenDatabase("localhost:/home/robi526/codes/e4d_dev/e4d_dev/tutorial_JupyterLab/Electrical Methods/Time-lapse inverse modeling/mode_ERT4/sinking_plume/forward/models.xmf", 0)
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
OpenDatabase("localhost:/home/robi526/codes/e4d_dev/e4d_dev/tutorial_JupyterLab/Electrical Methods/Time-lapse inverse modeling/mode_ERT4/sinking_plume/forward/models.xmf", 0)
OpenDatabase("localhost:/home/robi526/codes/e4d_dev/e4d_dev/tutorial_JupyterLab/Electrical Methods/Time-lapse inverse modeling/mode_ERT4/sinking_plume/forward/models.xmf", 0)
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
SetTimeSliderState(0)
OpenDatabase("localhost:/home/robi526/codes/e4d_dev/e4d_dev/tutorial_JupyterLab/Electrical Methods/Time-lapse inverse modeling/mode_ERT4/sinking_plume/forward/models.xmf", 0)
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
SetTimeSliderState(0)
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
OpenDatabase("localhost:/home/robi526/codes/e4d_dev/e4d_dev/tutorial_JupyterLab/Electrical Methods/Time-lapse inverse modeling/mode_ERT4/sinking_plume/forward/models.xmf", 0)
# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.
OpenDatabase("localhost:/home/robi526/codes/e4d_dev/e4d_dev/tutorial_JupyterLab/Electrical Methods/Time-lapse inverse modeling/mode_ERT4/sinking_plume/forward/models.xmf", 0)
OpenDatabase("localhost:/home/robi526/codes/e4d_dev/e4d_dev/tutorial_JupyterLab/Electrical Methods/Time-lapse inverse modeling/mode_ERT4/sinking_plume/forward/models.xmf", 0)
OpenDatabase("localhost:/home/robi526/codes/e4d_dev/e4d