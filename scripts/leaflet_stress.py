# state file generated using paraview version 5.11.1
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import os 

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1340, 764]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [-0.0004854832246049545, 0.00028971316678993997, 0.740667327319006]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-0.016843004181747997, -0.15296930865555933, 9.579856750768624]
renderView1.CameraFocalPoint = [-0.0004854832246049545, 0.00028971316678993997, 0.740667327319006]
renderView1.CameraViewUp = [0.007117357285673379, 0.9998241673588113, 0.017348705730629056]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.35 # 1.573239243304682
renderView1.CameraParallelProjection = 1
renderView1.UseColorPaletteForBackground = 0
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [1340, 762]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [-0.0004854832246049545, 0.00028971316678993997, 0.740667327319006]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [-0.15201620894692056, -5.369324944913632, 1.6342453866313187]
renderView2.CameraFocalPoint = [-0.00048548322460495317, 0.0002897131667899662, 0.7406673273190069]
renderView2.CameraViewUp = [-0.002916867980950387, 0.16423569881121142, 0.9864168120612985]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 2.0635267975416554
renderView2.UseColorPaletteForBackground = 0
renderView2.Background = [1.0, 1.0, 1.0]
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitVertical(0, 0.500000)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, renderView2)
# layout1.SetSize(1340, 1527)
x_width = 750

include_colorbar = True
if include_colorbar:
    x_width = 1340

layout1.SetSize(x_width, 1527)


# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML Unstructured Grid Reader'
aortic_no_partition_3840662_facesvtu = XMLUnstructuredGridReader(registrationName='aortic_no_partition_3840662_faces_mechanics.vtu', FileName=['aortic_no_partition_3840662_faces_mechanics.vtu'])
aortic_no_partition_3840662_facesvtu.PointArrayStatus = ['sigma_circ', 'sigma_rad', 'stress_circ', 'stress_rad', 'sigma_circ_ventricular', 'sigma_rad_ventricular', 'stress_circ_ventricular', 'stress_rad_ventricular']
aortic_no_partition_3840662_facesvtu.TimeArray = 'None'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from aortic_no_partition_3840662_facesvtu
aortic_no_partition_3840662_facesvtuDisplay = Show(aortic_no_partition_3840662_facesvtu, renderView1, 'UnstructuredGridRepresentation')

# get 2D transfer function for 'stress_circ_ventricular'
stress_circ_ventricularTF2D = GetTransferFunction2D('stress_circ_ventricular')
stress_circ_ventricularTF2D.ScalarRangeInitialized = 1
stress_circ_ventricularTF2D.Range = [0.0, 5000000.0, 0.0, 1.0]

# get color transfer function/color map for 'stress_circ_ventricular'
stress_circ_ventricularLUT = GetColorTransferFunction('stress_circ_ventricular')
stress_circ_ventricularLUT.AutomaticRescaleRangeMode = 'Never'
stress_circ_ventricularLUT.TransferFunction2D = stress_circ_ventricularTF2D
stress_circ_ventricularLUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 238094.99999999983, 1.0, 1.0, 0.8125, 555554.9999999998, 1.0, 1.0, 0.5625, 873015.0000000003, 1.0, 1.0, 0.3125, 1190474.9999999998, 1.0, 0.947392, 0.15625, 1507937.5, 1.0, 0.858805, 0.03125, 1825397.5, 1.0, 0.708333, 0.0, 2142857.5, 1.0, 0.541667, 0.0, 2460317.5, 1.0, 0.375, 0.0, 2777777.5, 1.0, 0.208333, 0.0, 3095237.5, 0.937488, 0.039062, 0.0, 3412697.5, 0.854137, 0.0, 0.0, 3730157.5, 0.708333, 0.0, 0.0, 4047620.0, 0.541667, 0.0, 0.0, 4365080.0, 0.375, 0.0, 0.0, 4682540.0, 0.208333, 0.0, 0.0, 5000000.0, 0.0416667, 0.0, 0.0]
stress_circ_ventricularLUT.ColorSpace = 'Lab'
stress_circ_ventricularLUT.NanColor = [0.25, 0.0, 0.0]
stress_circ_ventricularLUT.ScalarRangeInitialized = 1.0
stress_circ_ventricularLUT.Discretize = 0

# get opacity transfer function/opacity map for 'stress_circ_ventricular'
stress_circ_ventricularPWF = GetOpacityTransferFunction('stress_circ_ventricular')
stress_circ_ventricularPWF.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 5000000.0, 1.0, 0.5, 0.0]
stress_circ_ventricularPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
aortic_no_partition_3840662_facesvtuDisplay.Representation = 'Surface'
aortic_no_partition_3840662_facesvtuDisplay.ColorArrayName = ['POINTS', 'stress_circ_ventricular']
aortic_no_partition_3840662_facesvtuDisplay.LookupTable = stress_circ_ventricularLUT
aortic_no_partition_3840662_facesvtuDisplay.SelectTCoordArray = 'None'
aortic_no_partition_3840662_facesvtuDisplay.SelectNormalArray = 'None'
aortic_no_partition_3840662_facesvtuDisplay.SelectTangentArray = 'None'
aortic_no_partition_3840662_facesvtuDisplay.OSPRayScaleArray = 'sigma_circ'
aortic_no_partition_3840662_facesvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
aortic_no_partition_3840662_facesvtuDisplay.SelectOrientationVectors = 'None'
aortic_no_partition_3840662_facesvtuDisplay.ScaleFactor = 0.25115852355956997
aortic_no_partition_3840662_facesvtuDisplay.SelectScaleArray = 'sigma_circ'
aortic_no_partition_3840662_facesvtuDisplay.GlyphType = 'Arrow'
aortic_no_partition_3840662_facesvtuDisplay.GlyphTableIndexArray = 'sigma_circ'
aortic_no_partition_3840662_facesvtuDisplay.GaussianRadius = 0.012557926177978499
aortic_no_partition_3840662_facesvtuDisplay.SetScaleArray = ['POINTS', 'sigma_circ']
aortic_no_partition_3840662_facesvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
aortic_no_partition_3840662_facesvtuDisplay.OpacityArray = ['POINTS', 'sigma_circ']
aortic_no_partition_3840662_facesvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
aortic_no_partition_3840662_facesvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
aortic_no_partition_3840662_facesvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
aortic_no_partition_3840662_facesvtuDisplay.ScalarOpacityFunction = stress_circ_ventricularPWF
aortic_no_partition_3840662_facesvtuDisplay.ScalarOpacityUnitDistance = 0.11323792152563056
aortic_no_partition_3840662_facesvtuDisplay.OpacityArrayName = ['POINTS', 'sigma_circ']
aortic_no_partition_3840662_facesvtuDisplay.SelectInputVectors = [None, '']
aortic_no_partition_3840662_facesvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
aortic_no_partition_3840662_facesvtuDisplay.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
aortic_no_partition_3840662_facesvtuDisplay.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 9776207296.914888, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
aortic_no_partition_3840662_facesvtuDisplay.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 9776207296.914888, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for stress_circ_ventricularLUT in view renderView1
stress_circ_ventricularLUTColorBar = GetScalarBar(stress_circ_ventricularLUT, renderView1)
stress_circ_ventricularLUTColorBar.Position = [0.7622201492537314, 0.39659685863874344]
stress_circ_ventricularLUTColorBar.Title = 'circumferential stress\n(dynes/cm^2)'
stress_circ_ventricularLUTColorBar.ComponentTitle = ''
stress_circ_ventricularLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
stress_circ_ventricularLUTColorBar.TitleFontFamily = 'Times'
stress_circ_ventricularLUTColorBar.TitleBold = 0
stress_circ_ventricularLUTColorBar.TitleFontSize = 26
stress_circ_ventricularLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
stress_circ_ventricularLUTColorBar.LabelFontFamily = 'Times'
stress_circ_ventricularLUTColorBar.LabelBold = 0
stress_circ_ventricularLUTColorBar.LabelFontSize = 26
stress_circ_ventricularLUTColorBar.ScalarBarLength = 0.4
stress_circ_ventricularLUTColorBar.ScalarBarThickness = 24
stress_circ_ventricularLUTColorBar.AutomaticLabelFormat = 0
stress_circ_ventricularLUTColorBar.LabelFormat = '%-#6.1e'
stress_circ_ventricularLUTColorBar.UseCustomLabels = 1
stress_circ_ventricularLUTColorBar.CustomLabels = [5000000.0, 4000000.0, 3000000.0, 2000000.0, 1000000.0, 0.0]
stress_circ_ventricularLUTColorBar.AddRangeLabels = 0
stress_circ_ventricularLUTColorBar.RangeLabelFormat = '%-#6.2g'

# set color bar visibility
stress_circ_ventricularLUTColorBar.Visibility = 1

# show color legend
#aortic_no_partition_3840662_facesvtuDisplay.SetScalarBarVisibility(renderView1, False)
aortic_no_partition_3840662_facesvtuDisplay.SetScalarBarVisibility(renderView1, include_colorbar)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from aortic_no_partition_3840662_facesvtu
aortic_no_partition_3840662_facesvtuDisplay_1 = Show(aortic_no_partition_3840662_facesvtu, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
aortic_no_partition_3840662_facesvtuDisplay_1.Representation = 'Surface'
aortic_no_partition_3840662_facesvtuDisplay_1.ColorArrayName = ['POINTS', 'stress_circ_ventricular']
aortic_no_partition_3840662_facesvtuDisplay_1.LookupTable = stress_circ_ventricularLUT
aortic_no_partition_3840662_facesvtuDisplay_1.SelectTCoordArray = 'None'
aortic_no_partition_3840662_facesvtuDisplay_1.SelectNormalArray = 'None'
aortic_no_partition_3840662_facesvtuDisplay_1.SelectTangentArray = 'None'
aortic_no_partition_3840662_facesvtuDisplay_1.OSPRayScaleArray = 'sigma_circ'
aortic_no_partition_3840662_facesvtuDisplay_1.OSPRayScaleFunction = 'PiecewiseFunction'
aortic_no_partition_3840662_facesvtuDisplay_1.SelectOrientationVectors = 'None'
aortic_no_partition_3840662_facesvtuDisplay_1.ScaleFactor = 0.25882584727689
aortic_no_partition_3840662_facesvtuDisplay_1.SelectScaleArray = 'sigma_circ'
aortic_no_partition_3840662_facesvtuDisplay_1.GlyphType = 'Arrow'
aortic_no_partition_3840662_facesvtuDisplay_1.GlyphTableIndexArray = 'sigma_circ'
aortic_no_partition_3840662_facesvtuDisplay_1.GaussianRadius = 0.012941292363844501
aortic_no_partition_3840662_facesvtuDisplay_1.SetScaleArray = ['POINTS', 'sigma_circ']
aortic_no_partition_3840662_facesvtuDisplay_1.ScaleTransferFunction = 'PiecewiseFunction'
aortic_no_partition_3840662_facesvtuDisplay_1.OpacityArray = ['POINTS', 'sigma_circ']
aortic_no_partition_3840662_facesvtuDisplay_1.OpacityTransferFunction = 'PiecewiseFunction'
aortic_no_partition_3840662_facesvtuDisplay_1.DataAxesGrid = 'GridAxesRepresentation'
aortic_no_partition_3840662_facesvtuDisplay_1.PolarAxes = 'PolarAxesRepresentation'
aortic_no_partition_3840662_facesvtuDisplay_1.ScalarOpacityFunction = stress_circ_ventricularPWF
aortic_no_partition_3840662_facesvtuDisplay_1.ScalarOpacityUnitDistance = 0.11266604367451537
aortic_no_partition_3840662_facesvtuDisplay_1.OpacityArrayName = ['POINTS', 'sigma_circ']
aortic_no_partition_3840662_facesvtuDisplay_1.SelectInputVectors = [None, '']
aortic_no_partition_3840662_facesvtuDisplay_1.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
aortic_no_partition_3840662_facesvtuDisplay_1.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
aortic_no_partition_3840662_facesvtuDisplay_1.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 9776207296.917452, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
aortic_no_partition_3840662_facesvtuDisplay_1.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 9776207296.917452, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# restore active source
SetActiveSource(aortic_no_partition_3840662_facesvtu)
# ----------------------------------------------------------------


if __name__ == '__main__':

    if len(sys.argv) >= 2:
        basename = sys.argv[1]

    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')

    if len(sys.argv) >= 2:
        basename = sys.argv[1]
    else: 
        (head, tail) = os.path.split(os.getcwd())
        basename = tail + 'stress_662'

    # res = (5120, 4320)
    if include_colorbar:
        scaling = 4
        basename = 'colorbar_stress'
    else:
        scaling = 2

    res = (x_width * scaling, 1527 * scaling)

    animationScene1 = GetAnimationScene()
    # animationScene1.GoToFirst()

    Render()
    SaveScreenshot(basename + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)




