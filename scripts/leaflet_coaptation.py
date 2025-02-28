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
SetActiveView(renderView2)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML Unstructured Grid Reader'
aortic_no_partition_3840662_faces_mechanicsvtu = XMLUnstructuredGridReader(registrationName='aortic_no_partition_3840662_faces_mechanics.vtu', FileName=['aortic_no_partition_3840662_faces_mechanics.vtu'])
aortic_no_partition_3840662_faces_mechanicsvtu.PointArrayStatus = ['sigma_circ', 'sigma_rad', 'stress_circ', 'stress_rad', 'sigma_circ_ventricular', 'sigma_rad_ventricular', 'stress_circ_ventricular', 'stress_rad_ventricular', 'coapt_distances']
aortic_no_partition_3840662_faces_mechanicsvtu.TimeArray = 'None'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from aortic_no_partition_3840662_faces_mechanicsvtu
aortic_no_partition_3840662_faces_mechanicsvtuDisplay = Show(aortic_no_partition_3840662_faces_mechanicsvtu, renderView1, 'UnstructuredGridRepresentation')

# # get 2D transfer function for 'coapt_distances'
# coapt_distancesTF2D = GetTransferFunction2D('coapt_distances')
# coapt_distancesTF2D.ScalarRangeInitialized = 1
# coapt_distancesTF2D.Range = [0.0, 1.4, 0.0, 1.0]

# # get color transfer function/color map for 'coapt_distances'
# coapt_distancesLUT = GetColorTransferFunction('coapt_distances')
# coapt_distancesLUT.TransferFunction2D = coapt_distancesTF2D
# coapt_distancesLUT.RGBPoints = [0.0, 1.0, 1.0, 0.988235, 0.0028, 1.0, 1.0, 0.988235, 0.06999999999999999, 0.984314, 0.988235, 0.843137, 0.13999999999999999, 0.988235, 0.988235, 0.741176, 0.21, 0.980392, 0.968627, 0.654902, 0.27999999999999997, 0.980392, 0.945098, 0.576471, 0.35, 0.968627, 0.905882, 0.486275, 0.42, 0.968627, 0.862745, 0.388235, 0.48999999999999994, 0.960784, 0.803922, 0.286275, 0.5599999999999999, 0.94902, 0.741176, 0.219608, 0.63, 0.941176, 0.678431, 0.14902, 0.7, 0.929412, 0.607843, 0.094118, 0.77, 0.921569, 0.545098, 0.054902, 0.84, 0.909804, 0.486275, 0.035294, 0.9099999999999999, 0.890196, 0.411765, 0.019608, 0.9799999999999999, 0.8, 0.305882, 0.0, 1.0499999999999998, 0.760784, 0.239216, 0.0, 1.1199999999999999, 0.678431, 0.180392, 0.011765, 1.19, 0.6, 0.121569, 0.023529, 1.26, 0.501961, 0.054902, 0.031373, 1.3299999999999998, 0.4, 0.039216, 0.058824, 1.4, 0.301961, 0.047059, 0.090196]
# coapt_distancesLUT.ColorSpace = 'Lab'
# coapt_distancesLUT.NanColor = [0.25, 0.0, 0.0]
# coapt_distancesLUT.ScalarRangeInitialized = 1.0

# # get opacity transfer function/opacity map for 'coapt_distances'
# coapt_distancesPWF = GetOpacityTransferFunction('coapt_distances')
# coapt_distancesPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.4, 1.0, 0.5, 0.0]
# coapt_distancesPWF.ScalarRangeInitialized = 1

# get 2D transfer function for 'coapt_distances'
coapt_distancesTF2D = GetTransferFunction2D('coapt_distances')
coapt_distancesTF2D.ScalarRangeInitialized = 1

# get color transfer function/color map for 'coapt_distances'
coapt_distancesLUT = GetColorTransferFunction('coapt_distances')
coapt_distancesLUT.AutomaticRescaleRangeMode = 'Never'
coapt_distancesLUT.TransferFunction2D = coapt_distancesTF2D
coapt_distancesLUT.RGBPoints = [0.0, 0.0416667, 0.0, 0.0, 0.06349199999999999, 0.208333, 0.0, 0.0, 0.12698399999999999, 0.375, 0.0, 0.0, 0.19047599999999998, 0.541667, 0.0, 0.0, 0.25396850000000004, 0.708333, 0.0, 0.0, 0.31746050000000003, 0.854137, 0.0, 0.0, 0.3809525, 0.937488, 0.039062, 0.0, 0.4444445, 1.0, 0.208333, 0.0, 0.5079365, 1.0, 0.375, 0.0, 0.5714285, 1.0, 0.541667, 0.0, 0.6349205, 1.0, 0.708333, 0.0, 0.6984125, 1.0, 0.858805, 0.03125, 0.761905, 1.0, 0.947392, 0.15625, 0.8253969999999999, 1.0, 1.0, 0.3125, 0.888889, 1.0, 1.0, 0.5625, 0.9523809999999999, 1.0, 1.0, 0.8125, 1.0, 1.0, 1.0, 1.0]
coapt_distancesLUT.ColorSpace = 'Lab'
coapt_distancesLUT.NanColor = [0.25, 0.0, 0.0]
coapt_distancesLUT.ScalarRangeInitialized = 1.0
coapt_distancesLUT.Discretize = 0

# get opacity transfer function/opacity map for 'coapt_distances'
coapt_distancesPWF = GetOpacityTransferFunction('coapt_distances')
coapt_distancesPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.Representation = 'Surface'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.ColorArrayName = ['POINTS', 'coapt_distances']
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.LookupTable = coapt_distancesLUT
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.SelectTCoordArray = 'None'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.SelectNormalArray = 'None'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.SelectTangentArray = 'None'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.OSPRayScaleArray = 'sigma_circ'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.SelectOrientationVectors = 'None'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.ScaleFactor = 0.2588429466599
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.SelectScaleArray = 'sigma_circ'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.GlyphType = 'Arrow'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.GlyphTableIndexArray = 'sigma_circ'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.GaussianRadius = 0.012942147332995
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.SetScaleArray = ['POINTS', 'sigma_circ']
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.OpacityArray = ['POINTS', 'sigma_circ']
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.ScalarOpacityFunction = coapt_distancesPWF
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.ScalarOpacityUnitDistance = 0.11397292299352607
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.OpacityArrayName = ['POINTS', 'sigma_circ']
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.SelectInputVectors = [None, '']
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 6063351753.840244, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 6063351753.840244, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for coapt_distancesLUT in view renderView1
coapt_distancesLUTColorBar = GetScalarBar(coapt_distancesLUT, renderView1)
coapt_distancesLUTColorBar.Title = 'distance (cm)'
# coapt_distancesLUTColorBar.ComponentTitle = ''
# coapt_distancesLUTColorBar.AutomaticLabelFormat = 0
# coapt_distancesLUTColorBar.LabelFormat = '%-#6.1f'
# coapt_distancesLUTColorBar.UseCustomLabels = 1
# coapt_distancesLUTColorBar.CustomLabels = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
# coapt_distancesLUTColorBar.RangeLabelFormat = '%-#6.1f'

coapt_distancesLUTColorBar.Position = [0.7622201492537314, 0.39659685863874344]
coapt_distancesLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
coapt_distancesLUTColorBar.TitleFontFamily = 'Times'
coapt_distancesLUTColorBar.TitleBold = 0
coapt_distancesLUTColorBar.TitleFontSize = 26
coapt_distancesLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
coapt_distancesLUTColorBar.LabelFontFamily = 'Times'
coapt_distancesLUTColorBar.LabelBold = 0
coapt_distancesLUTColorBar.LabelFontSize = 26
coapt_distancesLUTColorBar.ScalarBarLength = 0.4
coapt_distancesLUTColorBar.ScalarBarThickness = 24
coapt_distancesLUTColorBar.AutomaticLabelFormat = 0
coapt_distancesLUTColorBar.LabelFormat = '%-#6.1f'
coapt_distancesLUTColorBar.UseCustomLabels = 1
coapt_distancesLUTColorBar.CustomLabels = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
coapt_distancesLUTColorBar.RangeLabelFormat = '%-#6.1f'

coapt_distancesLUTColorBar.ComponentTitle = ''
# coapt_distancesLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
# coapt_distancesLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
# coapt_distancesLUTColorBar.AutomaticLabelFormat = 0
# coapt_distancesLUTColorBar.LabelFormat = '%-#6.1f'
# coapt_distancesLUTColorBar.UseCustomLabels = 1
# coapt_distancesLUTColorBar.CustomLabels = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
# coapt_distancesLUTColorBar.RangeLabelFormat = '%-#6.1f'

# set color bar visibility
coapt_distancesLUTColorBar.Visibility = 1

# show color legend
aortic_no_partition_3840662_faces_mechanicsvtuDisplay.SetScalarBarVisibility(renderView1, include_colorbar)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from aortic_no_partition_3840662_faces_mechanicsvtu
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1 = Show(aortic_no_partition_3840662_faces_mechanicsvtu, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.Representation = 'Surface'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.ColorArrayName = ['POINTS', 'coapt_distances']
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.LookupTable = coapt_distancesLUT
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.SelectTCoordArray = 'None'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.SelectNormalArray = 'None'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.SelectTangentArray = 'None'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.OSPRayScaleArray = 'sigma_circ'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.OSPRayScaleFunction = 'PiecewiseFunction'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.SelectOrientationVectors = 'None'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.ScaleFactor = 0.2588429466599
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.SelectScaleArray = 'sigma_circ'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.GlyphType = 'Arrow'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.GlyphTableIndexArray = 'sigma_circ'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.GaussianRadius = 0.012942147332995
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.SetScaleArray = ['POINTS', 'sigma_circ']
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.ScaleTransferFunction = 'PiecewiseFunction'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.OpacityArray = ['POINTS', 'sigma_circ']
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.OpacityTransferFunction = 'PiecewiseFunction'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.DataAxesGrid = 'GridAxesRepresentation'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.PolarAxes = 'PolarAxesRepresentation'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.ScalarOpacityFunction = coapt_distancesPWF
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.ScalarOpacityUnitDistance = 0.11397292299352607
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.OpacityArrayName = ['POINTS', 'sigma_circ']
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.SelectInputVectors = [None, '']
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 6063351753.840244, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 6063351753.840244, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for coapt_distancesLUT in view renderView2
coapt_distancesLUTColorBar_1 = GetScalarBar(coapt_distancesLUT, renderView2)
coapt_distancesLUTColorBar_1.Title = 'coapt_distances'
coapt_distancesLUTColorBar_1.ComponentTitle = ''

# set color bar visibility
coapt_distancesLUTColorBar_1.Visibility = 1

# show color legend
aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.SetScalarBarVisibility(renderView2, False)
# aortic_no_partition_3840662_faces_mechanicsvtuDisplay_1.SetScalarBarVisibility(renderView2, include_colorbar)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# restore active source
SetActiveSource(aortic_no_partition_3840662_faces_mechanicsvtu)
# ----------------------------------------------------------------


if __name__ == '__main__':

    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')

    # if len(sys.argv) >= 2:
    #     basename = sys.argv[1]
    # else: 
    (head, tail) = os.path.split(os.getcwd())
    basename = tail
    basename += 'coaptation_662'

    # res = (5120, 4320)
    # res = (5120, 4320)
    if include_colorbar:
        scaling = 4
        basename = 'colorbar_coaptation'
    else:
        scaling = 2
    res = (x_width * scaling, 1527 * scaling)

    animationScene1 = GetAnimationScene()
    # animationScene1.GoToFirst()

    Render()
    SaveScreenshot(basename + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)
