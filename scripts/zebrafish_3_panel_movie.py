print("in script")

# state file generated using paraview version 5.11.1
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
import os 
import math 

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------


import re


version = paraview.simple.GetParaViewVersion()
print("version = ", version)
version_num_str_list = re.findall(r'\d+', str(version))

if int(version_num_str_list[0]) == 5: 
    version_major = 5
else: 
    raise ValueError("unsupported major version")

if int(version_num_str_list[1]) == 8:
    version_sub_num = 8
elif int(version_num_str_list[1]) == 11:
    version_sub_num = 11
else: 
    raise ValueError("unsupported minor version")

print("working with version ", version_major, ".", version_sub_num)

print("through imports")

if len(sys.argv) >= 5:
    vertical_paper = int(sys.argv[4])
else: 
    vertical_paper = 0


# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1786, 2828]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 'Crystal Eyes'
if vertical_paper:
    # renderView1.CenterOfRotation = [-0.0, -0.0, 0.0334999980404973]
    renderView1.CenterOfRotation = [-0.0, -0.0, 0.0]
    renderView1.CameraFocalPoint = [-0.0, -0.0, 0.033225393645198456]
    # renderView1.CameraPosition = [0.4405789377345504, -0.00892120153778747, 0.033225393645198456]
    renderView1.CameraPosition = [0.4405789377345504, -0.0, 0.033225393645198456]
else: 
    renderView1.CenterOfRotation = [-0.0031767683103680836, -0.0005000000819563866, 0.0334999980404973]
    renderView1.CameraFocalPoint = [-0.0031767683103680836, -0.00892120153778747, 0.033225393645198456]
    renderView1.CameraPosition = [0.4405789377345504, -0.00892120153778747, 0.033225393645198456]

renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.06471510249209315
renderView1.CameraParallelProjection = 1
# renderView1.UseColorPaletteForBackground = 0
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [1356, 1376]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [-4.214234650150653e-06, -3.632158039959532e-08, 0.014272129305027208]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [-4.214234650150653e-06, -3.632158039959532e-08, 0.12778617853365398]
renderView2.CameraFocalPoint = [-4.214234650150653e-06, -3.632158039959532e-08, 0.014272129305027208]
renderView2.CameraViewUp = [-1.0, 2.220446049250313e-16, 0.0]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 0.020429645795568864
renderView2.CameraParallelProjection = 1
# renderView2.UseColorPaletteForBackground = 0
renderView2.Background = [1.0, 1.0, 1.0]
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView3 = CreateView('RenderView')
renderView3.ViewSize = [1356, 1374]
renderView3.AxesGrid = 'GridAxes3DActor'
renderView3.OrientationAxesVisibility = 0
renderView3.CenterOfRotation = [-2.4391338228988563e-06, -3.725290294992467e-08, 0.01427239351687602]
renderView3.StereoType = 'Crystal Eyes'
renderView3.CameraPosition = [-0.08237521506951388, -0.06721252031715848, 0.05407495975593145]
renderView3.CameraFocalPoint = [-2.439133822899114e-06, -3.725290295026318e-08, 0.01427239351687602]
renderView3.CameraViewUp = [0.270401265009845, 0.22320349959021077, 0.9365166061804556]
renderView3.CameraFocalDisk = 1.0
renderView3.CameraParallelScale = 0.02008755131496027
renderView3.CameraParallelProjection = 1
# renderView3.UseColorPaletteForBackground = 0
renderView3.Background = [1.0, 1.0, 1.0]
renderView3.BackEnd = 'OSPRay raycaster'
renderView3.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

print("to setup view layouts")

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

if vertical_paper:

    # x_orig = 3143
    # y_orig = 2828

    # # compute the original widths and ratios 
    # panel_1_x = math.floor(0.4 * x_orig)
    # panel_1_y = y_orig

    # panel_2_x = 1356
    # panel_2_y = 1356 # math.floor(y_orig / 2)

    # panel_3_x = 1356
    # panel_3_y = 1356

    # width = x_orig
    # height = panel_1_y + panel_2_y + panel_3_y

    width = 1356 
    height = 1356 * 4

    # create new layout object 'Layout #1'
    layout1 = CreateLayout(name='Layout #1')
    layout1.SplitVertical(0, 0.5)
    layout1.AssignView(1, renderView1)
    layout1.SplitVertical(2, 0.5)
    layout1.AssignView(5, renderView2)
    layout1.AssignView(6, renderView3)
    layout1.SetSize((width, height))

else: 
    # default tiled 
    # create new layout object 'Layout #1'
    layout1 = CreateLayout(name='Layout #1')
    layout1.SplitHorizontal(0, 0.567738)
    layout1.AssignView(1, renderView1)
    layout1.SplitVertical(2, 0.500000)
    layout1.AssignView(5, renderView2)
    layout1.AssignView(6, renderView3)
    layout1.SetSize((3143, 2828))

print("through layout setup")

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

if '_96_' in os.getcwd(): 
    resolution_string = '96'
if '_192_' in os.getcwd(): 
    resolution_string = '192'
if '_384_' in os.getcwd(): 
    resolution_string = '384'
if '_768_' in os.getcwd(): 
    resolution_string = '768'

if version_sub_num == 8:
    font_size_colorbar = 48
    font_size_time = 8
elif version_sub_num == 11: 
    font_size_colorbar = 96
    font_size_time = 96

# create a new 'PVD Reader'
eulerian_varspvd = PVDReader(FileName='eulerian_vars.pvd')
eulerian_varspvd.CellArrays = ['P', 'U']

if not vertical_paper:
    # create a new 'Annotate Time Filter'
    annotateTimeFilter1 = AnnotateTimeFilter(registrationName='AnnotateTimeFilter1', Input=eulerian_varspvd)
    if version_sub_num == 8: 
        annotateTimeFilter1.Format = 't = %.3f s'
    elif version_sub_num == 11: 
        annotateTimeFilter1.Format = 't = {time:2.3f} s'

# create a new 'PVD Reader'
aortic_fish_192_cylinderpvd = PVDReader(FileName='aortic_fish_' + resolution_string + '_cylinder_faces.pvd')

# create a new 'PVD Reader'
aortic_fish_192_facespvd = PVDReader(FileName='aortic_fish_' + resolution_string + '_faces.pvd')

# create a new 'Slice'
slice4 = Slice(registrationName='Slice4', Input=eulerian_varspvd)
slice4.SliceType = 'Plane'
slice4.HyperTreeGridSlicer = 'Plane'
slice4.Triangulatetheslice = 0
slice4.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice4.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice4.HyperTreeGridSlicer.Origin = [-0.0005000000819563866, -0.0005000000819563866, 0.0334999980404973]

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=eulerian_varspvd)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.Triangulatetheslice = 0
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [-0.0005000000819563866, -0.0005000000819563866, 0.0334999980404973]

# create a new 'Clip'
clip2 = Clip(registrationName='Clip2', Input=aortic_fish_192_cylinderpvd)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', '']

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Origin = [2.9243528842926025e-07, 2.952292561531067e-07, 0.014446462504565716]

# added 
if version_sub_num == 8: 
    clip2.ClipType.Normal = [-1.0, 0.0, 0.0]
elif version_sub_num == 11:
    clip2.ClipType.Normal = [1.0, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [2.9243528842926025e-07, 2.952292561531067e-07, 0.014446462504565716]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from clip2
clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip2Display.Representation = 'Surface'
clip2Display.AmbientColor = [0.0, 0.0, 0.0]
clip2Display.ColorArrayName = [None, '']
clip2Display.DiffuseColor = [0.0, 0.0, 0.0]
# clip2Display.SelectTCoordArray = 'None'
# clip2Display.SelectNormalArray = 'None'
# clip2Display.SelectTangentArray = 'None'
clip2Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip2Display.SelectOrientationVectors = 'None'
clip2Display.ScaleFactor = 0.008889281935989857
clip2Display.SelectScaleArray = 'None'
clip2Display.GlyphType = 'Arrow'
clip2Display.GlyphTableIndexArray = 'None'
clip2Display.GaussianRadius = 0.00044446409679949284
clip2Display.SetScaleArray = [None, '']
clip2Display.ScaleTransferFunction = 'PiecewiseFunction'
clip2Display.OpacityArray = [None, '']
clip2Display.OpacityTransferFunction = 'PiecewiseFunction'
clip2Display.DataAxesGrid = 'GridAxesRepresentation'
clip2Display.PolarAxes = 'PolarAxesRepresentation'
clip2Display.ScalarOpacityUnitDistance = 0.001932391059992836
# clip2Display.OpacityArrayName = ['FIELD', 'avtOriginalBounds']
# clip2Display.SelectInputVectors = [None, '']
# clip2Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip2Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip2Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip2Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get 2D transfer function for 'U'
# uTF2D = GetTransferFunction2D('U')
# uTF2D.ScalarRangeInitialized = 1
# uTF2D.Range = [-12.0, 12.0, 0.0, 1.0]

# get color transfer function/color map for 'U'
uLUT = GetColorTransferFunction('U')
uLUT.AutomaticRescaleRangeMode = 'Never'
# uLUT.TransferFunction2D = uTF2D
uLUT.RGBPoints = [-12.0, 0.0862745098039216, 0.00392156862745098, 0.298039215686275, -11.271978588674958, 0.113725, 0.0235294, 0.45098, -10.667352128216567, 0.105882, 0.0509804, 0.509804, -10.24781457014383, 0.0392157, 0.0392157, 0.560784, -9.84061589806928, 0.0313725, 0.0980392, 0.6, -9.445757515095906, 0.0431373, 0.164706, 0.639216, -8.878149115735198, 0.054902, 0.243137, 0.678431, -8.125449007783285, 0.054902, 0.317647, 0.709804, -7.200000000000038, 0.0509804, 0.396078, 0.741176, -6.600000000000033, 0.0392157, 0.466667, 0.768627, -6.000000000000028, 0.0313725, 0.537255, 0.788235, -5.3737778098649835, 0.0313725, 0.615686, 0.811765, -4.732133288308978, 0.0235294, 0.709804, 0.831373, -4.09048736365005, 0.0509804, 0.8, 0.85098, -3.559897025387029, 0.0705882, 0.854902, 0.870588, -3.0663233451189242, 0.262745, 0.901961, 0.862745, -2.634446901047861, 0.423529, 0.941176, 0.87451, -1.9681232043923487, 0.572549, 0.964706, 0.835294, -1.5239064712201866, 0.658824, 0.980392, 0.843137, -1.203084210442178, 0.764706, 0.980392, 0.866667, -0.8575834761161332, 0.827451, 0.980392, 0.886275, -0.1789201919109793, 0.913725, 0.988235, 0.937255, 0.03084858712542271, 1.0, 1.0, 0.972549019607843, 0.24061736616183005, 0.988235, 0.980392, 0.870588, 0.549100740941503, 0.992156862745098, 0.972549019607843, 0.803921568627451, 0.7835479935259517, 0.992157, 0.964706, 0.713725, 1.1784063764993284, 0.988235, 0.956863, 0.643137, 1.8077120120570775, 0.980392, 0.917647, 0.509804, 2.338303753422945, 0.968627, 0.87451, 0.407843, 2.8812343807870704, 0.94902, 0.823529, 0.321569, 3.2760927637604293, 0.929412, 0.776471, 0.278431, 3.8560414522222324, 0.909804, 0.717647, 0.235294, 4.374292904486968, 0.890196, 0.658824, 0.196078, 4.799999999999997, 0.878431, 0.619608, 0.168627, 5.399999999999995, 0.870588, 0.54902, 0.156863, 5.9999999999999964, 0.85098, 0.47451, 0.145098, 6.599999999999994, 0.831373, 0.411765, 0.133333, 7.200000000000003, 0.811765, 0.345098, 0.113725, 7.799999999999997, 0.788235, 0.266667, 0.0941176, 8.399999999999999, 0.741176, 0.184314, 0.0745098, 9.0, 0.690196, 0.12549, 0.0627451, 9.600000000000001, 0.619608, 0.0627451, 0.0431373, 10.161439631544773, 0.54902, 0.027451, 0.0705882, 10.655012895266779, 0.470588, 0.0156863, 0.0901961, 11.210282781328374, 0.4, 0.00392157, 0.101961, 12.0, 0.188235294117647, 0.0, 0.0705882352941176]
uLUT.ColorSpace = 'Lab'
uLUT.ScalarRangeInitialized = 1.0
uLUT.VectorComponent = 2
uLUT.VectorMode = 'Component'

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'U']
slice1Display.LookupTable = uLUT
slice1Display.Ambient = 1.0
slice1Display.Diffuse = 0.0
# slice1Display.SelectTCoordArray = 'None'
# slice1Display.SelectNormalArray = 'None'
# slice1Display.SelectTangentArray = 'None'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = -2.0000000000000002e+298
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = -1e+297
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'
# slice1Display.SelectInputVectors = [None, '']
# slice1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from aortic_fish_192_facespvd
aortic_fish_192_facespvdDisplay = Show(aortic_fish_192_facespvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
aortic_fish_192_facespvdDisplay.Representation = 'Surface'
aortic_fish_192_facespvdDisplay.ColorArrayName = [None, '']
# aortic_fish_192_facespvdDisplay.SelectTCoordArray = 'None'
# aortic_fish_192_facespvdDisplay.SelectNormalArray = 'None'
# aortic_fish_192_facespvdDisplay.SelectTangentArray = 'None'
aortic_fish_192_facespvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
aortic_fish_192_facespvdDisplay.SelectOrientationVectors = 'None'
aortic_fish_192_facespvdDisplay.ScaleFactor = 0.00363107323646545
aortic_fish_192_facespvdDisplay.SelectScaleArray = 'None'
aortic_fish_192_facespvdDisplay.GlyphType = 'Arrow'
aortic_fish_192_facespvdDisplay.GlyphTableIndexArray = 'None'
aortic_fish_192_facespvdDisplay.GaussianRadius = 0.0001815536618232725
aortic_fish_192_facespvdDisplay.SetScaleArray = [None, '']
aortic_fish_192_facespvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
aortic_fish_192_facespvdDisplay.OpacityArray = [None, '']
aortic_fish_192_facespvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
aortic_fish_192_facespvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
aortic_fish_192_facespvdDisplay.PolarAxes = 'PolarAxesRepresentation'
aortic_fish_192_facespvdDisplay.ScalarOpacityUnitDistance = 0.002533878303323313
# aortic_fish_192_facespvdDisplay.OpacityArrayName = [None, '']
# aortic_fish_192_facespvdDisplay.SelectInputVectors = [None, '']
# aortic_fish_192_facespvdDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
aortic_fish_192_facespvdDisplay.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
aortic_fish_192_facespvdDisplay.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
aortic_fish_192_facespvdDisplay.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

if not vertical_paper:
    # show data from annotateTimeFilter1
    annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1, 'TextSourceRepresentation')

    # trace defaults for the display properties.
    if version_sub_num == 8:
        annotateTimeFilter1Display.WindowLocation = 'AnyLocation'
    elif version_sub_num == 11: 
        annotateTimeFilter1Display.WindowLocation = 'Any Location'
    annotateTimeFilter1Display.Position = [0.01, 0.66]
    annotateTimeFilter1Display.Color = [0.0, 0.0, 0.0]
    annotateTimeFilter1Display.Bold = 1
    annotateTimeFilter1Display.FontSize = font_size_time

if not vertical_paper:
    # setup the color legend parameters for each legend in this view
    # get color legend/bar for uLUT in view renderView1
    uLUTColorBar = GetScalarBar(uLUT, renderView1)
    uLUTColorBar.AutoOrient = 0
    if version_sub_num == 8:
        uLUTColorBar.WindowLocation = 'AnyLocation'
    elif version_sub_num == 11:
        uLUTColorBar.WindowLocation = 'Any Location'
    uLUTColorBar.Position = [0.056531851933813465, 0.7255508682325303]
    uLUTColorBar.Title = 'uz'
    uLUTColorBar.ComponentTitle = '(cm/s)'
    uLUTColorBar.HorizontalTitle = 1
    uLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
    uLUTColorBar.TitleBold = 1
    uLUTColorBar.TitleFontSize = font_size_colorbar
    uLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
    uLUTColorBar.LabelBold = 1
    uLUTColorBar.LabelFontSize = font_size_colorbar

    # added 
    uLUTColorBar.ScalarBarThickness = 24

    uLUTColorBar.ScalarBarLength = 0.19999999999999996
    uLUTColorBar.AutomaticLabelFormat = 0
    uLUTColorBar.LabelFormat = '%.0f'
    uLUTColorBar.UseCustomLabels = 1
    uLUTColorBar.CustomLabels = [-12.0, -8.0, -4.0, 0.0, 4.0, 8.0, 12.0]
    uLUTColorBar.AddRangeLabels = 0
    uLUTColorBar.RangeLabelFormat = '%-#6.0f'

    # set color bar visibility
    uLUTColorBar.Visibility = 1

    # show color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from aortic_fish_192_facespvd
aortic_fish_192_facespvdDisplay_1 = Show(aortic_fish_192_facespvd, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
aortic_fish_192_facespvdDisplay_1.Representation = 'Surface'
aortic_fish_192_facespvdDisplay_1.ColorArrayName = [None, '']
# aortic_fish_192_facespvdDisplay_1.SelectTCoordArray = 'None'
# aortic_fish_192_facespvdDisplay_1.SelectNormalArray = 'None'
# aortic_fish_192_facespvdDisplay_1.SelectTangentArray = 'None'
aortic_fish_192_facespvdDisplay_1.OSPRayScaleFunction = 'PiecewiseFunction'
aortic_fish_192_facespvdDisplay_1.SelectOrientationVectors = 'None'
aortic_fish_192_facespvdDisplay_1.ScaleFactor = 0.0036344625055790003
aortic_fish_192_facespvdDisplay_1.SelectScaleArray = 'None'
aortic_fish_192_facespvdDisplay_1.GlyphType = 'Arrow'
aortic_fish_192_facespvdDisplay_1.GlyphTableIndexArray = 'None'
aortic_fish_192_facespvdDisplay_1.GaussianRadius = 0.00018172312527895
aortic_fish_192_facespvdDisplay_1.SetScaleArray = [None, '']
aortic_fish_192_facespvdDisplay_1.ScaleTransferFunction = 'PiecewiseFunction'
aortic_fish_192_facespvdDisplay_1.OpacityArray = [None, '']
aortic_fish_192_facespvdDisplay_1.OpacityTransferFunction = 'PiecewiseFunction'
aortic_fish_192_facespvdDisplay_1.DataAxesGrid = 'GridAxesRepresentation'
aortic_fish_192_facespvdDisplay_1.PolarAxes = 'PolarAxesRepresentation'
aortic_fish_192_facespvdDisplay_1.ScalarOpacityUnitDistance = 0.002538281472507009
# aortic_fish_192_facespvdDisplay_1.OpacityArrayName = [None, '']
# aortic_fish_192_facespvdDisplay_1.SelectInputVectors = [None, '']
# aortic_fish_192_facespvdDisplay_1.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
aortic_fish_192_facespvdDisplay_1.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
aortic_fish_192_facespvdDisplay_1.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
aortic_fish_192_facespvdDisplay_1.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from slice4
slice4Display = Show(slice4, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
slice4Display.Representation = 'Surface'
slice4Display.ColorArrayName = ['CELLS', 'U']
slice4Display.LookupTable = uLUT
slice4Display.Ambient = 1.0
slice4Display.Diffuse = 0.0
# slice4Display.SelectTCoordArray = 'None'
# slice4Display.SelectNormalArray = 'None'
# slice4Display.SelectTangentArray = 'None'
slice4Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice4Display.SelectOrientationVectors = 'None'
slice4Display.ScaleFactor = 0.006300000287592412
slice4Display.SelectScaleArray = 'P'
slice4Display.GlyphType = 'Arrow'
slice4Display.GlyphTableIndexArray = 'P'
slice4Display.GaussianRadius = 0.00031500001437962055
slice4Display.SetScaleArray = [None, '']
slice4Display.ScaleTransferFunction = 'PiecewiseFunction'
slice4Display.OpacityArray = [None, '']
slice4Display.OpacityTransferFunction = 'PiecewiseFunction'
slice4Display.DataAxesGrid = 'GridAxesRepresentation'
slice4Display.PolarAxes = 'PolarAxesRepresentation'
# slice4Display.SelectInputVectors = [None, '']
# slice4Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice4Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice4Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice4Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from aortic_fish_192_cylinderpvd
aortic_fish_192_cylinderpvdDisplay = Show(aortic_fish_192_cylinderpvd, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
aortic_fish_192_cylinderpvdDisplay.Representation = 'Surface'
aortic_fish_192_cylinderpvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
aortic_fish_192_cylinderpvdDisplay.ColorArrayName = [None, '']
aortic_fish_192_cylinderpvdDisplay.DiffuseColor = [0.0, 0.0, 0.0]
# aortic_fish_192_cylinderpvdDisplay.SelectTCoordArray = 'None'
# aortic_fish_192_cylinderpvdDisplay.SelectNormalArray = 'None'
# aortic_fish_192_cylinderpvdDisplay.SelectTangentArray = 'None'
aortic_fish_192_cylinderpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
aortic_fish_192_cylinderpvdDisplay.SelectOrientationVectors = 'None'
aortic_fish_192_cylinderpvdDisplay.ScaleFactor = 0.01265398170799017
aortic_fish_192_cylinderpvdDisplay.SelectScaleArray = 'None'
aortic_fish_192_cylinderpvdDisplay.GlyphType = 'Arrow'
aortic_fish_192_cylinderpvdDisplay.GlyphTableIndexArray = 'None'
aortic_fish_192_cylinderpvdDisplay.GaussianRadius = 0.0006326990853995085
aortic_fish_192_cylinderpvdDisplay.SetScaleArray = [None, '']
aortic_fish_192_cylinderpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
aortic_fish_192_cylinderpvdDisplay.OpacityArray = [None, '']
aortic_fish_192_cylinderpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
aortic_fish_192_cylinderpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
aortic_fish_192_cylinderpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
aortic_fish_192_cylinderpvdDisplay.ScalarOpacityUnitDistance = 0.0033256435452326512
# aortic_fish_192_cylinderpvdDisplay.OpacityArrayName = [None, '']
# aortic_fish_192_cylinderpvdDisplay.SelectInputVectors = [None, '']
# aortic_fish_192_cylinderpvdDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
aortic_fish_192_cylinderpvdDisplay.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
aortic_fish_192_cylinderpvdDisplay.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
aortic_fish_192_cylinderpvdDisplay.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView3'
# ----------------------------------------------------------------

# show data from aortic_fish_192_facespvd
aortic_fish_192_facespvdDisplay_2 = Show(aortic_fish_192_facespvd, renderView3, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
aortic_fish_192_facespvdDisplay_2.Representation = 'Surface'
aortic_fish_192_facespvdDisplay_2.ColorArrayName = [None, '']
# aortic_fish_192_facespvdDisplay_2.SelectTCoordArray = 'None'
# aortic_fish_192_facespvdDisplay_2.SelectNormalArray = 'None'
# aortic_fish_192_facespvdDisplay_2.SelectTangentArray = 'None'
aortic_fish_192_facespvdDisplay_2.OSPRayScaleFunction = 'PiecewiseFunction'
aortic_fish_192_facespvdDisplay_2.SelectOrientationVectors = 'None'
aortic_fish_192_facespvdDisplay_2.ScaleFactor = 0.0036431187763810106
aortic_fish_192_facespvdDisplay_2.SelectScaleArray = 'None'
aortic_fish_192_facespvdDisplay_2.GlyphType = 'Arrow'
aortic_fish_192_facespvdDisplay_2.GlyphTableIndexArray = 'None'
aortic_fish_192_facespvdDisplay_2.GaussianRadius = 0.0001821559388190505
aortic_fish_192_facespvdDisplay_2.SetScaleArray = [None, '']
aortic_fish_192_facespvdDisplay_2.ScaleTransferFunction = 'PiecewiseFunction'
aortic_fish_192_facespvdDisplay_2.OpacityArray = [None, '']
aortic_fish_192_facespvdDisplay_2.OpacityTransferFunction = 'PiecewiseFunction'
aortic_fish_192_facespvdDisplay_2.DataAxesGrid = 'GridAxesRepresentation'
aortic_fish_192_facespvdDisplay_2.PolarAxes = 'PolarAxesRepresentation'
aortic_fish_192_facespvdDisplay_2.ScalarOpacityUnitDistance = 0.002546334423661121
# aortic_fish_192_facespvdDisplay_2.OpacityArrayName = [None, '']
# aortic_fish_192_facespvdDisplay_2.SelectInputVectors = [None, '']
# aortic_fish_192_facespvdDisplay_2.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
aortic_fish_192_facespvdDisplay_2.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
aortic_fish_192_facespvdDisplay_2.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
aortic_fish_192_facespvdDisplay_2.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'U'
uPWF = GetOpacityTransferFunction('U')
uPWF.Points = [-12.0, 0.2589285671710968, 0.5, 0.0, 12.0, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# restore active source
if not vertical_paper:
    SetActiveSource(annotateTimeFilter1)
# ----------------------------------------------------------------


tk = GetTimeKeeper()
timesteps = tk.TimestepValues
numTimesteps = len(timesteps)

if len(sys.argv) >= 2:
    basename = sys.argv[1]
else: 
    basename = 'frames'    

if vertical_paper: 
    print('processing vertical')
    basename += '_vertical'

if len(sys.argv) >= 4:
    nprocs = int(sys.argv[2])
    proc_num = int(sys.argv[3])
else: 
    print("using default proc_num 0, nprocs = 1")
    proc_num = 0
    nprocs = 1

# old default (5120, 4320)
if vertical_paper:
    res = (width, height) 
else: 
    res = (3143, 2828) 

animationScene1 = GetAnimationScene()
# animationScene1.GoToFirst()

# run first iteration, then stop and start over 
for frame in range(len(timesteps)):
    print("frame = ", frame)
    if (frame % nprocs) == proc_num:
        animationScene1.AnimationTime = timesteps[frame]

        Render()
        SaveScreenshot(basename + str(frame).zfill(4) + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)
        print("about to break loop 1")
        break 

for frame in range(len(timesteps)-1):
    print("frame = ", frame)
    if (frame % nprocs) == proc_num:
        animationScene1.AnimationTime = timesteps[frame]

        Render()
        SaveScreenshot(basename + str(frame).zfill(4) + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)
        print("saved frame ", frame)

