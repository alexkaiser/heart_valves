# state file generated using paraview version 5.8.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import pdb
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [700, 1172]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [1.5652042627334595, -18.485937118530273, -22.303841590881348]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.6438031846839976, -36.99706527434366, -23.08780017268454]
renderView1.CameraFocalPoint = [0.6438031846839976, 1.516905378849943, -23.08780017268454]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 8.372844950840651
renderView1.Background = [1.0, 0.9999694819562066, 0.9999847409781033]
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [700, 1172]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [1.565246820449829, -18.48591709136963, -22.37158489227295]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [6.9233121594997415, -29.053646876183148, -12.506256321311906]
renderView2.CameraFocalPoint = [-8.385050586005658, -3.271773237536895, -35.45672074468072]
renderView2.CameraViewUp = [-0.40854053478514424, 0.45946975997037115, 0.7886584629040875]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 9.972453281730647
renderView2.Background = [1.0, 0.9999694819562066, 0.9999847409781033]
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.500000)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, renderView2)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView2)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
aorta_384pvd = PVDReader(FileName='aorta_384.pvd')

# create a new 'Text'
# text1 = Text()
# text1.Text = 'vertical component of velocity'

# # create a new 'Text'
# text2 = Text()
# text2.Text = 'velocity normal to slice'

# create a new 'Annotate Time'
# annotateTime1 = AnnotateTime()
# annotateTime1.Format = 't = %.3f s'

# interpolated_mesh = False 
# if interpolated_mesh:
#     eulerian_name = 'eulerian_vars.pvd'
# else: 
#     eulerian_name = 'eulerian_vars_restricted_cells.pvd'

cell_data = False

if cell_data:
    # create a new 'PVD Reader'
    # eulerian_varspvd = PVDReader(FileName='eulerian_vars.pvd')
    eulerian_varspvd = PVDReader(FileName='eulerian_vars_restricted_cells.pvd')
    eulerian_varspvd.CellArrays = ['U']
else: 
    eulerian_varspvd = PVDReader(FileName='eulerian_vars_restricted_points.pvd')
    eulerian_varspvd.PointArrays = ['U']

# create a new 'Slice'
slice9 = Slice(Input=eulerian_varspvd)
slice9.SliceType = 'Plane'
slice9.HyperTreeGridSlicer = 'Plane'
slice9.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice9.SliceType.Origin = [1.0154999494552612, -18.4862003326416, -22.669400215148926]
slice9.SliceType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice9.HyperTreeGridSlicer.Origin = [1.0154999494552612, -18.4862003326416, -22.669400215148926]

# create a new 'XML Unstructured Grid Reader'
# aorta_384_volumetricmeshvtu = XMLUnstructuredGridReader(FileName=['aorta_384_volumetric.mesh.vtu'])
# aorta_384_volumetricmeshvtu.CellArrayStatus = ['ModelRegionID', 'GlobalElementID']
# aorta_384_volumetricmeshvtu.PointArrayStatus = ['GlobalNodeID']

# # create a new 'Extract Selection'
# extractSelection1 = ExtractSelection(Input=aorta_384_volumetricmeshvtu)

# create a new 'PVD Reader'
aortic_no_partition_384_facespvd = PVDReader(FileName='aortic_no_partition_384_faces.pvd')

# create a new 'Resample With Dataset'
# resampleWithDataset1 = ResampleWithDataset(SourceDataArrays=eulerian_varspvd,
#     DestinationMesh=aorta_384_volumetricmeshvtu)
# resampleWithDataset1.CellLocator = 'Static Cell Locator'

# pick the right mesh 
# if interpolated_mesh:
#     eulerian_source = resampleWithDataset1
# else: 
#     eulerian_source = eulerian_varspvd

# create a new 'Calculator'
calculator5 = Calculator(Input=eulerian_varspvd)
if cell_data:
    calculator5.AttributeType = 'Cell Data'
calculator5.ResultArrayName = 'contour_6'
calculator5.Function = '-0.0974948538117627 * U_X + -0.028568296191069 * U_Y + 0.994825917401111* U_Z'

# create a new 'Slice'
slice7 = Slice(Input=calculator5)
slice7.SliceType = 'Plane'
slice7.HyperTreeGridSlicer = 'Plane'
slice7.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice7.SliceType.Origin = [0.0488354041365377, -18.3884646217316, -20.8624742881165]
slice7.SliceType.Normal = [-0.0974948538117627, -0.028568296191069, 0.994825917401111]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice7.HyperTreeGridSlicer.Origin = [1.5904362797737122, -18.485923767089844, -22.346393585205078]

# create a new 'Calculator'
calculator3 = Calculator(Input=eulerian_varspvd)
if cell_data:
    calculator3.AttributeType = 'Cell Data'
calculator3.ResultArrayName = 'contour_15'
calculator3.Function = '-0.758899471490809 * U_X + 0.0231783966274574 * U_Y + 0.650795170618799 * U_Z'

# create a new 'Slice'
slice5 = Slice(Input=calculator3)
slice5.SliceType = 'Plane'
slice5.HyperTreeGridSlicer = 'Plane'
slice5.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice5.SliceType.Origin = [1.9166492846547107, -18.4721636302181, -24.7070353270031]
slice5.SliceType.Normal = [-0.758899471490809, 0.0231783966274574, 0.650795170618799]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice5.HyperTreeGridSlicer.Origin = [1.5904362797737122, -18.485923767089844, -22.346393585205078]

# create a new 'Calculator'
calculator1 = Calculator(Input=eulerian_varspvd)
if cell_data:
    calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'annulus_normal_projected'
calculator1.Function = '-0.644118902159037 * U_X + 0.007931209298900 * U_Y + 0.764884263010094 * U_Z'

# create a new 'Slice'
slice3 = Slice(Input=calculator1)
slice3.SliceType = 'Plane'
slice3.HyperTreeGridSlicer = 'Plane'
slice3.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice3.SliceType.Origin = [2.447601610336508, -18.376013668463226, -25.219513937752243]
slice3.SliceType.Normal = [-0.644118902159037, 0.0079312092989, 0.764884263010094]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice3.HyperTreeGridSlicer.Origin = [1.5904362797737122, -18.485923767089844, -22.346393585205078]

# create a new 'Calculator'
calculator6 = Calculator(Input=eulerian_varspvd)
if cell_data:
    calculator6.AttributeType = 'Cell Data'
calculator6.Function = '0.277454395233418 * U_X + 0.0616022174821306*U_Y + 0.958761818892963*U_Z'

# create a new 'Slice'
slice8 = Slice(Input=calculator6)
slice8.SliceType = 'Plane'
slice8.HyperTreeGridSlicer = 'Plane'
slice8.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice8.SliceType.Origin = [0.2803383547117727, -18.3479701228845, -18.8436930740923]
slice8.SliceType.Normal = [0.277454395233418, 0.0616022174821306, 0.958761818892963]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice8.HyperTreeGridSlicer.Origin = [1.5904362797737122, -18.485923767089844, -22.346393585205078]

# create a new 'Slice'
slice2 = Slice(Input=aorta_384pvd)
slice2.SliceType = 'Plane'
slice2.HyperTreeGridSlicer = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [1.565198004245758, -18.485950469970703, -22.37162494659424]
slice2.SliceType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice2.HyperTreeGridSlicer.Origin = [1.565198004245758, -18.485950469970703, -22.37162494659424]

# create a new 'Calculator'
calculator2 = Calculator(Input=eulerian_varspvd)
if cell_data:
    calculator2.AttributeType = 'Cell Data'
calculator2.ResultArrayName = 'coutour_12'
calculator2.Function = '-0.546829510850715 * U_X + 0.0721783996060898 * U_Y + 0.834126947588358 * U_Z'

# create a new 'Slice'
slice4 = Slice(Input=calculator2)
slice4.SliceType = 'Plane'
slice4.HyperTreeGridSlicer = 'Plane'
slice4.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice4.SliceType.Origin = [1.1135387247389723, -18.4115922051707, -23.8033429814558]
slice4.SliceType.Normal = [-0.546829510850715, 0.0721783996060898, 0.834126947588358]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice4.HyperTreeGridSlicer.Origin = [1.5904362797737122, -18.485923767089844, -22.346393585205078]

# create a new 'Calculator'
calculator4 = Calculator(Input=eulerian_varspvd)
if cell_data:
    calculator4.AttributeType = 'Cell Data'
calculator4.ResultArrayName = 'countour_9'
calculator4.Function = '-0.369605687499951 * U_X + 0.0302064524432408*U_Y + 0.928697585868771*U_Z'

# create a new 'Slice'
slice6 = Slice(Input=calculator4)
slice6.SliceType = 'Plane'
slice6.HyperTreeGridSlicer = 'Plane'
slice6.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice6.SliceType.Origin = [0.5826182201137723, -18.3351820981491, -22.7279221637119]
slice6.SliceType.Normal = [-0.369605687499951, 0.0302064524432408, 0.928697585868771]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice6.HyperTreeGridSlicer.Origin = [1.5904362797737122, -18.485923767089844, -22.346393585205078]

# create a new 'Slice'
slice1 = Slice(Input=eulerian_varspvd)
if cell_data:
    slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [1.590436339378357, -18.485923767089844, -22.346393585205078]
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [1.590436339378357, -18.485923767089844, -22.346393585205078]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'U'
uLUT = GetColorTransferFunction('U')
uLUT.AutomaticRescaleRangeMode = 'Never'
uLUT.RGBPoints = [-300.0, 0.0862745098039216, 0.00392156862745098, 0.298039215686275, -281.799464716874, 0.113725, 0.0235294, 0.45098, -266.68380320541416, 0.105882, 0.0509804, 0.509804, -256.1953642535958, 0.0392157, 0.0392157, 0.560784, -246.015397451732, 0.0313725, 0.0980392, 0.6, -236.1439378773976, 0.0431373, 0.164706, 0.639216, -221.95372789337992, 0.054902, 0.243137, 0.678431, -203.1362251945821, 0.054902, 0.317647, 0.709804, -180.00000000000094, 0.0509804, 0.396078, 0.741176, -165.00000000000082, 0.0392157, 0.466667, 0.768627, -150.00000000000074, 0.0313725, 0.537255, 0.788235, -134.3444452466246, 0.0313725, 0.615686, 0.811765, -118.30333220772445, 0.0235294, 0.709804, 0.831373, -102.26218409125124, 0.0509804, 0.8, 0.85098, -88.99742563467572, 0.0705882, 0.854902, 0.870588, -76.65808362797307, 0.262745, 0.901961, 0.862745, -65.8611725261965, 0.423529, 0.941176, 0.87451, -49.203080109808695, 0.572549, 0.964706, 0.835294, -38.09766178050461, 0.658824, 0.980392, 0.843137, -30.0771052610545, 0.764706, 0.980392, 0.866667, -21.43958690290333, 0.827451, 0.980392, 0.886275, -4.473004797774422, 0.913725, 0.988235, 0.937255, 0.7712146781356068, 1.0, 1.0, 0.972549019607843, 6.015434154045749, 0.988235, 0.980392, 0.870588, 13.727518523537583, 0.992156862745098, 0.972549019607843, 0.803921568627451, 19.588699838148784, 0.992157, 0.964706, 0.713725, 29.460159412483222, 0.988235, 0.956863, 0.643137, 45.192800301426985, 0.980392, 0.917647, 0.509804, 58.457593835573675, 0.968627, 0.87451, 0.407843, 72.03085951967677, 0.94902, 0.823529, 0.321569, 81.90231909401075, 0.929412, 0.776471, 0.278431, 96.4010363055558, 0.909804, 0.717647, 0.235294, 109.35732261217419, 0.890196, 0.658824, 0.196078, 119.99999999999989, 0.878431, 0.619608, 0.168627, 134.99999999999983, 0.870588, 0.54902, 0.156863, 149.99999999999994, 0.85098, 0.47451, 0.145098, 164.99999999999994, 0.831373, 0.411765, 0.133333, 180.0, 0.811765, 0.345098, 0.113725, 195.0, 0.788235, 0.266667, 0.0941176, 210.0, 0.741176, 0.184314, 0.0745098, 225.0, 0.690196, 0.12549, 0.0627451, 240.0, 0.619608, 0.0627451, 0.0431373, 254.03599078861942, 0.54902, 0.027451, 0.0705882, 266.37532238166943, 0.470588, 0.0156863, 0.0901961, 280.25706953320935, 0.4, 0.00392157, 0.101961, 300.0, 0.188235294117647, 0.0, 0.0705882352941176]
uLUT.ColorSpace = 'Lab'
uLUT.ScalarRangeInitialized = 1.0
uLUT.VectorComponent = 2
uLUT.VectorMode = 'Component'

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
if cell_data:
    slice1Display.ColorArrayName = ['CELLS', 'U']
else: 
    slice1Display.ColorArrayName = ['POINTS', 'U']
slice1Display.LookupTable = uLUT
slice1Display.Ambient = 1.0
slice1Display.Diffuse = 0.0
slice1Display.OSPRayScaleArray = 'U'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 1.0412282943725586
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.05206141471862793
slice1Display.SetScaleArray = ['CELLS', 'U']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['CELLS', 'U']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [-153.36936950683594, 0.2589285671710968, 0.5, 0.0, 96.7054443359375, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [-153.36936950683594, 0.2589285671710968, 0.5, 0.0, 96.7054443359375, 1.0, 0.5, 0.0]

# show data from slice2
slice2Display = Show(slice2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice2Display.Representation = 'Surface'
slice2Display.AmbientColor = [0.0, 0.0, 0.0]
slice2Display.ColorArrayName = [None, '']
slice2Display.DiffuseColor = [0.0, 0.0, 0.0]
slice2Display.PointSize = 4.0
slice2Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice2Display.SelectOrientationVectors = 'None'
slice2Display.ScaleFactor = 1.046976661682129
slice2Display.SelectScaleArray = 'None'
slice2Display.GlyphType = 'Arrow'
slice2Display.GlyphTableIndexArray = 'None'
slice2Display.GaussianRadius = 0.05234883308410645
slice2Display.SetScaleArray = [None, '']
slice2Display.ScaleTransferFunction = 'PiecewiseFunction'
slice2Display.OpacityArray = [None, '']
slice2Display.OpacityTransferFunction = 'PiecewiseFunction'
slice2Display.DataAxesGrid = 'GridAxesRepresentation'
slice2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice2Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice2Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice2Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from aortic_no_partition_384_facespvd
aortic_no_partition_384_facespvdDisplay = Show(aortic_no_partition_384_facespvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
aortic_no_partition_384_facespvdDisplay.Representation = 'Surface'
aortic_no_partition_384_facespvdDisplay.Opacity = 0.2
aortic_no_partition_384_facespvdDisplay.ColorArrayName = [None, '']
# aortic_no_partition_384_facespvdDisplay.SelectTCoordArray = 'None'
# aortic_no_partition_384_facespvdDisplay.SelectNormalArray = 'None'
# aortic_no_partition_384_facespvdDisplay.SelectTangentArray = 'None'
aortic_no_partition_384_facespvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
aortic_no_partition_384_facespvdDisplay.SelectOrientationVectors = 'None'
aortic_no_partition_384_facespvdDisplay.ScaleFactor = 0.30792751312255007
aortic_no_partition_384_facespvdDisplay.SelectScaleArray = 'None'
aortic_no_partition_384_facespvdDisplay.GlyphType = 'Arrow'
aortic_no_partition_384_facespvdDisplay.GlyphTableIndexArray = 'None'
aortic_no_partition_384_facespvdDisplay.GaussianRadius = 0.015396375656127503
aortic_no_partition_384_facespvdDisplay.SetScaleArray = [None, '']
aortic_no_partition_384_facespvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
aortic_no_partition_384_facespvdDisplay.OpacityArray = [None, '']
aortic_no_partition_384_facespvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
aortic_no_partition_384_facespvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
aortic_no_partition_384_facespvdDisplay.PolarAxes = 'PolarAxesRepresentation'
aortic_no_partition_384_facespvdDisplay.ScalarOpacityUnitDistance = 0.136635908798675
# aortic_no_partition_384_facespvdDisplay.OpacityArrayName = [None, '']
# aortic_no_partition_384_facespvdDisplay.SelectInputVectors = [None, '']
# aortic_no_partition_384_facespvdDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
aortic_no_partition_384_facespvdDisplay.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
aortic_no_partition_384_facespvdDisplay.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
aortic_no_partition_384_facespvdDisplay.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from annotateTime1
# annotateTime1Display = Show(annotateTime1, renderView1, 'TextSourceRepresentation')

# # trace defaults for the display properties.
# annotateTime1Display.Color = [0.0, 0.0, 0.0]
# annotateTime1Display.FontFamily = 'Times'
# annotateTime1Display.FontSize = 9
# annotateTime1Display.WindowLocation = 'AnyLocation'
# annotateTime1Display.Position = [0.02, 0.07]

# show data from text1
# text1Display = Show(text1, renderView1, 'TextSourceRepresentation')

# # trace defaults for the display properties.
# text1Display.Color = [0.0, 0.0, 0.0]
# text1Display.FontFamily = 'Times'
# text1Display.FontSize = 9
# text1Display.WindowLocation = 'AnyLocation'
# text1Display.Position = [0.02, 0.02]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uLUT in view renderView1
# uLUTColorBar = GetScalarBar(uLUT, renderView1)
# uLUTColorBar.WindowLocation = 'AnyLocation'
# uLUTColorBar.Position = [0.02, 0.15]
# uLUTColorBar.Title = 'velocity (cm/s)'
# uLUTColorBar.ComponentTitle = ''
# uLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
# uLUTColorBar.TitleFontFamily = 'Times'
# uLUTColorBar.TitleFontSize = 26
# uLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
# uLUTColorBar.LabelFontFamily = 'Times'
# uLUTColorBar.LabelFontSize = 18
# uLUTColorBar.AutomaticLabelFormat = 0
# uLUTColorBar.LabelFormat = '%-#6.0f'
# uLUTColorBar.UseCustomLabels = 1
# uLUTColorBar.CustomLabels = [-300.0, -150.0, 0.0, 150.0, 300.0]
# uLUTColorBar.RangeLabelFormat = '%-#6.0f'
# uLUTColorBar.ScalarBarThickness = 24
# uLUTColorBar.ScalarBarLength = 0.19999999999999996

# # set color bar visibility
# uLUTColorBar.Visibility = 1

# show color legend
# slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from aortic_no_partition_384pvd
aortic_no_partition_384_facespvdDisplay_1 = Show(aortic_no_partition_384_facespvd, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
aortic_no_partition_384_facespvdDisplay_1.Representation = 'Surface'
aortic_no_partition_384_facespvdDisplay_1.ColorArrayName = [None, '']
aortic_no_partition_384_facespvdDisplay_1.Opacity = 0.2
# aortic_no_partition_384_facespvdDisplay_1.SelectTCoordArray = 'None'
# aortic_no_partition_384_facespvdDisplay_1.SelectNormalArray = 'None'
# aortic_no_partition_384_facespvdDisplay_1.SelectTangentArray = 'None'
aortic_no_partition_384_facespvdDisplay_1.OSPRayScaleFunction = 'PiecewiseFunction'
aortic_no_partition_384_facespvdDisplay_1.SelectOrientationVectors = 'None'
aortic_no_partition_384_facespvdDisplay_1.ScaleFactor = 0.30792751312255007
aortic_no_partition_384_facespvdDisplay_1.SelectScaleArray = 'None'
aortic_no_partition_384_facespvdDisplay_1.GlyphType = 'Arrow'
aortic_no_partition_384_facespvdDisplay_1.GlyphTableIndexArray = 'None'
aortic_no_partition_384_facespvdDisplay_1.GaussianRadius = 0.015396375656127503
aortic_no_partition_384_facespvdDisplay_1.SetScaleArray = [None, '']
aortic_no_partition_384_facespvdDisplay_1.ScaleTransferFunction = 'PiecewiseFunction'
aortic_no_partition_384_facespvdDisplay_1.OpacityArray = [None, '']
aortic_no_partition_384_facespvdDisplay_1.OpacityTransferFunction = 'PiecewiseFunction'
aortic_no_partition_384_facespvdDisplay_1.DataAxesGrid = 'GridAxesRepresentation'
aortic_no_partition_384_facespvdDisplay_1.PolarAxes = 'PolarAxesRepresentation'
aortic_no_partition_384_facespvdDisplay_1.ScalarOpacityUnitDistance = 0.136635908798675
# aortic_no_partition_384_facespvdDisplay_1.OpacityArrayName = [None, '']
# aortic_no_partition_384_facespvdDisplay_1.SelectInputVectors = [None, '']
# aortic_no_partition_384_facespvdDisplay_1.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
aortic_no_partition_384_facespvdDisplay_1.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
aortic_no_partition_384_facespvdDisplay_1.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
aortic_no_partition_384_facespvdDisplay_1.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]


# show data from aorta_384pvd
aorta_384pvdDisplay = Show(aorta_384pvd, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
aorta_384pvdDisplay.Representation = 'Surface'
# aorta_384pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
aorta_384pvdDisplay.ColorArrayName = [None, '']
# aorta_384pvdDisplay.DiffuseColor = [0.0, 0.0, 0.0]
aorta_384pvdDisplay.Opacity = 0.02
aorta_384pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
aorta_384pvdDisplay.SelectOrientationVectors = 'None'
aorta_384pvdDisplay.ScaleFactor = 1.06043758392334
aorta_384pvdDisplay.SelectScaleArray = 'None'
aorta_384pvdDisplay.GlyphType = 'Arrow'
aorta_384pvdDisplay.GlyphTableIndexArray = 'None'
aorta_384pvdDisplay.GaussianRadius = 0.05302187919616699
aorta_384pvdDisplay.SetScaleArray = [None, '']
aorta_384pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
aorta_384pvdDisplay.OpacityArray = [None, '']
aorta_384pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
aorta_384pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
aorta_384pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
aorta_384pvdDisplay.ScalarOpacityUnitDistance = 0.11505942928416406

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
aorta_384pvdDisplay.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
aorta_384pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
aorta_384pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from slice3
slice3Display = Show(slice3, renderView2, 'GeometryRepresentation')

# get color transfer function/color map for 'annulus_normal_projected'
annulus_normal_projectedLUT = GetColorTransferFunction('annulus_normal_projected')
annulus_normal_projectedLUT.AutomaticRescaleRangeMode = 'Never'
annulus_normal_projectedLUT.RGBPoints = [-300.0, 0.0862745098039216, 0.00392156862745098, 0.298039215686275, -281.79946471687396, 0.113725, 0.0235294, 0.45098, -266.68380320541416, 0.105882, 0.0509804, 0.509804, -256.1953642535958, 0.0392157, 0.0392157, 0.560784, -246.015397451732, 0.0313725, 0.0980392, 0.6, -236.1439378773976, 0.0431373, 0.164706, 0.639216, -221.95372789337995, 0.054902, 0.243137, 0.678431, -203.1362251945821, 0.054902, 0.317647, 0.709804, -180.00000000000094, 0.0509804, 0.396078, 0.741176, -165.00000000000082, 0.0392157, 0.466667, 0.768627, -150.0000000000007, 0.0313725, 0.537255, 0.788235, -134.3444452466246, 0.0313725, 0.615686, 0.811765, -118.30333220772445, 0.0235294, 0.709804, 0.831373, -102.26218409125124, 0.0509804, 0.8, 0.85098, -88.99742563467572, 0.0705882, 0.854902, 0.870588, -76.65808362797307, 0.262745, 0.901961, 0.862745, -65.8611725261965, 0.423529, 0.941176, 0.87451, -49.203080109808724, 0.572549, 0.964706, 0.835294, -38.09766178050461, 0.658824, 0.980392, 0.843137, -30.0771052610545, 0.764706, 0.980392, 0.866667, -21.43958690290333, 0.827451, 0.980392, 0.886275, -4.473004797774422, 0.913725, 0.988235, 0.937255, 0.7712146781356068, 1.0, 1.0, 0.972549019607843, 6.015434154045749, 0.988235, 0.980392, 0.870588, 13.727518523537583, 0.992156862745098, 0.972549019607843, 0.803921568627451, 19.588699838148784, 0.992157, 0.964706, 0.713725, 29.460159412483222, 0.988235, 0.956863, 0.643137, 45.192800301426985, 0.980392, 0.917647, 0.509804, 58.457593835573675, 0.968627, 0.87451, 0.407843, 72.03085951967677, 0.94902, 0.823529, 0.321569, 81.90231909401075, 0.929412, 0.776471, 0.278431, 96.4010363055558, 0.909804, 0.717647, 0.235294, 109.35732261217419, 0.890196, 0.658824, 0.196078, 119.99999999999989, 0.878431, 0.619608, 0.168627, 134.99999999999983, 0.870588, 0.54902, 0.156863, 149.99999999999994, 0.85098, 0.47451, 0.145098, 164.9999999999999, 0.831373, 0.411765, 0.133333, 180.0, 0.811765, 0.345098, 0.113725, 195.0, 0.788235, 0.266667, 0.0941176, 210.0, 0.741176, 0.184314, 0.0745098, 225.0, 0.690196, 0.12549, 0.0627451, 240.0, 0.619608, 0.0627451, 0.0431373, 254.03599078861942, 0.54902, 0.027451, 0.0705882, 266.37532238166943, 0.470588, 0.0156863, 0.0901961, 280.25706953320935, 0.4, 0.00392157, 0.101961, 300.0, 0.188235294117647, 0.0, 0.0705882352941176]
annulus_normal_projectedLUT.ColorSpace = 'Lab'
annulus_normal_projectedLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice3Display.Representation = 'Surface'
if cell_data:
    slice3Display.ColorArrayName = ['CELLS', 'annulus_normal_projected']
else:
    slice3Display.ColorArrayName = ['POINTS', 'annulus_normal_projected']
slice3Display.LookupTable = annulus_normal_projectedLUT
slice3Display.Ambient = 1.0
slice3Display.Diffuse = 0.0
slice3Display.OSPRayScaleArray = 'annulus_normal_projected'
slice3Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice3Display.SelectOrientationVectors = 'None'
slice3Display.ScaleFactor = 0.2277305603027344
slice3Display.SelectScaleArray = 'annulus_normal_projected'
slice3Display.GlyphType = 'Arrow'
slice3Display.GlyphTableIndexArray = 'annulus_normal_projected'
slice3Display.GaussianRadius = 0.011386528015136718
slice3Display.SetScaleArray = ['CELLS', 'annulus_normal_projected']
slice3Display.ScaleTransferFunction = 'PiecewiseFunction'
slice3Display.OpacityArray = ['CELLS', 'annulus_normal_projected']
slice3Display.OpacityTransferFunction = 'PiecewiseFunction'
slice3Display.DataAxesGrid = 'GridAxesRepresentation'
slice3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice3Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice3Display.ScaleTransferFunction.Points = [-8.173152896346384, 0.2589285671710968, 0.5, 0.0, 141.10649334402115, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice3Display.OpacityTransferFunction.Points = [-8.173152896346384, 0.2589285671710968, 0.5, 0.0, 141.10649334402115, 1.0, 0.5, 0.0]

# show data from slice4
slice4Display = Show(slice4, renderView2, 'GeometryRepresentation')

# get color transfer function/color map for 'coutour_12'
coutour_12LUT = GetColorTransferFunction('coutour_12')
coutour_12LUT.AutomaticRescaleRangeMode = 'Never'
coutour_12LUT.RGBPoints = [-300.0, 0.0862745098039216, 0.00392156862745098, 0.298039215686275, -281.79946471687396, 0.113725, 0.0235294, 0.45098, -266.68380320541416, 0.105882, 0.0509804, 0.509804, -256.1953642535958, 0.0392157, 0.0392157, 0.560784, -246.015397451732, 0.0313725, 0.0980392, 0.6, -236.1439378773976, 0.0431373, 0.164706, 0.639216, -221.95372789337995, 0.054902, 0.243137, 0.678431, -203.1362251945821, 0.054902, 0.317647, 0.709804, -180.00000000000094, 0.0509804, 0.396078, 0.741176, -165.00000000000082, 0.0392157, 0.466667, 0.768627, -150.0000000000007, 0.0313725, 0.537255, 0.788235, -134.3444452466246, 0.0313725, 0.615686, 0.811765, -118.30333220772445, 0.0235294, 0.709804, 0.831373, -102.26218409125124, 0.0509804, 0.8, 0.85098, -88.99742563467572, 0.0705882, 0.854902, 0.870588, -76.65808362797307, 0.262745, 0.901961, 0.862745, -65.8611725261965, 0.423529, 0.941176, 0.87451, -49.203080109808724, 0.572549, 0.964706, 0.835294, -38.09766178050461, 0.658824, 0.980392, 0.843137, -30.0771052610545, 0.764706, 0.980392, 0.866667, -21.43958690290333, 0.827451, 0.980392, 0.886275, -4.473004797774422, 0.913725, 0.988235, 0.937255, 0.7712146781356068, 1.0, 1.0, 0.972549019607843, 6.015434154045749, 0.988235, 0.980392, 0.870588, 13.727518523537583, 0.992156862745098, 0.972549019607843, 0.803921568627451, 19.588699838148784, 0.992157, 0.964706, 0.713725, 29.460159412483222, 0.988235, 0.956863, 0.643137, 45.192800301426985, 0.980392, 0.917647, 0.509804, 58.457593835573675, 0.968627, 0.87451, 0.407843, 72.03085951967677, 0.94902, 0.823529, 0.321569, 81.90231909401075, 0.929412, 0.776471, 0.278431, 96.4010363055558, 0.909804, 0.717647, 0.235294, 109.35732261217419, 0.890196, 0.658824, 0.196078, 119.99999999999989, 0.878431, 0.619608, 0.168627, 134.99999999999983, 0.870588, 0.54902, 0.156863, 149.99999999999994, 0.85098, 0.47451, 0.145098, 164.9999999999999, 0.831373, 0.411765, 0.133333, 180.0, 0.811765, 0.345098, 0.113725, 195.0, 0.788235, 0.266667, 0.0941176, 210.0, 0.741176, 0.184314, 0.0745098, 225.0, 0.690196, 0.12549, 0.0627451, 240.0, 0.619608, 0.0627451, 0.0431373, 254.03599078861942, 0.54902, 0.027451, 0.0705882, 266.37532238166943, 0.470588, 0.0156863, 0.0901961, 280.25706953320935, 0.4, 0.00392157, 0.101961, 300.0, 0.188235294117647, 0.0, 0.0705882352941176]
coutour_12LUT.ColorSpace = 'Lab'
coutour_12LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice4Display.Representation = 'Surface'
if cell_data:
    slice4Display.ColorArrayName = ['CELLS', 'coutour_12']
else:
    slice4Display.ColorArrayName = ['POINTS', 'coutour_12']
slice4Display.LookupTable = coutour_12LUT
slice4Display.Ambient = 1.0
slice4Display.Diffuse = 0.0
slice4Display.OSPRayScaleArray = 'coutour_12'
slice4Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice4Display.SelectOrientationVectors = 'None'
slice4Display.ScaleFactor = 0.27580108642578127
slice4Display.SelectScaleArray = 'coutour_12'
slice4Display.GlyphType = 'Arrow'
slice4Display.GlyphTableIndexArray = 'coutour_12'
slice4Display.GaussianRadius = 0.013790054321289063
slice4Display.SetScaleArray = ['POINTS', 'coutour_12']
slice4Display.ScaleTransferFunction = 'PiecewiseFunction'
slice4Display.OpacityArray = ['POINTS', 'coutour_12']
slice4Display.OpacityTransferFunction = 'PiecewiseFunction'
slice4Display.DataAxesGrid = 'GridAxesRepresentation'
slice4Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice4Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice4Display.ScaleTransferFunction.Points = [-165.53967271753004, 0.2589285671710968, 0.5, 0.0, 48.76916323898445, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice4Display.OpacityTransferFunction.Points = [-165.53967271753004, 0.2589285671710968, 0.5, 0.0, 48.76916323898445, 1.0, 0.5, 0.0]

# show data from slice6
slice6Display = Show(slice6, renderView2, 'GeometryRepresentation')

# get color transfer function/color map for 'countour_9'
countour_9LUT = GetColorTransferFunction('countour_9')
countour_9LUT.AutomaticRescaleRangeMode = 'Never'
countour_9LUT.RGBPoints = [-300.0, 0.0862745098039216, 0.00392156862745098, 0.298039215686275, -281.79946471687396, 0.113725, 0.0235294, 0.45098, -266.68380320541416, 0.105882, 0.0509804, 0.509804, -256.1953642535958, 0.0392157, 0.0392157, 0.560784, -246.015397451732, 0.0313725, 0.0980392, 0.6, -236.1439378773976, 0.0431373, 0.164706, 0.639216, -221.95372789337995, 0.054902, 0.243137, 0.678431, -203.1362251945821, 0.054902, 0.317647, 0.709804, -180.00000000000094, 0.0509804, 0.396078, 0.741176, -165.00000000000082, 0.0392157, 0.466667, 0.768627, -150.0000000000007, 0.0313725, 0.537255, 0.788235, -134.3444452466246, 0.0313725, 0.615686, 0.811765, -118.30333220772445, 0.0235294, 0.709804, 0.831373, -102.26218409125124, 0.0509804, 0.8, 0.85098, -88.99742563467572, 0.0705882, 0.854902, 0.870588, -76.65808362797307, 0.262745, 0.901961, 0.862745, -65.8611725261965, 0.423529, 0.941176, 0.87451, -49.203080109808724, 0.572549, 0.964706, 0.835294, -38.09766178050461, 0.658824, 0.980392, 0.843137, -30.0771052610545, 0.764706, 0.980392, 0.866667, -21.43958690290333, 0.827451, 0.980392, 0.886275, -4.473004797774422, 0.913725, 0.988235, 0.937255, 0.7712146781356068, 1.0, 1.0, 0.972549019607843, 6.015434154045749, 0.988235, 0.980392, 0.870588, 13.727518523537583, 0.992156862745098, 0.972549019607843, 0.803921568627451, 19.588699838148784, 0.992157, 0.964706, 0.713725, 29.460159412483222, 0.988235, 0.956863, 0.643137, 45.192800301426985, 0.980392, 0.917647, 0.509804, 58.457593835573675, 0.968627, 0.87451, 0.407843, 72.03085951967677, 0.94902, 0.823529, 0.321569, 81.90231909401075, 0.929412, 0.776471, 0.278431, 96.4010363055558, 0.909804, 0.717647, 0.235294, 109.35732261217419, 0.890196, 0.658824, 0.196078, 119.99999999999989, 0.878431, 0.619608, 0.168627, 134.99999999999983, 0.870588, 0.54902, 0.156863, 149.99999999999994, 0.85098, 0.47451, 0.145098, 164.9999999999999, 0.831373, 0.411765, 0.133333, 180.0, 0.811765, 0.345098, 0.113725, 195.0, 0.788235, 0.266667, 0.0941176, 210.0, 0.741176, 0.184314, 0.0745098, 225.0, 0.690196, 0.12549, 0.0627451, 240.0, 0.619608, 0.0627451, 0.0431373, 254.03599078861942, 0.54902, 0.027451, 0.0705882, 266.37532238166943, 0.470588, 0.0156863, 0.0901961, 280.25706953320935, 0.4, 0.00392157, 0.101961, 300.0, 0.188235294117647, 0.0, 0.0705882352941176]
countour_9LUT.ColorSpace = 'Lab'
countour_9LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice6Display.Representation = 'Surface'
if cell_data:
    slice6Display.ColorArrayName = ['CELLS', 'countour_9']
else: 
    slice6Display.ColorArrayName = ['POINTS', 'countour_9']
slice6Display.LookupTable = countour_9LUT
slice6Display.Ambient = 1.0
slice6Display.Diffuse = 0.0
slice6Display.OSPRayScaleArray = 'countour_9'
slice6Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice6Display.SelectOrientationVectors = 'None'
slice6Display.ScaleFactor = 0.28134746551513673
slice6Display.SelectScaleArray = 'countour_9'
slice6Display.GlyphType = 'Arrow'
slice6Display.GlyphTableIndexArray = 'countour_9'
slice6Display.GaussianRadius = 0.014067373275756837
slice6Display.SetScaleArray = ['POINTS', 'countour_9']
slice6Display.ScaleTransferFunction = 'PiecewiseFunction'
slice6Display.OpacityArray = ['POINTS', 'countour_9']
slice6Display.OpacityTransferFunction = 'PiecewiseFunction'
slice6Display.DataAxesGrid = 'GridAxesRepresentation'
slice6Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice6Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice6Display.ScaleTransferFunction.Points = [-8.056463866290857, 0.2589285671710968, 0.5, 0.0, 113.75935885662665, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice6Display.OpacityTransferFunction.Points = [-8.056463866290857, 0.2589285671710968, 0.5, 0.0, 113.75935885662665, 1.0, 0.5, 0.0]

# show data from slice7
slice7Display = Show(slice7, renderView2, 'GeometryRepresentation')

# get color transfer function/color map for 'contour_6'
contour_6LUT = GetColorTransferFunction('contour_6')
contour_6LUT.AutomaticRescaleRangeMode = 'Never'
contour_6LUT.RGBPoints = [-300.0, 0.0862745098039216, 0.00392156862745098, 0.298039215686275, -281.79946471687396, 0.113725, 0.0235294, 0.45098, -266.68380320541416, 0.105882, 0.0509804, 0.509804, -256.1953642535958, 0.0392157, 0.0392157, 0.560784, -246.015397451732, 0.0313725, 0.0980392, 0.6, -236.1439378773976, 0.0431373, 0.164706, 0.639216, -221.95372789337995, 0.054902, 0.243137, 0.678431, -203.1362251945821, 0.054902, 0.317647, 0.709804, -180.00000000000094, 0.0509804, 0.396078, 0.741176, -165.00000000000082, 0.0392157, 0.466667, 0.768627, -150.0000000000007, 0.0313725, 0.537255, 0.788235, -134.3444452466246, 0.0313725, 0.615686, 0.811765, -118.30333220772445, 0.0235294, 0.709804, 0.831373, -102.26218409125124, 0.0509804, 0.8, 0.85098, -88.99742563467572, 0.0705882, 0.854902, 0.870588, -76.65808362797307, 0.262745, 0.901961, 0.862745, -65.8611725261965, 0.423529, 0.941176, 0.87451, -49.203080109808724, 0.572549, 0.964706, 0.835294, -38.09766178050461, 0.658824, 0.980392, 0.843137, -30.0771052610545, 0.764706, 0.980392, 0.866667, -21.43958690290333, 0.827451, 0.980392, 0.886275, -4.473004797774422, 0.913725, 0.988235, 0.937255, 0.7712146781356068, 1.0, 1.0, 0.972549019607843, 6.015434154045749, 0.988235, 0.980392, 0.870588, 13.727518523537583, 0.992156862745098, 0.972549019607843, 0.803921568627451, 19.588699838148784, 0.992157, 0.964706, 0.713725, 29.460159412483222, 0.988235, 0.956863, 0.643137, 45.192800301426985, 0.980392, 0.917647, 0.509804, 58.457593835573675, 0.968627, 0.87451, 0.407843, 72.03085951967677, 0.94902, 0.823529, 0.321569, 81.90231909401075, 0.929412, 0.776471, 0.278431, 96.4010363055558, 0.909804, 0.717647, 0.235294, 109.35732261217419, 0.890196, 0.658824, 0.196078, 119.99999999999989, 0.878431, 0.619608, 0.168627, 134.99999999999983, 0.870588, 0.54902, 0.156863, 149.99999999999994, 0.85098, 0.47451, 0.145098, 164.9999999999999, 0.831373, 0.411765, 0.133333, 180.0, 0.811765, 0.345098, 0.113725, 195.0, 0.788235, 0.266667, 0.0941176, 210.0, 0.741176, 0.184314, 0.0745098, 225.0, 0.690196, 0.12549, 0.0627451, 240.0, 0.619608, 0.0627451, 0.0431373, 254.03599078861942, 0.54902, 0.027451, 0.0705882, 266.37532238166943, 0.470588, 0.0156863, 0.0901961, 280.25706953320935, 0.4, 0.00392157, 0.101961, 300.0, 0.188235294117647, 0.0, 0.0705882352941176]
contour_6LUT.ColorSpace = 'Lab'
contour_6LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice7Display.Representation = 'Surface'
if cell_data:
    slice7Display.ColorArrayName = ['CELLS', 'contour_6']
else:
    slice7Display.ColorArrayName = ['POINTS', 'contour_6']
slice7Display.LookupTable = contour_6LUT
slice7Display.Ambient = 1.0
slice7Display.Diffuse = 0.0
slice7Display.OSPRayScaleArray = 'contour_6'
slice7Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice7Display.SelectOrientationVectors = 'None'
slice7Display.ScaleFactor = 0.28354568481445314
slice7Display.SelectScaleArray = 'contour_6'
slice7Display.GlyphType = 'Arrow'
slice7Display.GlyphTableIndexArray = 'contour_6'
slice7Display.GaussianRadius = 0.014177284240722657
slice7Display.SetScaleArray = ['POINTS', 'contour_6']
slice7Display.ScaleTransferFunction = 'PiecewiseFunction'
slice7Display.OpacityArray = ['POINTS', 'contour_6']
slice7Display.OpacityTransferFunction = 'PiecewiseFunction'
slice7Display.DataAxesGrid = 'GridAxesRepresentation'
slice7Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice7Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice7Display.ScaleTransferFunction.Points = [-11.532170355951907, 0.2589285671710968, 0.5, 0.0, 95.53556180917279, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice7Display.OpacityTransferFunction.Points = [-11.532170355951907, 0.2589285671710968, 0.5, 0.0, 95.53556180917279, 1.0, 0.5, 0.0]

# show data from slice8
slice8Display = Show(slice8, renderView2, 'GeometryRepresentation')

# get color transfer function/color map for 'Result'
resultLUT = GetColorTransferFunction('Result')
resultLUT.AutomaticRescaleRangeMode = 'Never'
resultLUT.RGBPoints = [-300.0, 0.0862745098039216, 0.00392156862745098, 0.298039215686275, -281.79946471687396, 0.113725, 0.0235294, 0.45098, -266.68380320541416, 0.105882, 0.0509804, 0.509804, -256.1953642535958, 0.0392157, 0.0392157, 0.560784, -246.015397451732, 0.0313725, 0.0980392, 0.6, -236.1439378773976, 0.0431373, 0.164706, 0.639216, -221.95372789337995, 0.054902, 0.243137, 0.678431, -203.1362251945821, 0.054902, 0.317647, 0.709804, -180.00000000000094, 0.0509804, 0.396078, 0.741176, -165.00000000000082, 0.0392157, 0.466667, 0.768627, -150.0000000000007, 0.0313725, 0.537255, 0.788235, -134.3444452466246, 0.0313725, 0.615686, 0.811765, -118.30333220772445, 0.0235294, 0.709804, 0.831373, -102.26218409125124, 0.0509804, 0.8, 0.85098, -88.99742563467572, 0.0705882, 0.854902, 0.870588, -76.65808362797307, 0.262745, 0.901961, 0.862745, -65.8611725261965, 0.423529, 0.941176, 0.87451, -49.203080109808724, 0.572549, 0.964706, 0.835294, -38.09766178050461, 0.658824, 0.980392, 0.843137, -30.0771052610545, 0.764706, 0.980392, 0.866667, -21.43958690290333, 0.827451, 0.980392, 0.886275, -4.473004797774422, 0.913725, 0.988235, 0.937255, 0.7712146781356068, 1.0, 1.0, 0.972549019607843, 6.015434154045749, 0.988235, 0.980392, 0.870588, 13.727518523537583, 0.992156862745098, 0.972549019607843, 0.803921568627451, 19.588699838148784, 0.992157, 0.964706, 0.713725, 29.460159412483222, 0.988235, 0.956863, 0.643137, 45.192800301426985, 0.980392, 0.917647, 0.509804, 58.457593835573675, 0.968627, 0.87451, 0.407843, 72.03085951967677, 0.94902, 0.823529, 0.321569, 81.90231909401075, 0.929412, 0.776471, 0.278431, 96.4010363055558, 0.909804, 0.717647, 0.235294, 109.35732261217419, 0.890196, 0.658824, 0.196078, 119.99999999999989, 0.878431, 0.619608, 0.168627, 134.99999999999983, 0.870588, 0.54902, 0.156863, 149.99999999999994, 0.85098, 0.47451, 0.145098, 164.9999999999999, 0.831373, 0.411765, 0.133333, 180.0, 0.811765, 0.345098, 0.113725, 195.0, 0.788235, 0.266667, 0.0941176, 210.0, 0.741176, 0.184314, 0.0745098, 225.0, 0.690196, 0.12549, 0.0627451, 240.0, 0.619608, 0.0627451, 0.0431373, 254.03599078861942, 0.54902, 0.027451, 0.0705882, 266.37532238166943, 0.470588, 0.0156863, 0.0901961, 280.25706953320935, 0.4, 0.00392157, 0.101961, 300.0, 0.188235294117647, 0.0, 0.0705882352941176]
resultLUT.ColorSpace = 'Lab'
resultLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice8Display.Representation = 'Surface'
if cell_data:
    slice8Display.ColorArrayName = ['CELLS', 'Result']
else:
    slice8Display.ColorArrayName = ['POINTS', 'Result']
slice8Display.LookupTable = resultLUT
slice8Display.Ambient = 1.0
slice8Display.Diffuse = 0.0
slice8Display.OSPRayScaleArray = 'Result'
slice8Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice8Display.SelectOrientationVectors = 'None'
slice8Display.ScaleFactor = 0.28310604095458985
slice8Display.SelectScaleArray = 'Result'
slice8Display.GlyphType = 'Arrow'
slice8Display.GlyphTableIndexArray = 'Result'
slice8Display.GaussianRadius = 0.014155302047729492
slice8Display.SetScaleArray = ['POINTS', 'Result']
slice8Display.ScaleTransferFunction = 'PiecewiseFunction'
slice8Display.OpacityArray = ['POINTS', 'Result']
slice8Display.OpacityTransferFunction = 'PiecewiseFunction'
slice8Display.DataAxesGrid = 'GridAxesRepresentation'
slice8Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice8Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice8Display.ScaleTransferFunction.Points = [-7.05830921300197, 0.2589285671710968, 0.5, 0.0, 85.67541241027644, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice8Display.OpacityTransferFunction.Points = [-7.05830921300197, 0.2589285671710968, 0.5, 0.0, 85.67541241027644, 1.0, 0.5, 0.0]

# show data from text2
# text2Display = Show(text2, renderView2, 'TextSourceRepresentation')

# trace defaults for the display properties.
# text2Display.Color = [0.0, 0.0, 0.0]
# text2Display.FontFamily = 'Times'
# text2Display.FontSize = 9
# text2Display.WindowLocation = 'AnyLocation'
# text2Display.Position = [0.02, 0.02]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'annulus_normal_projected'
annulus_normal_projectedPWF = GetOpacityTransferFunction('annulus_normal_projected')
annulus_normal_projectedPWF.Points = [-300.0, 0.2589285671710968, 0.5, 0.0, 300.0, 1.0, 0.5, 0.0]
annulus_normal_projectedPWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'coutour_12'
coutour_12PWF = GetOpacityTransferFunction('coutour_12')
coutour_12PWF.Points = [-300.0, 0.2589285671710968, 0.5, 0.0, 300.0, 1.0, 0.5, 0.0]
coutour_12PWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'U'
uPWF = GetOpacityTransferFunction('U')
uPWF.Points = [-300.0, 0.2589285671710968, 0.5, 0.0, 300.0, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'contour_6'
contour_6PWF = GetOpacityTransferFunction('contour_6')
contour_6PWF.Points = [-300.0, 0.2589285671710968, 0.5, 0.0, 300.0, 1.0, 0.5, 0.0]
contour_6PWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'Result'
resultPWF = GetOpacityTransferFunction('Result')
resultPWF.Points = [-300.0, 0.2589285671710968, 0.5, 0.0, 300.0, 1.0, 0.5, 0.0]
resultPWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'countour_9'
countour_9PWF = GetOpacityTransferFunction('countour_9')
countour_9PWF.Points = [-300.0, 0.2589285671710968, 0.5, 0.0, 300.0, 1.0, 0.5, 0.0]
countour_9PWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(slice1)
# ----------------------------------------------------------------


tk = GetTimeKeeper()
timesteps = tk.TimestepValues
numTimesteps = len(timesteps)

print("timesteps = ", timesteps)
print("numTimesteps = ", numTimesteps)

if len(sys.argv) >= 2:
    basename = sys.argv[1]
else: 
    basename = 'frames'    
    if not cell_data:
        basename += '_points'
    basename += '_paper'

if 'paper' not in basename:
    basename += '_paper'

# take these from passed in 
# if not cell_data:
#     basename += '_points'

# basename += '_paper'

if len(sys.argv) >= 4:
    nprocs = int(sys.argv[2])
    proc_num = int(sys.argv[3])
else: 
    #print "using default proc_num 0, nprocs = 1"
    proc_num = 0
    nprocs = 1

scaling = 2 
res = (scaling*1280, scaling*1080)

frame_range = [1309,1341,1365,1427]
# frame_range = [1365]
# frame_range = [350]


animationScene1 = GetAnimationScene()
# animationScene1.GoToFirst()

# # run first iteration, then stop and start over 
for frame in range(len(timesteps)):
    if (frame % nprocs) == proc_num:
        animationScene1.AnimationTime = timesteps[frame]

        Render()
        SaveScreenshot(basename + str(frame).zfill(4) + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)
        break 

for frame in range(len(timesteps)-1):
    if (frame % nprocs) == proc_num:
        animationScene1.AnimationTime = timesteps[frame]

        Render()
        SaveScreenshot(basename + str(frame).zfill(4) + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)


# serial version just on frame_range

# for frame in frame_range: # range(len(timesteps)-1):
    
#     animationScene1.AnimationTime = timesteps[frame]

#     print("frame = ", frame)
#     print("animationScene1.AnimationTime = ", animationScene1.AnimationTime)
#     Render()
#     SaveScreenshot(basename + str(frame).zfill(4) + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)
