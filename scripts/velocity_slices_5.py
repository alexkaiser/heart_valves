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
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [640, 1080]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.0, 0.0, 0.75]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-0.0118654110947045, 0.06453941407452532, 25.008722343099873]
renderView1.CameraFocalPoint = [-0.0118654110947045, 0.06453941407452532, 0.75]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 2.1650035488494987
renderView1.CameraParallelProjection = 1
renderView1.Background = [1.0, 0.9999694819562066, 0.9999847409781033]
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [640, 1080]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [0.0, -0.00016164779663085938, 0.5]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [0.0, -18.143050647253336, 0.5]
renderView2.CameraFocalPoint = [0.0, 17.856788040611463, 0.5]
renderView2.CameraViewUp = [0.0, 0.0, 1.0]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 9.532591314548378
renderView2.Background = [1.0, 0.9999694819562066, 0.9999847409781033]
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.5)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, renderView2)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
# pVDReader1 = PVDReader(FileName='/Users/alex/data_to_remove/aortic_65595790_384_4495be5_circ_pt15_rad_pt54_2mm_radial_4mm_circ_circ_model_basic_updated_output_semifinal/exported_viz/particles_vel.pvd')
# pVDReader1.PointArrays = ['velocity']

# create a new 'PVD Reader'
pVDReader2 = PVDReader(FileName='aortic_384_cylinder.pvd')

# create a new 'PVD Reader'
eulerian_varspvd = PVDReader(FileName='eulerian_vars.pvd')
eulerian_varspvd.CellArrays = ['U']

# create a new 'Annotate Time'
annotateTime1 = AnnotateTime()
annotateTime1.Format = 't = %.4f s'

# create a new 'PVD Reader'
eulerian_varspvd_1 = PVDReader(FileName='eulerian_vars.pvd')
eulerian_varspvd_1.CellArrays = ['U']

# create a new 'Temporal Particles To Pathlines'
# temporalParticlesToPathlines1 = TemporalParticlesToPathlines(Input=pVDReader1,
#     Selection=None)
# temporalParticlesToPathlines1.MaskPoints = 1

# create a new 'Slice'
slice1 = Slice(Input=eulerian_varspvd_1)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0.0, 0.0, 0.5]

# create a new 'PVD Reader'
pVDReader3 = PVDReader(FileName='aortic_384.pvd')

# create a new 'Clip'
clip1 = Clip(Input=pVDReader3)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', '']

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [-1.3, -0.00016176700592041016, 0.8519301824271679]
clip1.ClipType.Normal = [-1.0, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [-0.00016355514526367188, -0.00016176700592041016, 0.8519301824271679]

# create a new 'Clip'
clip2 = Clip(Input=clip1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', '']

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Origin = [1.3, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [0.46972882747650146, -0.00016176700592041016, 0.8519301824271679]

# create a new 'Clip'
clip3 = Clip(Input=clip2)
clip3.ClipType = 'Plane'
clip3.HyperTreeGridClipper = 'Plane'
clip3.Scalars = ['POINTS', '']

# init the 'Plane' selected for 'ClipType'
clip3.ClipType.Origin = [0.0, -1.3, 0.0]
clip3.ClipType.Normal = [0.0, -1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip3.HyperTreeGridClipper.Origin = [0.0, -0.00016176700592041016, 0.8519301824271679]

# create a new 'Clip'
clip4 = Clip(Input=clip3)
clip4.ClipType = 'Plane'
clip4.HyperTreeGridClipper = 'Plane'
clip4.Scalars = ['POINTS', '']

# init the 'Plane' selected for 'ClipType'
clip4.ClipType.Origin = [0.0, 1.3, 0.0]
clip4.ClipType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip4.HyperTreeGridClipper.Origin = [0.0, 0.46972858905792236, 0.8519301824271679]

# create a new 'Slice'
slice2 = Slice(Input=eulerian_varspvd)
slice2.SliceType = 'Plane'
slice2.HyperTreeGridSlicer = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice2.HyperTreeGridSlicer.Origin = [0.0, 0.0, 0.5]

# create a new 'Clip'
clip5 = Clip(Input=slice2)
clip5.ClipType = 'Plane'
clip5.HyperTreeGridClipper = 'Plane'
clip5.Scalars = ['POINTS', '']

# init the 'Plane' selected for 'ClipType'
clip5.ClipType.Origin = [-1.3, 0.0, 0.0]
clip5.ClipType.Normal = [-1.0, 0.0, 0.0]

# create a new 'Clip'
clip6 = Clip(Input=clip5)
clip6.ClipType = 'Plane'
clip6.HyperTreeGridClipper = 'Plane'
clip6.Scalars = ['POINTS', '']

# init the 'Plane' selected for 'ClipType'
clip6.ClipType.Origin = [1.3, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip6.HyperTreeGridClipper.Origin = [0.4750000238418579, 0.0, 0.0]

# create a new 'Clip'
clip7 = Clip(Input=clip6)
clip7.ClipType = 'Plane'
clip7.HyperTreeGridClipper = 'Plane'
clip7.Scalars = ['POINTS', '']

# init the 'Plane' selected for 'ClipType'
clip7.ClipType.Origin = [0.0, -1.3, 0.0]
clip7.ClipType.Normal = [0.0, -1.0, 0.0]

# create a new 'Clip'
clip8 = Clip(Input=clip7)
clip8.ClipType = 'Plane'
clip8.HyperTreeGridClipper = 'Plane'
clip8.Scalars = ['POINTS', '']

# init the 'Plane' selected for 'ClipType'
clip8.ClipType.Origin = [0.0, 1.3, 0.0]
clip8.ClipType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip8.HyperTreeGridClipper.Origin = [0.0, 0.4750000238418579, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from pVDReader2
pVDReader2Display = Show(pVDReader2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
pVDReader2Display.Representation = 'Surface'
pVDReader2Display.AmbientColor = [0.0, 0.0, 0.0]
pVDReader2Display.ColorArrayName = ['POINTS', '']
pVDReader2Display.DiffuseColor = [0.0, 0.0, 0.0]
pVDReader2Display.OSPRayScaleFunction = 'PiecewiseFunction'
pVDReader2Display.SelectOrientationVectors = 'None'
pVDReader2Display.ScaleFactor = 0.25863996744155887
pVDReader2Display.SelectScaleArray = 'None'
pVDReader2Display.GlyphType = 'Arrow'
pVDReader2Display.GlyphTableIndexArray = 'None'
pVDReader2Display.GaussianRadius = 0.012931998372077941
pVDReader2Display.SetScaleArray = ['POINTS', '']
pVDReader2Display.ScaleTransferFunction = 'PiecewiseFunction'
pVDReader2Display.OpacityArray = ['POINTS', '']
pVDReader2Display.OpacityTransferFunction = 'PiecewiseFunction'
pVDReader2Display.DataAxesGrid = 'GridAxesRepresentation'
pVDReader2Display.PolarAxes = 'PolarAxesRepresentation'
pVDReader2Display.ScalarOpacityUnitDistance = 0.06717646765003486

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
pVDReader2Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
pVDReader2Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
pVDReader2Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from annotateTime1
annotateTime1Display = Show(annotateTime1, renderView1, 'TextSourceRepresentation')

# trace defaults for the display properties.
annotateTime1Display.Color = [0.0, 0.0, 0.0]
annotateTime1Display.FontFamily = 'Times'
annotateTime1Display.FontSize = 10
annotateTime1Display.WindowLocation = 'AnyLocation'
# annotateTime1Display.Position = [0.05556511864406785, 0.8133102846441947]

annotateTime1Display.Position = [0.64, 0.79]

# show data from clip4
clip4Display = Show(clip4, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip4Display.Representation = 'Surface'
clip4Display.AmbientColor = [0.0, 0.0, 0.0]
clip4Display.ColorArrayName = ['POINTS', '']
clip4Display.DiffuseColor = [0.0, 0.0, 0.0]
clip4Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip4Display.SelectOrientationVectors = 'None'
clip4Display.ScaleFactor = 0.25999999046325684
clip4Display.SelectScaleArray = 'None'
clip4Display.GlyphType = 'Arrow'
clip4Display.GlyphTableIndexArray = 'None'
clip4Display.GaussianRadius = 0.012999999523162843
clip4Display.SetScaleArray = ['POINTS', '']
clip4Display.ScaleTransferFunction = 'PiecewiseFunction'
clip4Display.OpacityArray = ['POINTS', '']
clip4Display.OpacityTransferFunction = 'PiecewiseFunction'
clip4Display.DataAxesGrid = 'GridAxesRepresentation'
clip4Display.PolarAxes = 'PolarAxesRepresentation'
clip4Display.ScalarOpacityUnitDistance = 0.06449917355947005

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip4Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip4Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip4Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from clip8
clip8Display = Show(clip8, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'U'
uLUT = GetColorTransferFunction('U')
uLUT.AutomaticRescaleRangeMode = 'Never'
uLUT.RGBPoints = [-200.0, 0.0862745098039216, 0.00392156862745098, 0.298039215686275, -187.86630981124935, 0.113725, 0.0235294, 0.45098, -177.78920213694278, 0.105882, 0.0509804, 0.509804, -170.7969095023972, 0.0392157, 0.0392157, 0.560784, -164.01026496782134, 0.0313725, 0.0980392, 0.6, -157.42929191826505, 0.0431373, 0.164706, 0.639216, -147.96915192891996, 0.054902, 0.243137, 0.678431, -135.4241501297214, 0.054902, 0.317647, 0.709804, -120.00000000000063, 0.0509804, 0.396078, 0.741176, -110.00000000000055, 0.0392157, 0.466667, 0.768627, -100.00000000000048, 0.0313725, 0.537255, 0.788235, -89.56296349774972, 0.0313725, 0.615686, 0.811765, -78.86888813848296, 0.0235294, 0.709804, 0.831373, -68.1747893941675, 0.0509804, 0.8, 0.85098, -59.331617089783805, 0.0705882, 0.854902, 0.870588, -51.10538908531538, 0.262745, 0.901961, 0.862745, -43.907448350797665, 0.423529, 0.941176, 0.87451, -32.80205340653913, 0.572549, 0.964706, 0.835294, -25.398441187003073, 0.658824, 0.980392, 0.843137, -20.051403507369656, 0.764706, 0.980392, 0.866667, -14.293057935268905, 0.827451, 0.980392, 0.886275, -2.9820031985162814, 0.913725, 0.988235, 0.937255, 0.5141431187570618, 1.0, 1.0, 0.972549019607843, 4.01028943603049, 0.988235, 0.980392, 0.870588, 9.151679015691712, 0.992156862745098, 0.972549019607843, 0.803921568627451, 13.05913322543256, 0.992157, 0.964706, 0.713725, 19.640106274988824, 0.988235, 0.956863, 0.643137, 30.128533534284628, 0.980392, 0.917647, 0.509804, 38.971729223715755, 0.968627, 0.87451, 0.407843, 48.020573013117826, 0.94902, 0.823529, 0.321569, 54.60154606267383, 0.929412, 0.776471, 0.278431, 64.26735753703719, 0.909804, 0.717647, 0.235294, 72.90488174144946, 0.890196, 0.658824, 0.196078, 79.99999999999994, 0.878431, 0.619608, 0.168627, 89.99999999999989, 0.870588, 0.54902, 0.156863, 99.99999999999994, 0.85098, 0.47451, 0.145098, 109.99999999999994, 0.831373, 0.411765, 0.133333, 120.0, 0.811765, 0.345098, 0.113725, 130.0, 0.788235, 0.266667, 0.0941176, 140.0, 0.741176, 0.184314, 0.0745098, 150.0, 0.690196, 0.12549, 0.0627451, 160.0, 0.619608, 0.0627451, 0.0431373, 169.35732719241292, 0.54902, 0.027451, 0.0705882, 177.58354825444627, 0.470588, 0.0156863, 0.0901961, 186.8380463554729, 0.4, 0.00392157, 0.101961, 200.0, 0.188235294117647, 0.0, 0.0705882352941176]
uLUT.ColorSpace = 'Lab'
uLUT.ScalarRangeInitialized = 1.0
uLUT.VectorComponent = 2
uLUT.VectorMode = 'Component'

# get opacity transfer function/opacity map for 'U'
uPWF = GetOpacityTransferFunction('U')
uPWF.Points = [-200.0, 0.2589285671710968, 0.5, 0.0, 200.0, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
clip8Display.Representation = 'Surface'
clip8Display.ColorArrayName = ['CELLS', 'U']
clip8Display.LookupTable = uLUT
clip8Display.SpecularPower = 0.0
clip8Display.Ambient = 1.0
clip8Display.Diffuse = 0.0
clip8Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip8Display.SelectOrientationVectors = 'None'
clip8Display.ScaleFactor = 0.25999999046325684
clip8Display.SelectScaleArray = 'None'
clip8Display.GlyphType = 'Arrow'
clip8Display.GlyphTableIndexArray = 'None'
clip8Display.GaussianRadius = 0.012999999523162843
clip8Display.SetScaleArray = [None, '']
clip8Display.ScaleTransferFunction = 'PiecewiseFunction'
clip8Display.OpacityArray = [None, '']
clip8Display.OpacityTransferFunction = 'PiecewiseFunction'
clip8Display.DataAxesGrid = 'GridAxesRepresentation'
clip8Display.PolarAxes = 'PolarAxesRepresentation'
clip8Display.ScalarOpacityFunction = uPWF
clip8Display.ScalarOpacityUnitDistance = 0.199360979343879
clip8Display.ExtractedBlockIndex = 2

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip8Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip8Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip8Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)
uLUTColorBar.AutoOrient = 0
uLUTColorBar.WindowLocation = 'AnyLocation'
# uLUTColorBar.Position = [0.719774011299435, 0.8058614232209736]
uLUTColorBar.Position = [0.05, 0.8058614232209736]
uLUTColorBar.Title = 'u_z (cm/s)'
uLUTColorBar.ComponentTitle = ''
uLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
uLUTColorBar.TitleFontFamily = 'Times'
uLUTColorBar.TitleFontSize = 24
uLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
uLUTColorBar.LabelFontFamily = 'Times'
uLUTColorBar.LabelFontSize = 20
uLUTColorBar.LabelFormat = '%-#1.0f'
uLUTColorBar.UseCustomLabels = 1
uLUTColorBar.CustomLabels = [-200.0, -100.0, 0.0, 100.0, 200.0]
uLUTColorBar.AddRangeLabels = 0
uLUTColorBar.DrawAnnotations = 0
uLUTColorBar.ScalarBarThickness = 20
uLUTColorBar.ScalarBarLength = 0.17000000000000004

# set color bar visibility
uLUTColorBar.Visibility = 1

# show color legend
clip8Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from pVDReader3
pVDReader3Display = Show(pVDReader3, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
pVDReader3Display.Representation = 'Surface'
pVDReader3Display.AmbientColor = [0.0, 0.0, 0.0]
pVDReader3Display.ColorArrayName = ['POINTS', '']
pVDReader3Display.DiffuseColor = [0.0, 0.0, 0.0]
pVDReader3Display.Opacity = 0.2
pVDReader3Display.OSPRayScaleFunction = 'PiecewiseFunction'
pVDReader3Display.SelectOrientationVectors = 'None'
pVDReader3Display.ScaleFactor = 0.4479223966598511
pVDReader3Display.SelectScaleArray = 'None'
pVDReader3Display.GlyphType = 'Arrow'
pVDReader3Display.GlyphTableIndexArray = 'None'
pVDReader3Display.GaussianRadius = 0.022396119832992552
pVDReader3Display.SetScaleArray = ['POINTS', '']
pVDReader3Display.ScaleTransferFunction = 'PiecewiseFunction'
pVDReader3Display.OpacityArray = ['POINTS', '']
pVDReader3Display.OpacityTransferFunction = 'PiecewiseFunction'
pVDReader3Display.DataAxesGrid = 'GridAxesRepresentation'
pVDReader3Display.PolarAxes = 'PolarAxesRepresentation'
pVDReader3Display.ScalarOpacityUnitDistance = 0.08164692686596023

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
pVDReader3Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
pVDReader3Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
pVDReader3Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from pVDReader2
pVDReader2Display_1 = Show(pVDReader2, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
pVDReader2Display_1.Representation = 'Surface'
pVDReader2Display_1.AmbientColor = [0.0, 0.0, 0.0]
pVDReader2Display_1.ColorArrayName = ['POINTS', '']
pVDReader2Display_1.DiffuseColor = [0.0, 0.0, 0.0]
pVDReader2Display_1.Opacity = 0.04
pVDReader2Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
pVDReader2Display_1.SelectOrientationVectors = 'None'
pVDReader2Display_1.ScaleFactor = 0.2581812381744385
pVDReader2Display_1.SelectScaleArray = 'None'
pVDReader2Display_1.GlyphType = 'Arrow'
pVDReader2Display_1.GlyphTableIndexArray = 'None'
pVDReader2Display_1.GaussianRadius = 0.012909061908721924
pVDReader2Display_1.SetScaleArray = ['POINTS', '']
pVDReader2Display_1.ScaleTransferFunction = 'PiecewiseFunction'
pVDReader2Display_1.OpacityArray = ['POINTS', '']
pVDReader2Display_1.OpacityTransferFunction = 'PiecewiseFunction'
pVDReader2Display_1.DataAxesGrid = 'GridAxesRepresentation'
pVDReader2Display_1.PolarAxes = 'PolarAxesRepresentation'
pVDReader2Display_1.ScalarOpacityUnitDistance = 0.06706963278681245

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
pVDReader2Display_1.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
pVDReader2Display_1.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
pVDReader2Display_1.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from slice1
slice1Display = Show(slice1, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'U']
slice1Display.LookupTable = uLUT
slice1Display.InterpolateScalarsBeforeMapping = 0
slice1Display.Interpolation = 'Flat'
slice1Display.SpecularColor = [1.0, 0.9999694819562066, 0.9999847409781033]
slice1Display.SpecularPower = 0.0
slice1Display.Ambient = 1.0
slice1Display.Diffuse = 0.0
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.9
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.045
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(clip8)
# ----------------------------------------------------------------


tk = GetTimeKeeper()
timesteps = tk.TimestepValues
numTimesteps = len(timesteps)

if len(sys.argv) >= 2:
    basename = sys.argv[1]
else: 
    basename = 'frames'    

if len(sys.argv) >= 4:
    nprocs = int(sys.argv[2])
    proc_num = int(sys.argv[3])
else: 
    print "using default proc_num 0, nprocs = 1"
    proc_num = 0
    nprocs = 1

res = (5120, 4320)

animationScene1 = GetAnimationScene()
# animationScene1.GoToFirst()

for frame in range(len(timesteps)):
    if (frame % nprocs) == proc_num:
        animationScene1.AnimationTime = timesteps[frame]

        Render()
        SaveScreenshot(basename + str(frame).zfill(4) + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)

