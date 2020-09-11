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
renderView1.ViewSize = [640*4, 1080*4] # was [640, 1080]
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
renderView2.ViewSize = [640*4, 1080*4] # [640, 1080]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [-0.0001615285873413086, -0.0001615285873413086, 0.75]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [-0.0001615285873413086, -18.953078232283797, 0.75]
renderView2.CameraFocalPoint = [-0.0001615285873413086, 9.865960517829109, 0.75]
renderView2.CameraViewUp = [0.0, 0.0, 1.0]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 7.591561086876913
renderView2.Background = [1.0, 0.9999694819562066, 0.9999847409781033]
renderView2.EnableRayTracing = 1
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
SetActiveView(renderView2)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Annotate Time'
annotateTime1 = AnnotateTime()
annotateTime1.Format = 't = %.4f s'

# create a new 'PVD Reader'
# particles_velpvd = PVDReader(FileName='/Users/alex/data_to_remove/aortic_65595790_384_4495be5_circ_pt15_rad_pt54_2mm_radial_4mm_circ_circ_model_basic_updated_output_semifinal/exported_viz/particles_vel.pvd')
particles_velpvd = PVDReader(FileName='aortic_384_particles_vel.pvd')
particles_velpvd.PointArrays = ['velocity']

# create a new 'Temporal Particles To Pathlines'
temporalParticlesToPathlines1 = TemporalParticlesToPathlines(Input=particles_velpvd,
    Selection=None)
temporalParticlesToPathlines1.MaskPoints = 1
temporalParticlesToPathlines1.MaxTrackLength = 25 

# create a new 'PVD Reader'
# cylinderpvd = PVDReader(FileName='/Users/alex/data_to_remove/aortic_65595790_384_4495be5_circ_pt15_rad_pt54_2mm_radial_4mm_circ_circ_model_basic_updated_output_semifinal/exported_viz/cylinder.pvd')
cylinderpvd = PVDReader(FileName='aortic_384_cylinder.pvd')

# create a new 'PVD Reader'
# aortic_valvepvd = PVDReader(FileName='/Users/alex/data_to_remove/aortic_65595790_384_4495be5_circ_pt15_rad_pt54_2mm_radial_4mm_circ_circ_model_basic_updated_output_semifinal/exported_viz/aortic_valve.pvd')
aortic_valvepvd = PVDReader(FileName='aortic_384.pvd')


# create a new 'Clip'
clip1 = Clip(Input=aortic_valvepvd)
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

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from temporalParticlesToPathlines1
temporalParticlesToPathlines1Display = Show(temporalParticlesToPathlines1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'velocity'
velocityLUT = GetColorTransferFunction('velocity')
velocityLUT.AutomaticRescaleRangeMode = 'Never'
velocityLUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 50.0, 0.901960784314, 0.901960784314, 0.0, 100.0, 0.901960784314, 0.0, 0.0, 150.0, 0.3333333333333333, 0.0, 0.0, 200.0, 0.0, 0.0, 0.0]
velocityLUT.ColorSpace = 'RGB'
velocityLUT.NanColor = [0.0, 0.498039215686, 1.0]
velocityLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
temporalParticlesToPathlines1Display.Representation = 'Surface'
temporalParticlesToPathlines1Display.ColorArrayName = ['POINTS', 'velocity']
temporalParticlesToPathlines1Display.LookupTable = velocityLUT
temporalParticlesToPathlines1Display.Opacity = 0.0
temporalParticlesToPathlines1Display.OSPRayScaleArray = 'TrailId'
temporalParticlesToPathlines1Display.OSPRayScaleFunction = 'PiecewiseFunction'
temporalParticlesToPathlines1Display.SelectOrientationVectors = 'None'
temporalParticlesToPathlines1Display.ScaleFactor = 0.675
temporalParticlesToPathlines1Display.SelectScaleArray = 'TrailId'
temporalParticlesToPathlines1Display.GlyphType = 'Arrow'
temporalParticlesToPathlines1Display.GlyphTableIndexArray = 'TrailId'
temporalParticlesToPathlines1Display.GaussianRadius = 0.03375
temporalParticlesToPathlines1Display.SetScaleArray = ['POINTS', 'TrailId']
temporalParticlesToPathlines1Display.ScaleTransferFunction = 'PiecewiseFunction'
temporalParticlesToPathlines1Display.OpacityArray = ['POINTS', 'TrailId']
temporalParticlesToPathlines1Display.OpacityTransferFunction = 'PiecewiseFunction'
temporalParticlesToPathlines1Display.DataAxesGrid = 'GridAxesRepresentation'
temporalParticlesToPathlines1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
temporalParticlesToPathlines1Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
temporalParticlesToPathlines1Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1000.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
temporalParticlesToPathlines1Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1000.0, 1.0, 0.5, 0.0]

# show data from cylinderpvd
cylinderpvdDisplay = Show(cylinderpvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cylinderpvdDisplay.Representation = 'Surface'
cylinderpvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
cylinderpvdDisplay.ColorArrayName = [None, '']
cylinderpvdDisplay.DiffuseColor = [0.0, 0.0, 0.0]
cylinderpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
cylinderpvdDisplay.SelectOrientationVectors = 'None'
cylinderpvdDisplay.ScaleFactor = 0.25863996744155887
cylinderpvdDisplay.SelectScaleArray = 'None'
cylinderpvdDisplay.GlyphType = 'Arrow'
cylinderpvdDisplay.GlyphTableIndexArray = 'None'
cylinderpvdDisplay.GaussianRadius = 0.012931998372077941
cylinderpvdDisplay.SetScaleArray = [None, '']
cylinderpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
cylinderpvdDisplay.OpacityArray = [None, '']
cylinderpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
cylinderpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
cylinderpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
cylinderpvdDisplay.ScalarOpacityUnitDistance = 0.06717646765003486

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
cylinderpvdDisplay.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
cylinderpvdDisplay.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cylinderpvdDisplay.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

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
clip4Display.ColorArrayName = [None, '']
clip4Display.DiffuseColor = [0.0, 0.0, 0.0]
clip4Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip4Display.SelectOrientationVectors = 'None'
clip4Display.ScaleFactor = 0.25999999046325684
clip4Display.SelectScaleArray = 'None'
clip4Display.GlyphType = 'Arrow'
clip4Display.GlyphTableIndexArray = 'None'
clip4Display.GaussianRadius = 0.012999999523162843
clip4Display.SetScaleArray = [None, '']
clip4Display.ScaleTransferFunction = 'PiecewiseFunction'
clip4Display.OpacityArray = [None, '']
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

# setup the color legend parameters for each legend in this view

# get color legend/bar for velocityLUT in view renderView1
velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)
velocityLUTColorBar.AutoOrient = 0
velocityLUTColorBar.Orientation = 'Vertical'
velocityLUTColorBar.WindowLocation = 'AnyLocation'
velocityLUTColorBar.Position = [0.05, 0.8058614232209736]
velocityLUTColorBar.Title = '|u| (cm/s)'
velocityLUTColorBar.ComponentTitle = ''
velocityLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
velocityLUTColorBar.TitleFontFamily = 'Times'
velocityLUTColorBar.TitleFontSize = 24*4
velocityLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
velocityLUTColorBar.LabelFontFamily = 'Times'
velocityLUTColorBar.LabelFontSize = 20*4
velocityLUTColorBar.LabelFormat = '%-#1.0f'
velocityLUTColorBar.UseCustomLabels = 1
velocityLUTColorBar.CustomLabels = [0.0, 50.0, 100.0, 150.0, 200.0]
velocityLUTColorBar.AddRangeLabels = 0
velocityLUTColorBar.RangeLabelFormat = '%-#1.0f'
velocityLUTColorBar.DrawAnnotations = 0
velocityLUTColorBar.ScalarBarThickness = 20*4
velocityLUTColorBar.ScalarBarLength = 0.17000000000000004

# set color bar visibility
velocityLUTColorBar.Visibility = 1

# show color legend
temporalParticlesToPathlines1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from temporalParticlesToPathlines1
temporalParticlesToPathlines1Display_1 = Show(temporalParticlesToPathlines1, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
temporalParticlesToPathlines1Display_1.Representation = 'Surface'
temporalParticlesToPathlines1Display_1.ColorArrayName = ['POINTS', 'velocity']
temporalParticlesToPathlines1Display_1.LookupTable = velocityLUT
temporalParticlesToPathlines1Display_1.OSPRayScaleArray = 'TrailId'
temporalParticlesToPathlines1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
temporalParticlesToPathlines1Display_1.SelectOrientationVectors = 'None'
temporalParticlesToPathlines1Display_1.ScaleFactor = 0.84375
temporalParticlesToPathlines1Display_1.SelectScaleArray = 'TrailId'
temporalParticlesToPathlines1Display_1.GlyphType = 'Arrow'
temporalParticlesToPathlines1Display_1.GlyphTableIndexArray = 'TrailId'
temporalParticlesToPathlines1Display_1.GaussianRadius = 0.0421875
temporalParticlesToPathlines1Display_1.SetScaleArray = ['POINTS', 'TrailId']
temporalParticlesToPathlines1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
temporalParticlesToPathlines1Display_1.OpacityArray = ['POINTS', 'TrailId']
temporalParticlesToPathlines1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
temporalParticlesToPathlines1Display_1.DataAxesGrid = 'GridAxesRepresentation'
temporalParticlesToPathlines1Display_1.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
temporalParticlesToPathlines1Display_1.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
temporalParticlesToPathlines1Display_1.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1023.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
temporalParticlesToPathlines1Display_1.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1023.0, 1.0, 0.5, 0.0]

# show data from aortic_valvepvd
aortic_valvepvdDisplay = Show(aortic_valvepvd, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
aortic_valvepvdDisplay.Representation = 'Surface'
aortic_valvepvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
aortic_valvepvdDisplay.ColorArrayName = [None, '']
aortic_valvepvdDisplay.DiffuseColor = [0.0, 0.0, 0.0]
aortic_valvepvdDisplay.Opacity = 0.2
aortic_valvepvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
aortic_valvepvdDisplay.SelectOrientationVectors = 'None'
aortic_valvepvdDisplay.ScaleFactor = 0.4479223966598511
aortic_valvepvdDisplay.SelectScaleArray = 'None'
aortic_valvepvdDisplay.GlyphType = 'Arrow'
aortic_valvepvdDisplay.GlyphTableIndexArray = 'None'
aortic_valvepvdDisplay.GaussianRadius = 0.022396119832992552
aortic_valvepvdDisplay.SetScaleArray = [None, '']
aortic_valvepvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
aortic_valvepvdDisplay.OpacityArray = [None, '']
aortic_valvepvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
aortic_valvepvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
aortic_valvepvdDisplay.PolarAxes = 'PolarAxesRepresentation'
aortic_valvepvdDisplay.ScalarOpacityUnitDistance = 0.08164692686596023

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
aortic_valvepvdDisplay.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
aortic_valvepvdDisplay.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
aortic_valvepvdDisplay.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from cylinderpvd
cylinderpvdDisplay_1 = Show(cylinderpvd, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cylinderpvdDisplay_1.Representation = 'Surface'
cylinderpvdDisplay_1.AmbientColor = [0.0, 0.0, 0.0]
cylinderpvdDisplay_1.ColorArrayName = [None, '']
cylinderpvdDisplay_1.DiffuseColor = [0.0, 0.0, 0.0]
cylinderpvdDisplay_1.Opacity = 0.04
cylinderpvdDisplay_1.OSPRayScaleFunction = 'PiecewiseFunction'
cylinderpvdDisplay_1.SelectOrientationVectors = 'None'
cylinderpvdDisplay_1.ScaleFactor = 0.2581812381744385
cylinderpvdDisplay_1.SelectScaleArray = 'None'
cylinderpvdDisplay_1.GlyphType = 'Arrow'
cylinderpvdDisplay_1.GlyphTableIndexArray = 'None'
cylinderpvdDisplay_1.GaussianRadius = 0.012909061908721924
cylinderpvdDisplay_1.SetScaleArray = [None, '']
cylinderpvdDisplay_1.ScaleTransferFunction = 'PiecewiseFunction'
cylinderpvdDisplay_1.OpacityArray = [None, '']
cylinderpvdDisplay_1.OpacityTransferFunction = 'PiecewiseFunction'
cylinderpvdDisplay_1.DataAxesGrid = 'GridAxesRepresentation'
cylinderpvdDisplay_1.PolarAxes = 'PolarAxesRepresentation'
cylinderpvdDisplay_1.ScalarOpacityUnitDistance = 0.06706963278681245

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
cylinderpvdDisplay_1.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
cylinderpvdDisplay_1.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cylinderpvdDisplay_1.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# find source
temporalParticlesToPathlines1_1 = FindSource('TemporalParticlesToPathlines1')

# show data from temporalParticlesToPathlines1_1
temporalParticlesToPathlines1_1Display = Show(OutputPort(temporalParticlesToPathlines1_1, 1), renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
temporalParticlesToPathlines1_1Display.Representation = 'Surface'
temporalParticlesToPathlines1_1Display.ColorArrayName = ['POINTS', 'velocity']
temporalParticlesToPathlines1_1Display.LookupTable = velocityLUT
temporalParticlesToPathlines1_1Display.PointSize = 1.0
temporalParticlesToPathlines1_1Display.OSPRayScaleArray = 'velocity'
temporalParticlesToPathlines1_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
temporalParticlesToPathlines1_1Display.SelectOrientationVectors = 'None'
temporalParticlesToPathlines1_1Display.ScaleFactor = 0.84375
temporalParticlesToPathlines1_1Display.SelectScaleArray = 'None'
temporalParticlesToPathlines1_1Display.GlyphType = 'Arrow'
temporalParticlesToPathlines1_1Display.GlyphTableIndexArray = 'None'
temporalParticlesToPathlines1_1Display.GaussianRadius = 0.0421875
temporalParticlesToPathlines1_1Display.SetScaleArray = ['POINTS', 'velocity']
temporalParticlesToPathlines1_1Display.ScaleTransferFunction = 'PiecewiseFunction'
temporalParticlesToPathlines1_1Display.OpacityArray = ['POINTS', 'velocity']
temporalParticlesToPathlines1_1Display.OpacityTransferFunction = 'PiecewiseFunction'
temporalParticlesToPathlines1_1Display.DataAxesGrid = 'GridAxesRepresentation'
temporalParticlesToPathlines1_1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
temporalParticlesToPathlines1_1Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
temporalParticlesToPathlines1_1Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
temporalParticlesToPathlines1_1Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'velocity'
velocityPWF = GetOpacityTransferFunction('velocity')
velocityPWF.Points = [0.0, 0.0, 0.5, 0.0, 200.0, 1.0, 0.5, 0.0]
velocityPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(annotateTime1)
# ----------------------------------------------------------------

tk = GetTimeKeeper()
timesteps = tk.TimestepValues
numTimesteps = len(timesteps)

# print "timesteps = ", timesteps
# print "numTimesteps = ", numTimesteps

if not ((numTimesteps == 1442) or (numTimesteps == 1443)): 
    raise ValueError('incorrect numTimesteps')

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

import time 

save_start = 0
previous = temporalParticlesToPathlines1.MaxTrackLength + 1 

no_save_start = save_start - previous
if no_save_start < 0:
    no_save_start = 0

# run first iteration, then stop and start over 
# hack to remove first 
# for frame in range(len(timesteps)):
#     animationScene1.AnimationTime = timesteps[frame]
#     Render()
#     SaveScreenshot(basename + str(frame).zfill(4) + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)
#     break 

# for frame in range(no_save_start,save_start):
#     animationScene1.AnimationTime = timesteps[frame]
#     Render()

for frame in range(save_start, len(timesteps)-1):
    if (frame % nprocs) == proc_num:
        render_min = max(frame - previous,0)

        t_start = time.time()

        # run all the previous frames for making this tail 
        for frame_temp in range(render_min, frame):
            animationScene1.AnimationTime = timesteps[frame_temp]
            Render()

        t_mid = time.time()
        t_preproc = t_mid - t_start

        # finally the current frame 
        animationScene1.AnimationTime = timesteps[frame]
        Render()
        SaveScreenshot(basename + str(frame).zfill(4) + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)

        t_final_and_render = time.time() - t_mid

        print "On proc ", proc_num, " of ", nprocs, " writing frame ", frame, " preprocess time ",  t_preproc, "s, final frame ", t_final_and_render, "s"

    # time.sleep(1)

    # if (frame % 10) == 0:
    #     time.sleep(10)

