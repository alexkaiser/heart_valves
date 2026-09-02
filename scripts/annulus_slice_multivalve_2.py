# state file generated using paraview version 5.13.3
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [842, 1162]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.39124321937561035, -6.339656591415405, -20.356863498687744]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-0.993780237321503, -46.47044463371684, -22.28741158396091]
renderView1.CameraFocalPoint = [-0.993780237321503, -6.339656591415405, -22.28741158396091]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 4.500776581113036
renderView1.CameraParallelProjection = 1
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.PolarGrid = 'Polar Grid Actor'
renderView1.UseColorPaletteForBackground = 0
renderView1.Background = [0.0, 0.0, 0.0]
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [778, 1162]
renderView2.AxesGrid = 'Grid Axes 3D Actor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [0.3915112018585205, -6.544241189956665, -20.362027168273926]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [-44.50829767092782, -5.791111557572403, -0.644453544273079]
renderView2.CameraFocalPoint = [-0.84005473912766, -6.497987764191308, -22.791799814083713]
renderView2.CameraViewUp = [-0.3927101709802329, 0.4718777615585791, -0.7893732322261386]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 2.238357017053534
renderView2.CameraParallelProjection = 1
renderView2.LegendGrid = 'Legend Grid Actor'
renderView2.PolarGrid = 'Polar Grid Actor'
renderView2.UseColorPaletteForBackground = 0
renderView2.Background = [0.0, 0.0, 0.0]
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView3 = CreateView('RenderView')
renderView3.ViewSize = [525, 562]
renderView3.AxesGrid = 'Grid Axes 3D Actor'
renderView3.OrientationAxesVisibility = 0
renderView3.CenterOfRotation = [0.3915412425994873, -6.544245004653931, -20.362072467803955]
renderView3.StereoType = 'Crystal Eyes'
renderView3.CameraPosition = [5.645375544416632, -43.48900223551617, -4.832205917175775]
renderView3.CameraFocalPoint = [0.25312823717469873, -7.792207816769519, -22.799888996372225]
renderView3.CameraViewUp = [-0.9006381858866369, 0.07904092602676187, 0.4273211791330489]
renderView3.CameraFocalDisk = 1.0
renderView3.CameraParallelScale = 1.8916643966537983
renderView3.CameraParallelProjection = 1
renderView3.LegendGrid = 'Legend Grid Actor'
renderView3.PolarGrid = 'Polar Grid Actor'
renderView3.UseColorPaletteForBackground = 0
renderView3.Background = [0.0, 0.0, 0.0]
renderView3.BackEnd = 'OSPRay raycaster'
renderView3.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView4 = CreateView('RenderView')
renderView4.ViewSize = [525, 561]
renderView4.AxesGrid = 'Grid Axes 3D Actor'
renderView4.OrientationAxesVisibility = 0
renderView4.CenterOfRotation = [0.39190173149108887, -6.492579698562622, -20.36062526702881]
renderView4.StereoType = 'Crystal Eyes'
renderView4.CameraPosition = [-10.808290570186946, -16.749785357164996, -48.27489554456502]
renderView4.CameraFocalPoint = [1.5624652168179471, -5.354844642778539, -21.381128589039797]
renderView4.CameraViewUp = [-0.9192404224395572, 0.09806314036647111, 0.38128816694805845]
renderView4.CameraFocalDisk = 1.0
renderView4.CameraParallelScale = 1.7941561895411071
renderView4.CameraParallelProjection = 1
renderView4.LegendGrid = 'Legend Grid Actor'
renderView4.PolarGrid = 'Polar Grid Actor'
renderView4.UseColorPaletteForBackground = 0
renderView4.Background = [0.0, 0.0, 0.0]
renderView4.BackEnd = 'OSPRay raycaster'
renderView4.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.390886)
layout1.AssignView(1, renderView1)
layout1.SplitHorizontal(2, 0.595598)
layout1.AssignView(5, renderView2)
layout1.SplitVertical(6, 0.500000)
layout1.AssignView(13, renderView3)
layout1.AssignView(14, renderView4)
layout1.SetSize(2147, 1162)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView4)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
aortic_no_partition_384_facespvd = PVDReader(registrationName='aortic_no_partition_384_faces.pvd', FileName='aortic_no_partition_384_faces.pvd')

# create a new 'Connectivity'
connectivity1 = Connectivity(registrationName='Connectivity1', Input=aortic_no_partition_384_facespvd)

# create a new 'Threshold'
threshold2 = Threshold(registrationName='Threshold2', Input=connectivity1)
threshold2.Scalars = ['POINTS', 'RegionId']
threshold2.LowerThreshold = 0.24
threshold2.UpperThreshold = 1.0

# create a new 'Threshold'
threshold1 = Threshold(registrationName='Threshold1', Input=connectivity1)
threshold1.Scalars = ['POINTS', 'RegionId']
threshold1.UpperThreshold = 0.76

# create a new 'PVD Reader'
aorta_native_hist3_pt25mm_384_facespvd = PVDReader(registrationName='aorta_native_hist3_pt25mm_384_faces.pvd', FileName='aorta_native_hist3_pt25mm_384_faces.pvd')
aorta_native_hist3_pt25mm_384_facespvd.CellArrays = ['Normals']
aorta_native_hist3_pt25mm_384_facespvd.PointArrays = ['Normals', 'pyvistaOriginalPointIds']

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=aorta_native_hist3_pt25mm_384_facespvd)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'pyvistaOriginalPointIds']
clip1.Value = 79450.5
clip1.Invert = 0

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [-0.6721283644003664, -6.384510840593238, -23.11101180938841]
clip1.ClipType.Normal = [0.7934753815452537, 0.007764814857901565, -0.608552813264253]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [0.39405035972595215, -6.545432806015015, -20.360831260681152]

# create a new 'Clip'
clip3 = Clip(registrationName='Clip3', Input=aorta_native_hist3_pt25mm_384_facespvd)
clip3.ClipType = 'Plane'
clip3.HyperTreeGridClipper = 'Plane'
clip3.Scalars = ['POINTS', 'pyvistaOriginalPointIds']
clip3.Value = 79460.5
clip3.Invert = 0

# init the 'Plane' selected for 'ClipType'
clip3.ClipType.Origin = [0.2994667519130252, -6.367706087613032, -23.62526988386739]
clip3.ClipType.Normal = [0.3899992701891607, 0.3592358170329242, 0.847850338805505]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip3.HyperTreeGridClipper.Origin = [0.2925170660018921, -6.325026035308838, -23.640156745910645]

# create a new 'Clip'
clip5 = Clip(registrationName='Clip5', Input=threshold2)
clip5.ClipType = 'Plane'
clip5.HyperTreeGridClipper = 'Plane'
clip5.Scalars = ['POINTS', 'pyvistaOriginalPointIds']
clip5.Value = 79460.5
clip5.Invert = 0

# init the 'Plane' selected for 'ClipType'
clip5.ClipType.Origin = [0.2994667519130252, -6.367706087613032, -23.62526988386739]
clip5.ClipType.Normal = [0.3899992701891607, 0.3592358170329242, 0.847850338805505]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip5.HyperTreeGridClipper.Origin = [0.2925170660018921, -6.325026035308838, -23.640156745910645]

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=aorta_native_hist3_pt25mm_384_facespvd)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.PointMergeMethod = 'Uniform Binning'

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.3940403461456299, -6.545432090759277, -20.360841274261475]
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0.3940403461456299, -6.545432090759277, -20.360841274261475]

# create a new 'PVD Reader'
eulerian_vars_restricted_pointspvd = PVDReader(registrationName='eulerian_vars_restricted_points.pvd', FileName='eulerian_vars_restricted_points.pvd')
eulerian_vars_restricted_pointspvd.PointArrays = ['U', 'P', 'SelectedPoints']

# create a new 'Clip'
clip2 = Clip(registrationName='Clip2', Input=clip1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'pyvistaOriginalPointIds']
clip2.Value = 79450.0

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Origin = [1.5941076708404216, -6.43388460103102, -23.79526330724851]
clip2.ClipType.Normal = [0.8903918008353289, -0.08987440323827121, -0.4462342800007467]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [1.8817338347434998, -6.463640451431274, -23.942702293395996]

# create a new 'Clip'
clip4 = Clip(registrationName='Clip4', Input=threshold1)
clip4.ClipType = 'Plane'
clip4.HyperTreeGridClipper = 'Plane'
clip4.Scalars = ['POINTS', 'pyvistaOriginalPointIds']
clip4.Value = 79460.5
clip4.Invert = 0

# init the 'Plane' selected for 'ClipType'
clip4.ClipType.Origin = [0.2994667519130252, -6.367706087613032, -23.62526988386739]
clip4.ClipType.Normal = [0.3899992701891607, 0.3592358170329242, 0.847850338805505]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip4.HyperTreeGridClipper.Origin = [0.2925170660018921, -6.325026035308838, -23.640156745910645]

# create a new 'Slice'
slice2 = Slice(registrationName='Slice2', Input=eulerian_vars_restricted_pointspvd)
slice2.SliceType = 'Plane'
slice2.HyperTreeGridSlicer = 'Plane'
slice2.SliceOffsetValues = [0.0]
slice2.PointMergeMethod = 'Uniform Binning'

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [0.6097414493560791, -6.470432758331299, -20.085942268371582]
slice2.SliceType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice2.HyperTreeGridSlicer.Origin = [0.6097414493560791, -6.470432758331299, -20.085942268371582]

# create a new 'Clip'
clip6 = Clip(registrationName='Clip6', Input=aorta_native_hist3_pt25mm_384_facespvd)
clip6.ClipType = 'Plane'
clip6.HyperTreeGridClipper = 'Plane'
clip6.Scalars = ['POINTS', 'pyvistaOriginalPointIds']
clip6.Value = 79450.5
clip6.Invert = 0

# init the 'Plane' selected for 'ClipType'
clip6.ClipType.Origin = [0.979886285797685, -6.451693189558228, -23.893054984589142]
clip6.ClipType.Normal = [-0.10344927702281514, 0.9383240770417535, -0.32991843465802057]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip6.HyperTreeGridClipper.Origin = [0.3915412425994873, -6.544245004653931, -20.362072467803955]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = [None, '']
slice1Display.LineWidth = 2.0
slice1Display.SelectNormalArray = 'Normals'
slice1Display.SelectTangentArray = 'None'
slice1Display.SelectTCoordArray = 'None'
slice1Display.TextureTransform = 'Transform2'
slice1Display.OSPRayScaleArray = 'Normals'
slice1Display.OSPRayScaleFunction = 'Piecewise Function'
slice1Display.Assembly = ''
slice1Display.SelectedBlockSelectors = ['']
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 1.1179544448852539
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.0558977222442627
slice1Display.SetScaleArray = ['POINTS', 'Normals']
slice1Display.ScaleTransferFunction = 'Piecewise Function'
slice1Display.OpacityArray = ['POINTS', 'Normals']
slice1Display.OpacityTransferFunction = 'Piecewise Function'
slice1Display.DataAxesGrid = 'Grid Axes Representation'
slice1Display.PolarAxes = 'Polar Axes Representation'
slice1Display.SelectInputVectors = ['POINTS', 'Normals']
slice1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [-0.9999951124191284, 0.2589285671710968, 0.5, 0.0, 0.9961395859718323, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [-0.9999951124191284, 0.2589285671710968, 0.5, 0.0, 0.9961395859718323, 1.0, 0.5, 0.0]

# show data from slice2
slice2Display = Show(slice2, renderView1, 'GeometryRepresentation')

# get 2D transfer function for 'U'
uTF2D = GetTransferFunction2D('U')
uTF2D.ScalarRangeInitialized = 1
uTF2D.Range = [0.0, 300.0, 0.0, 1.0]

# get color transfer function/color map for 'U'
uLUT = GetColorTransferFunction('U')
uLUT.AutomaticRescaleRangeMode = 'Never'
uLUT.TransferFunction2D = uTF2D
uLUT.RGBPoints = [0.0, 0.0, 0.0, 0.0, 120.0, 0.901960784314, 0.0, 0.0, 240.0, 0.901960784314, 0.901960784314, 0.0, 300.0, 1.0, 1.0, 1.0]
uLUT.ColorSpace = 'RGB'
uLUT.NanColor = [0.0, 0.498039215686, 1.0]
uLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice2Display.Representation = 'Surface'
slice2Display.ColorArrayName = ['POINTS', 'U']
slice2Display.LookupTable = uLUT
slice2Display.SelectNormalArray = 'None'
slice2Display.SelectTangentArray = 'None'
slice2Display.SelectTCoordArray = 'None'
slice2Display.TextureTransform = 'Transform2'
slice2Display.OSPRayScaleArray = 'U'
slice2Display.OSPRayScaleFunction = 'Piecewise Function'
slice2Display.Assembly = ''
slice2Display.SelectedBlockSelectors = ['']
slice2Display.SelectOrientationVectors = 'None'
slice2Display.ScaleFactor = 1.0600000381469727
slice2Display.SelectScaleArray = 'U'
slice2Display.GlyphType = 'Arrow'
slice2Display.GlyphTableIndexArray = 'U'
slice2Display.GaussianRadius = 0.05300000190734863
slice2Display.SetScaleArray = ['POINTS', 'U']
slice2Display.ScaleTransferFunction = 'Piecewise Function'
slice2Display.OpacityArray = ['POINTS', 'U']
slice2Display.OpacityTransferFunction = 'Piecewise Function'
slice2Display.DataAxesGrid = 'Grid Axes Representation'
slice2Display.PolarAxes = 'Polar Axes Representation'
slice2Display.SelectInputVectors = ['POINTS', 'U']
slice2Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
slice2Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice2Display.ScaleTransferFunction.Points = [-43.12060950706473, 0.2589285671710968, 0.5, 0.0, 84.73390347060922, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice2Display.OpacityTransferFunction.Points = [-43.12060950706473, 0.2589285671710968, 0.5, 0.0, 84.73390347060922, 1.0, 0.5, 0.0]

# show data from threshold1
threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.AmbientColor = [1.0, 0.0, 0.0]
threshold1Display.ColorArrayName = ['POINTS', '']
threshold1Display.DiffuseColor = [1.0, 0.0, 0.0]
threshold1Display.SelectNormalArray = 'None'
threshold1Display.SelectTangentArray = 'None'
threshold1Display.SelectTCoordArray = 'None'
threshold1Display.TextureTransform = 'Transform2'
threshold1Display.OSPRayScaleArray = 'RegionId'
threshold1Display.OSPRayScaleFunction = 'Piecewise Function'
threshold1Display.Assembly = ''
threshold1Display.SelectedBlockSelectors = ['']
threshold1Display.SelectOrientationVectors = 'None'
threshold1Display.ScaleFactor = 0.28522924184799203
threshold1Display.SelectScaleArray = 'RegionId'
threshold1Display.GlyphType = 'Arrow'
threshold1Display.GlyphTableIndexArray = 'RegionId'
threshold1Display.GaussianRadius = 0.014261462092399603
threshold1Display.SetScaleArray = ['POINTS', 'RegionId']
threshold1Display.ScaleTransferFunction = 'Piecewise Function'
threshold1Display.OpacityArray = ['POINTS', 'RegionId']
threshold1Display.OpacityTransferFunction = 'Piecewise Function'
threshold1Display.DataAxesGrid = 'Grid Axes Representation'
threshold1Display.PolarAxes = 'Polar Axes Representation'
threshold1Display.ScalarOpacityUnitDistance = 0.143736451135407
threshold1Display.OpacityArrayName = ['POINTS', 'RegionId']
threshold1Display.SelectInputVectors = [None, '']
threshold1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
threshold1Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
threshold1Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
threshold1Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# show data from threshold2
threshold2Display = Show(threshold2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold2Display.Representation = 'Surface'
threshold2Display.AmbientColor = [0.0, 0.0, 1.0]
threshold2Display.ColorArrayName = ['POINTS', '']
threshold2Display.DiffuseColor = [0.0, 0.0, 1.0]
threshold2Display.SelectNormalArray = 'None'
threshold2Display.SelectTangentArray = 'None'
threshold2Display.SelectTCoordArray = 'None'
threshold2Display.TextureTransform = 'Transform2'
threshold2Display.OSPRayScaleArray = 'RegionId'
threshold2Display.OSPRayScaleFunction = 'Piecewise Function'
threshold2Display.Assembly = ''
threshold2Display.SelectedBlockSelectors = ['']
threshold2Display.SelectOrientationVectors = 'None'
threshold2Display.ScaleFactor = 0.2903535842895501
threshold2Display.SelectScaleArray = 'RegionId'
threshold2Display.GlyphType = 'Arrow'
threshold2Display.GlyphTableIndexArray = 'RegionId'
threshold2Display.GaussianRadius = 0.014517679214477504
threshold2Display.SetScaleArray = ['POINTS', 'RegionId']
threshold2Display.ScaleTransferFunction = 'Piecewise Function'
threshold2Display.OpacityArray = ['POINTS', 'RegionId']
threshold2Display.OpacityTransferFunction = 'Piecewise Function'
threshold2Display.DataAxesGrid = 'Grid Axes Representation'
threshold2Display.PolarAxes = 'Polar Axes Representation'
threshold2Display.ScalarOpacityUnitDistance = 0.15079223836623884
threshold2Display.OpacityArrayName = ['POINTS', 'RegionId']
threshold2Display.SelectInputVectors = [None, '']
threshold2Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
threshold2Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
threshold2Display.ScaleTransferFunction.Points = [1.0, 0.2589285671710968, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
threshold2Display.OpacityTransferFunction.Points = [1.0, 0.2589285671710968, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)
uLUTColorBar.WindowLocation = 'Any Location'
uLUTColorBar.Position = [0.82, 0.65]
uLUTColorBar.Title = 'U'
uLUTColorBar.ComponentTitle = 'Magnitude'
uLUTColorBar.TitleColor = [0.9999, 0.9999, 0.9999]
uLUTColorBar.TitleFontSize = 18
uLUTColorBar.LabelColor = [0.9999, 0.9999, 0.9999]
uLUTColorBar.LabelFontSize = 18
uLUTColorBar.UseCustomLabels = 1
uLUTColorBar.CustomLabels = [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0]
uLUTColorBar.AddRangeLabels = 0
uLUTColorBar.RangeLabelFormat = '%-#6.3g'

# set color bar visibility
uLUTColorBar.Visibility = 1

# show color legend
slice2Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from threshold1
threshold1Display_1 = Show(threshold1, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold1Display_1.Representation = 'Surface'
threshold1Display_1.AmbientColor = [1.0, 0.0, 0.0]
threshold1Display_1.ColorArrayName = ['POINTS', '']
threshold1Display_1.DiffuseColor = [1.0, 0.0, 0.0]
threshold1Display_1.SelectNormalArray = 'None'
threshold1Display_1.SelectTangentArray = 'None'
threshold1Display_1.SelectTCoordArray = 'None'
threshold1Display_1.TextureTransform = 'Transform2'
threshold1Display_1.OSPRayScaleArray = 'RegionId'
threshold1Display_1.OSPRayScaleFunction = 'Piecewise Function'
threshold1Display_1.Assembly = ''
threshold1Display_1.SelectedBlockSelectors = ['']
threshold1Display_1.SelectOrientationVectors = 'None'
threshold1Display_1.ScaleFactor = 0.28522924184799203
threshold1Display_1.SelectScaleArray = 'RegionId'
threshold1Display_1.GlyphType = 'Arrow'
threshold1Display_1.GlyphTableIndexArray = 'RegionId'
threshold1Display_1.GaussianRadius = 0.014261462092399603
threshold1Display_1.SetScaleArray = ['POINTS', 'RegionId']
threshold1Display_1.ScaleTransferFunction = 'Piecewise Function'
threshold1Display_1.OpacityArray = ['POINTS', 'RegionId']
threshold1Display_1.OpacityTransferFunction = 'Piecewise Function'
threshold1Display_1.DataAxesGrid = 'Grid Axes Representation'
threshold1Display_1.PolarAxes = 'Polar Axes Representation'
threshold1Display_1.ScalarOpacityUnitDistance = 0.143736451135407
threshold1Display_1.OpacityArrayName = ['POINTS', 'RegionId']
threshold1Display_1.SelectInputVectors = [None, '']
threshold1Display_1.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
threshold1Display_1.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
threshold1Display_1.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
threshold1Display_1.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# show data from threshold2
threshold2Display_1 = Show(threshold2, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold2Display_1.Representation = 'Surface'
threshold2Display_1.AmbientColor = [0.0, 0.0, 1.0]
threshold2Display_1.ColorArrayName = ['POINTS', '']
threshold2Display_1.DiffuseColor = [0.0, 0.0, 1.0]
threshold2Display_1.SelectNormalArray = 'None'
threshold2Display_1.SelectTangentArray = 'None'
threshold2Display_1.SelectTCoordArray = 'None'
threshold2Display_1.TextureTransform = 'Transform2'
threshold2Display_1.OSPRayScaleArray = 'RegionId'
threshold2Display_1.OSPRayScaleFunction = 'Piecewise Function'
threshold2Display_1.Assembly = ''
threshold2Display_1.SelectedBlockSelectors = ['']
threshold2Display_1.SelectOrientationVectors = 'None'
threshold2Display_1.ScaleFactor = 0.2903535842895501
threshold2Display_1.SelectScaleArray = 'RegionId'
threshold2Display_1.GlyphType = 'Arrow'
threshold2Display_1.GlyphTableIndexArray = 'RegionId'
threshold2Display_1.GaussianRadius = 0.014517679214477504
threshold2Display_1.SetScaleArray = ['POINTS', 'RegionId']
threshold2Display_1.ScaleTransferFunction = 'Piecewise Function'
threshold2Display_1.OpacityArray = ['POINTS', 'RegionId']
threshold2Display_1.OpacityTransferFunction = 'Piecewise Function'
threshold2Display_1.DataAxesGrid = 'Grid Axes Representation'
threshold2Display_1.PolarAxes = 'Polar Axes Representation'
threshold2Display_1.ScalarOpacityUnitDistance = 0.15079223836623884
threshold2Display_1.OpacityArrayName = ['POINTS', 'RegionId']
threshold2Display_1.SelectInputVectors = [None, '']
threshold2Display_1.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
threshold2Display_1.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
threshold2Display_1.ScaleTransferFunction.Points = [1.0, 0.2589285671710968, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
threshold2Display_1.OpacityTransferFunction.Points = [1.0, 0.2589285671710968, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# show data from clip2
clip2Display = Show(clip2, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip2Display.Representation = 'Surface'
clip2Display.ColorArrayName = [None, '']
clip2Display.SelectNormalArray = 'Normals'
clip2Display.SelectTangentArray = 'None'
clip2Display.SelectTCoordArray = 'None'
clip2Display.TextureTransform = 'Transform2'
clip2Display.OSPRayScaleArray = 'Normals'
clip2Display.OSPRayScaleFunction = 'Piecewise Function'
clip2Display.Assembly = ''
clip2Display.SelectedBlockSelectors = ['']
clip2Display.SelectOrientationVectors = 'None'
clip2Display.ScaleFactor = 0.3613387703895569
clip2Display.SelectScaleArray = 'None'
clip2Display.GlyphType = 'Arrow'
clip2Display.GlyphTableIndexArray = 'None'
clip2Display.GaussianRadius = 0.018066938519477844
clip2Display.SetScaleArray = ['POINTS', 'Normals']
clip2Display.ScaleTransferFunction = 'Piecewise Function'
clip2Display.OpacityArray = ['POINTS', 'Normals']
clip2Display.OpacityTransferFunction = 'Piecewise Function'
clip2Display.DataAxesGrid = 'Grid Axes Representation'
clip2Display.PolarAxes = 'Polar Axes Representation'
clip2Display.ScalarOpacityUnitDistance = 0.10408644332096584
clip2Display.OpacityArrayName = ['POINTS', 'Normals']
clip2Display.SelectInputVectors = ['POINTS', 'Normals']
clip2Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
clip2Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
clip2Display.ScaleTransferFunction.Points = [-0.7214022278785706, 0.2589285671710968, 0.5, 0.0, 0.7480586767196655, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
clip2Display.OpacityTransferFunction.Points = [-0.7214022278785706, 0.2589285671710968, 0.5, 0.0, 0.7480586767196655, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView3'
# ----------------------------------------------------------------

# show data from threshold1
threshold1Display_2 = Show(threshold1, renderView3, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold1Display_2.Representation = 'Surface'
threshold1Display_2.AmbientColor = [1.0, 0.0, 0.0]
threshold1Display_2.ColorArrayName = ['POINTS', '']
threshold1Display_2.DiffuseColor = [1.0, 0.0, 0.0]
threshold1Display_2.SelectNormalArray = 'None'
threshold1Display_2.SelectTangentArray = 'None'
threshold1Display_2.SelectTCoordArray = 'None'
threshold1Display_2.TextureTransform = 'Transform2'
threshold1Display_2.OSPRayScaleArray = 'RegionId'
threshold1Display_2.OSPRayScaleFunction = 'Piecewise Function'
threshold1Display_2.Assembly = ''
threshold1Display_2.SelectedBlockSelectors = ['']
threshold1Display_2.SelectOrientationVectors = 'None'
threshold1Display_2.ScaleFactor = 0.288239109516144
threshold1Display_2.SelectScaleArray = 'RegionId'
threshold1Display_2.GlyphType = 'Arrow'
threshold1Display_2.GlyphTableIndexArray = 'RegionId'
threshold1Display_2.GaussianRadius = 0.0144119554758072
threshold1Display_2.SetScaleArray = ['POINTS', 'RegionId']
threshold1Display_2.ScaleTransferFunction = 'Piecewise Function'
threshold1Display_2.OpacityArray = ['POINTS', 'RegionId']
threshold1Display_2.OpacityTransferFunction = 'Piecewise Function'
threshold1Display_2.DataAxesGrid = 'Grid Axes Representation'
threshold1Display_2.PolarAxes = 'Polar Axes Representation'
threshold1Display_2.ScalarOpacityUnitDistance = 0.14711442215262313
threshold1Display_2.OpacityArrayName = ['POINTS', 'RegionId']
threshold1Display_2.SelectInputVectors = [None, '']
threshold1Display_2.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
threshold1Display_2.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
threshold1Display_2.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
threshold1Display_2.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# show data from threshold2
threshold2Display_2 = Show(threshold2, renderView3, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold2Display_2.Representation = 'Surface'
threshold2Display_2.AmbientColor = [0.0, 0.0, 1.0]
threshold2Display_2.ColorArrayName = ['POINTS', '']
threshold2Display_2.DiffuseColor = [0.0, 0.0, 1.0]
threshold2Display_2.SelectNormalArray = 'None'
threshold2Display_2.SelectTangentArray = 'None'
threshold2Display_2.SelectTCoordArray = 'None'
threshold2Display_2.TextureTransform = 'Transform2'
threshold2Display_2.OSPRayScaleArray = 'RegionId'
threshold2Display_2.OSPRayScaleFunction = 'Piecewise Function'
threshold2Display_2.Assembly = ''
threshold2Display_2.SelectedBlockSelectors = ['']
threshold2Display_2.SelectOrientationVectors = 'None'
threshold2Display_2.ScaleFactor = 0.29037780761719
threshold2Display_2.SelectScaleArray = 'RegionId'
threshold2Display_2.GlyphType = 'Arrow'
threshold2Display_2.GlyphTableIndexArray = 'RegionId'
threshold2Display_2.GaussianRadius = 0.0145188903808595
threshold2Display_2.SetScaleArray = ['POINTS', 'RegionId']
threshold2Display_2.ScaleTransferFunction = 'Piecewise Function'
threshold2Display_2.OpacityArray = ['POINTS', 'RegionId']
threshold2Display_2.OpacityTransferFunction = 'Piecewise Function'
threshold2Display_2.DataAxesGrid = 'Grid Axes Representation'
threshold2Display_2.PolarAxes = 'Polar Axes Representation'
threshold2Display_2.ScalarOpacityUnitDistance = 0.1498385119816704
threshold2Display_2.OpacityArrayName = ['POINTS', 'RegionId']
threshold2Display_2.SelectInputVectors = [None, '']
threshold2Display_2.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
threshold2Display_2.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
threshold2Display_2.ScaleTransferFunction.Points = [1.0, 0.2589285671710968, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
threshold2Display_2.OpacityTransferFunction.Points = [1.0, 0.2589285671710968, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# show data from clip6
clip6Display = Show(clip6, renderView3, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip6Display.Representation = 'Surface'
clip6Display.ColorArrayName = [None, '']
clip6Display.SelectNormalArray = 'Normals'
clip6Display.SelectTangentArray = 'None'
clip6Display.SelectTCoordArray = 'None'
clip6Display.TextureTransform = 'Transform2'
clip6Display.OSPRayScaleArray = 'Normals'
clip6Display.OSPRayScaleFunction = 'Piecewise Function'
clip6Display.Assembly = ''
clip6Display.SelectedBlockSelectors = ['']
clip6Display.SelectOrientationVectors = 'None'
clip6Display.ScaleFactor = 1.0090550422668458
clip6Display.SelectScaleArray = 'None'
clip6Display.GlyphType = 'Arrow'
clip6Display.GlyphTableIndexArray = 'None'
clip6Display.GaussianRadius = 0.050452752113342283
clip6Display.SetScaleArray = ['POINTS', 'Normals']
clip6Display.ScaleTransferFunction = 'Piecewise Function'
clip6Display.OpacityArray = ['POINTS', 'Normals']
clip6Display.OpacityTransferFunction = 'Piecewise Function'
clip6Display.DataAxesGrid = 'Grid Axes Representation'
clip6Display.PolarAxes = 'Polar Axes Representation'
clip6Display.ScalarOpacityUnitDistance = 0.18971218892335076
clip6Display.OpacityArrayName = ['POINTS', 'Normals']
clip6Display.SelectInputVectors = ['POINTS', 'Normals']
clip6Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
clip6Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
clip6Display.ScaleTransferFunction.Points = [-1.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
clip6Display.OpacityTransferFunction.Points = [-1.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView4'
# ----------------------------------------------------------------

# show data from clip3
clip3Display = Show(clip3, renderView4, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip3Display.Representation = 'Surface'
clip3Display.ColorArrayName = [None, '']
clip3Display.SelectNormalArray = 'Normals'
clip3Display.SelectTangentArray = 'None'
clip3Display.SelectTCoordArray = 'None'
clip3Display.TextureTransform = 'Transform2'
clip3Display.OSPRayScaleArray = 'Normals'
clip3Display.OSPRayScaleFunction = 'Piecewise Function'
clip3Display.Assembly = ''
clip3Display.SelectedBlockSelectors = ['']
clip3Display.SelectOrientationVectors = 'None'
clip3Display.ScaleFactor = 0.3051143050193787
clip3Display.SelectScaleArray = 'None'
clip3Display.GlyphType = 'Arrow'
clip3Display.GlyphTableIndexArray = 'None'
clip3Display.GaussianRadius = 0.015255715250968933
clip3Display.SetScaleArray = ['POINTS', 'Normals']
clip3Display.ScaleTransferFunction = 'Piecewise Function'
clip3Display.OpacityArray = ['POINTS', 'Normals']
clip3Display.OpacityTransferFunction = 'Piecewise Function'
clip3Display.DataAxesGrid = 'Grid Axes Representation'
clip3Display.PolarAxes = 'Polar Axes Representation'
clip3Display.ScalarOpacityUnitDistance = 0.11484031791506064
clip3Display.OpacityArrayName = ['POINTS', 'Normals']
clip3Display.SelectInputVectors = ['POINTS', 'Normals']
clip3Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
clip3Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
clip3Display.ScaleTransferFunction.Points = [-0.46199992299079895, 0.2589285671710968, 0.5, 0.0, 0.7480586767196655, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
clip3Display.OpacityTransferFunction.Points = [-0.46199992299079895, 0.2589285671710968, 0.5, 0.0, 0.7480586767196655, 1.0, 0.5, 0.0]

# show data from clip4
clip4Display = Show(clip4, renderView4, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip4Display.Representation = 'Surface'
clip4Display.AmbientColor = [1.0, 0.0, 0.0]
clip4Display.ColorArrayName = ['POINTS', '']
clip4Display.DiffuseColor = [1.0, 0.0, 0.0]
clip4Display.SelectNormalArray = 'None'
clip4Display.SelectTangentArray = 'None'
clip4Display.SelectTCoordArray = 'None'
clip4Display.TextureTransform = 'Transform2'
clip4Display.OSPRayScaleArray = 'RegionId'
clip4Display.OSPRayScaleFunction = 'Piecewise Function'
clip4Display.Assembly = ''
clip4Display.SelectedBlockSelectors = ['']
clip4Display.SelectOrientationVectors = 'None'
clip4Display.ScaleFactor = 0.18832161970199454
clip4Display.SelectScaleArray = 'RegionId'
clip4Display.GlyphType = 'Arrow'
clip4Display.GlyphTableIndexArray = 'RegionId'
clip4Display.GaussianRadius = 0.009416080985099726
clip4Display.SetScaleArray = ['POINTS', 'RegionId']
clip4Display.ScaleTransferFunction = 'Piecewise Function'
clip4Display.OpacityArray = ['POINTS', 'RegionId']
clip4Display.OpacityTransferFunction = 'Piecewise Function'
clip4Display.DataAxesGrid = 'Grid Axes Representation'
clip4Display.PolarAxes = 'Polar Axes Representation'
clip4Display.ScalarOpacityUnitDistance = 0.13106112356992897
clip4Display.OpacityArrayName = ['POINTS', 'RegionId']
clip4Display.SelectInputVectors = [None, '']
clip4Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
clip4Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
clip4Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
clip4Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# show data from clip5
clip5Display = Show(clip5, renderView4, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip5Display.Representation = 'Surface'
clip5Display.AmbientColor = [0.0, 0.0, 1.0]
clip5Display.ColorArrayName = ['POINTS', '']
clip5Display.DiffuseColor = [0.0, 0.0, 1.0]
clip5Display.SelectNormalArray = 'None'
clip5Display.SelectTangentArray = 'None'
clip5Display.SelectTCoordArray = 'None'
clip5Display.TextureTransform = 'Transform2'
clip5Display.OSPRayScaleArray = 'RegionId'
clip5Display.OSPRayScaleFunction = 'Piecewise Function'
clip5Display.Assembly = ''
clip5Display.SelectedBlockSelectors = ['']
clip5Display.SelectOrientationVectors = 'None'
clip5Display.ScaleFactor = 0.2237316082766121
clip5Display.SelectScaleArray = 'RegionId'
clip5Display.GlyphType = 'Arrow'
clip5Display.GlyphTableIndexArray = 'RegionId'
clip5Display.GaussianRadius = 0.011186580413830605
clip5Display.SetScaleArray = ['POINTS', 'RegionId']
clip5Display.ScaleTransferFunction = 'Piecewise Function'
clip5Display.OpacityArray = ['POINTS', 'RegionId']
clip5Display.OpacityTransferFunction = 'Piecewise Function'
clip5Display.DataAxesGrid = 'Grid Axes Representation'
clip5Display.PolarAxes = 'Polar Axes Representation'
clip5Display.ScalarOpacityUnitDistance = 0.13414990950581338
clip5Display.OpacityArrayName = ['POINTS', 'RegionId']
clip5Display.SelectInputVectors = [None, '']
clip5Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
clip5Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
clip5Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
clip5Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'U'
uPWF = GetOpacityTransferFunction('U')
uPWF.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 300.0, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation scene

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = [renderView1, renderView2, renderView3, renderView4]
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime = 0.0
# hardcoded by the state export at only ~47% of the run, so times past it were
# clamped and the tail of the movie froze. leave commented to let the scene
# take its range from the loaded data.
# animationScene1.EndTime = 0.661005
animationScene1.PlayMode = 'Snap To TimeSteps'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(eulerian_vars_restricted_pointspvd)
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
    print("using default proc_num 0, nprocs = 1")
    proc_num = 0
    nprocs = 1

res = (4734, 2536)

animationScene1 = GetAnimationScene()
# animationScene1.GoToFirst()

# run first iteration, then stop and start over 
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
