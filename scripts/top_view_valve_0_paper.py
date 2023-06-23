# state file generated using paraview version 5.9.0

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
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [892, 804]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [1.9946428835391998, -18.485949516296387, -24.677393913269043]
renderView2.StereoType = 'Crystal Eyes'

# position 1
# renderView2.CameraPosition = [-1.530326693564913, -18.44633308984853, -20.4966926916046]
# renderView2.CameraFocalPoint = [9.33116904822627, -18.580073608467885, -33.394601971708234]
# renderView2.CameraViewUp = [-0.4158018371368334, 0.8356787190797338, -0.3588173779388197]

# position 2
renderView2.CameraPosition = [-1.7111433708949269, -18.444106645464068, -20.281974833588187]
renderView2.CameraFocalPoint = [9.15035237089623, -18.577847164083423, -33.17988411369181]
renderView2.CameraViewUp = [-0.4158018371368334, 0.8356787190797338, -0.3588173779388197]


renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 4.364352523816178
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView2)
layout1.SetSize((1180, 1080))

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView2)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
aorta_384pvd = PVDReader(registrationName='aorta_384.pvd', FileName='aorta_384.pvd')

# create a new 'XML Unstructured Grid Reader'
# aorta_384_volumetricmeshvtu = XMLUnstructuredGridReader(registrationName='aorta_384_volumetric.mesh.vtu', FileName=['/Users/alex/Dropbox/stanford/research_stanford/aortic_bicuspid_2020/aorta_384_volumetric.mesh.vtu'])
# aorta_384_volumetricmeshvtu.CellArrayStatus = ['ModelRegionID', 'GlobalElementID']
# aorta_384_volumetricmeshvtu.PointArrayStatus = ['GlobalNodeID']

# create a new 'Extract Selection'
# extractSelection1 = ExtractSelection(registrationName='ExtractSelection1', Input=aorta_384_volumetricmeshvtu)

# create a new 'PVD Reader'
aortic_no_partition_384_facespvd = PVDReader(registrationName='aortic_no_partition_384_faces.pvd', FileName='aortic_no_partition_384_faces.pvd')

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=aorta_384pvd)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', '']

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [1.3493156091274934, -18.499512458863364, -21.992646633131386]
clip1.ClipType.Normal = [-0.6572517378107009, -0.026955799993162867, 0.7531889125521846]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [1.5651809573173523, -18.485949516296387, -22.371617317199707]

# create a new 'PVD Reader'
# eulerian_varspvd = PVDReader(registrationName='eulerian_vars.pvd', FileName='eulerian_vars.pvd')
# eulerian_varspvd.CellArrays = ['U']

eulerian_varspvd = PVDReader(FileName='eulerian_vars_restricted_points.pvd')
eulerian_varspvd.PointArrays = ['U']

# create a new 'Resample With Dataset'
# resampleWithDataset1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=eulerian_varspvd,
#     DestinationMesh=aorta_384_volumetricmeshvtu)
# resampleWithDataset1.CellLocator = 'Static Cell Locator'

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=eulerian_varspvd)
calculator1.ResultArrayName = 'annulus_normal_projected'
calculator1.Function = '-0.644118902159037 * U_X + 0.007931209298900 * U_Y + 0.764884263010094 * U_Z'

# create a new 'Slice'
slice3 = Slice(registrationName='Slice3', Input=calculator1)
slice3.SliceType = 'Plane'
slice3.HyperTreeGridSlicer = 'Plane'
slice3.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice3.SliceType.Origin = [2.447601610336508, -18.376013668463226, -25.219513937752243]
slice3.SliceType.Normal = [0.644118902159037, -0.0079312092989, -0.764884263010094]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice3.HyperTreeGridSlicer.Origin = [1.5904362797737122, -18.485923767089844, -22.346393585205078]

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=eulerian_varspvd)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [1.590436339378357, -18.485923767089844, -22.346393585205078]
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [1.590436339378357, -18.485923767089844, -22.346393585205078]

# create a new 'Slice'
slice9 = Slice(registrationName='Slice9', Input=eulerian_varspvd)
slice9.SliceType = 'Plane'
slice9.HyperTreeGridSlicer = 'Plane'
slice9.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice9.SliceType.Origin = [1.0154999494552612, -18.4862003326416, -22.669400215148926]
slice9.SliceType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice9.HyperTreeGridSlicer.Origin = [1.0154999494552612, -18.4862003326416, -22.669400215148926]

# create a new 'Slice'
slice2 = Slice(registrationName='Slice2', Input=aorta_384pvd)
slice2.SliceType = 'Plane'
slice2.HyperTreeGridSlicer = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [1.565198004245758, -18.485950469970703, -22.37162494659424]
slice2.SliceType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice2.HyperTreeGridSlicer.Origin = [1.565198004245758, -18.485950469970703, -22.37162494659424]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from aortic_no_partition_384_facespvd
aortic_no_partition_384_facespvdDisplay = Show(aortic_no_partition_384_facespvd, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
aortic_no_partition_384_facespvdDisplay.Representation = 'Surface'
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
slice3Display.SetScaleArray = ['POINTS', 'annulus_normal_projected']
slice3Display.ScaleTransferFunction = 'PiecewiseFunction'
slice3Display.OpacityArray = ['POINTS', 'annulus_normal_projected']
slice3Display.OpacityTransferFunction = 'PiecewiseFunction'
slice3Display.DataAxesGrid = 'GridAxesRepresentation'
slice3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice3Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice3Display.ScaleTransferFunction.Points = [-8.173152896346384, 0.2589285671710968, 0.5, 0.0, 141.10649334402115, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice3Display.OpacityTransferFunction.Points = [-8.173152896346384, 0.2589285671710968, 0.5, 0.0, 141.10649334402115, 1.0, 0.5, 0.0]

# show data from clip1
clip1Display = Show(clip1, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
# clip1Display.SelectTCoordArray = 'None'
# clip1Display.SelectNormalArray = 'None'
# clip1Display.SelectTangentArray = 'None'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 0.7066963195800782
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 0.035334815979003904
clip1Display.SetScaleArray = [None, '']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = [None, '']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 0.10680600835959692
# clip1Display.OpacityArrayName = [None, '']

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip1Display.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [0.0, 0.2589285671710968, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'annulus_normal_projected'
annulus_normal_projectedPWF = GetOpacityTransferFunction('annulus_normal_projected')
annulus_normal_projectedPWF.Points = [-300.0, 0.2589285671710968, 0.5, 0.0, 300.0, 1.0, 0.5, 0.0]
annulus_normal_projectedPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# restore active source
SetActiveSource(slice3)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    
    tk = GetTimeKeeper()
    timesteps = tk.TimestepValues
    numTimesteps = len(timesteps)

    if len(sys.argv) >= 2:
        basename = sys.argv[1]
    else: 
        basename = 'frames'    
        basename += '_top_view'

    if 'top_view' not in basename:
        basename += '_top_view'        

    if len(sys.argv) >= 4:
        nprocs = int(sys.argv[2])
        proc_num = int(sys.argv[3])
    else: 
        print("using default proc_num 0, nprocs = 1")
        proc_num = 0
        nprocs = 1

    scaling = 1 
    res = (scaling*1180, scaling*1080)

    animationScene1 = GetAnimationScene()
    # animationScene1.GoToFirst()


    #frame_range = [1309,1341,1365,1427]

    # frame_range = [1365]
    run_specific = False 
    if run_specific:
        frame_range = [241,406, 723, 888]
    else: 
        frame_range = range(len(timesteps)-1)

    # run first iteration, then stop and start over 
    for frame in range(len(timesteps)):
        if (frame % nprocs) == proc_num:
            animationScene1.AnimationTime = timesteps[frame]

            Render()
            SaveScreenshot(basename + str(frame).zfill(4) + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)
            break 

    # for frame in range(len(timesteps)-1):
    #     if (frame % nprocs) == proc_num:
    #         animationScene1.AnimationTime = timesteps[frame]

    #         Render()
    #         SaveScreenshot(basename + str(frame).zfill(4) + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)

    for frame in frame_range:
        if (frame % nprocs) == proc_num:
            animationScene1.AnimationTime = timesteps[frame]

            Render()
            #SaveScreenshot(basename + str(frame).zfill(4) + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)
            SaveScreenshot(basename + str(frame).zfill(4) + '.jpeg', viewOrLayout=layout1, ImageResolution=res, Quality=100, Progressive=1)




