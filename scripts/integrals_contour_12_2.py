from __future__ import print_function
import os 
import math 
# state file generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'


if __name__ == '__main__':

    if len(sys.argv) >= 2:
        contour_name = sys.argv[1]

        valid_names = ['annulus_normal_projected', 'contour_3', 'contour_6', 'contour_9', 'contour_12', 'contour_21'] 

        # valid contour names 
        if contour_name not in valid_names:
            raise ValueError('must provide countour_name, valid names: annulus_normal_projected, contour_3, contour_6, contour_9, contour_12, contour_21')

        if contour_name == 'annulus_normal_projected':
            origin = [2.447601610336508, -18.376013668463226, -25.219513937752243]
            normal = [-0.644118902159037, 0.0079312092989, 0.764884263010094]
        elif contour_name == 'contour_3':
            origin = [0.2803383547117727, -18.3479701228845, -18.8436930740923]
            normal = [0.277454395233418, 0.0616022174821306, 0.958761818892963]
        elif contour_name == 'contour_6':
            origin = [0.0488354041365377, -18.3884646217316, -20.8624742881165]
            normal = [-0.0974948538117627, -0.028568296191069, 0.994825917401111]
        elif contour_name == 'contour_9':
            origin = [0.5826182201137723, -18.3351820981491, -22.7279221637119]
            normal = [-0.369605687499951, 0.0302064524432408, 0.928697585868771]
        elif contour_name == 'contour_12':
            origin = [1.1135387247389723, -18.4115922051707, -23.8033429814558]
            normal = [-0.546829510850715, 0.0721783996060898, 0.834126947588358]
        elif contour_name == 'contour_21':
            # LVOT, sub annular 
            origin = [3.044823118615524, -18.4669463809972, -25.4758572325839]
            normal = [-0.871486988343708, -0.0231008494055719, 0.489874249072508]

    else: 
        raise ValueError('must provide countour_name, valid names: annulus_normal_projected, contour_3, contour_6, contour_9, contour_12')
    



    paraview.simple._DisableFirstRenderCameraReset()

    # ----------------------------------------------------------------
    # setup views used in the visualization
    # ----------------------------------------------------------------

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # Create a new 'Render View'
    renderView2 = CreateView('RenderView')
    renderView2.ViewSize = [776, 1081]
    renderView2.AxesGrid = 'GridAxes3DActor'
    renderView2.OrientationAxesVisibility = 0
    renderView2.CenterOfRotation = [1.1908643022179604, -18.412199020385742, -23.76491069793701]
    renderView2.StereoType = 'Crystal Eyes'
    renderView2.CameraPosition = [8.37062857270826, -17.5255319461868, -31.232425539351397]
    renderView2.CameraFocalPoint = [1.1908643022179604, -18.412199020385742, -23.76491069793701]
    renderView2.CameraViewUp = [-0.6908210277780568, -0.21636844906573194, -0.6898920218617977]
    renderView2.CameraFocalDisk = 1.0
    renderView2.CameraParallelScale = 2.735736701526508
    renderView2.BackEnd = 'OSPRay raycaster'
    renderView2.OSPRayMaterialLibrary = materialLibrary1

    # Create a new 'SpreadSheet View'
    spreadSheetView1 = CreateView('SpreadSheetView')
    spreadSheetView1.ColumnToSort = ''
    spreadSheetView1.BlockSize = 1024

    SetActiveView(None)

    # ----------------------------------------------------------------
    # setup view layouts
    # ----------------------------------------------------------------

    # create new layout object 'Layout #1'
    layout1 = CreateLayout(name='Layout #1')
    layout1.SplitHorizontal(0, 0.472511)
    layout1.AssignView(1, renderView2)
    layout1.AssignView(2, spreadSheetView1)
    layout1.SetSize(1177, 1081)

    # ----------------------------------------------------------------
    # restore active view
    SetActiveView(spreadSheetView1)
    # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    # setup the data processing pipelines
    # ----------------------------------------------------------------

    # create a new 'Annotate Time'
    annotateTime1 = AnnotateTime(registrationName='AnnotateTime1')
    annotateTime1.Format = 't = %.3f s'

    # create a new 'PVD Reader'
    # eulerian_varspvd = PVDReader(registrationName='eulerian_vars.pvd', FileName='eulerian_vars.pvd')
    # eulerian_varspvd.CellArrays = ['U']
    eulerian_varspvd = PVDReader(FileName='eulerian_vars_restricted_points.pvd')
    eulerian_varspvd.PointArrays = ['U','P']

    # create a new 'PVD Reader'
    # aorta_384pvd = PVDReader(registrationName='aorta_384.pvd', FileName='aorta_384.pvd')


    # create a new 'Text'
    # text1 = Text(registrationName='Text1')
    # text1.Text = 'vertical component of velocity'

    # # create a new 'PVD Reader'
    # aortic_no_partition_384pvd = PVDReader(registrationName='aortic_no_partition_384.pvd', FileName='aortic_no_partition_384.pvd')



    # create a new 'Calculator'
    calculator2 = Calculator(registrationName='Calculator2', Input=eulerian_varspvd)
    calculator2.ResultArrayName = 'u_dot_n'
    calculator2.Function = str(normal[0]) + ' * U_X + ' + str(normal[1]) + ' * U_Y + ' + str(normal[2]) + ' * U_Z'
    # calculator2.Function = '-0.546829510850715 * U_X + 0.0721783996060898 * U_Y + 0.834126947588358 * U_Z'

    # create a new 'Slice'
    slice4 = Slice(registrationName='Slice4', Input=calculator2)
    slice4.SliceType = 'Plane'
    slice4.Triangulatetheslice = 0
    slice4.HyperTreeGridSlicer = 'Plane'
    slice4.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice4.SliceType.Origin = origin #[1.1135387247389723, -18.4115922051707, -23.8033429814558]
    slice4.SliceType.Normal = normal #[-0.546829510850715, 0.0721783996060898, 0.834126947588358]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice4.HyperTreeGridSlicer.Origin = origin #[1.5904362797737122, -18.485923767089844, -22.346393585205078]

    # create a new 'Calculator'
    calculator1_Q_A = Calculator(registrationName='Calculator1_Q_A', Input=slice4)
    calculator1_Q_A.ResultArrayName = 'One_scalar'
    calculator1_Q_A.Function = '1'

    # create a new 'Slice'
    # slice9 = Slice(registrationName='Slice9', Input=eulerian_varspvd)
    # slice9.SliceType = 'Plane'
    # slice9.HyperTreeGridSlicer = 'Plane'
    # slice9.SliceOffsetValues = [0.0]

    # # init the 'Plane' selected for 'SliceType'
    # slice9.SliceType.Origin = [1.0154999494552612, -18.4862003326416, -22.669400215148926]
    # slice9.SliceType.Normal = [0.0, 1.0, 0.0]

    # # init the 'Plane' selected for 'HyperTreeGridSlicer'
    # slice9.HyperTreeGridSlicer.Origin = [1.0154999494552612, -18.4862003326416, -22.669400215148926]

    # create a new 'Calculator'
    calculator1_normal = Calculator(registrationName='Calculator1_normal', Input=slice4)
    calculator1_normal.ResultArrayName = 'normal'
    calculator1_normal.Function = str(normal[0]) + ' * iHat + ' + str(normal[1]) + ' * jHat + ' + str(normal[2]) + ' * kHat'
    # calculator1_normal.Function = '-0.546829510850715 * iHat + 0.0721783996060898 * jHat + 0.834126947588358 * kHat'

    # create a new 'Calculator'
    calculator1_u_tangential = Calculator(registrationName='Calculator1_u_tangential', Input=calculator1_normal)
    calculator1_u_tangential.ResultArrayName = 'u_tangential'
    calculator1_u_tangential.Function = 'U - u_dot_n*normal'

    # create a new 'Glyph'
    glyph1 = Glyph(registrationName='Glyph1', Input=calculator1_u_tangential,
        GlyphType='Arrow')
    glyph1.OrientationArray = ['POINTS', 'u_tangential']
    glyph1.ScaleArray = ['POINTS', 'u_tangential']
    glyph1.ScaleFactor = 0.01
    glyph1.GlyphTransform = 'Transform2'
    glyph1.GlyphMode = 'All Points'

    # create a new 'Calculator'
    calculator1_norm_sq_u_tangential = Calculator(registrationName='Calculator1_norm_sq_u_tangential', Input=calculator1_u_tangential)
    calculator1_norm_sq_u_tangential.ResultArrayName = 'norm_sq_u_tangential'
    calculator1_norm_sq_u_tangential.Function = 'mag(u_tangential)^2'

    # create a new 'Integrate Variables'
    integrateVariables3 = IntegrateVariables(registrationName='IntegrateVariables3', Input=calculator1_norm_sq_u_tangential)

    # create a new 'Slice'
    # slice2 = Slice(registrationName='Slice2', Input=aorta_384pvd)
    # slice2.SliceType = 'Plane'
    # slice2.Triangulatetheslice = 0
    # slice2.HyperTreeGridSlicer = 'Plane'
    # slice2.SliceOffsetValues = [0.0]

    # # init the 'Plane' selected for 'SliceType'
    # slice2.SliceType.Origin = [1.565198004245758, -18.485950469970703, -22.37162494659424]
    # slice2.SliceType.Normal = [0.0, 1.0, 0.0]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    # slice2.HyperTreeGridSlicer.Origin = [1.565198004245758, -18.485950469970703, -22.37162494659424]

    # create a new 'Calculator'
    calculator1_u_dot_n_sq = Calculator(registrationName='Calculator1_u_dot_n_sq', Input=slice4)
    calculator1_u_dot_n_sq.ResultArrayName = 'u_dot_n_sq'
    calculator1_u_dot_n_sq.Function = 'u_dot_n^2'

    # create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=calculator1_Q_A)

    # create a new 'Integrate Variables'
    integrateVariables2 = IntegrateVariables(registrationName='IntegrateVariables2', Input=calculator1_u_dot_n_sq)

    calculator1_positive_flow = Calculator(registrationName='Calculator1_positive_flow', Input=slice4)
    calculator1_positive_flow.ResultArrayName = 'indicator_u_dot_n_positive'
    calculator1_positive_flow.Function = '0.5 * (u_dot_n/abs(u_dot_n) + 1.0)'

    # create a new 'Integrate Variables'
    integrateVariables4 = IntegrateVariables(registrationName='IntegrateVariables4', Input=calculator1_positive_flow)


    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView2'
    # ----------------------------------------------------------------

    # show data from calculator1_norm_sq_u_tangential
    calculator1_norm_sq_u_tangentialDisplay = Show(calculator1_norm_sq_u_tangential, renderView2, 'GeometryRepresentation')

    # get color transfer function/color map for 'norm_sq_u_tangential'
    norm_sq_u_tangentialLUT = GetColorTransferFunction('norm_sq_u_tangential')
    norm_sq_u_tangentialLUT.RGBPoints = [0.002455722592898142, 0.231373, 0.298039, 0.752941, 5476.125002854154, 0.865003, 0.865003, 0.865003, 10952.247549985717, 0.705882, 0.0156863, 0.14902]
    norm_sq_u_tangentialLUT.ScalarRangeInitialized = 1.0

    # trace defaults for the display properties.
    calculator1_norm_sq_u_tangentialDisplay.Representation = 'Surface'
    calculator1_norm_sq_u_tangentialDisplay.ColorArrayName = ['POINTS', 'norm_sq_u_tangential']
    calculator1_norm_sq_u_tangentialDisplay.LookupTable = norm_sq_u_tangentialLUT
    calculator1_norm_sq_u_tangentialDisplay.SelectTCoordArray = 'None'
    calculator1_norm_sq_u_tangentialDisplay.SelectNormalArray = 'None'
    calculator1_norm_sq_u_tangentialDisplay.SelectTangentArray = 'None'
    calculator1_norm_sq_u_tangentialDisplay.OSPRayScaleArray = 'norm_sq_u_tangential'
    calculator1_norm_sq_u_tangentialDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1_norm_sq_u_tangentialDisplay.SelectOrientationVectors = 'u_tangential'
    calculator1_norm_sq_u_tangentialDisplay.ScaleFactor = 0.27580108642578127
    calculator1_norm_sq_u_tangentialDisplay.SelectScaleArray = 'norm_sq_u_tangential'
    calculator1_norm_sq_u_tangentialDisplay.GlyphType = 'Arrow'
    calculator1_norm_sq_u_tangentialDisplay.GlyphTableIndexArray = 'norm_sq_u_tangential'
    calculator1_norm_sq_u_tangentialDisplay.GaussianRadius = 0.013790054321289063
    calculator1_norm_sq_u_tangentialDisplay.SetScaleArray = ['POINTS', 'norm_sq_u_tangential']
    calculator1_norm_sq_u_tangentialDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1_norm_sq_u_tangentialDisplay.OpacityArray = ['POINTS', 'norm_sq_u_tangential']
    calculator1_norm_sq_u_tangentialDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1_norm_sq_u_tangentialDisplay.DataAxesGrid = 'GridAxesRepresentation'
    calculator1_norm_sq_u_tangentialDisplay.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    calculator1_norm_sq_u_tangentialDisplay.OSPRayScaleFunction.Points = [-256.96728515625, 0.2589285671710968, 0.5, 0.0, 102.52591705322266, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1_norm_sq_u_tangentialDisplay.ScaleTransferFunction.Points = [0.002455722592898142, 0.2589285671710968, 0.5, 0.0, 10952.247549985717, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1_norm_sq_u_tangentialDisplay.OpacityTransferFunction.Points = [0.002455722592898142, 0.2589285671710968, 0.5, 0.0, 10952.247549985717, 1.0, 0.5, 0.0]

    # setup the color legend parameters for each legend in this view

    # get color legend/bar for norm_sq_u_tangentialLUT in view renderView2
    norm_sq_u_tangentialLUTColorBar = GetScalarBar(norm_sq_u_tangentialLUT, renderView2)
    norm_sq_u_tangentialLUTColorBar.Title = 'norm_sq_u_tangential'
    norm_sq_u_tangentialLUTColorBar.ComponentTitle = ''

    # set color bar visibility
    norm_sq_u_tangentialLUTColorBar.Visibility = 1

    # show color legend
    calculator1_norm_sq_u_tangentialDisplay.SetScalarBarVisibility(renderView2, True)

    # ----------------------------------------------------------------
    # setup the visualization in view 'spreadSheetView1'
    # ----------------------------------------------------------------

    # show data from integrateVariables3
    integrateVariables3Display = Show(integrateVariables3, spreadSheetView1, 'SpreadSheetRepresentation')

    # ----------------------------------------------------------------
    # setup color maps and opacity mapes used in the visualization
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------

    # get opacity transfer function/opacity map for 'norm_sq_u_tangential'
    norm_sq_u_tangentialPWF = GetOpacityTransferFunction('norm_sq_u_tangential')
    norm_sq_u_tangentialPWF.Points = [0.002455722592898142, 0.2589285671710968, 0.5, 0.0, 10952.247549985717, 1.0, 0.5, 0.0]
    norm_sq_u_tangentialPWF.ScalarRangeInitialized = 1

    # ----------------------------------------------------------------
    # restore active source
    SetActiveSource(integrateVariables3)
    # ----------------------------------------------------------------


    # basic example of getting the value out of the integral 
    # integrated_filter_1 = paraview.servermanager.Fetch(integrateVariables1)
    
    # num_points = integrated_filter_1.GetNumberOfPoints()

    # print("num_points = ", num_points)

    # print("integrated_filter_1 = ", integrated_filter_1) 
    # print("integrated_filter_1.GetPointData() = ", integrated_filter_1.GetPointData()) 
    # print("integrated_filter_1.GetPointData().GetArray('u_dot_n') = ", integrated_filter_1.GetPointData().GetArray('u_dot_n')) 

    # for pt in range(num_points): 
    #     print("integrated_filter_1.GetPointData().GetArray('u_dot_n').GetValue(pt) = ", integrated_filter_1.GetPointData().GetArray('u_dot_n').GetValue(pt))         

    # # written assuming one point each as result of integrant 
    # assert num_points == 1

    # integrated_filter_2 = paraview.servermanager.Fetch(integrateVariables2)    
    # num_points_2 = integrated_filter_2.GetNumberOfPoints()

    # print("num_points_2 = ", num_points_2)
    # print("integrated_filter_2.GetPointData() = ", integrated_filter_2.GetPointData()) 

    integrated_filter_3 = paraview.servermanager.Fetch(integrateVariables3)    
    print("integrated_filter_3.GetPointData() = ", integrated_filter_3.GetPointData()) 

    integrated_filter_4 = paraview.servermanager.Fetch(integrateVariables4)    
    print("integrated_filter_4.GetPointData() = ", integrated_filter_4.GetPointData()) 

    tk = GetTimeKeeper()
    timesteps = tk.TimestepValues

    #timesteps = tk.TimestepValues[0:10]

    numTimesteps = len(timesteps)    

    filename = "integral_metric_data_" + contour_name + ".m"
    f = open(filename, 'a')

    # write times 
    print('t = [', file=f, end='')
    for frame in range(len(timesteps)-1):
        print(timesteps[frame] , ', ', file=f, end='')
    print('];\n', file=f, end='')

    animationScene1 = GetAnimationScene()

    # # compute area 
    animationScene1.AnimationTime = timesteps[0]
    Show(integrateVariables1)
    integrated_filter_1 = paraview.servermanager.Fetch(integrateVariables1)
    A = integrated_filter_1.GetPointData().GetArray('One_scalar').GetValue(0)
    print("A = ",  A, ";\n", file=f, end='')
    f.flush()
    os.fsync(f)
    print('area passed')


    # # flux 
    print("Q = [", file=f, end='')
    for frame in range(len(timesteps)-1):
        animationScene1.AnimationTime = timesteps[frame]

        Show(integrateVariables1)
        integrated_filter_1 = paraview.servermanager.Fetch(integrateVariables1)

        Q = integrated_filter_1.GetPointData().GetArray('u_dot_n').GetValue(0)
        print(Q, ", ", file=f, end='')

        if (frame % math.floor(len(timesteps)/10)) == 0:
            print("frame ", frame, "compute Q")

    print('];\n', file=f, end='')
    f.flush()
    os.fsync(f)

    # # normal component, two norm squared 
    # print("int_u_dot_n_squared = [", file=f, end='')
    # for frame in range(len(timesteps)-1):
    #     animationScene1.AnimationTime = timesteps[frame]
    #     Show(integrateVariables2)
    #     integrated_filter_2 = paraview.servermanager.Fetch(integrateVariables2)

    #     int_u_dot_n_squared = integrated_filter_2.GetPointData().GetArray('u_dot_n_sq').GetValue(0)

    #     print(int_u_dot_n_squared, ", ", file=f, end='')

    #     if (frame % math.floor(len(timesteps)/10)) == 0:
    #         print("frame ", frame, "compute u_dot_n_squared")

    # print('];\n', file=f, end='')
    # f.flush()
    # os.fsync(f)


    # # norm_sq_u_tangential
    # print("int_norm_sq_u_tangential = [", file=f, end='')
    # for frame in range(len(timesteps)-1):
    #     animationScene1.AnimationTime = timesteps[frame]
    #     Show(integrateVariables3)
    #     integrated_filter_3 = paraview.servermanager.Fetch(integrateVariables3)

    #     int_norm_sq_u_tangential = integrated_filter_3.GetPointData().GetArray('norm_sq_u_tangential').GetValue(0)

    #     print(int_norm_sq_u_tangential, ", ", file=f, end='')

    #     if (frame % math.floor(len(timesteps)/10)) == 0:
    #         print("frame ", frame, "compute int_norm_sq_u_tangential")


    # print('];\n', file=f, end='')
    # f.flush()
    # os.fsync(f)



    # # indicator_u_dot_n_positive
    # print("int_indicator_u_dot_n_positive = [", file=f, end='')
    # for frame in range(len(timesteps)-1):
    #     animationScene1.AnimationTime = timesteps[frame]
    #     Show(integrateVariables4)
    #     integrated_filter_4 = paraview.servermanager.Fetch(integrateVariables4)

    #     int_indicator_u_dot_n_positive = integrated_filter_4.GetPointData().GetArray('indicator_u_dot_n_positive').GetValue(0)

    #     print(int_indicator_u_dot_n_positive, ", ", file=f, end='')

    #     if (frame % math.floor(len(timesteps)/10)) == 0:
    #         print("frame ", frame, "compute int_indicator_u_dot_n_positive")

    # print('];\n', file=f, end='')
    # f.flush()
    # os.fsync(f)

    # # flux 
    print("P = [", file=f, end='')
    for frame in range(len(timesteps)-1):
        animationScene1.AnimationTime = timesteps[frame]

        Show(integrateVariables1)
        integrated_filter_1 = paraview.servermanager.Fetch(integrateVariables1)

        P = integrated_filter_1.GetPointData().GetArray('P').GetValue(0)
        print(P, ", ", file=f, end='')

        if (frame % math.floor(len(timesteps)/10)) == 0:
            print("frame ", frame, "compute P")

    print('];\n', file=f, end='')
    f.flush()
    os.fsync(f)

    f.close()

    quit()
