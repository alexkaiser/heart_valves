function skeleton = get_skeleton_aortic_hist3()

% scan is in mm, con
ring_pts_raw_LR_cusp = 0.1 * [ 
-16.983776092529297, -161.2555389404297, -163.5882568359375,
-13.077391624450684, -157.9479522705078, -169.11514282226563,
-7.8202805519104, -158.03750610351563, -173.54421997070313,
-2.486410140991211, -160.36691284179688, -175.14462280273438,
1.3799179792404175, -165.36549377441406, -175.52175903320313,
-1.3226488828659058, -174.33387756347656, -178.4099578857422,
-9.562684544740621, -181.61004635596518, -180.93456165010247,
-15.33113956451416, -183.05392456054688, -179.898681640625,
-19.65570831298828, -183.8021697998047, -177.125,
-23.952800750732422, -184.1968994140625, -170.5049285888672]';

% these happen to be in right to left order from segmentation 
ring_pts_LR_cusp = fliplr(ring_pts_raw_LR_cusp);
[~, n_pts_LR] = size(ring_pts_LR_cusp); 

ring_pts_raw_Non_cusp = 0.1 * [
-17.160755157470703, -161.1587371826172, -163.78045654296875,
-15.720479965209961, -157.78114318847656, -170.81082153320313,
-15.656978607177734, -157.5762939453125, -177.89515686035156,
-18.72540855407715, -159.20443725585938, -184.58387756347656,
-23.016414642333984, -163.69322204589844, -188.88502502441406,
-22.19114112854004, -170.0157928466797, -189.7144317626953,
-20.328929901123047, -175.75787353515625, -188.66978454589844,
-19.336400985717773, -180.0279998779297, -183.9769744873047,
-21.443647384643555, -183.0331573486328, -178.16384887695313,
-23.968223571777344, -184.16403198242188, -170.601806640625]';

ring_pts_Non_cusp = ring_pts_raw_Non_cusp; 
[~, n_pts_Non] = size(ring_pts_Non_cusp);

% sets the commissure points to be exactly equal on two leafelets 
patch_commissure = true; 
if patch_commissure
    ring_pts_Non_cusp(:,1) = ring_pts_raw_LR_cusp(:,1);
    ring_pts_Non_cusp(:,n_pts_Non) = ring_pts_raw_LR_cusp(:,n_pts_LR);
end 



% rotations made in paraview for axis alignment 
% scan coordinates to simulation coordinates 
rot_1_x = deg2rad( 34);
rot_2_y = deg2rad(-15);
rot_3_z = deg2rad(-45);

R = rotation_matrix_z(rot_3_z) * rotation_matrix_y(rot_2_y) * rotation_matrix_x(rot_1_x); 

% simulation coordinates 
ring_pts_LR_cusp_sim_coords = R * ring_pts_LR_cusp; 
ring_pts_Non_cusp_sim_coords = R * ring_pts_Non_cusp; 

% convert from sim to model construction coordinates 

% put midpoint of the commissures at origin 
comm_1 = ring_pts_LR_cusp_sim_coords(:,1); 
comm_2 = ring_pts_LR_cusp_sim_coords(:,n_pts_LR); 
midpoint_comms = 0.5 * (comm_1 + comm_2);

% approx midpoint of nadir of attachment 
% use this to find a third point in the cusp plane 
nadir_LR = ring_pts_LR_cusp_sim_coords(:,floor(n_pts_LR/2));
nadir_Non = ring_pts_Non_cusp_sim_coords(:,floor(n_pts_Non/2));
midpoint_nadir = 0.5 * (nadir_LR + nadir_Non); 

% vevtor from midpoint of vbr to midpoint of commissural plane 
nadir_LR_to_commissure_plane = midpoint_comms - midpoint_nadir; 

% project a vector 
nadir_LR_translated = nadir_LR + nadir_LR_to_commissure_plane; 


basis_nadir = nadir_LR_translated - midpoint_comms; 
basis_commissure = ring_pts_LR_cusp_sim_coords(:,n_pts_LR) - midpoint_comms; 

normal = cross(basis_nadir, basis_commissure); 
normal = normal / norm(normal);

R_0_temp = eye(3);
T_0_temp = zeros(3,1);
% fwd transformation goes 
apply_inverse = true; 
[ring_pts_LR_cusp_model_coords, normal, centroid, R] = coordinate_transformation_vertices(ring_pts_LR_cusp_sim_coords, [], R_0_temp, T_0_temp, apply_inverse, normal, midpoint_comms);
[ring_pts_Non_cusp_model_coords, normal, centroid, R] = coordinate_transformation_vertices(ring_pts_Non_cusp_sim_coords, [], R_0_temp, T_0_temp, apply_inverse, normal, midpoint_comms);

% put comm_2  to the y axis     
R_extra_z = rotation_matrix_z(pi - atan2(ring_pts_LR_cusp_model_coords(2,end), ring_pts_LR_cusp_model_coords(1,end)));

ring_pts_LR_cusp_model_coords  = R_extra_z * ring_pts_LR_cusp_model_coords; 
ring_pts_Non_cusp_model_coords = R_extra_z * ring_pts_Non_cusp_model_coords; 


skeleton.centroid = midpoint_comms; 
skeleton.r = norm(comm_1 - comm_2)/2;
skeleton.r_commissure = skeleton.r; 
skeleton.valve_ring_pts = ring_pts_LR_cusp_model_coords; 

% scale comm height proportional to STJ radius 
r_stj = skeleton.r; % 1.67 / 2;
% r_temp = 2.3 / 2; % vbr radius
hc = 0.5 * r_stj; 
h1 = 0.87 * 2 * r_stj - hc; % 1.4 * r_stj - hc; 

skeleton.height_min_comm = h1; 


debug_plots = true;
if debug_plots
    plot3(ring_pts_raw_LR_cusp(1,:), ring_pts_raw_LR_cusp(2,:), ring_pts_raw_LR_cusp(3,:),'*-');
    hold on 
    plot3(ring_pts_raw_Non_cusp(1,:), ring_pts_raw_Non_cusp(2,:), ring_pts_raw_Non_cusp(3,:),'*-');
    plot3(ring_pts_LR_cusp_sim_coords(1,:), ring_pts_LR_cusp_sim_coords(2,:), ring_pts_LR_cusp_sim_coords(3,:),'*-');
    plot3(ring_pts_Non_cusp_sim_coords(1,:), ring_pts_Non_cusp_sim_coords(2,:), ring_pts_Non_cusp_sim_coords(3,:),'*-');

    plot3(ring_pts_LR_cusp_model_coords(1,:), ring_pts_LR_cusp_model_coords(2,:), ring_pts_LR_cusp_model_coords(3,:),'*-');
    plot3(ring_pts_Non_cusp_model_coords(1,:), ring_pts_Non_cusp_model_coords(2,:), ring_pts_Non_cusp_model_coords(3,:),'*-');
    axis equal
end 








