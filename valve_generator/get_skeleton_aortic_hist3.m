function skeleton = get_skeleton_aortic_hist3()

low_comms = true; 

if low_comms
    % scan is in mm
    ring_pts_raw_LR_cusp = 0.1 * [ 
        -14.31479556638554, -158.5828612844842, -166.14133402684905,
        -12.81344985961914, -158.18797302246094, -167.7631072998047,
        -10.745355606079102, -158.41268920898438, -169.69651794433594,
        -7.073822021484375, -158.87242126464844, -171.82298278808594,
        -3.2711124420166016, -160.6879425048828, -172.90689086914063,
        0.8725491166114807, -164.6865692138672, -173.50421142578125,
        0.19005458056926727, -170.52297973632813, -175.7110137939453,
        -2.0231025218963623, -175.45574951171875, -178.04913330078125,
        -5.795762647116474, -178.78905465817175, -180.17210184751303,
        -9.575104713439941, -181.5930938720703, -180.76571655273438,
        -13.996484756469727, -182.8255157470703, -180.2509002685547,
        -19.65570831298828, -183.8021697998047, -177.125,
        -21.88850212097168, -183.95953369140625, -174.61090087890625,]';

    ring_pts_raw_Non_cusp = 0.1 * [
        -14.392179489135742, -158.59210205078125, -166.1156463623047,
        -13.665631294250488, -157.94525146484375, -168.35780334472656,
        -13.461244583129883, -157.824951171875, -171.11607360839844,
        -13.600553512573242, -157.50714111328125, -174.67413330078125,
        -14.667784690856934, -157.43319702148438, -179.31419372558594,
        -18.72540855407715, -159.20443725585938, -184.58387756347656,
        -23.097063064575195, -163.47999572753906, -188.67919921875,
        -22.47825050354004, -169.8924102783203, -189.45077514648438,
        -20.889745712280273, -175.38388061523438, -188.02960205078125,
        -19.336400985717773, -180.0279998779297, -183.9769744873047,
        -20.845043182373047, -183.25579833984375, -178.1044464111328,
        -21.88850212097168, -183.95953369140625, -174.61090087890625, 
        ]';

else 
    % scan is in mm
    ring_pts_raw_LR_cusp = 0.1 * [ 
    -16.983776092529297, -161.2555389404297, -163.5882568359375,
    -14.98015308380127, -159.37615966796875, -164.9232940673828,
    -13.159804344177246, -158.35263061523438, -166.87225341796875,
    -10.745355606079102, -158.41268920898438, -169.69651794433594,
    -6.853701114654541, -159.12757873535156, -171.39987182617188,
    -2.7679173946380615, -161.17489624023438, -172.51446533203125,
    0.8725491166114807, -164.6865692138672, -173.50421142578125,
    0.19005458056926727, -170.52297973632813, -175.7110137939453,
    -2.0231025218963623, -175.45574951171875, -178.04913330078125,
    -5.795762647116474, -178.78905465817175, -180.17210184751306,
    -9.56268454474062, -181.6100463559652, -180.93456165010247,
    -15.33113956451416, -183.05392456054688, -179.898681640625,
    -19.65570831298828, -183.8021697998047, -177.125,
    -23.952800750732425, -184.1968994140625, -170.5049285888672,]';
    
    ring_pts_raw_Non_cusp = 0.1 * [
    -17.160755157470703, -161.1587371826172, -163.78045654296875,
    -15.77912712097168, -158.7215576171875, -166.4163818359375,
    -14.807823181152344, -157.79409790039063, -170.8694610595703,
    -15.324409484863281, -157.5122833251953, -177.5084991455078,
    -18.72540855407715, -159.20443725585938, -184.58387756347656,
    -22.71842384338379, -163.48788452148438, -189.0548858642578,
    -22.19114112854004, -170.0157928466797, -189.7144317626953,
    -20.328929901123047, -175.75787353515625, -188.66978454589844,
    -19.336400985717773, -180.0279998779297, -183.9769744873047,
    -21.443647384643555, -183.0331573486328, -178.16384887695313,
    -23.968223571777344, -184.16403198242188, -170.601806640625,
    ]';
end 

% these happen to be in right to left order from segmentation 
ring_pts_LR_cusp = fliplr(ring_pts_raw_LR_cusp);
[~, n_pts_LR] = size(ring_pts_LR_cusp); 

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
% rot_1_x = deg2rad( 34);
% rot_2_y = deg2rad(-15);
% rot_3_z = deg2rad(-45);
% 
% R = rotation_matrix_z(rot_3_z) * rotation_matrix_y(rot_2_y) * rotation_matrix_x(rot_1_x); 

rot_1_x = deg2rad(34);
rot_2_y = deg2rad(15); % required sign swap for this component only 
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
[ring_pts_LR_cusp_model_coords, ~, ~, R] = coordinate_transformation_vertices(ring_pts_LR_cusp_sim_coords, [], R_0_temp, T_0_temp, apply_inverse, normal, midpoint_comms);
[ring_pts_Non_cusp_model_coords, ~, ~, R] = coordinate_transformation_vertices(ring_pts_Non_cusp_sim_coords, [], R_0_temp, T_0_temp, apply_inverse, normal, midpoint_comms);

% put comm_2  to the y axis     

R_extra_z = rotation_matrix_z(pi - atan2(ring_pts_LR_cusp_model_coords(2,end), ring_pts_LR_cusp_model_coords(1,end)));

ring_pts_LR_cusp_model_coords  = R_extra_z * ring_pts_LR_cusp_model_coords; 
ring_pts_Non_cusp_model_coords = R_extra_z * ring_pts_Non_cusp_model_coords; 

% skeleton.normal = midpoint_comms; 
% skeleton.centroid = midpoint_comms;
% 
% skeleton.initial_rotation_aortic = R_extra_z; 
% skeleton.initial_translation_aortic = T_0_temp; 

% this function applied at end 
% inverse of sim to model coordinate map 
skeleton.inverse_transformation_initial_condition = @(x) R * R_extra_z' * x + midpoint_comms; 

ring_pts_LR_cusp_sim_coords_transformed = skeleton.inverse_transformation_initial_condition(ring_pts_LR_cusp_model_coords);
ring_pts_Non_cusp_sim_coords_transformed = skeleton.inverse_transformation_initial_condition(ring_pts_Non_cusp_model_coords);

skeleton.r = norm(comm_1 - comm_2)/2;
skeleton.r_commissure = skeleton.r; 
skeleton.ring_pts_LR_cusp_model_coords = ring_pts_LR_cusp_model_coords; 
skeleton.ring_pts_Non_cusp_model_coords = ring_pts_Non_cusp_model_coords; 

% scale comm height proportional to STJ radius 
r_stj = skeleton.r; % 1.67 / 2;
% r_temp = 2.3 / 2; % vbr radius
hc = 1.0 * r_stj; % 0.5 * r_stj; 
h1 = 0.87 * 2 * r_stj - hc; % 1.4 * r_stj - hc; 

skeleton.height_min_comm = h1; 


debug_plots = true;
if debug_plots
    
    figure; 
    hold on; 

    % plot3(ring_pts_raw_LR_cusp(1,:), ring_pts_raw_LR_cusp(2,:), ring_pts_raw_LR_cusp(3,:),'*-');
    % plot3(ring_pts_raw_Non_cusp(1,:), ring_pts_raw_Non_cusp(2,:), ring_pts_raw_Non_cusp(3,:),'*-');

    plot3(ring_pts_LR_cusp_sim_coords(1,:), ring_pts_LR_cusp_sim_coords(2,:), ring_pts_LR_cusp_sim_coords(3,:),'*-');
    plot3(ring_pts_Non_cusp_sim_coords(1,:), ring_pts_Non_cusp_sim_coords(2,:), ring_pts_Non_cusp_sim_coords(3,:),'*-');

    % plot3(ring_pts_LR_cusp_model_coords(1,:), ring_pts_LR_cusp_model_coords(2,:), ring_pts_LR_cusp_model_coords(3,:),'*-');
    % plot3(ring_pts_Non_cusp_model_coords(1,:), ring_pts_Non_cusp_model_coords(2,:), ring_pts_Non_cusp_model_coords(3,:),'*-');

    plot3(ring_pts_LR_cusp_sim_coords_transformed(1,:), ring_pts_LR_cusp_sim_coords_transformed(2,:), ring_pts_LR_cusp_sim_coords_transformed(3,:),'o-');
    plot3(ring_pts_Non_cusp_sim_coords_transformed(1,:), ring_pts_Non_cusp_sim_coords_transformed(2,:), ring_pts_Non_cusp_sim_coords_transformed(3,:),'o-');
    axis equal
end 








