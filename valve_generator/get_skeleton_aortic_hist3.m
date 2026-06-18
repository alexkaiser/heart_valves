function skeleton = get_skeleton_aortic_hist3()


from_mat_file = true; 


if from_mat_file

    % note that these are already in sim coordinates 

    fname = '/Users/alex/Dropbox/stanford/research_stanford/annulus_morph_2025/meshes/ring_points_hist_3_processed.mat'
    load(fname, "ring_pts_LR_cusp_sim_coords_mm", "ring_pts_Non_cusp_sim_coords_mm")
    
    ring_pts_LR_cusp_sim_coords  = 0.1 * ring_pts_LR_cusp_sim_coords_mm';
    ring_pts_Non_cusp_sim_coords = 0.1 * ring_pts_Non_cusp_sim_coords_mm';

    % these are coming in opposite order
    ring_pts_Non_cusp_sim_coords = fliplr(ring_pts_Non_cusp_sim_coords);

    [~, n_pts_LR] = size(ring_pts_LR_cusp_sim_coords); 
    [~, n_pts_Non] = size(ring_pts_Non_cusp_sim_coords);


else 

    % raw data here
    % high_comm_neo_adjust
    
    % scan is in mm
    ring_pts_raw_LR_cusp = 0.1 * [ 
        -18.561134338378903, -162.1835479736328, -161.6475372314453,
        -16.709476470947266, -161.3505859375, -163.36756896972656,
        -14.98015308380127, -159.37615966796875, -164.9232940673828,
        -13.159804344177246, -158.35263061523438, -166.87225341796875,
        -10.745355606079102, -158.41268920898438, -169.69651794433594,
        -6.853701114654541, -159.12757873535156, -171.39987182617188,
        -2.7679173946380615, -161.17489624023438, -172.51446533203125,
        0.8725491166114807, -164.6865692138672, -173.50421142578125,
        0.19005458056926727, -170.52297973632813, -175.7110137939453,
        -2.0231025218963623, -175.45574951171875, -178.04913330078125,
        -5.795762647116474, -178.78905465817175, -180.17210184751303,
        -9.56268454474062, -181.6100463559652, -180.93456165010247,
        -15.33113956451416, -183.05392456054688, -179.898681640625,
        -19.65570831298828, -183.8021697998047, -177.125,
        -22.135957717895508, -184.1526641845703, -173.285888671875,
        -24.248830795288086, -184.28268432617188, -169.50759887695313,
        ]';

    ring_pts_raw_Non_cusp = 0.1 * [
        -18.561134338378903, -162.1835479736328, -161.6475372314453,
        -17.20952606201172, -161.2138671875, -163.71502685546875,
        -16.027801513671875, -159.2170867919922, -165.73770141601565,
        -15.225613581057996, -158.20417933741552, -168.30629406741485,
        -14.807823181152344, -157.79409790039065, -170.8694610595703,
        -15.32440948486328, -157.5122833251953, -177.5084991455078,
        -18.72540855407715, -159.20443725585938, -184.58387756347656,
        -22.71842384338379, -163.48788452148438, -189.0548858642578,
        -22.19114112854004, -170.0157928466797, -189.7144317626953,
        -20.328929901123047, -175.75787353515625, -188.66978454589844,
        -19.336400985717773, -180.0279998779297, -183.9769744873047,
        -21.324996948242188, -183.1256103515625, -178.061767578125,
        -22.989686965942383, -183.84304809570313, -173.56512451171875,
        -24.250965118408203, -184.2740478515625, -169.54635620117188,
        ]';    
    

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
    
    save_raw_pts = false; 
    if save_raw_pts
        
        ring_pts_LR_cusp_sim_coords_mm = 10 * ring_pts_LR_cusp_sim_coords'; 
        ring_pts_Non_cusp_sim_coords_mm = 10 * ring_pts_Non_cusp_sim_coords'; 
    
        save('ring_points_hist_3.mat', 'ring_pts_LR_cusp_sim_coords_mm', 'ring_pts_Non_cusp_sim_coords_mm')
    
    end 

end 


% % ruler test sim coordinates back to scan 
% paraview_ruler_sim_coords_mm = [-4.326762676239014, -59.78771209716797, -220.92550659179688
%                                 -14.859077453613281, -66.8626937866211, -240.3558349609375]';
% 
% % R is orthogonal so R' is inverse 
% paraview_ruler_scan_coords_mm = [R' * paraview_ruler_sim_coords_mm]'; 
% 
% paraview_ruler_sim_coords_mm_redone = [-4.353950023651123, -59.81804656982422, -220.897216796875
%                                        -14.537799835205078, -68.7021255493164, -239.96475219726562]';
% 
% paraview_ruler_scan_coords_mm_redone = [R' * paraview_ruler_sim_coords_mm]';
% 
% 
% paraview_ruler_scan_coords_mm_ruler2 = [-4.145022392272949, -59.74376678466797, -221.0463104248047
%                                        -15.036727905273438, -68.81126403808594, -239.31927490234375]'; 
% 
% paraview_ruler_sim_coords_mm_ruler2 = [R' * paraview_ruler_scan_coords_mm_ruler2]';



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
% artifact of other scaling 
r_stj_temp = skeleton.r; % 1.67 / 2;
% r_temp = 2.3 / 2; % vbr radius
hc = 1.0 * r_stj_temp; % 0.5 * r_stj; 
h1 = 0.87 * 2 * r_stj_temp - hc; % 1.4 * r_stj - hc; 

skeleton.height_min_comm = h1; 

skeleton.r_stj = 2.31/2; 


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








