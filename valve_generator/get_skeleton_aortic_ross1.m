function skeleton = get_skeleton_aortic_ross1()

aortic_annulus_right_raw = 0.1 * [
-0.46831756830215454, -169.6324462890625, -667.995361328125
1.655497670173645, -169.13372802734375, -672.6422119140625
4.604877948760986, -168.5965118408203, -677.3262939453125
9.914689064025879, -169.38131713867188, -680.3304443359375
15.234301567077637, -169.55740356445312, -679.2272338867188
19.462570190429688, -165.55258178710938, -676.4866333007812
20.77568817138672, -163.998046875, -672.1533203125
20.454574584960938, -162.8600311279297, -667.1326904296875
% 19.78349494934082, -162.34991455078125, -659.9765625
19.774280548095703, -162.3553466796875, -659.881591796875
]';

aortic_annulus_non_raw = 0.1 * [
5.236440658569336, -149.076416015625, -663.696044921875
5.8963446617126465, -149.24378967285156, -670.3046264648438
5.979079246520996, -149.17974853515625, -675.0145263671875
5.497542381286621, -150.6417694091797, -679.2730712890625
3.034792184829712, -154.01437377929688, -681.9052124023438
0.5804076194763184, -158.937744140625, -682.6820068359375
2.0220649242401123, -163.97756958007812, -681.3392944335938
2.1807055473327637, -166.66575622558594, -677.3714599609375
1.0360711812973022, -168.52215576171875, -672.9295654296875
-0.35125765204429626, -169.6897430419922, -668.0812377929688
]';

aortic_annulus_left_raw = 0.1 * [
19.857906341552734, -162.0392608642578, -659.920166015625
20.641239166259766, -162.15647888183594, -667.0173950195312
21.807863235473633, -160.075927734375, -672.738525390625
21.760543823242188, -154.13853454589844, -674.685546875
18.094959259033203, -149.1846923828125, -675.5745849609375
12.451701164245605, -147.70068359375, -675.5194091796875
8.696045875549316, -148.48374938964844, -673.3262939453125
7.0716753005981445, -148.86367797851562, -669.8989868164062
5.306053638458252, -149.06787109375, -663.6882934570312
]';


[~, n_pts_left ]  = size(aortic_annulus_left_raw); 
[~, n_pts_non  ]  = size(aortic_annulus_non_raw); 
[~, n_pts_right] = size(aortic_annulus_right_raw); 


aortic_annulus_left  = aortic_annulus_left_raw; 
aortic_annulus_non   = aortic_annulus_non_raw; 
aortic_annulus_right = aortic_annulus_right_raw; 


% sets the commissure points to be exactly equal on two leafelets 
patch_commissure = true; 
if patch_commissure    
    aortic_annulus_left (:,n_pts_left ) = aortic_annulus_non(:,1);
    aortic_annulus_non  (:,n_pts_non  ) = aortic_annulus_right(:,1);
    aortic_annulus_right(:,n_pts_right) = aortic_annulus_left(:,1);
end 


rot_1_x = deg2rad(22);
rot_2_y = deg2rad(31); % required sign swap for this component only 
rot_3_z = deg2rad(-38);
R = rotation_matrix_z(rot_3_z) * rotation_matrix_y(rot_2_y) * rotation_matrix_x(rot_1_x); 


% simulation coordinates 
aortic_annulus_left_sim_coords = R * aortic_annulus_left; 
aortic_annulus_non_sim_coords = R * aortic_annulus_non; 
aortic_annulus_right_sim_coords = R * aortic_annulus_right; 

comm_LN = aortic_annulus_left_sim_coords(:,n_pts_left); 
comm_NR = aortic_annulus_non_sim_coords (:,n_pts_non); 
comm_LR = aortic_annulus_right_sim_coords(:,n_pts_right); 
comms_all = [comm_LN, comm_NR, comm_LR]; 

midpoint_comms = mean(comms_all, 2); 

basis_1 = comm_LR - comm_NR; 
basis_2 = comm_LN - comm_NR; 

normal = cross(basis_1, basis_2); 
normal = normal / norm(normal);


R_0_temp = eye(3);
T_0_temp = zeros(3,1);
% fwd transformation goes 
apply_inverse = true; 
[aortic_annulus_left_model_coords , ~, ~, R] = coordinate_transformation_vertices(aortic_annulus_left_sim_coords,  [], R_0_temp, T_0_temp, apply_inverse, normal, midpoint_comms);
[aortic_annulus_non_model_coords  , ~, ~, R] = coordinate_transformation_vertices(aortic_annulus_non_sim_coords,   [], R_0_temp, T_0_temp, apply_inverse, normal, midpoint_comms);
[aortic_annulus_right_model_coords, ~, ~, R] = coordinate_transformation_vertices(aortic_annulus_right_sim_coords, [], R_0_temp, T_0_temp, apply_inverse, normal, midpoint_comms);



% this function applied at end 
% inverse of sim to model coordinate map 
skeleton.inverse_transformation_initial_condition = @(x) R * x + midpoint_comms; 

aortic_annulus_left_sim_coords_transformed  = skeleton.inverse_transformation_initial_condition(aortic_annulus_left_model_coords );
aortic_annulus_non_sim_coords_transformed   = skeleton.inverse_transformation_initial_condition(aortic_annulus_non_model_coords  );
aortic_annulus_right_sim_coords_transformed = skeleton.inverse_transformation_initial_condition(aortic_annulus_right_model_coords);


distances = sqrt(sum((comms_all - midpoint_comms).^2,2));

dist_LN = norm(comm_LN - midpoint_comms); 
dist_NR = norm(comm_NR - midpoint_comms); 
dist_LR = norm(comm_LR - midpoint_comms);

mean_radius = mean(distances);

skeleton.r = mean_radius; 
skeleton.r_commissure = skeleton.r; 

skeleton.aortic_annulus_left_model_coords = aortic_annulus_left_model_coords;
skeleton.aortic_annulus_non_model_coords = aortic_annulus_non_model_coords;
skeleton.aortic_annulus_right_model_coords = aortic_annulus_right_model_coords;


% scale comm height proportional to STJ radius 
% artifact of other scaling 
r_stj_temp = skeleton.r; % 1.67 / 2;
% r_temp = 2.3 / 2; % vbr radius
hc = 1.0 * r_stj_temp; % 0.5 * r_stj; 
h1 = 0.87 * 2 * r_stj_temp - hc; % 1.4 * r_stj - hc; 

skeleton.height_min_comm = h1; 


debug_plots = true;
if debug_plots
    
    figure; 
    hold on; 

    % plot3(aortic_annulus_left_raw(1,:), aortic_annulus_left_raw(2,:), aortic_annulus_left_raw(3,:),'*-');


    % plot3(aortic_annulus_left_sim_coords(1,:), aortic_annulus_left_sim_coords(2,:), aortic_annulus_left_sim_coords(3,:),'*-');


    plot3(aortic_annulus_left_model_coords(1,:), aortic_annulus_left_model_coords(2,:), aortic_annulus_left_model_coords(3,:),'*-');
    plot3(aortic_annulus_non_model_coords(1,:), aortic_annulus_non_model_coords(2,:), aortic_annulus_non_model_coords(3,:),'*-');
    plot3(aortic_annulus_right_model_coords(1,:), aortic_annulus_right_model_coords(2,:), aortic_annulus_right_model_coords(3,:),'*-');

    % plot3(aortic_annulus_left_sim_coords_transformed(1,:), aortic_annulus_left_sim_coords_transformed(2,:), aortic_annulus_left_sim_coords_transformed(3,:),'o-');
    
    axis equal
end 





