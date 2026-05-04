function skeleton = get_skeleton_aortic_ross1(hc_coefficient)

aortic_annulus_right = 0.1 * [
-0.46831756830215454, -169.6324462890625, -667.995361328125
1.655497670173645, -169.13372802734375, -672.6422119140625
4.604877948760986, -168.5965118408203, -677.3262939453125
6.932901859283447, -168.97323608398438, -679.2528686523438
9.91468906402588, -169.38131713867188, -680.3304443359375
13.179856107740964, -169.96281031777673, -680.0252232064211
15.645309448242188, -169.53533935546875, -678.9284057617188
17.707929611206055, -167.88890075683594, -677.8658447265625
20.345012664794922, -166.10545349121094, -675.6703491210938
20.968385696411133, -163.62380981445312, -672.0158081054688
20.454574584960938, -162.8600311279297, -667.1326904296875
19.432260513305664, -163.73670959472656, -660.2597045898438
]';

aortic_annulus_non = 0.1 * [
4.243202209472656, -149.2470245361328, -664.0732421875
5.8963446617126465, -149.24378967285156, -670.3046264648438
5.979079246520996, -149.17974853515625, -675.0145263671875
5.555403232574463, -150.14662170410156, -678.0017700195312
4.8892316818237305, -151.9224395751953, -680.411376953125
3.4961860179901123, -154.25082397460938, -681.8854370117188
1.2636120877285006, -156.74777883724894, -682.4976935964681
0.5603188872337341, -159.28330993652344, -682.5777587890625
1.216586947441101, -161.9753875732422, -682.1973266601562
1.936771273612976, -164.267333984375, -680.6436157226562
2.1807055473327637, -166.66575622558594, -677.3714599609375
1.0360711812973022, -168.52215576171875, -672.9295654296875
-0.35125765204429626, -169.6897430419922, -668.0812377929688
]';

aortic_annulus_left = 0.1 * [
19.453561782836914, -163.67344665527344, -660.2803955078125
20.64123916625977, -162.15647888183594, -667.0173950195312
21.966136932373047, -160.14491271972656, -672.1704711914062
22.104199020134214, -156.85719660646612, -674.1333463077202
21.760543823242188, -154.13853454589844, -674.685546875
20.2633665438459, -151.49329726433314, -675.154004442689
18.094959259033203, -149.1846923828125, -675.5745849609375
15.137526793271627, -148.11981032675791, -675.7646523293416
12.451701164245604, -147.70068359375, -675.5194091796875
8.696045875549316, -148.48374938964844, -673.3262939453125
6.875917434692383, -148.8879852294922, -670.0118408203125
4.268559455871582, -149.24269104003906, -664.0699462890625
]';



[~, n_pts_left ]  = size(aortic_annulus_left); 
[~, n_pts_non  ]  = size(aortic_annulus_non); 
[~, n_pts_right] = size(aortic_annulus_right); 


% aortic_annulus_left  = aortic_annulus_left_raw; 
% aortic_annulus_non   = aortic_annulus_non_raw; 
% aortic_annulus_right = aortic_annulus_right_raw; 


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

dist_LN = norm(comm_LN - midpoint_comms) 
dist_NR = norm(comm_NR - midpoint_comms) 
dist_LR = norm(comm_LR - midpoint_comms)

mean_radius = mean(distances);

skeleton.r = mean_radius; 
skeleton.r_commissure = skeleton.r; 

skeleton.aortic_annulus_left_model_coords = aortic_annulus_left_model_coords;
skeleton.aortic_annulus_non_model_coords = aortic_annulus_non_model_coords;
skeleton.aortic_annulus_right_model_coords = aortic_annulus_right_model_coords;


% scale comm height proportional to STJ radius 
% artifact of other scaling 
skeleton.r_stj = skeleton.r; % 1.67 / 2;


if ~exist("hc_coefficient", "var") 
    hc_coefficient = 0.5; 
end 

% default 
% hc = 0.5 * r_stj; 
hc = hc_coefficient * skeleton.r_stj; 
h1 = 1.4 * skeleton.r_stj - hc; 
skeleton.height_min_comm = h1; 


debug_plots = true;
if debug_plots
    
    figure; 
    hold on; 

    % plot3(aortic_annulus_left_raw(1,:), aortic_annulus_left_raw(2,:), aortic_annulus_left_raw(3,:),'*-');


    % plot3(aortic_annulus_left_sim_coords(1,:), aortic_annulus_left_sim_coords(2,:), aortic_annulus_left_sim_coords(3,:),'*-');


    plot3(aortic_annulus_left_model_coords(1,:), aortic_annulus_left_model_coords(2,:), aortic_annulus_left_model_coords(3,:),'k*-');
    plot3(aortic_annulus_non_model_coords(1,:), aortic_annulus_non_model_coords(2,:), aortic_annulus_non_model_coords(3,:),'r*-');
    plot3(aortic_annulus_right_model_coords(1,:), aortic_annulus_right_model_coords(2,:), aortic_annulus_right_model_coords(3,:),'g*-');
    legend('left','non','right')

    % plot3(aortic_annulus_left_sim_coords_transformed(1,:), aortic_annulus_left_sim_coords_transformed(2,:), aortic_annulus_left_sim_coords_transformed(3,:),'o-');
    
    axis equal
end 





