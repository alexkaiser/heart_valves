function skeleton = get_skeleton_aortic_ross1(hc_coefficient)

aortic_annulus_right = 0.1 * [
19.433538,-163.73706,-660.2596
19.459713,-163.69975,-660.4658
19.53419,-163.60135,-661.0337
19.65783,-163.46451,-661.88666
19.836311,-163.31181,-662.9451
20.038195,-163.15636,-664.1356
20.212803,-163.00824,-665.3876
20.364815,-162.8894,-666.62177
20.514816,-162.82658,-667.7628
20.641315,-162.83022,-668.7905
20.730425,-162.90749,-669.73016
20.8121,-163.07538,-670.60956
20.924078,-163.35997,-671.4347
21.046417,-163.77632,-672.2216
21.091318,-164.28946,-673.02875
21.06606,-164.86896,-673.7943
20.971754,-165.5098,-674.5053
20.708475,-166.13551,-675.22144
20.305454,-166.70213,-675.8549
19.750938,-167.24696,-676.4133
19.11944,-167.79134,-676.8684
18.41706,-168.33383,-677.39465
17.659946,-168.78087,-677.8553
16.90435,-169.14621,-678.3293
16.0989,-169.41524,-678.71234
15.309302,-169.58751,-679.09375
14.482385,-169.66243,-679.4933
13.6198635,-169.66151,-679.8598
12.718301,-169.6093,-680.13965
11.789796,-169.55338,-680.3167
10.855456,-169.50697,-680.3893
9.924171,-169.4258,-680.3373
9.002627,-169.30168,-680.16034
8.120754,-169.15364,-679.85754
7.27881,-169.0023,-679.44086
6.475516,-168.85934,-678.9335
5.721984,-168.72424,-678.35376
5.0271006,-168.61725,-677.7295
4.410381,-168.55533,-677.0815
3.86845,-168.55576,-676.4166
3.3723733,-168.61877,-675.72095
2.8981633,-168.73242,-674.9725
2.4298093,-168.87648,-674.1428
1.9493697,-169.02963,-673.2028
1.434102,-169.17903,-672.1317
0.8975259,-169.31488,-670.97485
0.39007354,-169.43018,-669.8544
-0.04023284,-169.52454,-668.8993
-0.34271625,-169.59325,-668.23645
-0.4578438,-169.62035,-667.9891
]';

aortic_annulus_non = 0.1 * [
4.2439857,-149.25241,-664.07306
4.314662,-149.25061,-664.2935
4.503331,-149.25139,-664.8955
4.773851,-149.25787,-665.7898
5.089114,-149.26509,-666.8875
5.4122543,-149.27104,-668.10016
5.707207,-149.27333,-669.3385
5.9336095,-149.24777,-670.51575
6.068756,-149.19376,-671.5767
6.123294,-149.13579,-672.5423
6.1161256,-149.11397,-673.44434
6.0613947,-149.14476,-674.3173
5.9700823,-149.24452,-675.1956
5.853616,-149.42674,-676.09753
5.7227206,-149.71355,-677.0009
5.5747957,-150.08385,-677.88275
5.3831067,-150.49184,-678.73114
5.1289763,-150.95189,-679.52423
4.820482,-151.4881,-680.23755
4.404068,-152.06158,-680.8246
3.9050379,-152.63853,-681.3333
3.3123765,-153.30002,-681.7515
2.7079368,-154.04572,-682.0684
2.1659806,-154.84073,-682.29663
1.6796508,-155.64888,-682.4329
1.2632793,-156.4981,-682.50964
0.9425846,-157.3917,-682.54224
0.7402282,-158.31801,-682.55804
0.66799897,-159.25545,-682.54987
0.7435207,-160.18974,-682.52203
0.9711184,-161.10165,-682.4094
1.2564453,-161.975,-682.18005
1.5355598,-162.78595,-681.7894
1.7711624,-163.54433,-681.266
1.9595852,-164.24788,-680.6367
2.1086397,-164.89275,-679.9257
2.2025962,-165.48553,-679.1636
2.2349586,-166.03069,-678.3804
2.215378,-166.51384,-677.5932
2.1353502,-166.94072,-676.825
1.9907551,-167.32277,-676.0655
1.789898,-167.67145,-675.27905
1.5357909,-168.00218,-674.4283
1.2370191,-168.3207,-673.4734
0.9061394,-168.63277,-672.37427
0.5548462,-168.94173,-671.1812
0.22416002,-169.2186,-670.0185
-0.059677925,-169.45004,-669.02405
-0.26185608,-169.61395,-668.3327
-0.34058297,-169.67737,-668.07495
]';

aortic_annulus_left = 0.1 * [
19.452919,-163.67326,-660.28046
19.481085,-163.62697,-660.4805
19.562323,-163.50002,-661.0321
19.695978,-163.31151,-661.8617
19.886208,-163.08096,-662.8945
20.106058,-162.82114,-664.0609
20.313852,-162.54182,-665.2933
20.510635,-162.26152,-666.51495
20.705908,-162.00078,-667.6565
20.886736,-161.75081,-668.6885
21.048382,-161.48747,-669.6229
21.23581,-161.19316,-670.4879
21.50972,-160.84715,-671.2619
21.852413,-160.40514,-671.8972
22.172497,-159.80908,-672.4082
22.434351,-159.05983,-672.843
22.633417,-158.18399,-673.23785
22.72887,-157.27174,-673.5633
22.717665,-156.2854,-673.8407
22.51607,-155.26636,-674.10803
22.195688,-154.38336,-674.3824
21.803247,-153.60162,-674.6071
21.509481,-152.74823,-674.7395
21.131615,-151.86287,-674.76013
20.594097,-150.9778,-674.7453
19.902098,-150.1164,-674.8394
19.191925,-149.37547,-675.0451
18.404894,-148.81644,-675.2426
17.503258,-148.3981,-675.4526
16.566925,-148.15741,-675.59406
15.620263,-147.99873,-675.698
14.668692,-147.90027,-675.7439
13.728168,-147.83374,-675.7234
12.800973,-147.80486,-675.6073
11.892329,-147.83218,-675.37256
11.00966,-147.94379,-675.0171
10.175829,-148.15492,-674.5591
9.409611,-148.35803,-674.0014
8.746337,-148.5129,-673.37604
8.210282,-148.65967,-672.6967
7.7727213,-148.81798,-671.9489
7.3734207,-148.97159,-671.10004
6.9546733,-149.08414,-670.126
6.4754066,-149.14583,-669.02454
5.9678164,-149.1774,-667.86255
5.4594884,-149.20537,-666.7234
4.9973335,-149.22227,-665.6897
4.6185565,-149.23265,-664.8465
4.3630123,-149.24312,-664.27826
4.269419,-149.2486,-664.06976
]';

% left is already correct 
% aortic_annulus_left  = fliplr(aortic_annulus_left);
% aortic_annulus_non   = fliplr(aortic_annulus_non);

% these happen to be in right to left order from segmentation 
aortic_annulus_right = fliplr(aortic_annulus_right);


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





