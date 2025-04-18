
MMHG_TO_CGS = 1333.22368;


line_width = 4; 
font = 24; 
cycle_duration = 0.8; 

min_plot_time = 0.4 + cycle_duration; 
max_plot_time = 0.8 + cycle_duration; 

plot_width = 1000; 
plot_height = 750; 


% colors = [1 0 0; 
%           0.8500 0.3250 0.0980;
%           0.9290 0.6940 0.1250;
%           0.4660 0.6740 0.1880;
%           0 0 0;
%           0.3010 0.7450 0.9330;
%           0 0 1;
%           0.4940 0.1840 0.5560];
      
purple = [0.4940 0.1840 0.5560] * 1.1;
purple = purple / norm(purple); 
      
orange = [255,140,0]; 
orange = orange / max(orange); 

yellow = [0.9290 .75 0.1250]; 

colors = [1 0 0;
          orange;
          yellow; 
          0.4660 0.6740 0.1880;
          0.3010 0.7450 0.9330;
          0 0 1;
          purple;
          0 0 0;];      
      

basedir = '~/Dropbox/stanford/research_stanford/aortic_graft_2024/'
      
      
names_struct(1).dir =  'aortic_57411843_384_f9882c3_bicusp_c1pt57d_r1pt4_sv0d_lv_stj_vbr_25_mesh_4e03466_3cm_extender'; 
names_struct(1).circ =  '3.0'; 
names_struct(1).circ_over_d =  '1.57d'; 
names_struct(1).rad =  '1.4'; 
names_struct(1).comment = 'f9882c3_bicusp_c1pt57d_r1pt4';
names_struct(1).frame_sys = 837; 
names_struct(1).frame_dia = 662; 

names_struct(2).dir =  'aortic_57414099_384_3dd0ce9_bicusp_c1pt57d_r0pt8_sv0d_lv_stj_vbr_25_mesh_4e03466_3cm_extender'; 
names_struct(2).circ =  '3.0'; 
names_struct(2).circ_over_d =  '1.57d'; 
names_struct(2).rad =  '0.8'; 
names_struct(2).comment = '3dd0ce9_bicusp_c1pt57d_r0pt8';
names_struct(2).frame_sys = 837; 
names_struct(2).frame_dia = 662; 

names_struct(3).dir =  'aortic_57414559_384_40d7682_bicusp_c1pt4d_r1pt4_sv0d_lv_stj_vbr_25_mesh_4e03466_3cm_extender'; 
names_struct(3).circ =  '3.0'; 
names_struct(3).circ_over_d =  '1.4d'; 
names_struct(3).rad =  '1.4'; 
names_struct(3).comment = '40d7682_bicusp_c1pt4d_r1pt4';
names_struct(3).frame_sys = 837; 
names_struct(3).frame_dia = 662; 

names_struct(4).dir =  'aortic_59135540_384_b086bc_c1p57stj_r_1p4_stj25_vbr25_ACTUALLY_C_1pt2'; 
names_struct(4).circ =  '3.0'; 
names_struct(4).circ_over_d =  '1.2d'; 
names_struct(4).rad =  '1.4'; 
names_struct(4).comment = 'f9882c3_bicusp_c1pt2d_r1pt4';
names_struct(4).frame_sys = 837; 
names_struct(4).frame_dia = 662; 

names_struct(5).dir =  'aortic_59342289_384_a4fa8f9_c1p0stj_r_1pt4_stj25_vbr25'; 
names_struct(5).circ =  '3.0'; 
names_struct(5).circ_over_d =  '1.0d'; 
names_struct(5).rad =  '1.4'; 
names_struct(5).comment = 'a4fa8f9_c1p0stj_r_1pt4';
names_struct(5).frame_sys = 837; 
names_struct(5).frame_dia = 662; 



y_max_pressure = 200; 
y_min_flow = -400;
y_max_flow = 800;

t_min_plot = cycle_duration; 
t_max_plot = 2*cycle_duration;

shift_time_to_zero = true; 
if shift_time_to_zero
    t_shift = t_min_plot; 
else 
    t_shift = 0; 
end 

for data_idx = 1:length(names_struct) 
        
    data_dir = strcat(basedir, '/', names_struct(data_idx).dir); 

    load(fullfile(data_dir, 'bc_data.mat'), 'times', 'q_aorta', 'p_aorta', 'p_lv', 'V_ventricle'); 


    run(fullfile(data_dir, 'integral_metric_data_contour_21'));
    times_paraview = t;
    p_lvot = P / (A * MMHG_TO_CGS); 
    q_lvot = Q;
    clear P A Q t 
    
    run(fullfile(data_dir, 'integral_metric_data_contour_12'));
    p_stj = P / (A * MMHG_TO_CGS); 
    q_stj = Q;
    clear P A Q t 


    fig = figure;
    plot(times - t_shift, p_aorta, 'color', 'g', 'LineWidth', line_width)
    hold on
    plot(times_paraview - t_shift, p_stj, 'r', 'LineWidth', line_width)
    plot(times_paraview - t_shift, p_lvot, 'b', 'LineWidth', line_width)
    plot(times, p_lv, 'k', 'LineWidth', line_width)
    
    xlim([t_min_plot  - t_shift, t_max_plot - t_shift])
    ylim([0 y_max_pressure])
    ax = gca; 
    ax.FontSize = font; 
    legend('Aorta', 'STJ', 'LVOT', 'Left Ventricle', 'Location','NorthEast');
    xlabel('Time (s)');
    ylabel('Pressure (mmHg)');
    set(fig, 'Position', [100, 100, plot_width, plot_height])
    set(fig,'PaperPositionMode','auto')
    file_name_pressure = sprintf('pressure_cycle2_%s.eps', names_struct(data_idx).comment);                 
    printfig(fig, fullfile(data_dir, file_name_pressure))


    fig = figure; 
    plot(times - t_shift, q_aorta, 'k', 'LineWidth', line_width)
    xlim([t_min_plot - t_shift, t_max_plot - t_shift])
    ylim([y_min_flow y_max_flow])
    ax = gca; 
    ax.FontSize = font; 
    hold on
    plot(times, 0*q_aorta, ':k')
    xlabel('Time (s)')
    ylabel('Flow (ml/s)')
    % legend('Flow', 'Location','NorthEast')
    set(fig, 'Position', [100, 100, plot_width, plot_height])
    set(fig,'PaperPositionMode','auto')
    file_name_flow = sprintf('flow_cycle2_%s.eps', names_struct(data_idx).comment);                 
    printfig(fig, fullfile(data_dir, file_name_flow))
    
    min_time_idx_cycle = find(times > t_min_plot,1);
    
    min_time_idx_cycle_paraview = find(times_paraview > t_min_plot,1);
    
    p_min = 0; 
    p_max = 200; 
    v_min = 0; 
    v_max = 240;
    
    V_ventricle_paraview_times = interp1(times, V_ventricle, times_paraview);
    
    
    fig = figure; 
    plot(V_ventricle(min_time_idx_cycle:end), p_lv(min_time_idx_cycle:end), 'k', 'LineWidth', line_width)    
    ax = gca; 
    ax.FontSize = font; 
    hold on
    
    plot(V_ventricle_paraview_times(min_time_idx_cycle_paraview:end), p_lvot(min_time_idx_cycle_paraview:end), 'b', 'LineWidth', line_width)    
    
    legend('LV inlet', 'LVOT mean')
    
    xlabel('Volume (ml)')
    ylabel('Pressure (mmHg)')
    axis([v_min v_max p_min p_max])
    set(fig, 'Position', [100, 100, plot_width, plot_height])
    set(fig,'PaperPositionMode','auto')
    
    file_name_pv_loop = sprintf('pv_loop_cycle2_%s.eps', names_struct(data_idx).comment);                 
    printfig(fig, fullfile(data_dir, file_name_pv_loop))
    
    clear times q_aorta p_aorta p_lv V_ventricle
    clear p_lvot q_lvot p_stj q_stj
end 
    
    


% % colors = parula(8); 
%       
% 
% cd aortic_57411843_384_f9882c3_bicusp_c1pt57d_r1pt4_sv0d_lv_stj_vbr_25_mesh_4e03466_3cm_extender
% 
% integral_metric_data_contour_21; 
% p_lvot = P / (A * MMHG_TO_CGS); 
% 
% 
% integral_metric_data_contour_12; 
% p_st_junction = P / (A * MMHG_TO_CGS); 
% 
% p_diff_mmHg_1pt4 = p_lvot - p_st_junction; 
% 
% plot(t, p_diff_mmHg_1pt4, 'color', 'k', 'LineWidth',line_width)
% 
% cd .. 
% 
% 
% cd aortic_57414099_384_3dd0ce9_bicusp_c1pt57d_r0pt8_sv0d_lv_stj_vbr_25_mesh_4e03466_3cm_extender
% 
% integral_metric_data_contour_21; 
% p_lvot = P / (A * MMHG_TO_CGS); 
% 
% 
% integral_metric_data_contour_12; 
% p_st_junction = P / (A * MMHG_TO_CGS); 
% 
% p_diff_mmHg_0pt8 = p_lvot - p_st_junction; 
% 
% plot(t, p_diff_mmHg_0pt8, 'color', colors(1,:), 'LineWidth',line_width)
% 
% cd .. 
% 
% 
% cd aortic_59135540_384_b086bc_c1p57stj_r_1p4_stj25_vbr25_ACTUALLY_C_1pt2
% 
% integral_metric_data_contour_21; 
% p_lvot = P / (A * MMHG_TO_CGS); 
% 
% integral_metric_data_contour_12; 
% p_st_junction = P / (A * MMHG_TO_CGS); 
% 
% p_diff_mmHg_1pt2d_1pt4 = p_lvot - p_st_junction; 
% plot(t, p_diff_mmHg_1pt2d_1pt4, 'color', colors(2,:), 'LineWidth',line_width)
% cd .. 
% 
% 
% cd aortic_59342289_384_a4fa8f9_c1p0stj_r_1pt4_stj25_vbr25
% 
% integral_metric_data_contour_21; 
% p_lvot = P / (A * MMHG_TO_CGS); 
% 
% integral_metric_data_contour_12; 
% p_st_junction = P / (A * MMHG_TO_CGS); 
% 
% p_diff_mmHg_1pt0d_1pt4 = p_lvot - p_st_junction; 
% plot(t, p_diff_mmHg_1pt0d_1pt4, 'color', colors(3,:), 'LineWidth',line_width)
% cd .. 
% 
% 
% 
% cd aortic_57411843_384_f9882c3_bicusp_c1pt57d_r1pt4_sv0d_lv_stj_vbr_25_mesh_4e03466_3cm_extender 
% load('bc_data.mat', 'times', 'q_aorta', 'p_aorta', 'p_lv', 'V_ventricle'); 
% times_1pt4 = times; 
% q_aorta_1pt4 = q_aorta; 
% p_aorta_1pt4 = p_aorta; 
% V_ventricle_1pt4 = V_ventricle;
% p_lv_1pt4 = p_lv; 
% clear times q_aorta p_aorta p_lv V_ventricle
% 
% run('integral_metric_data_contour_21')
% times_paraview_0pt8 = t;
% p_lvot_1pt4  = P / (A * MMHG_TO_CGS); 
% q_lvot_1pt4  = Q;
% clear P A Q t 
% run('integral_metric_data_contour_12')
% p_stj_1pt4  = P / (A * MMHG_TO_CGS); 
% q_stj_1pt4  = Q;
% clear P A Q t 
% 
% cd .. 
% 
% cd aortic_57414099_384_3dd0ce9_bicusp_c1pt57d_r0pt8_sv0d_lv_stj_vbr_25_mesh_4e03466_3cm_extender
% load('bc_data.mat', 'times', 'q_aorta', 'p_aorta', 'p_lv', 'V_ventricle'); 
% times_0pt8 = times; 
% q_aorta_0pt8 = q_aorta; 
% p_aorta_0pt8 = p_aorta; 
% V_ventricle_0pt8 = V_ventricle;
% p_lv_0pt8 = p_lv; 
% clear times q_aorta p_aorta p_lv V_ventricle
% 
% run('integral_metric_data_contour_21')
% times_paraview_1pt4  = t;
% p_lvot_0pt8 = P / (A * MMHG_TO_CGS); 
% q_lvot_0pt8 = Q;
% clear P A Q t 
% run('integral_metric_data_contour_12')
% p_stj_0pt8 = P / (A * MMHG_TO_CGS); 
% q_stj_0pt8 = Q;
% clear P A Q t 
% 
% cd .. 
% 
% cd aortic_59135540_384_b086bc_c1p57stj_r_1p4_stj25_vbr25_ACTUALLY_C_1pt2
% load('bc_data.mat', 'times', 'q_aorta', 'p_aorta', 'p_lv', 'V_ventricle'); 
% times_1pt2d_1pt4= times; 
% q_aorta_1pt2d_1pt4 = q_aorta; 
% p_aorta_1pt2d_1pt4 = p_aorta; 
% V_ventricle_1pt2d_1pt4 = V_ventricle;
% p_lv_1pt2d_1pt4 = p_lv; 
% clear times q_aorta p_aorta p_lv V_ventricle
% 
% run('integral_metric_data_contour_21')
% times_paraview_1pt2d_1pt4 = t;
% p_lvot_1pt2d_1pt4 = P / (A * MMHG_TO_CGS); 
% q_lvot_1pt2d_1pt4 = Q;
% clear P A Q t 
% run('integral_metric_data_contour_12')
% p_stj_1pt2d_1pt4 = P / (A * MMHG_TO_CGS); 
% q_stj_1pt2d_1pt4 = Q;
% clear P A Q t 
% cd .. 
% 
% cd aortic_59342289_384_a4fa8f9_c1p0stj_r_1pt4_stj25_vbr25
% load('bc_data.mat', 'times', 'q_aorta', 'p_aorta', 'p_lv', 'V_ventricle'); 
% times_1pt0d_1pt4= times; 
% q_aorta_1pt0d_1pt4 = q_aorta; 
% p_aorta_1pt0d_1pt4 = p_aorta; 
% V_ventricle_1pt0d_1pt4 = V_ventricle;
% p_lv_1pt0d_1pt4 = p_lv; 
% clear times q_aorta p_aorta p_lv V_ventricle
% 
% run('integral_metric_data_contour_21')
% times_paraview_1pt0d_1pt4 = t;
% p_lvot_1pt0d_1pt4 = P / (A * MMHG_TO_CGS); 
% q_lvot_1pt0d_1pt4 = Q;
% clear P A Q t 
% run('integral_metric_data_contour_12')
% p_stj_1pt0d_1pt4 = P / (A * MMHG_TO_CGS); 
% q_stj_1pt0d_1pt4 = Q;
% clear P A Q t 
% cd .. 
% 
% 
% 
% y_max_pressure = 80; 
% 
% xlim([min_plot_time, max_plot_time])
% ylim([0 y_max_pressure])
% 
% ax = gca; 
% ax.FontSize = font; 
% 
% ax.YTick = 0:10:y_max_pressure; 
% 
% xlabel('Time (s)')
% % ylabel("\Delta" + sprintf("P (mmHg)"))
% ylabel("Pressure gradient (mmHg)")
% 
% % title('Aortic pressure gradient, various free edge rest lengths')
% 
% legend('1.57d, 1.4', '1.57d, 0.8', '1pt2d, 1pt4', '1pt0d, 1pt4')
% 
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% 
% 
% printfig(fig, 'aortic_pressure_gradients')
% 
% 
% 
% y_max_pressure = 200; 
% y_min_flow = -400
% y_max_flow = 800;
% 
% 
% 
% fig = figure;
% plot(times_1pt4, p_aorta_1pt4, 'k', 'LineWidth', line_width)
% hold on
% plot(times_1pt4, p_lv_1pt4, '--k', 'LineWidth', line_width)
% plot(times_paraview_1pt4, p_stj_1pt4, 'b', 'LineWidth', line_width)
% plot(times_paraview_1pt4, p_lvot_1pt4, 'r', 'LineWidth', line_width)
% xlim([cycle_duration, 2*cycle_duration])
% ylim([0 y_max_pressure])
% ax = gca; 
% ax.FontSize = font; 
% legend('Aorta', 'Left Ventricle', 'stj', 'lvot', 'Location','NorthEast');
% xlabel('Time (s)');
% ylabel('Pressure (mmHg)');
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, 'pressure_cycle2_r1pt4')
% 
% 
% fig = figure; 
% plot(times_1pt4, q_aorta_1pt4, 'k', 'LineWidth', line_width)
% xlim([cycle_duration, 2*cycle_duration])
% ylim([y_min_flow y_max_flow])
% ax = gca; 
% ax.FontSize = font; 
% hold on
% plot(times_1pt4, 0*q_aorta_1pt4, ':k')
% xlabel('Time (s)')
% ylabel('Flow (ml/s)')
% % legend('Flow', 'Location','NorthEast')
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, 'flow_cycle2_r1pt4')
% 
% 
% fig = figure;
% plot(times_0pt8, p_aorta_0pt8, 'k', 'LineWidth', line_width)
% hold on
% plot(times_0pt8, p_lv_0pt8, '--k', 'LineWidth', line_width)
% plot(times_paraview_0pt8, p_stj_0pt8, 'b', 'LineWidth', line_width)
% plot(times_paraview_0pt8, p_lvot_0pt8, 'r', 'LineWidth', line_width)
% xlim([cycle_duration, 2*cycle_duration])
% ylim([0 y_max_pressure])
% ax = gca; 
% ax.FontSize = font; 
% legend('Aorta', 'Left Ventricle', 'Location','NorthEast');
% xlabel('Time (s)');
% ylabel('Pressure (mmHg)');
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, 'pressure_cycle2_r0pt8')
% 
% 
% fig = figure; 
% plot(times_0pt8, q_aorta_0pt8, 'k', 'LineWidth', line_width)
% xlim([cycle_duration, 2*cycle_duration])
% ylim([y_min_flow y_max_flow])
% ax = gca; 
% ax.FontSize = font; 
% hold on
% plot(times_0pt8, 0*q_aorta_0pt8, ':k')
% xlabel('Time (s)')
% ylabel('Flow (ml/s)')
% % legend('Flow', 'Location','NorthEast')
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, 'flow_cycle2_r0pt8')
% 
% 
% fig = figure;
% plot(times_1pt2d_1pt4, p_aorta_1pt2d_1pt4, 'k', 'LineWidth', line_width)
% hold on
% plot(times_1pt2d_1pt4, p_lv_1pt2d_1pt4, '--k', 'LineWidth', line_width)
% plot(times_paraview_1pt2d_1pt4, p_stj_1pt2d_1pt4, 'b', 'LineWidth', line_width)
% plot(times_paraview_1pt2d_1pt4, p_lvot_1pt2d_1pt4, 'r', 'LineWidth', line_width)
% xlim([cycle_duration, 2*cycle_duration])
% ylim([0 y_max_pressure])
% ax = gca; 
% ax.FontSize = font; 
% legend('Aorta', 'Left Ventricle', 'Location','NorthEast');
% xlabel('Time (s)');
% ylabel('Pressure (mmHg)');
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, 'pressure_cycle2_1pt2d_1pt4')
% 
% 
% fig = figure; 
% plot(times_1pt2d_1pt4, q_aorta_1pt2d_1pt4, 'k', 'LineWidth', line_width)
% xlim([cycle_duration, 2*cycle_duration])
% ylim([y_min_flow y_max_flow])
% ax = gca; 
% ax.FontSize = font; 
% hold on
% plot(times_1pt2d_1pt4, 0*q_aorta_1pt2d_1pt4, ':k')
% xlabel('Time (s)')
% ylabel('Flow (ml/s)')
% % legend('Flow', 'Location','NorthEast')
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, 'flow_cycle2_1pt2d_1pt4')
% 
% 
% fig = figure;
% plot(times_1pt0d_1pt4, p_aorta_1pt0d_1pt4, 'k', 'LineWidth', line_width)
% hold on
% plot(times_1pt0d_1pt4, p_lv_1pt0d_1pt4, '--k', 'LineWidth', line_width)
% plot(times_paraview_1pt0d_1pt4, p_stj_1pt0d_1pt4, 'b', 'LineWidth', line_width)
% plot(times_paraview_1pt0d_1pt4, p_lvot_1pt0d_1pt4, 'r', 'LineWidth', line_width)
% xlim([cycle_duration, 2*cycle_duration])
% ylim([0 y_max_pressure])
% ax = gca; 
% ax.FontSize = font; 
% legend('Aorta', 'Left Ventricle', 'Location','NorthEast');
% xlabel('Time (s)');
% ylabel('Pressure (mmHg)');
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, 'pressure_cycle2_1pt0d_1pt4')
% 
% 
% fig = figure; 
% plot(times_1pt0d_1pt4, q_aorta_1pt0d_1pt4, 'k', 'LineWidth', line_width)
% xlim([cycle_duration, 2*cycle_duration])
% ylim([y_min_flow y_max_flow])
% ax = gca; 
% ax.FontSize = font; 
% hold on
% plot(times_1pt0d_1pt4, 0*q_aorta_1pt0d_1pt4, ':k')
% xlabel('Time (s)')
% ylabel('Flow (ml/s)')
% % legend('Flow', 'Location','NorthEast')
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, 'flow_cycle2_1pt0d_1pt4')
% 
% 
% 
% min_time_idx_cycle = find(times_0pt8 > cycle_duration,1);
% 
% p_min = 0; 
% p_max = 200; 
% v_min = 0; 
% v_max = 240;
% 
% fig = figure; 
% plot(V_ventricle_1pt4(min_time_idx_cycle:end), p_lv_1pt4(min_time_idx_cycle:end), 'k', 'LineWidth', line_width)
% ax = gca; 
% ax.FontSize = font; 
% hold on
% xlabel('volume (ml)')
% ylabel('p (mmHg)')
% axis([v_min v_max p_min p_max])
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, 'pv_loop_r1pt4')
% 
% fig = figure; 
% plot(V_ventricle_0pt8(min_time_idx_cycle:end), p_lv_0pt8(min_time_idx_cycle:end), 'k', 'LineWidth', line_width)
% ax = gca; 
% ax.FontSize = font; 
% hold on
% xlabel('volume (ml)')
% ylabel('p (mmHg)')
% axis([v_min v_max p_min p_max])
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, 'pv_loop_r0pt8')
% 
% fig = figure; 
% plot(V_ventricle_1pt2d_1pt4(min_time_idx_cycle:end), p_lv_1pt2d_1pt4(min_time_idx_cycle:end), 'k', 'LineWidth', line_width)
% ax = gca; 
% ax.FontSize = font; 
% hold on
% xlabel('volume (ml)')
% ylabel('p (mmHg)')
% axis([v_min v_max p_min p_max])
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, 'pv_loop_1pt2d_1pt4')
% 
% fig = figure; 
% plot(V_ventricle_1pt0d_1pt4(min_time_idx_cycle:end), p_lv_1pt0d_1pt4(min_time_idx_cycle:end), 'k', 'LineWidth', line_width)
% ax = gca; 
% ax.FontSize = font; 
% hold on
% xlabel('volume (ml)')
% ylabel('p (mmHg)')
% axis([v_min v_max p_min p_max])
% set(fig, 'Position', [100, 100, plot_width, plot_height])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, 'pv_loop_1pt0d_1pt4')
% 
% 
% 
% 
