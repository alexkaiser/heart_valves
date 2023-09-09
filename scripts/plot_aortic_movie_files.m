function plot_aortic_movie_files(data_file, stride, nframes, outdir, basename, image_size, axis_vec_pressure, axis_vec_flow)

% read file 
load(data_file, 'bc_vals')

MMHG_TO_CGS = 1333.22368;

times = bc_vals(:,1);
p_aorta = bc_vals(:,2)/MMHG_TO_CGS;
q_aorta = bc_vals(:,3);
% p_wk = bc_vals(:,4)/MMHG_TO_CGS;
p_lv = bc_vals(:,5)/MMHG_TO_CGS;

bicuspid_paper = false; 
if bicuspid_paper
    p_aorta          = bc_vals(:,3); 
    q_aorta          = bc_vals(:,5);
    p_extender_point = bc_vals(:,9);
    p_lv             = p_extender_point; 
end 

bicuspidization_paper = true; 
if bicuspidization_paper
    times            =  bc_vals(:,1);
    p_lv             =  bc_vals(:,2);
    p_aorta          =  bc_vals(:,3); 
    q_aorta          =  bc_vals(:,5);
end 


width = 3; 
fontsize = 44; 
font = 'Times New Roman'; 

x_shift_position = 0.04; 
y_border_bottom = 0.07; 
y_height = .4; 

dt = times(2) - times(1);
net_flux = dt*cumsum(q_aorta);

for step = 0:(nframes-1)  % floor(nframes/2)
        
    range = 1:(step * stride); 
    if step == 0
        range = 1; 
    end 
    
    fig = figure('visible','off');
    fig.Units = 'points';
    fig.Position = (1/2) * [0 0 image_size(1), image_size(2)];

    outname = sprintf('%s%s_%04d.jpeg', outdir, basename, step); 
    
    hAxis(1) = subplot(2,1,1);
    plot(times(range), p_aorta(range), 'k', 'LineWidth', width)
    axis(axis_vec_pressure)
    
    ax = gca; 
    ax.FontSize = fontsize; 
    ax.FontName = font;  

    hold on
    % plot(times, p_wk, ':k')
    plot(times(range), p_lv(range), '--k', 'LineWidth', width)
    legend('Aorta', 'Left Ventricle', 'Location','NorthEast', 'AutoUpdate','off');
    xlabel('Time (s)');
    ylabel('Pressure (mmHg)');
    
    plot(times(max(range)), p_aorta(max(range)), 'ko', 'MarkerSize', 4*width, 'MarkerFaceColor', 'k')
    plot(times(max(range)), p_lv(max(range)), 'ko', 'MarkerSize', 4*width, 'MarkerFaceColor', 'k')
        
    hAxis(2) = subplot(2,1,2); 
    plot(times(range), q_aorta(range), 'k', 'LineWidth', width)
    axis(axis_vec_flow)
    ax = gca; 
    ax.FontSize = fontsize; 
    ax.FontName = font;  
        
    hold on
    plot(times(range), net_flux(range), '--k', 'LineWidth', width)
    legend('Flow', 'Cumulative Flow', 'Location', 'NorthEast', 'AutoUpdate','off')
    plot(times, 0*net_flux,':k', 'LineWidth', width/2); 
    xlabel('Time (s)')
    ylabel('Flow (ml/s), Cumulative Flow (ml)'); 
    
    plot(times(max(range)), q_aorta(max(range)), 'ko', 'MarkerSize', 4*width, 'MarkerFaceColor', 'k')
    plot(times(max(range)), net_flux(max(range)), 'ko', 'MarkerSize', 4*width, 'MarkerFaceColor', 'k')
    
    % little more rightward placement 
    hAxis(1).Position(1) = hAxis(1).Position(1) + x_shift_position; 
    hAxis(2).Position(1) = hAxis(2).Position(1) + x_shift_position; 
    
    hAxis(1).Position(2) = y_border_bottom; 
    hAxis(2).Position(2) = 0.5 + y_border_bottom; 
    
    hAxis(1).Position(4) = y_height; 
    hAxis(2).Position(4) = y_height; 

    fig.PaperPositionMode = 'manual'; 
    fig.PaperUnits = 'points';
    fig.PaperPosition = [0 0 image_size(1), image_size(2)];

    print(fig, '-djpeg', outname);
    close(fig)
    
end 
    