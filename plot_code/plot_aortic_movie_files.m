function plot_aortic_movie_files(data_file, stride, nframes, outdir, basename, image_size)

% read file 
load(data_file, 'bc_vals')

MMHG_TO_CGS = 1333.22368;

times = bc_vals(:,1);
p_aorta = bc_vals(:,2)/MMHG_TO_CGS;
q_aorta = bc_vals(:,3);
% p_wk = bc_vals(:,4)/MMHG_TO_CGS;
p_lv = bc_vals(:,5)/MMHG_TO_CGS;

width = 6; 
fontsize = 96; 
font = 'Times New Roman'; 

x_shift_position = 0.04; 
y_border_bottom = 0.07; 
y_height = .4; 

dt = times(2) - times(1);
net_flux = dt*cumsum(q_aorta);

for step = floor(nframes/2) % 0:(nframes-1)  
        
    range = 1:(step * stride); 
    if step == 0
        range = 1; 
    end 
    
    fig = figure;
    fig.Units = 'points';
    fig.Position = [0 0 image_size(1)/4, image_size(2)/4];

    outname = sprintf('%s%s_%04d.jpeg', outdir, basename, step); 
    
    hAxis(1) = subplot(2,1,1);
    plot(times(range), p_aorta(range), 'k', 'LineWidth', width)
    axis([0 2.4 -5 160])
    
    ax = gca; 
    ax.FontSize = fontsize; 
    ax.FontName = font;  

    hold on
    % plot(times, p_wk, ':k')
    plot(times(range), p_lv(range), '--k', 'LineWidth', width)
    legend('Ao', 'LV', 'Location','NorthEast', 'AutoUpdate','off');
    xlabel('Time (s)');
    ylabel('P (mmHg)');
    
    plot(times(max(range)), p_aorta(max(range)), 'ko', 'MarkerSize', 4*width, 'MarkerFaceColor', 'k')
    plot(times(max(range)), p_lv(max(range)), 'ko', 'MarkerSize', 4*width, 'MarkerFaceColor', 'k')
        
    hAxis(2) = subplot(2,1,2); 
    plot(times(range), q_aorta(range), 'k', 'LineWidth', width)
    axis([0 2.4 -750 750])
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
    