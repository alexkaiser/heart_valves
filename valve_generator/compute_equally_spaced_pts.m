function [pts, fig] = compute_equally_spaced_pts(X, n_pts_query, fig)

[dim, n_pts_orig] = size(X); 

if ~((dim == 2) || (dim == 3))
    warning('dimension %f points being respaced\n', dim); 
end 

if ~ismatrix(X)
    error('Must pass R^{dimension x n} matrix to interpolate')
end 


arc_len_cumulative = zeros(1, n_pts_orig);
for j=2:n_pts_orig
    arc_len_cumulative(j) = arc_len_cumulative(j-1) + norm(X(:,j-1) - X(:,j)); 
end 

arc_len_total = arc_len_cumulative(end); 

query_pts = arc_len_total * linspace(0,1,n_pts_query);

% two transposes for convention that first index is component 
pts = interp1(arc_len_cumulative, X', query_pts)'; 


debug = false; 
% error checking 
if debug 
    
    x_component = squeeze(pts(1,:,1)); 
    y_component = squeeze(pts(2,:,1)); 
    z_component = squeeze(pts(3,:,1)); 

    if ~exist('fig', 'var')
        fig = figure; 
    end        

    hold on 
    width = 1.0; 
    plot3(x_component, y_component, z_component, 'o', 'LineWidth', width);

    x_component_X = squeeze(X(1,:)); 
    y_component_X = squeeze(X(2,:)); 
    z_component_X = squeeze(X(3,:)); 

    plot3(x_component_X, y_component_X, z_component_X, 'LineWidth',width);

    axis equal 
    axis auto 

%         xlabel('x'); 
%         ylabel('y'); 
    title('spacing ')

    % quick length consistency check 
    tol_lengths = 1e-3; 
    max_diff_absolute = 0; 


    % comm point below 
    arc_len_first = norm(pts(:,1,1) - pts(:,2,1)); 

    for j=2:n_pts_query
        j
        arc_len_current = norm(pts(:,j-1,1) - pts(:,j,1)) 
        max_diff_absolute = max(max_diff_absolute, abs((arc_len_current - arc_len_first)));
        rel_err = abs((arc_len_current - arc_len_first)/arc_len_first)
        if rel_err > tol_lengths
            warning('inconsistent arc lengths after equalizing');                     
        end 
    end 

    max_diff_absolute
    
    
    
end 


