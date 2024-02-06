function valve_with_reference = aortic_dialate_annulus(valve_with_reference)

    % dialation with preserved arc length 

    if isfield(valve_with_reference, 'dilation_dist')
        dilation_dist = valve_with_reference.dilation_dist; 
    else 
        error('must provide dilation_dist if valve.dilate_graft is true'); 
    end 

    leaflet = valve_with_reference.leaflets(1); 

    % compute annular length initial 
    [~, len_annulus_min_initial] = build_initial_fibers_aortic(leaflet, valve_with_reference); 

    % adjust radius 
    valve_with_reference.skeleton.r            = valve_with_reference.skeleton.r            + dilation_dist; 
    valve_with_reference.skeleton.r_commissure = valve_with_reference.skeleton.r_commissure + dilation_dist; 

    height_min_comm_override_initial_guess = valve_with_reference.skeleton.height_min_comm; 

    options = optimset('Display','off','TolFun',1e-16);

    % this has two return values 
    len_annulus_tmp = @(height_min_comm) build_initial_fibers_aortic(leaflet, valve_with_reference, height_min_comm); 

    % hack to return only the next one 
    len_annulus = @(height_min_comm) Out2(len_annulus_tmp, height_min_comm); 
    
    len_annulus_minus_initial = @(height_min_comm) abs(len_annulus(height_min_comm) - len_annulus_min_initial); 

    height_min_comm_new = fsolve(len_annulus_minus_initial,height_min_comm_override_initial_guess,options);     
    
    X = build_initial_fibers_aortic(leaflet, valve_with_reference, height_min_comm_new); 
    
    valve_with_reference.leaflets(1).X = X;  

    % for j=1:j_max 
    %     for k=1:k_max 
    %        [th, r, z] = cart2pol(X(1,j,k), X(2,j,k), X(3,j,k));
    %        [x_tmp, y_tmp, z_tmp] = pol2cart(th,r + dilation_dist,z); 
    %        X(:,j,k) = [x_tmp; y_tmp; z_tmp]; 
    %     end 
    % end 
    % 
    % valve_with_reference.leaflets(1).X = X; 
    % 
    % valve_with_reference.skeleton.r = valve_with_reference.skeleton.r + dilation_dist; 
    % valve_with_reference.skeleton.r_of_z = @(z) valve_with_reference.skeleton.r_of_z(z) + dilation_dist; 

    debug = false; 
    if debug 
    %     figure; 
    % 
    %     x_component = squeeze(X(1,:,:)); 
    %     y_component = squeeze(X(2,:,:)); 
    %     z_component = squeeze(X(3,:,:)); 
    % 
    %     width = 1.0; 
    %     surf(x_component, y_component, z_component, 'LineWidth',width);
    %     axis equal 
    %     axis auto 

        valve_plot(valve_with_reference); 

    end 

end 


function Z = Out2(FUN,varargin)
    % Z = Out2(FUN,VARARGIN);
    %
    %	Provides the second output from the function
    [~,Z] = FUN(varargin{:});

end 