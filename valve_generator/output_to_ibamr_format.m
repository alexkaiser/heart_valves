function [] = output_to_ibamr_format(valve)
    % 
    % Outputs the current configuration of the leaflets to IBAMR format
    % Spring constants are computed in dimensional form 
    % 
    %
    % Input: 
    %    base_name                  File base name
    %    L                          Outputs extra mesh to use a [-L,L]^3 cube
    %    ratio                      Ratio of pressure to nondimensionalized spring constant   
    %    params_posterior           Parameters for various leaflets 
    %    filter_params_posterior
    %    params_anterior
    %    p_physical                 Pressure in mmHg
    %    target_multiplier          Target spring strength is target_multiplier * p_physical/ratio
    %    refinement                 N/32. Mesh is this many times finer than original debug width 
    %    n_lagrangian_tracers       Places a 3d uniform mesh of tracers, this many (plus one) per side
    %    X_config_is_reference      Replaces the reference configuration with the current configuration
    %                               This attempts to remove initial transients entirely 
    % 
    % 
    % Output: 
    %    Files written in IBAMR format 
    % 
    
    % Copyright (c) 2019, Alexander D. Kaiser
    % All rights reserved.
    % 
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are met:
    % 
    % 1. Redistributions of source code must retain the above copyright notice, this
    %    list of conditions and the following disclaimer.
    % 
    % 2. Redistributions in binary form must reproduce the above copyright notice,
    %    this list of conditions and the following disclaimer in the documentation
    %    and/or other materials provided with the distribution.
    % 
    % 3. Neither the name of the copyright holder nor the names of its
    %    contributors may be used to endorse or promote products derived from
    %    this software without specific prior written permission.
    % 
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    % DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    % DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    % SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    % CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    % OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    % OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    N                           = valve.N; 
    base_name                   = valve.base_name; 
    L                           = valve.L; 
    tension_base                = valve.tension_base; 
    target_net                  = valve.target_net; 
    target_papillary            = valve.target_papillary; 
    n_lagrangian_tracers        = valve.n_lagrangian_tracers; 
    
    if isfield(valve, 'collagen_constitutive')
        collagen_constitutive       = valve.collagen_constitutive; 
    else
        collagen_constitutive       = false; 
    end 
    
    % if this is true, all partition gets ignored except right at valve
    % ring 
    if ~isfield(valve, 'in_heart')
        in_heart = false; 
    else 
        in_heart = valve.in_heart; 
    end 
    
    if in_heart 
        n_lagrangian_tracers = 0; 
    end 
    
    if ~isfield(valve, 'name')
        valve.name = 'default_mitral'; 
    end 
    
    if ~strcmp(valve.name, 'aortic') && ~valve.split_papillary
        error('Must have split papillary locations in current implementation.'); 
    end

    if strcmp(valve.name, 'aortic') 
        params.type = 'aortic'; 
    else 
        params.type = 'default_mitral'; 
    end 
    
    if strcmp(valve.name, 'aortic') 
        if length(valve.leaflets) ~= 1
            error('only one aortic leaflet supported'); 
        end 
    end 
    
    params.vertex        = fopen(strcat(base_name, '.vertex'), 'w'); 
    params.spring        = fopen(strcat(base_name, '.spring'), 'w'); 
    params.target        = fopen(strcat(base_name, '.target'), 'w'); 
    params.inst          = fopen(strcat(base_name, '.inst'), 'w'); 
    if ~strcmp(params.type, 'aortic') 
        params.papillary     = fopen(strcat(base_name, '.papillary'), 'w'); 
    end 
    
    if isfield(valve, 'k_bend_radial')           || isfield(valve, 'k_bend_circ') || ... 
       isfield(valve, 'k_bend_radial_free_edge') || isfield(valve, 'k_bend_cross_layer') || ...
       isfield(valve, 'k_bend_nodule') 
        
        if isfield(valve, 'k_bend_radial')
            k_bend_radial = valve.k_bend_radial; 
        else 
            k_bend_radial = 0.0; 
        end
        
        if isfield(valve, 'k_bend_circ')
            k_bend_circ = valve.k_bend_circ; 
        else
            k_bend_circ = 0.0; 
        end
        
        if isfield(valve, 'k_bend_radial_free_edge')
            k_bend_radial_free_edge = valve.k_bend_radial_free_edge; 
        else 
            k_bend_radial_free_edge = 0.0; 
        end
        
        if isfield(valve, 'k_bend_cross_layer')
            k_bend_cross_layer = valve.k_bend_cross_layer; 
        else 
            k_bend_cross_layer = 0.0; 
        end
        
        if isfield(valve, 'k_bend_nodule')
            k_bend_nodule = valve.k_bend_nodule; 
            if ~isfield(valve, 'k_bend_nodule_length')
                error('must provide nodule length if adding nodule bending'); 
            end 
        else 
            k_bend_nodule = 0.0; 
        end
        
        if any(k_bend_radial > 0) || (k_bend_circ > 0) || (k_bend_radial_free_edge > 0) || (k_bend_cross_layer > 0) || (k_bend_nodule)
            params.beam_on = true; 
        else 
            params.beam_on = false;
        end 
    else
        params.beam_on = false; 
    end     
    if params.beam_on
        params.beam = fopen(strcat(base_name, '.beam'), 'w');
    end 
    
    % just make this ridiculously big for now 
    % would be better to implement some resizing but that will also clutter things up 
    params.vertices_capacity = 100 * N^2; 
    params.vertices = zeros(3,params.vertices_capacity);        % this has the initial points for target points 
    
    if valve.targets_for_bcs
        params.targets_for_bcs = true; 
        params.vertices_target = nan * zeros(3,params.vertices_capacity); % this has the target locations for target points 
        params.vertex_target_pos_file = fopen(strcat(base_name, '_target_position.vertex'), 'w'); 
        params.total_vertices_target = 0; 
    else 
        params.targets_for_bcs = false; 
    end 
    
    if isfield(valve.skeleton, 'inverse_transformation_initial_condition')
        params.inverse_transformation_initial_condition = valve.skeleton.inverse_transformation_initial_condition; 
    end 
    
    % keep one global index through the whole thing 
    % every time a vertex is placed this is incremented 
    % and params.vertics(:,i) contains the coordinates 
    params.global_idx = 0;
    
    % also keep an index for placing after known indices 
    params.max_idx_after_reserved_indices = 0; 
    
    % just count the number of vertices and strings throughout 
    params.total_vertices  = 0; 
    params.total_springs   = 0; 
    params.total_targets   = 0; 
    if ~strcmp(params.type, 'aortic') 
        params.total_papillary = 0; 
    end 
    if params.beam_on
        params.total_beams = 0; 
    end 
    
    if strcmp(params.type, 'aortic') 
        ray_springs = true;
    else 
        ray_springs = false;
    end 
    
    
    
    % keep a single parameter for outputting copies 
    params.z_offset = 0; 

    % box sizes 
    params.x_min = -L; 
    params.x_max =  L; 
    params.y_min = -L; 
    params.y_max =  L; 
    params.z_min =  3.0 - 4*L; 
    params.z_max =  3.0; 
    
    if isfield(valve.skeleton, 'ring_center')
        ring_center = valve.skeleton.ring_center
        params.x_min = params.x_min + valve.skeleton.ring_center(1); 
        params.x_max = params.x_max + valve.skeleton.ring_center(1); 
        params.y_min = params.y_min + valve.skeleton.ring_center(2); 
        params.y_max = params.y_max + valve.skeleton.ring_center(2); 
        params.z_min = params.z_min + valve.skeleton.ring_center(3); 
        params.z_max = params.z_max + valve.skeleton.ring_center(3); 
    else
        valve.skeleton.ring_center = zeros(3,1); 
    end 
    params.ring_center = valve.skeleton.ring_center; 
    
    % parameters for scaling of other constants 
    params.num_copies = valve.num_copies; 
    if isfield(valve, 'copy_spring_weights')
        params.copy_spring_weights = valve.copy_spring_weights; 
    end 
    params.eta_multiplier_linear   = valve.eta_multiplier_linear; 
    params.eta_multiplier_collagen = valve.eta_multiplier_collagen; 
    
    % parameters for output flags 
    params.output = valve.output; 
    
    % Spring constant base for targets and 
    % Approximate force is tension_base multiplied by a length element 
    du = 1/N; 
    k_rel = tension_base * du; 
        

    % papillary target constants 
    % this does not scale when the mesh is changed 
    k_target_papillary = target_papillary; 
    eta_papillary      = valve.eta_papillary; 
    
    k_target_net       = target_net; 
    eta_net            = valve.eta_net; 

    
    % Lagrangian mesh spacing 
    ds = valve.ds
    params.ds = ds; 
    
    if isfield(valve, 'kappa_cross_layer_multipler') && (valve.kappa_cross_layer_multipler ~= 0) && (params.num_copies > 1)
        params.cross_layer_on          = true; 
        params.kappa_cross_layer       = valve.kappa_cross_layer_multipler * tension_base; 
        params.rest_len_cross_layer    = ds; 
        params.total_per_layer         = nan; 
        params.min_idx_for_cross_layer = nan; 
        params.max_idx_for_cross_layer = nan; 
        
        if isfield(valve, 'k_bend_cross_layer') && (valve.k_bend_cross_layer > 0)
            params.k_bend_cross_layer = valve.k_bend_cross_layer; 
        end 
    else 
        params.cross_layer_on          = false; 
    end 
    
    if isfield(valve, 'kappa_radial_free_edge_compressive') && ...
            isfield(valve, 'kappa_radial_free_edge_compressive_percentage') && ...
            isfield(valve, 'kappa_radial_free_edge_compressive_stretch') 
        params.kappa_radial_free_edge_compressive = valve.kappa_radial_free_edge_compressive_unscaled  * tension_base; 
        params.kappa_radial_free_edge_compressive_percentage = valve.kappa_radial_free_edge_compressive_percentage; 
        params.kappa_radial_free_edge_compressive_stretch = valve.kappa_radial_free_edge_compressive_stretch;
        params.kappa_radial_free_edge_compressive_fn_idx = valve.kappa_radial_free_edge_compressive_fn_idx; 
        params.compressive_extra_springs_on = true; 
    else 
        params.compressive_extra_springs_on = false; 
    end 
    
    if isfield(valve, 'normal_thicken') && valve.normal_thicken
        params.normal_thicken = valve.normal_thicken; 
        params.ds_extrude = valve.normal_thickness / (valve.num_copies-1); 
        params.rest_len_cross_layer = params.ds_extrude; 
        
        if isfield(valve, 'center_extrusion') && valve.center_extrusion
            params.center_extrusion = valve.center_extrusion;
            
            if isfield(valve, 'pre_extrude') && valve.pre_extrude
                % already handled 
            else
                % remove outer layer at comms everywhere 
                valve.leaflets(1) = add_bc_layer_at_commmissure_aortic(valve.leaflets(1)); 
                % remove next layer too 
                % valve.leaflets(1) = add_bc_layer_at_commmissure_aortic(valve.leaflets(1)); 
                % valve.leaflets(1) = add_bc_layer_at_commmissure_aortic(valve.leaflets(1)); 
            end 
        end
        
    end 
    
    % copies, if needed, will be placed this far down 
    if params.num_copies > 1
        z_offset_vals = -ds*[0:(params.num_copies-1)]
    else 
        z_offset_vals = 0; 
    end  
    
    % print the increment for systolic motion
    % this is the negative of the motion that occurred 
    % from the systolic to diastolic before 
    if ~strcmp(valve.name, 'aortic') 
        diastolic_increment = valve.diastolic_increment; 
        fprintf(params.papillary, '%.14f\t %.14f\t %.14f\n', diastolic_increment(1), diastolic_increment(2), diastolic_increment(3)); 
        times = valve.papillary_movement_times; 
        if length(times) ~= 5
            error('Must have five times for movement'); 
        end 
        fprintf(params.papillary, '%.14f\t %.14f\t %.14f %.14f %.14f\n', times(1), times(2), times(3), times(4), times(5)); 
    end 
    
    
    % main loop over copies 
    for copy = 1:params.num_copies
        
        params.copy = copy; 
        
        if strcmp(valve.name, 'aortic') && ... 
           isfield(valve, 'pre_extrude') && valve.pre_extrude && ... 
           isfield(valve.leaflets(1), 'pre_extrude') && valve.leaflets(1).pre_extrude
       
              % extra leaflet copies already extruded  
              valve.leaflets(1).X           = valve.leaflets_pre_extruded(copy).X; 
              valve.leaflets(1).is_bc       = valve.leaflets_pre_extruded(copy).is_bc; 
              valve.leaflets(1).is_internal = valve.leaflets_pre_extruded(copy).is_internal; 
            
        end 
            
        if params.cross_layer_on
            params.min_idx_for_cross_layer = params.global_idx; 
            % fprintf('min_idx_for_cross_layer = %d\n', params.min_idx_for_cross_layer); 
        end 
        
        params.z_offset = z_offset_vals(copy); 
        fprintf('params.z_offset = %f\n', params.z_offset); 
        first_idx = params.global_idx + 1; 
        
        for i=1:length(valve.leaflets)
            j_max = valve.leaflets(i).j_max; 
            k_max = valve.leaflets(i).k_max; 
            valve.leaflets(i).indices_global = nan * zeros(j_max, k_max);  
        
            if ~strcmp(params.type, 'aortic') 
                % allocate chordae indexing arrays here 
                % to stop dissimilar struct assignment error 
                for tree_idx = 1:valve.leaflets(i).num_trees
                    C = valve.leaflets(i).chordae(tree_idx).C; 
                    [m N_chordae] = size(C);
                    valve.leaflets(i).chordae(tree_idx).idx_root = nan;                 
                    valve.leaflets(i).chordae(tree_idx).indices_global = nan * zeros(N_chordae, 1);
                end 
            end 
            
        end 
        
        for i=1:length(valve.leaflets)
            [params valve.leaflets(i)] = assign_indices_vertex_target(params, valve.leaflets(i), k_target_net, k_target_papillary, eta_net, eta_papillary); 
        end 
        
        if strcmp(valve.name, 'aortic') 
            params.layer_indices(copy).indices_global = valve.leaflets(1).indices_global; 
        end 
        
        for i=1:length(valve.leaflets)
            params = add_springs(params, valve.leaflets(i), ds, collagen_constitutive); 
        end 
        
        for i=1:length(valve.leaflets)
            if params.beam_on
                if isfield(valve, 'k_bend_radial_free_edge') && ... 
                   isfield(valve, 'k_bend_radial_free_edge_percentage') 
                    k_bend_radial_free_edge = valve.k_bend_radial_free_edge; 
                    k_bend_radial_free_edge_percentage = valve.k_bend_radial_free_edge_percentage; 
                else 
                    k_bend_radial_free_edge = 0; 
                    k_bend_radial_free_edge_percentage = 0; 
                end
                
                if isfield(valve, 'k_bend_circ_free_edge') && ... 
                   isfield(valve, 'k_bend_circ_free_edge_percentage') 
                    k_bend_circ_free_edge = valve.k_bend_circ_free_edge; 
                    k_bend_circ_free_edge_percentage = valve.k_bend_circ_free_edge_percentage; 
                else 
                    k_bend_circ_free_edge = 0; 
                    k_bend_circ_free_edge_percentage = 0; 
                end
                
                if isfield(valve, 'k_bend_radial_interp_pts')
                    params = add_beams(params, valve.leaflets(i), k_bend_radial, k_bend_circ, ...
                                                                  k_bend_radial_free_edge, k_bend_radial_free_edge_percentage, ...
                                                                  k_bend_circ_free_edge,   k_bend_circ_free_edge_percentage, ...
                                                                  valve.k_bend_radial_interp_pts);
                else 
                    params = add_beams(params, valve.leaflets(i), k_bend_radial, k_bend_circ, ...
                                                                  k_bend_radial_free_edge, k_bend_radial_free_edge_percentage, ...
                                                                  k_bend_circ_free_edge,   k_bend_circ_free_edge_percentage);
                end 
                
                
            end 
        end 

        if strcmp(params.type, 'aortic') 
            if isfield(params, 'normal_thicken') && params.normal_thicken
                % skip the bottom bc layer for cross layers, this is placed first 
                % params.min_idx_for_cross_layer = params.min_idx_for_cross_layer + valve.leaflets(1).j_max; 
                % these don't get adjusted, first idx is the minimum to adjust 
                first_idx = params.global_idx + 1; 
            end 
        end 
        
        if params.cross_layer_on
            params.max_idx_for_cross_layer = params.global_idx - 1; 
        end 
        
        % flat part of mesh 
        r = valve.r; 
        ref_frac_net = 1.0;

        if ~in_heart
            for i=1:length(valve.leaflets)
                if length(valve.leaflets) ~= 1 
                    error('Only one leaflet version currently supported'); 
                end 
                hoop_springs = true; 
                params = place_net(params, valve.leaflets(i), ds, r, L, k_rel, k_target_net, ref_frac_net, eta_net, hoop_springs, ray_springs); 
            end 

            if ~strcmp(params.type, 'aortic') 
                % approximate geodesic continutations of fibers 
                for i=1:length(valve.leaflets)
                    params = place_rays(params, valve.leaflets(i), ds, L, k_rel, k_target_net, ref_frac_net, eta_net);
                end 
            end 

            % flat part of mesh with Cartesian coordinates
            % inner radius, stop mesh here 
            r_extra = 4*ds; 
            for i=1:length(valve.leaflets)
                if length(valve.leaflets) ~= 1 
                    error('Only one leaflet version currently supported'); 
                end 
                params = place_cartesian_net(params, valve.leaflets(i), r_extra, L, ds, k_rel, k_target_net, ref_frac_net, eta_net); 
            end 
            
        else 
            
            
            if isfield(params, 'targets_for_bcs') && params.targets_for_bcs 
                if copy == 1
                    params = write_inst_for_targets_as_bcs(params, valve.leaflets(1));                 
                end 
            else             
                % pass L=r to get only one ring placed 
                hoop_springs = true; 
                
                if isfield(valve, 'extra_radius_hoops')
                    extra_radius_hoops = valve.extra_radius_hoops; 
                else 
                    extra_radius_hoops = 0; 
                end 
                
                if strcmp(params.type, 'aortic') 
                    % only place on the first copy here 
                    if copy == 1 
                        params = place_net(params, valve.leaflets(i), ds, r, r + extra_radius_hoops, k_rel, k_target_net, ref_frac_net, eta_net, hoop_springs); 
                    end 
                else
                    params = place_net(params, valve.leaflets(i), ds, r, r + extra_radius_hoops, k_rel, k_target_net, ref_frac_net, eta_net, hoop_springs); 
                end 
                
                if extra_radius_hoops > 0.0
                    max_to_place = max(0,floor(extra_radius_hoops / ds) - 1);
                    params = place_rays(params, valve.leaflets(i), ds, r + extra_radius_hoops, k_rel, k_target_net, ref_frac_net, eta_net, max_to_place);
                end 
                
            end 
        end 
        

        % first time through, count all included indices 
        if params.cross_layer_on && (copy == 1)
            params.total_per_layer = params.global_idx; 
            fprintf('params.total_per_layer = %d\n', params.total_per_layer); 
        end 
        
        % adjust for offset 
        last_idx = params.global_idx; 
        
        fprintf('this layer range, matlab indices, inclusive                          = %d:%d\n', first_idx, last_idx); 
        
        params.vertices(3,first_idx:last_idx) = params.vertices(3,first_idx:last_idx) + params.z_offset; 
        
        if params.targets_for_bcs
            params.vertices_target(3,first_idx:last_idx) = params.vertices_target(3,first_idx:last_idx) + params.z_offset; 
        end 
        
        if strcmp(params.type, 'aortic') && (copy == 1) && isfield(params, 'normal_thicken') && params.normal_thicken 
            % set end links in circ fibers to be bc in the normal thicken since they cross each other 
            % with the discontinuous normal at the commissure 
            
            if isfield(params, 'center_extrusion') && params.center_extrusion
                % nothing, already set  
            elseif isfield(valve, 'pre_extrude') && valve.pre_extrude 
                % nothing, already set      
            else
                valve.leaflets(1) = add_bc_layer_at_commmissure_aortic(valve.leaflets(1)); 
                % remove second layer too 
                valve.leaflets(1) = add_bc_layer_at_commmissure_aortic(valve.leaflets(1)); 
            end 
        end 
            
        if (copy > 1) && ~strcmp(params.type, 'aortic')  
            params = place_cross_layer_springs(params); 
        end 
        
    end 

    
    if strcmp(params.type, 'aortic') 
        
        params = place_cross_layer_springs_aortic(params, valve.leaflets(1)); 
        
        if isfield(params, 'k_bend_cross_layer') && params.cross_layer_on
            params = place_cross_layer_beams_aortic(params, valve.leaflets(1)); 
        end 
        
        if isfield(valve, 'k_bend_nodule') && isfield(valve, 'k_bend_nodule_length') && params.beam_on
            params = place_nodule_beams(params, valve.leaflets(1), valve.k_bend_nodule_length, valve.k_bend_nodule);
        end 
        
        if isfield(valve, 'place_cylinder') && valve.place_cylinder
            
            params_cylinder.vertex        = fopen(strcat(base_name, '_cylinder.vertex'), 'w'); 
            params_cylinder.spring        = fopen(strcat(base_name, '_cylinder.spring'), 'w'); 
            params_cylinder.target        = fopen(strcat(base_name, '_cylinder.target'), 'w'); 
            params_cylinder.total_vertices  = 0; 
            params_cylinder.total_springs   = 0; 
            params_cylinder.total_targets   = 0; 
            params_cylinder.global_idx = 0; 
            params_cylinder.num_copies = params.num_copies; 
            
            if ~isfield(valve, 'z_min_cylinder')
                valve.z_min_cylinder = params.z_min;                
            end 
            if ~isfield(valve, 'z_max_cylinder')
                valve.z_max_cylinder = params.z_max;                
            end 
            
            if ~isfield(valve, 'n_layers_cylinder')
                valve.n_layers_cylinder = 3;                
            end
            
            if isfield(valve, 'tight_cylinder')
                tight_cylinder = valve.tight_cylinder; 
            else 
                tight_cylinder = false; 
            end
            
            if isfield(valve, 'z_extra_cylinder')
                z_extra_cylinder = valve.z_extra_cylinder; 
            else 
                z_extra_cylinder = 0; 
            end 
            
            if isfield(valve, 'r_max_cylinder')
                r_max_cylinder = valve.r_max_cylinder; 
            else 
                r_max_cylinder = inf; 
            end 
            
            if isfield(valve.skeleton, 'r_of_z')
                % r that follows z 
                params_cylinder = place_cylinder(params_cylinder, valve.leaflets(1), r, ds, valve.z_min_cylinder, valve.z_max_cylinder, valve.n_layers_cylinder, k_rel, k_target_net, eta_net, tight_cylinder, z_extra_cylinder, r_max_cylinder, valve.skeleton.r_of_z); 
            else 
                % constant r 
                params_cylinder = place_cylinder(params_cylinder, valve.leaflets(1), r, ds, valve.z_min_cylinder, valve.z_max_cylinder, valve.n_layers_cylinder, k_rel, k_target_net, eta_net, tight_cylinder, z_extra_cylinder, r_max_cylinder); 
            end 
            
            if valve.in_heart && isfield(valve, 'transformation_vertex_file')
                if isfield(valve, 'initial_rotation_aortic')
                    R_0 = valve.initial_rotation_aortic; 
                else
                    R_0 = eye(3);
                end 
                if isfield(valve, 'initial_translation_aortic')
                    T_0 = valve.initial_translation_aortic; 
                else
                    T_0 = zeros(3,1);
                end 
                
                params_cylinder.vertices = coordinate_transformation_vertices(params_cylinder.vertices, valve.transformation_vertex_file, R_0, T_0);
            elseif valve.in_heart && isfield(valve.skeleton, 'inverse_transformation_initial_condition')
                params_cylinder.vertices = valve.skeleton.inverse_transformation_initial_condition(params_cylinder.vertices) 
            else 
                
                if isfield(valve, 'initial_translation_aortic')
                    % translation
                    params_cylinder.vertices = params_cylinder.vertices + valve.initial_translation_aortic;
                end 
                
                if isfield(valve, 'initial_rotation_aortic')
                    % rotation alone 
                    R_0 = valve.initial_rotation_aortic; 
                    params_cylinder.vertices = R_0 * params_cylinder.vertices;
                end 
            end 

            % finally, write all vertices 
            params_cylinder = write_all_vertices(params_cylinder); 

            % and clean up files with totals 
            fclose(params_cylinder.vertex   ); 
            fclose(params_cylinder.spring   ); 
            fclose(params_cylinder.target   ); 

            prepend_line_with_int(strcat(base_name, '_cylinder.vertex'), params_cylinder.total_vertices); 
            prepend_line_with_int(strcat(base_name, '_cylinder.spring'), params_cylinder.total_springs); 
            prepend_line_with_int(strcat(base_name, '_cylinder.target'), params_cylinder.total_targets);

        end
        
    end 
    
    
    
    
    if n_lagrangian_tracers > 0
        
        double_z = true; 
        
        params_particles.vertex = fopen(strcat(base_name, '_particles.vertex'), 'w'); 
        params_particles.total_vertices  = 0; 
        params_particles.global_idx = 0; 
        params_particles.num_copies = params.num_copies; 
        
        params_particles.x_min = params.x_min; 
        params_particles.x_max = params.x_max; 
        params_particles.y_min = params.y_min; 
        params_particles.y_max = params.y_max; 
        params_particles.z_min = params.z_min; 
        params_particles.z_max = params.z_max;
        
        
        [params_particles, total_lagrangian_placed] = place_lagrangian_tracers(params_particles, n_lagrangian_tracers, double_z); 
        
        if valve.in_heart && isfield(valve, 'transformation_vertex_file')
            if isfield(valve, 'initial_rotation_aortic')
                R_0 = valve.initial_rotation_aortic; 
            else
                R_0 = eye(3);
            end 
            if isfield(valve, 'initial_translation_aortic')
                T_0 = valve.initial_translation_aortic; 
            else
                T_0 = zeros(3,1);
            end 
            params_particles.vertices = coordinate_transformation_vertices(params_particles.vertices, valve.transformation_vertex_file, R_0, T_0);
        elseif valve.in_heart && isfield(valve.skeleton, 'inverse_transformation_initial_condition')
            params_particles.vertices = valve.skeleton.inverse_transformation_initial_condition(params_particles.vertices); 
        else 
                
                if isfield(valve, 'initial_translation_aortic')
                    % translation
                    params_particles.vertices = params_particles.vertices + valve.initial_translation_aortic;
                end 
                
                if isfield(valve, 'initial_rotation_aortic')
                    % rotation alone 
                    R_0 = valve.initial_rotation_aortic; 
                    params_particles.vertices = R_0 * params_particles.vertices;
                end
        end 
        
        params_particles = write_all_vertices(params_particles); 
        
        fclose(params_particles.vertex); 
        prepend_line_with_int(strcat(base_name, '_particles.vertex'), params_particles.total_vertices); 

    end 
    
    if valve.in_heart && isfield(valve, 'transformation_vertex_file')
        if isfield(valve, 'initial_rotation_aortic')
            R_0 = valve.initial_rotation_aortic; 
        else
            R_0 = eye(3);
        end 
        if isfield(valve, 'initial_translation_aortic')
            T_0 = valve.initial_translation_aortic; 
        else
            T_0 = zeros(3,1);
        end 
        
        params.vertices = coordinate_transformation_vertices(params.vertices, valve.transformation_vertex_file, R_0, T_0); 
    elseif valve.in_heart && isfield(valve.skeleton, 'inverse_transformation_initial_condition')
        params.vertices = valve.skeleton.inverse_transformation_initial_condition(params.vertices) 
    else 

        if isfield(valve, 'initial_translation_aortic')
            % translation
            params.vertices = params.vertices + valve.initial_translation_aortic;
        end 

        if isfield(valve, 'initial_rotation_aortic')
            % rotation alone 
            R_0 = valve.initial_rotation_aortic; 
            params.vertices = R_0 * params.vertices;
        end
    end 
    
    % finally, write all vertices 
    params = write_all_vertices(params); 

    % and clean up files with totals 
    fclose(params.vertex   ); 
    fclose(params.spring   ); 
    fclose(params.target   ); 
    fclose(params.inst     ); 
    if ~strcmp(params.type, 'aortic') 
        fclose(params.papillary); 
    end 
    if params.beam_on 
        fclose(params.beam);     
    end 
    
    prepend_line_with_int(strcat(base_name, '.vertex'), params.total_vertices); 
    prepend_line_with_int(strcat(base_name, '.spring'), params.total_springs); 
    prepend_line_with_int(strcat(base_name, '.target'), params.total_targets); 
    if ~strcmp(params.type, 'aortic') 
        prepend_line_with_int(strcat(base_name, '.papillary'), params.total_papillary); 
    end 
    if params.beam_on
        prepend_line_with_int(strcat(base_name, '.beam'), params.total_beams); 
    end 
    
    if params.targets_for_bcs
        fclose(params.vertex_target_pos_file); 
        prepend_line_with_int(strcat(base_name, '_target_position.vertex'), params.total_vertices_target); 
        
        if params.total_vertices_target ~= params.total_vertices
            error('total vertices in target position file and basic position file are not consistent')
        end 
        
    end
    

end 


function vertex_string(coords, file)
    % prints formatted string for current vertex to vertex file   
    fprintf(file, '%.14f\t %.14f\t %.14f\n', coords(1), coords(2), coords(3)); 
end

function params = write_all_vertices(params)
    % writes all vertices to file 

    max_idx = params.global_idx;
    fprintf('max_idx = %d\n', max_idx); 
    
    debug = false; 
    if debug 
        min_x = min(params.vertices(1,1:max_idx)) 
        max_x = max(params.vertices(1,1:max_idx)) 
        min_y = min(params.vertices(2,1:max_idx)) 
        max_y = max(params.vertices(2,1:max_idx)) 
        min_z = min(params.vertices(3,1:max_idx)) 
        max_z = max(params.vertices(3,1:max_idx)) 
    end 
    
    for i=1:max_idx
        vertex_string(params.vertices(:,i), params.vertex); 
        params.total_vertices = params.total_vertices + 1; 
    end 
    
    if isfield(params, 'targets_for_bcs') && params.targets_for_bcs
        
        nan_indices = isnan(params.vertices_target); 
        params.vertices_target(nan_indices) = params.vertices(nan_indices); 
        
        for i=1:max_idx
            vertex_string(params.vertices_target(:,i), params.vertex_target_pos_file); 
            params.total_vertices_target = params.total_vertices_target + 1; 
        end 
    end 
    
end 


function params = spring_string(params, idx, nbr, kappa, rest_len, function_idx, output)
    % prints a spring format string to string file 
    if nbr <= idx
        error('By convention, only place springs with the second index larger to prevent duplicates'); 
    end 
    
    if isfield(params, 'copy_spring_weights')
        weight = params.copy_spring_weights(params.copy); 
        fprintf(params.spring, '%d\t %d\t %.14f\t %.14f', idx, nbr, kappa * weight, rest_len); 
    else
        % default even division 
        fprintf(params.spring, '%d\t %d\t %.14f\t %.14f', idx, nbr, kappa/params.num_copies, rest_len); 
    end 
    % index for custom spring functions 
%     if ~exist('function_idx', 'var') 
%         function_idx = 0; 
%     end 
   
    fprintf(params.spring, '\t %d', function_idx); 
    
%     if function_idx == 0
%         eta = params.eta_multiplier_linear   * kappa / params.num_copies; 
%     elseif function_idx == 1
%         eta = params.eta_multiplier_collagen * kappa / params.num_copies; 
%     else 
%         error('Only linear (default) and collagen function indices implemented'); 
%     end 
%    
%     fprintf(params.spring, '\t %.14f', eta); 


    fprintf(params.spring, '\t # %d', output); 

    fprintf(params.spring, '\n'); 

    params.total_springs = params.total_springs + 1; 
end 


function params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, collagen_spring, output)
    % 
    % Add one or more springs 
    % If the rest length is more than 2 times the specificed mesh width
    % Then it is split into multiple springs 
    % 
    % Note that the constant passed in must be the relative spring constant 
    % 

    if nbr_idx <= idx
        error('By convention, only place springs with the second index larger to prevent duplicates'); 
    end 
    
    if strcmp(collagen_spring, 'aortic_rad')
        function_idx = 3;
    elseif strcmp(collagen_spring, 'aortic_circ')
        function_idx = 2;
    elseif collagen_spring
        % default mitral 
        function_idx = 1;
    else 
        function_idx = 0;
    end 
    
    if k_rel == 0 
        warning(sprintf('Zero strength spring on idx,nbr = %d,%d, not placed.', idx, nbr_idx)); 
        return; 
    end 

    % max_strain = .01; 
    
    
    X     = params.vertices(:,idx + 1); 
    X_nbr = params.vertices(:,nbr_idx + 1); 
    L     = norm(X_nbr - X); 
    
    % Can use either rest length or current length to determine number of springs 
    % N_springs = floor(rest_len / ds); 
    N_springs = floor(L / ds); 
    
    strain = (L - rest_len) / rest_len; 
    
    % fprintf('strain = %e, idx = %d, nbr = %d\n', strain, idx, nbr_idx)
    
    % 2x expected loaded strains 
    max_strain_circ = .15; 
    max_strain_rad = 2*.54; 
    
%     [j_idx, k_idx] = find(params.layer_indices(params.copy).indices_global == idx); 
%     [j_nbr_idx, k_nbr_idx] = find(params.layer_indices(params.copy).indices_global == nbr_idx); 
%     if k_idx == 65
%         fprintf('copy %d, %s, strain = %f, idx = %d, j = %d, k = %d, nbr = %d, j_nbr = %d, k_nbr = %d\n', params.copy, collagen_spring, strain, idx, j_idx, k_idx, nbr_idx, j_nbr_idx, k_nbr_idx);
%     end 
    
    
    if (strcmp(collagen_spring, 'aortic_circ') && (strain > max_strain_circ)) || ...
       (strcmp(collagen_spring, 'aortic_rad') && (strain > max_strain_rad)) 
        
        [j_idx, k_idx] = find(params.layer_indices(params.copy).indices_global == idx); 
        [j_nbr_idx, k_nbr_idx] = find(params.layer_indices(params.copy).indices_global == nbr_idx); 
        
        fprintf('copy %d, %s, strain = %f, idx = %d, j = %d, k = %d, nbr = %d, j_nbr = %d, k_nbr = %d\n', params.copy, collagen_spring, strain, idx, j_idx, k_idx, nbr_idx, j_nbr_idx, k_nbr_idx);
    end 
    
    % Just one spring placed here 
    if N_springs <= 1 
        if collagen_spring
            % Scaling constant here does not change with rest lengths 
            k_col = k_rel; 

            % Finally, write the spring string 
            params = spring_string(params, idx, nbr_idx, k_col, rest_len, function_idx, output); 
        else 
            k_abs = k_rel / rest_len; 
            params = spring_string(params, idx, nbr_idx, k_abs, rest_len, function_idx, output); 
        end 
    else 
        
        % increment of new points 
        step = (X_nbr - X)/N_springs; 
        
        X_vertices_new = zeros(3,N_springs + 1); 
        
        % set new coordinates 
        for i=0:N_springs
            X_vertices_new(:,i+1) = X + i*step; 
        end 
        
        if X ~= X_vertices_new(:,1)
            error('First vertex must be equal to X'); 
        end 
        
        if X_nbr ~= X_vertices_new(:,N_springs + 1)
            error('Last vertex must be equal to X_nbr'); 
        end
        
        
        % first index is always first index suppied 
        idx_tmp = idx; 
        
        for i=1:N_springs 
            
            % place the upper new vertex if it is not placed already 
            if i<N_springs 
                params.vertices(:,params.global_idx + 1) = X_vertices_new(:,i+1); 
                
                % new vertex's index is current (zero indexed) global index 
                nbr_idx_tmp = params.global_idx; 

                params.global_idx = params.global_idx + 1;
                
            else 
                % last spring upper limit is the previous neighbor 
                nbr_idx_tmp = nbr_idx; 
            end 
            
            % global indices placed in order by convention 
            min_idx = min(idx_tmp, nbr_idx_tmp); 
            max_idx = max(idx_tmp, nbr_idx_tmp); 
        
            % Current length 
            L = norm(X_vertices_new(:,i+1) - X_vertices_new(:,i)); 
            
            % Rest length determined by strain 
            R = L / (strain + 1); 
            
            if collagen_spring
                % Scaling constant here does not change with rest lengths 
                k_col = k_rel; 

                % Finally, write the spring string 
                params = spring_string(params, min_idx, max_idx, k_col, R, function_idx, output);
            
            else 
                % Absolute spring constant must be used in spring file 
                k_abs = k_rel / R; 

                % Finally, write the spring string 
                params = spring_string(params, min_idx, max_idx, k_abs, R, function_idx, output);
            end 
            
            % lower index is always previous upper index 
            idx_tmp = nbr_idx_tmp; 
            
        end 
    
    end 

end 


function params = target_string(params, idx, kappa, eta)
    % prints a target format string to target file 
    
    if exist('eta', 'var') && (eta > 0.0)
        fprintf(params.target, '%d\t %.14f\t %.14f\n', idx, kappa/params.num_copies, eta/params.num_copies);
    else
        fprintf(params.target, '%d\t %.14f\n', idx, kappa/params.num_copies);
    end 
    params.total_targets = params.total_targets + 1; 
end 


function params = papillary_string(params, idx, coords)
    % prints a papillary format string to target file 
    
    if isfield(params, 'inverse_transformation_initial_condition')
        coords = params.inverse_transformation_initial_condition(coords); 
    end 
    
    fprintf(params.papillary, '%d\t %.14f\t %.14f %.14f\n', idx, coords(1), coords(2), coords(3) + params.z_offset);
    params.total_papillary = params.total_papillary + 1; 
end 

function params = beam_string(params, idx_minux, idx, idx_plus, k_bend)
    % prints a beam format string to target file 
    
    fprintf(params.beam, '%d\t %d\t %d\t %.14f\n', idx_minux, idx, idx_plus, k_bend);
    params.total_beams = params.total_beams + 1; 
end 

function [] = prepend_line_with_int(file_name, val)
    % Adds a single line to the file with given name
    % at the beginning with the integer val 
    % writes a temp file then calls cat 

    
    write = fopen('temp.txt', 'w'); 
    fprintf(write, '%d\n', val); 
    fclose(write); 
    
    file_temp = strcat(file_name, '.tmp'); 
    
    system(sprintf('cat temp.txt %s > %s', file_name, file_temp));
    movefile(file_temp, file_name); 
    system('rm temp.txt'); 

end 



function [params leaflet] = assign_indices_vertex_target(params, leaflet, k_target_net, k_target_papillary, eta_net, eta_papillary)
    % 
    % Assigns global indices to the leaflet and chordae 
    % Includes all boundary condition points and internal 
    % 
    % 
    % Places all main data into IBAMR format for this leaflet
    % Updates running totals on the way 

    
    % Unpack needed data 
    X                 = leaflet.X; 
    j_max             = leaflet.j_max; 
    k_max             = leaflet.k_max; 
    is_internal       = leaflet.is_internal;
    is_bc             = leaflet.is_bc; 
    if ~strcmp(params.type, 'aortic')
        num_trees         = leaflet.num_trees; 
        chordae           = leaflet.chordae; 
    end 
    
    if strcmp(params.type, 'aortic')
        if ~(isfield(leaflet, 'pre_extrude') && leaflet.pre_extrude)
            if isfield(params, 'normal_thicken') && params.normal_thicken 

                if isfield(params, 'center_extrusion') && params.center_extrusion

                    if params.num_copies ~= 3
                        error('cneter extrusion only implemted with 3 layers') 
                    end 

                    extrude_length = (params.copy - 2) * params.ds_extrude; 
                else 
                    extrude_length = (params.copy - 1) * params.ds_extrude; 
                end 

                X = normal_extrude_aortic(leaflet, extrude_length); 
            end 
        end 
    end 
    
    % Keep track of vector index in 1d array 
   
    for k=1:k_max
        for j=1:j_max
            
            if leaflet.targets_for_bcs

                % every internal point written to the file 
                % targets written in separate array
                if is_internal(j,k)

                    % check all neighbors for boundary conditions 
                    found_target = false; 
                    target_j = []; 
                    target_k = []; 
                    for j_nbr_tmp = [j-1,j+1]
                        k_nbr_tmp = k;                 
                        [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                        if target_spring
                            if found_target
                                error('Found two targets attached to same location, should be impossible'); 
                            end 
                            
                            found_target = true; 
                            target_j = j_nbr; 
                            target_k = k_nbr; 
                        end 
                    end 
                        
                    for k_nbr_tmp = [k-1,k+1]
                        j_nbr_tmp = j; 
                        [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                    
                        if target_spring
                            if found_target
                                error('Found two targets attached to same location, should be impossible'); 
                            end 
                            
                            found_target = true; 
                            target_j = j_nbr; 
                            target_k = k_nbr; 
                        end 
                    end 
                    
                    params.vertices(:,params.global_idx + 1) = X(:,j,k); 
                    
                    if found_target
                        params.vertices_target(:,params.global_idx + 1) = X(:,target_j,target_k); 
                    end 
                    leaflet.indices_global(j,k) = params.global_idx; 

                    % if flag set, this is a target point 
                    if found_target
                        if exist('eta_net', 'var')
                            params = target_string(params, params.global_idx, k_target_net, eta_net);     
                        else
                            params = target_string(params, params.global_idx, k_target_net);     
                        end 
                    end 

                    params.global_idx = params.global_idx + 1;
                end 
                
            else % ~leaflet.targets_for_bcs 
                
                % default behavior, bcs are targets with current position 
                
                % every internal and boundary point written to the file 
                if is_internal(j,k) || is_bc(j,k)

                    params.vertices(:,params.global_idx + 1) = X(:,j,k); 
                    leaflet.indices_global(j,k) = params.global_idx; 

                    % if on boundary, this is a target point 
                    if is_bc(j,k)
                        if exist('eta_net', 'var')
                            params = target_string(params, params.global_idx, k_target_net, eta_net);     
                        else
                            params = target_string(params, params.global_idx, k_target_net);     
                        end 
                    end 

                    params.global_idx = params.global_idx + 1;
                end 
                
            end 
        end 
    end
    
    if ~strcmp(params.type, 'aortic')
        % chordae internal terms 
        for tree_idx = 1:num_trees

            C = chordae(tree_idx).C; 
            [m N_chordae] = size(C);         

            % always place root index first 

            if leaflet.targets_for_bcs
                % split location of targets and their current position   
                params.vertices       (:,params.global_idx + 1) = leaflet.chordae(tree_idx).root;
                params.vertices_target(:,params.global_idx + 1) = leaflet.chordae(tree_idx).root_target;
            else
                % default behavior with only one position 
                params.vertices(:,params.global_idx + 1) = leaflet.chordae(tree_idx).root;
            end 

            leaflet.chordae(tree_idx).idx_root = params.global_idx;                 

            % root is always a boundary condition 
            if exist('eta_papillary', 'var')
                params = target_string(params, params.global_idx, k_target_papillary, eta_papillary);     
            else
                params = target_string(params, params.global_idx, k_target_papillary);     
            end 

            % write papillary file 
            params = papillary_string(params, params.global_idx, leaflet.chordae(tree_idx).root); 

            params.global_idx = params.global_idx + 1;

            for i=1:N_chordae
                params.vertices(:,params.global_idx + 1) = C(:,i); 
                leaflet.chordae(tree_idx).indices_global(i) = params.global_idx;                 
                params.global_idx = params.global_idx + 1;
            end 

        end 
    end 

end 
 

function params = add_springs(params, leaflet, ds, collagen_spring)

    params = add_leaflet_springs(params, leaflet, ds, collagen_spring); 
    if ~strcmp(params.type, 'aortic')
        params = add_chordae_tree_springs(params, leaflet, ds, collagen_spring); 
    end 
    
end 


function params = add_leaflet_springs(params, leaflet, ds, collagen_spring)
                      
    % Places all main data into IBAMR format for this leaflet
    % Updates running totals on the way 

    % Unpack needed data 
    j_max             = leaflet.j_max; 
    k_max             = leaflet.k_max; 
    is_internal       = leaflet.is_internal;
    is_bc             = leaflet.is_bc; 
    
    if ~strcmp(params.type, 'aortic')
        chordae           = leaflet.chordae;
        chordae_idx       = leaflet.chordae_idx; 
        collagen_constitutive_circ = collagen_spring; 
        collagen_constitutive_rad  = collagen_spring; 
    else 
        collagen_constitutive_circ = leaflet.collagen_constitutive_circ; 
        collagen_constitutive_rad  = leaflet.collagen_constitutive_rad; 
        
        if strcmp(collagen_constitutive_circ, 'linear_no_compressive')
            collagen_constitutive_circ = false; 
        end 
        if strcmp(collagen_constitutive_rad, 'linear_no_compressive')
            collagen_constitutive_rad = false; 
        end 
    end 
    
    
    
    R_u               = leaflet.R_u;
    k_u               = leaflet.k_u;
    R_v               = leaflet.R_v;
    k_v               = leaflet.k_v;

       
    if isfield(leaflet, 'periodic_j')
        periodic_j = leaflet.periodic_j; 
    else
        periodic_j = zeros(k_max,1); 
    end 
   
    % output flag information 
    copy = params.copy; 
    if params.output.leaflets(copy)
        output        = true; 
        output_stride = params.output.stride_leaflet; 
    else 
        output        = false; 
    end 
    
    if ~strcmp(params.type, 'aortic')
        if params.output.chordae(copy)
            output_tmp_chordae = true; 
        end 
    end 
    
    for k=1:k_max
        for j=1:j_max
            
            % every internal and boundary point may have springs connected to it 
            if is_internal(j,k) || (is_bc(j,k) && ~leaflet.targets_for_bcs)
                
                if output && ((mod(j,output_stride) == 1) || (output_stride == 1))
                    output_tmp_k = true; 
                else 
                    output_tmp_k = false; 
                end 
                
                if output && ((mod(k,output_stride) == 1) || (output_stride == 1))
                    output_tmp_j = true; 
                else 
                    output_tmp_j = false; 
                end                 
                
                % global index of current point 
                idx = leaflet.indices_global(j,k); 
                
                % current node has a chordae connection
                if ~strcmp(params.type, 'aortic') && chordae_idx(j,k).tree_idx
                    
                    tree_idx = chordae_idx(j,k).tree_idx; 
                    
                    [m N_chordae] = size(chordae(tree_idx).C);

                    % index in current free edge array 
                    i = chordae_idx(j,k).leaf_idx;
                    
                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet
                    leaf_idx = chordae_idx(j,k).leaf_idx + N_chordae;

                    % then take the parent index of that number in chordae variables
                    idx_chordae = floor(leaf_idx/2);  
                    
                    nbr_idx = chordae(tree_idx).indices_global(idx_chordae); 

                    rest_len = chordae(tree_idx).R_free_edge(i); 

                    k_rel = chordae(tree_idx).k_free_edge(i); 
                    
                    params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, collagen_spring, output_tmp_chordae); 
                    
                end 
                                
                % springs in leaflet, only go in up direction 
                j_nbr_tmp = j + 1; 
                k_nbr_tmp = k; 
                [valid j_nbr k_nbr j_spr k_spr target_spring target_k_no_j_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                if valid && (~target_spring) && (~target_k_no_j_spring)
                    
                    % no bc to bc springs 
                    if ~(is_bc(j, k) && is_bc(j_nbr, k_nbr))
                    
                        rest_len = R_u(j_spr, k_spr); 
                        k_rel    = k_u(j_spr, k_spr); 

                        nbr_idx = leaflet.indices_global(j_nbr,k_nbr);
                        
                        
                        if j_nbr_tmp ~= j_nbr 
                            % periodic wrapping requires oppositite order 
                            params = place_spring_and_split(params, nbr_idx, idx, k_rel, rest_len, ds, collagen_constitutive_circ, output_tmp_j);
                        else 
                           % standard order 
                            params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, collagen_constitutive_circ, output_tmp_j);
                        end 

                    end 

                end 
                
                % springs in leaflet, only go in up direction 
                j_nbr_tmp = j; 
                k_nbr_tmp = k + 1; 
                [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                if valid && (~target_spring)
                    
                    % no bc to bc springs 
                    if ~(is_bc(j, k) && is_bc(j_nbr, k_nbr))
                    
                        rest_len = R_v(j_spr, k_spr); 
                        k_rel    = k_v(j_spr, k_spr); 

                        nbr_idx = leaflet.indices_global(j_nbr,k_nbr);

                        params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, collagen_constitutive_rad, output_tmp_k);

                        if isfield(params, 'compressive_extra_springs_on') && params.compressive_extra_springs_on
                            if (k/k_max) >= (1 - params.kappa_radial_free_edge_compressive_percentage)
                                k_compressive = params.kappa_radial_free_edge_compressive; 
                                rest_len_compressive = params.kappa_radial_free_edge_compressive_stretch * rest_len; 
                                compressive_fn_idx = params.kappa_radial_free_edge_compressive_fn_idx;  
                                output_temp_compressive = false; 
                                params = spring_string(params, idx, nbr_idx, k_compressive, rest_len_compressive, compressive_fn_idx, output_temp_compressive); 
                            end 
                        end 
                                                
                    end 

                end 
                
            end 
        end 
    end
end 


function params = add_chordae_tree_springs(params, leaflet, ds, collagen_spring)
    % 
    % Adds chordae tree to IBAMR format files 
    % 
    
    num_trees = leaflet.num_trees; 
    chordae   = leaflet.chordae; 
        
    copy = params.copy; 
    if params.output.chordae(copy)
        output_tmp = true;
    else 
        output_tmp = false;
    end 
    
    % chordae internal terms 
    for tree_idx = 1:num_trees
        
        C = chordae(tree_idx).C; 
        [m N_chordae] = size(C);   
        
        indices_global = chordae(tree_idx).indices_global; 
        
        for i=1:N_chordae

            idx = indices_global(i); 
            
            % place the spring which goes to the parent 
            parent = floor(i/2); 
            
            if parent == 0
                nbr_idx = chordae(tree_idx).idx_root; 
            else
                nbr_idx = indices_global(parent); 
            end 
            
            % get the neighbors coordinates, reference coordinate and spring constants
            [nbr rest_len k_rel] = get_nbr_chordae(leaflet, i, parent, tree_idx); 
            
            % list nbr index first because nbr is parent and has lower index
            params = place_spring_and_split(params, nbr_idx, idx, k_rel, rest_len, ds, collagen_spring, output_tmp);
                
        end 
        
    end 

end 


function params = place_cross_layer_springs(params)

    function_idx = 0; 
    kappa        = params.kappa_cross_layer; 
    rest_len     = params.rest_len_cross_layer;
    % don't include these for now 
    output_tmp = false; 
        
    debug = true; 
    lenient_tol = true; 
    if lenient_tol
        tol = 1e-3; 
    else
        tol = eps; 
    end
    err_count = 0; 
    
    fprintf('Cross layer spring range     (inclusive, zero based for ibamr files) = %d:%d\n', params.min_idx_for_cross_layer, params.max_idx_for_cross_layer); 
    fprintf('connect to range in previous (inclusive, zero based for ibamr files) = %d:%d\n', params.min_idx_for_cross_layer -  params.total_per_layer, params.max_idx_for_cross_layer -  params.total_per_layer); 
    
    n_to_place = length((params.min_idx_for_cross_layer):(params.max_idx_for_cross_layer)); 
    
    if isfield(params, 'num_cross_layer_placed')
        if params.num_cross_layer_placed ~= n_to_place
            error('number of cross links must be the same on each layer')
        end 
    else         
        params.num_cross_layer_placed = n_to_place; 
    end 
    
    for i=(params.min_idx_for_cross_layer):(params.max_idx_for_cross_layer)
        idx          = i -  params.total_per_layer; 
        nbr_idx      = i;  

        if debug
           coords = params.vertices(:,idx+1); 
           coords_nbr = params.vertices(:,nbr_idx+1);
           
           if abs(norm(coords - coords_nbr) - rest_len) > tol
               if err_count < 5
                   fprintf('nonrest len spring, rest_len = %f, current_len = %f\n', rest_len, norm(coords - coords_nbr)); 
                   fprintf('    idx = %d, coords = %f %f %f\n', idx, coords(1), coords(2), coords(3)); 
                   fprintf('nbr_idx = %d, coords = %f %f %f\n',nbr_idx, coords_nbr(1), coords_nbr(2), coords_nbr(3)); 
                   fprintf('\n'); 
                   err_count = err_count + 1; 
               end 
           end 
           
        end 
        
        params = spring_string(params, idx, nbr_idx, kappa, rest_len, function_idx, output_tmp); 
    end 

end 

function params = place_cross_layer_springs_aortic(params, leaflet)

    if params.num_copies <= 1 
        return
    end 
    
    function_idx = 0; 
    kappa        = params.kappa_cross_layer; 
    rest_len     = params.rest_len_cross_layer;
    % don't include these for now 
    output_tmp = false; 
        
    debug = true; 
    lenient_tol = true; 
    if lenient_tol
        tol = 1e-3; 
    else
        tol = eps; 
    end
    err_count = 0; 
    
    j_max = leaflet.j_max; 
    k_max = leaflet.k_max; 
    
    for copy = 1:(params.num_copies-1)
        for j=1:j_max
            for k=1:k_max
                idx          = params.layer_indices(copy  ).indices_global(j,k); 
                nbr_idx      = params.layer_indices(copy+1).indices_global(j,k);  

                % when comms are attached and normal is undefined
                % do not place these springs 
                if isfield(leaflet, 'N_to_pinch') && (leaflet.N_to_pinch > 0)
                    if (mod(j, leaflet.N_each) == 0)
                        continue; 
                    end 
                end 
                
                % override rest length to the initial length 
                if isfield(leaflet, 'pre_extrude') && leaflet.pre_extrude
                    coords = params.vertices(:,idx+1); 
                    coords_nbr = params.vertices(:,nbr_idx+1);
                    rest_len = norm(coords - coords_nbr); 
                end 
                
                if debug
                    coords = params.vertices(:,idx+1); 
                    coords_nbr = params.vertices(:,nbr_idx+1);

                    if abs(norm(coords - coords_nbr) - rest_len) > tol
                        if err_count < 5
                            fprintf('nonrest len spring, rest_len = %f, current_len = %f\n', rest_len, norm(coords - coords_nbr)); 
                            fprintf('    idx = %d, coords = %f %f %f\n', idx, coords(1), coords(2), coords(3)); 
                            fprintf('nbr_idx = %d, coords = %f %f %f\n',nbr_idx, coords_nbr(1), coords_nbr(2), coords_nbr(3)); 
                            fprintf('\n'); 
                            err_count = err_count + 1; 
                        end 
                    end 

                end 

                params = spring_string(params, idx, nbr_idx, kappa, rest_len, function_idx, output_tmp);
            end
        end 
    end 
end 


function params = place_cross_layer_beams_aortic(params, leaflet)

    if params.num_copies ~= 3
        warning('must use three copies for cross layer beams'); 
        return
    end 
    
    k_bend_cross_layer = params.k_bend_cross_layer; 
    
    j_max = leaflet.j_max; 
    k_max = leaflet.k_max; 
    
    for j=1:j_max
        for k=1:k_max
            idx_minus = params.layer_indices(1).indices_global(j,k); 
            idx       = params.layer_indices(2).indices_global(j,k);  
            idx_plus  = params.layer_indices(3).indices_global(j,k); 

            params = beam_string(params, idx_minus, idx, idx_plus, k_bend_cross_layer);
            
        end
    end 
    
end 

function params = place_nodule_beams(params, leaflet, k_bend_nodule_length, k_bend_nodule)

    if params.num_copies ~= 3
        warning('must use three copies for cross layer beams'); 
    end 

    j_max              = leaflet.j_max; 
    k_max              = leaflet.k_max; 
    is_internal        = leaflet.is_internal;
    R_v                = leaflet.R_v; 
    R_u                = leaflet.R_u; 
    N_each             = leaflet.N_each; 
    center_leaflet_idx = N_each/2; 
    
    k_min_bending_set = false; 
    j = center_leaflet_idx; 
    height = 0; 
    for k=(k_max:-1:1)
        j_nbr_tmp = j; 
        k_nbr_tmp = k-1; 
        [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
        
        height = height + R_v(j_spr,k_spr); 
        
        if height > k_bend_nodule_length
            k_min_bending_set = true; 
            k_min_bending = k; 
            k_center = 0.5 * (k_min_bending + k_max); 
            break 
        end 
    end 
    
    if ~k_min_bending_set
        error('boundary for bending not set'); 
    end 
    
    j_range_bending_set = false; 
    length = 0; 
    k=k_max;
    for j=(center_leaflet_idx+1):N_each
        j_nbr_tmp = j-1; 
        k_nbr_tmp = k; 
        [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

        length = length + R_u(j_spr,k_spr);

        if length > (k_bend_nodule_length/2)
            j_max_one_leaflet_bending = j; 
            j_min_one_leaflet_bending = N_each - j; % symmetric over center index  
            j_center = center_leaflet_idx; 
            j_range_bending_set = true; 
            break 
        end
    end 
    
    if ~j_range_bending_set
        error('boundary for bending not set'); 
    end 
    
    masks = zeros(size(is_internal)); 
    
    for j=1:j_max
        for k=1:k_max
            j_this_leaflet = mod(j,N_each); 
            
            if k==k_max
                'pause'; 
            end 
            
            if (j_min_one_leaflet_bending <= j_this_leaflet) && (j_this_leaflet <= j_max_one_leaflet_bending) && (k_min_bending <= k) 
                
                % fraction of distance from center of region
                j_ref_distance_normalized = abs(j_this_leaflet - j_center) / (j_max_one_leaflet_bending - j_center); 
                k_ref_distance_normalized = abs(k              - k_center) / (k_max                     - k_center); 
                
                masks(j,k) = cos(j_ref_distance_normalized * pi/2) * cos(k_ref_distance_normalized * pi/2); 
            end 
        end 
    end 
    
    for k=1:k_max
        for j=1:j_max
            % beams only centered at internal points 
            if is_internal(j,k) && (masks(j,k) ~= 0)

                k_bend_temp = masks(j,k) * k_bend_nodule; 
                
                if k_bend_temp > 1e-10
                    % skip effectively zero springs 
                
                    for layer_idx = 1:params.num_copies
                        % global index of current point 
                        idx = params.layer_indices(layer_idx).indices_global(j,k); 

                        % circ direction 
                        j_minus_tmp = j - 1; 
                        k_minus_tmp = k; 
                        [valid_minus j_minus k_minus] = get_indices(leaflet, j, k, j_minus_tmp, k_minus_tmp); 
                        if valid_minus
                            idx_minus = params.layer_indices(layer_idx).indices_global(j_minus,k_minus); 
                        end 

                        j_plus_tmp  = j + 1;
                        k_plus_tmp  = k; 
                        [valid_plus j_plus k_plus] = get_indices(leaflet, j, k, j_plus_tmp, k_plus_tmp); 
                        if valid_plus
                            idx_plus = params.layer_indices(layer_idx).indices_global(j_plus,k_plus); 
                        end 

                        % both neighbors must be valid 
                        if valid_minus && valid_plus                        
                            if (k_minus ~= k) || (k ~= k_plus)
                                error('k indices shuold not change when placing circ (j) beam');                            
                            end              
                            params = beam_string(params, idx_minus, idx, idx_plus, k_bend_temp);
                        end

                        % radial direction 
                        j_minus_tmp = j; 
                        k_minus_tmp = k - 1; 
                        [valid_minus j_minus k_minus] = get_indices(leaflet, j, k, j_minus_tmp, k_minus_tmp); 
                        if valid_minus
                            idx_minus = params.layer_indices(layer_idx).indices_global(j_minus,k_minus); 
                        end 

                        j_plus_tmp  = j;
                        k_plus_tmp  = k + 1; 
                        [valid_plus j_plus k_plus] = get_indices(leaflet, j, k, j_plus_tmp, k_plus_tmp); 
                        if valid_plus
                            idx_plus = params.layer_indices(layer_idx).indices_global(j_plus,k_plus); 
                        end 

                        % both neighbors must be valid 
                        if valid_minus && valid_plus                        
                            if (j_minus ~= j) || (j ~= j_plus)
                                error('j indices shuold not change when placing radial (k) beam');                            
                            end    
                            params = beam_string(params, idx_minus, idx, idx_plus, k_bend_temp); 
                        end
                    end 

                    if params.num_copies == 3
                        % cross layer beam 
                        idx_minus = params.layer_indices(1).indices_global(j,k); 
                        idx       = params.layer_indices(2).indices_global(j,k);  
                        idx_plus  = params.layer_indices(3).indices_global(j,k); 
                        params = beam_string(params, idx_minus, idx, idx_plus, k_bend_temp);
                    end
                    
                end 
            end 
        end 
    end
end 

function params = add_beams(params, leaflet, k_bend_radial, k_bend_circ, ...
                                             k_bend_radial_free_edge, k_bend_radial_free_edge_percentage, ...
                                             k_bend_circ_free_edge, k_bend_circ_free_edge_percentage, ...
                                             k_bend_radial_interp_pts)

    j_max             = leaflet.j_max; 
    k_max             = leaflet.k_max; 
    is_internal       = leaflet.is_internal;
    interpolating     = false; 
    
    if ~exist('k_bend_radial_free_edge_percentage', 'var')
        k_bend_radial_free_edge_percentage = 0; 
    elseif ~exist('k_bend_radial_free_edge_percentage', 'var')
        k_bend_radial_free_edge_percentage = 0; 
    end 
    
    if exist('k_bend_radial_interp_pts', 'var')
        if length(k_bend_radial_interp_pts) ~= length(k_bend_radial)
            error('Must use same number of interp points if interpolating'); 
        end 
        
        if k_bend_radial_free_edge_percentage > 0
            warning('Some points in interpolant will be overwritten by free edge value'); 
        end 
        
        interpolating = true; 
        
    end 
    
    debug = false; 
    if debug 
        k_vals_debug = zeros(k_max,1); 
        bend_coeff_debug = zeros(k_max,1); 
    end 
    
    for k=1:k_max
        for j=1:j_max
            % beams only centered at internal points 
            if is_internal(j,k)
                               
                % global index of current point 
                idx = leaflet.indices_global(j,k); 

                if (k_bend_circ > 0) || (k_bend_circ_free_edge > 0)
                    j_minus_tmp = j - 1; 
                    k_minus_tmp = k; 
                    [valid_minus j_minus k_minus] = get_indices(leaflet, j, k, j_minus_tmp, k_minus_tmp); 
                    if valid_minus
                        idx_minus = leaflet.indices_global(j_minus,k_minus); 
                    end 
                    
                    j_plus_tmp  = j + 1;
                    k_plus_tmp  = k; 
                    [valid_plus j_plus k_plus] = get_indices(leaflet, j, k, j_plus_tmp, k_plus_tmp); 
                    if valid_plus
                        idx_plus = leaflet.indices_global(j_plus,k_plus); 
                    end 
                    
                    % both neighbors must be valid 
                    if valid_minus && valid_plus                        
                        if (k_minus ~= k) || (k ~= k_plus)
                            error('k indices shuold not change when placing circ (j) beam');                            
                        end              
                        
                        if (k_bend_circ_free_edge_percentage > 0) && ... 
                           ((k/k_max) >= (1 - k_bend_circ_free_edge_percentage))
                            params = beam_string(params, idx_minus, idx, idx_plus, k_bend_circ_free_edge);
                        elseif k_bend_circ > 0 
                            params = beam_string(params, idx_minus, idx, idx_plus, k_bend_circ);
                        else 
                        end 
                    end
                end 
                
                if any(k_bend_radial > 0)
                    j_minus_tmp = j; 
                    k_minus_tmp = k - 1; 
                    [valid_minus j_minus k_minus] = get_indices(leaflet, j, k, j_minus_tmp, k_minus_tmp); 
                    if valid_minus
                        idx_minus = leaflet.indices_global(j_minus,k_minus); 
                    end 
                    
                    j_plus_tmp  = j;
                    k_plus_tmp  = k + 1; 
                    [valid_plus j_plus k_plus] = get_indices(leaflet, j, k, j_plus_tmp, k_plus_tmp); 
                    if valid_plus
                        idx_plus = leaflet.indices_global(j_plus,k_plus); 
                    end 
                    
                    % both neighbors must be valid 
                    if valid_minus && valid_plus                        
                        if (j_minus ~= j) || (j ~= j_plus)
                            error('j indices shuold not change when placing radial (k) beam');                            
                        end    
                        
                        if (k_bend_radial_free_edge_percentage > 0) && ... 
                           ((k/k_max) >= (1 - k_bend_radial_free_edge_percentage))
                            params = beam_string(params, idx_minus, idx, idx_plus, k_bend_radial_free_edge); 
                            
                            if debug 
                                k_vals_debug(k) = k; 
                                bend_coeff_debug(k) = k_bend_radial_free_edge; 
                            end 
                            
                        elseif interpolating
                            
                            % C1 shape preserving interpolant 
                            k_bend_tmp = interp1(k_bend_radial_interp_pts, k_bend_radial, k/k_max, 'pchip');  
                            if k_bend_tmp > 0
                                % only place nonzero beams 
                                params = beam_string(params, idx_minus, idx, idx_plus, k_bend_tmp); 
                            end 
                            
                            if debug 
                                k_vals_debug(k) = k; 
                                bend_coeff_debug(k) = k_bend_tmp; 
                            end 
                            
                        else
                            params = beam_string(params, idx_minus, idx, idx_plus, k_bend_radial); 
                        end 
                        
                    end
                end 
            end 
        end 
    end
    
    if debug 
        figure; 
        k_vals_debug
        bend_coeff_debug
        semilogy(k_vals_debug, bend_coeff_debug); 
        title('bending stiffness as function of k')
    end 
    
end 

                        
function params = place_net(params, leaflet, ds, r, L, k_rel, k_target, ref_frac, eta, hoop_springs, ray_springs)
    % 
    % Places a polar coordinate mesh in a box 
    % Starts with N points evenly spaced at radius R
    % Spacing ds is determined by this value 
    % Evenly spaced in R, theta spacing increases as radius increases 
    %
    % Mesh is placed with a rectangular topology
    % Springs everywhere internal 
    % All points are made targets
    % 
    % Input
    %     r                Internal radius 
    %     h                Height of valve ring 
    %     L                Box width -- point must have max norm less than L to be included 
    %     N                Number of points placed on the initial circle (valve ring)
    %     spring           Files for writing, must be open already 
    %     vertex 
    %     target
    %     inst 
    %     params.global_idx       Running totals  
    %     total_vertices 
    %     total_springs 
    %     total_targets 
    % 

    if ~exist('hoop_springs', 'var')
        hoop_springs = true; 
    end 
    
    if ~exist('ray_springs', 'var')
        ray_springs = false; 
    end 
    
    X     = leaflet.X; 
    j_max = leaflet.j_max; 
    k_max = leaflet.k_max; 
    
    % output flag information 
    copy = params.copy; 
    if params.output.mesh(copy)
        output        = true; 
        output_stride = params.output.stride_mesh; 
    else 
        output        = false; 
    end 
   
    function_idx = 0; 
    
    % only include those full rings which fit in the domain 
    % always place at least one, which is the ring itself 
    k_max_rings = max(1, floor( (L-r) / ds)); 
        
    points = zeros(3,j_max,k_max_rings);

    % just keep a list of valid indices, mark NAN if out of physical bounds 
    indices_global = zeros(j_max,k_max_rings); 
    
    instrument_idx = 0; 
    
    % compute vertex positions and add to array 
    for k=1:k_max_rings
        for j=1:j_max
            
            % valve ring points from leaflet 
            if strcmp(params.type, 'aortic')
                % ring_pt = X(1:2,j,1);
                % literal circle for aortic partition placement 
                ring_pt = [r*cos(2*pi*(j-1)/j_max); r*sin(2*pi*(j-1)/j_max)];                 
                ring_z  = params.ring_center(3); 
            else 
                ring_pt = X(1:2,j,k_max);
                ring_z  = X(3,j,k_max);
            end
            
            increment = ring_pt - params.ring_center(1:2);
            increment = ds * increment / norm(increment); 
            
            % expand in the direction of a vector from venter of ring to
            % current point 
            coords_horiz = ring_pt + (k-1)*increment; 
              
            % alternative, just scalar mutiply the point 
            % coords_horiz = (1 + (k-1)*ds) * ring_pt; 
                        
            % these might be the same... 
            
            % if one norm is less than L, then the point is within the domain  
            %if norm(coords_horiz, inf) < L 
            if (params.x_min    <= coords_horiz(1)) && ...   
               (coords_horiz(1) <= params.x_max   ) && ...   
               (params.y_min    <= coords_horiz(2)) && ...   
               (coords_horiz(2) <= params.y_max   ) 
                
                % points(:,j,k) = [coords_horiz; params.ring_center(3)]; 
                points(:,j,k) = [coords_horiz; ring_z]; 
                indices_global(j,k) = params.global_idx; 
                params.vertices(:,params.global_idx + 1) = points(:,j,k); 
                
                % every valid vertex is a target point here 
                if exist('eta', 'var')
                    params = target_string(params, params.global_idx, k_target, eta);     
                else
                    params = target_string(params, params.global_idx, k_target);     
                end 

                params.global_idx = params.global_idx + 1;                   

            else 
                points(:,j,k) = NaN * ones(3,1);
                indices_global(j,k) = NaN; 
            end     
            
        end 
        
    end 

    
    % write the instrument file header here 
    fprintf(params.inst, '1   # num meters in file\n'); 
    fprintf(params.inst, 'meter_0   # name\n'); 
    fprintf(params.inst, '%d  # number of meter points\n', j_max); 
    
    % below the first possible point 
    idx = -1; 
    
    for k=1:k_max_rings
        for j=1:j_max
            
            % take mod one be 
            if output && ((mod(k,output_stride) == 1) || (output_stride == 1))
                output_tmp = true;
            else 
                output_tmp = false;
            end 
            
            % just ignore the nan 
            if ~isnan(indices_global(j,k))
 
                last_idx = idx; 
                idx = indices_global(j,k); 
                
                % instrument file on 
                if k == 1
                    fprintf(params.inst, '%d \t0 \t %d\n', idx, instrument_idx); 
                    instrument_idx = instrument_idx + 1; 
                end 
                
                if last_idx >= idx
                    error('should always be placing points in order, something wrong'); 
                end 

                if hoop_springs
                    % check up directions for springs 
                    if j < j_max
                        j_nbr = j+1; 
                        if ~isnan(indices_global(j_nbr,k))
                            rest_len = ref_frac * norm(points(:,j,k) - points(:,j_nbr,k)); 
                            k_abs = k_rel / rest_len;
                            nbr_idx = indices_global(j_nbr,k); 
                            params = spring_string(params, idx, nbr_idx, k_abs, rest_len, function_idx, output_tmp); 
                        end 
                    end 

                    % don't forget the periodic direction in j
                    if j == j_max
                       j_nbr = 1; 
                       % need to make sure that the 1,k point is also not a NaN  
                       if ~isnan(indices_global(j_nbr,k)) 
                           rest_len = ref_frac * norm(points(:,j,k) - points(:,j_nbr,k)); 
                           k_abs = k_rel / rest_len;
                           nbr_idx = indices_global(j_nbr,k); 
                           params = spring_string(params, nbr_idx, idx, k_abs, rest_len, function_idx, output_tmp); 
                       end 
                    end 
                end 
                
                if ray_springs
                    % check up directions for springs 
                    if k < k_max_rings
                        k_nbr = k+1; 
                        if ~isnan(indices_global(j,k_nbr))
                            rest_len = ref_frac * norm(points(:,j,k) - points(:,j,k_nbr)); 
                            k_abs = k_rel / rest_len;
                            nbr_idx = indices_global(j,k_nbr); 
                            params = spring_string(params, idx, nbr_idx, k_abs, rest_len, function_idx, output_tmp); 
                        end 
                    end 
                    
                end
                
            end 
           
        end 
    end  

end 

function params = place_cylinder(params, leaflet, r, ds, z_min, z_max, n_layers, k_rel, k_target, eta, tight_cylinder, z_extra, r_max, r_of_z)
    % Places n_layers cylinders
    % cylinders have axial radial and circumferential fibers 
    % 
            
    if ~exist('tight_cylinder', 'var')
        tight_cylinder = false; 
    end 
    if ~exist('z_extra', 'var')
        z_extra = 0; 
    end 
    
    if ~isa(z_min, 'function_handle')
        if ~isscalar(z_min)
            error('z_min must be function handle or scalar')
        end 
        z_min = @(theta) z_min * ones(size(theta)); 
    end 
    % brute force minmum of function 
    z_min_scalar = min(z_min(linspace(0,2*pi,10000))); 
    
    if isa(z_max, 'function_handle')
        % brute force maximum of function 
        z_max_scalar = max(z_max(linspace(0,2*pi,10000))); 
    else 
        z_max_scalar = z_max; 
    end 
    
    
    
    N_theta = floor(2*pi*r / ds);  
    N_z     = ceil((z_max_scalar - z_min_scalar)/ds); 
    N_r     = n_layers; 
    
    dtheta = 1/N_theta; 
    dz     = ds; %1/N_z; 
    dr     = ds; 
    
    % bottom leaflet ring
    X = leaflet.X; 
    j_max = leaflet.j_max; 
    z_min_leaflet = X(3,:,1);     
    
    theta_leaflet = atan2(X(2,:,1), X(1,:,1));
    theta_leaflet = mod(theta_leaflet, 2*pi); % argument from 0 to 2 pi 
    
    debug = false; 
    if debug 
        figure;  
        plot(theta_leaflet, z_min_leaflet); 
        axis equal 
    end 
    
    points = zeros(3,N_theta,N_r,N_z);
    
    % just keep a list of valid indices, mark NAN if out of physical bounds 
    indices_global = zeros(N_theta,N_r,N_z);     
    
    % always linear springs 
    function_idx = 0; 
    
    % This loop should be one indexed, 
    % because we want to start not at the edge but ds in 
    for z_idx=1:N_z
        for r_idx=1:N_r
            for theta_idx=1:N_theta
                       
                theta_tmp = 2*pi*dtheta*(theta_idx-1); 
                
                z_coord = z_min(theta_tmp) + z_max_scalar * (z_idx - 1) * dz; 
                
                if exist('r_of_z', 'var')
                    r_tmp = r_of_z(z_coord) + (r_idx - 1)*dr; 
                else 
                    r_tmp = r + (r_idx - 1)*dr; 
                end
               
                
                % first check for exact mesh alignment
                theta_leaflet_idx = find((theta_tmp == theta_leaflet), 1); 
                if ~isempty(theta_leaflet_idx)
                    z_annulus = max(z_min_leaflet(theta_leaflet_idx));
                else                    
                    % otherwise search for correct angle                     
                    theta_leaflet_idx = find((theta_tmp < theta_leaflet), 1); 
                    if isempty(theta_leaflet_idx)
                        theta_leaflet_idx = 1; 
                    end 

                    % get minimum relative
                    theta_leaflet_idx_prev = theta_leaflet_idx - 1; 
                    if theta_leaflet_idx_prev == 0 
                        theta_leaflet_idx_prev = j_max; 
                    end
    %                theta_leaflet_idx = theta_idx;

                    if theta_leaflet_idx < j_max
                        theta_leaflet_idx_next = theta_leaflet_idx + 1; 
                    else 
                        theta_leaflet_idx_next = 1; 
                    end 

                    z_annulus = max(z_min_leaflet(theta_leaflet_idx_prev), max(z_min_leaflet(theta_leaflet_idx), z_min_leaflet(theta_leaflet_idx_next))); 
                end 
                
%                 if r_idx==1 && z_idx==1
%                     fprintf('theta_idx = %d,\ttheta_leaflet_idx = %d,\t z_annulus = %f\n', theta_idx, theta_leaflet_idx, z_annulus); 
%                 end 
                
                x_coord = r_tmp * cos(2*pi*dtheta*(theta_idx-1)); 
                y_coord = r_tmp * sin(2*pi*dtheta*(theta_idx-1)); 
                
                if tight_cylinder
                    if isa(z_max, 'function_handle')
                        valid = (z_coord <= z_max(theta_tmp)); 
                    else 
                        valid = (z_coord < (z_annulus + z_extra));   
                    end 
                else 
                    valid = true; 
                end 
                
                % no upper springs
                if exist('r_max', 'var') && (r_tmp > r_max)
                    valid = false; 
                end 
                                
                if valid
                    points(:,theta_idx,r_idx,z_idx) = [x_coord, y_coord, z_coord]; 
                    indices_global(theta_idx,r_idx,z_idx) = params.global_idx; 
                    params.vertices(:,params.global_idx + 1) = points(:,theta_idx,r_idx,z_idx);

                    % every valid vertex is a target point here 
                    if exist('eta', 'var') && eta > 0 
                        params = target_string(params, params.global_idx, k_target, eta);     
                    else
                        params = target_string(params, params.global_idx, k_target);     
                    end 

                    params.global_idx = params.global_idx + 1; 
                else 
                    indices_global(theta_idx,r_idx,z_idx) = nan; 
                end 
            end        
        end 
    end 


% springs     
    % below the first possible point 
    idx = -1; 
    
    for z_idx=1:N_z
        for r_idx=1:N_r
            for theta_idx=1:N_theta  
                if ~isnan(indices_global(theta_idx,r_idx,z_idx))

                    output = true; 

                    last_idx = idx; 
                    idx = indices_global(theta_idx,r_idx,z_idx); 
                
                    if last_idx >= idx
                        error('should always be placing points in order, something wrong'); 
                    end

                    % check up directions for springs 
                    if theta_idx < N_theta
                        rest_len = norm(points(:,theta_idx,r_idx,z_idx) - points(:,theta_idx+1,r_idx,z_idx)); 
                        k_abs = k_rel / rest_len;
                        nbr_idx = indices_global(theta_idx+1,r_idx,z_idx); 
                        if ~isnan(nbr_idx)
                            min_idx = min(idx,nbr_idx); 
                            max_idx = max(idx,nbr_idx); 
                            params = spring_string(params, min_idx, max_idx, k_abs, rest_len, function_idx, output);
                        end 
                    end

                    % periodic wrap 
                    if theta_idx == N_theta
                        rest_len = norm(points(:,theta_idx,r_idx,z_idx) - points(:,1,r_idx,z_idx)); 
                        k_abs = k_rel / rest_len;
                        nbr_idx = indices_global(1,r_idx,z_idx); 
                        if ~isnan(nbr_idx)
                            min_idx = min(idx,nbr_idx); 
                            max_idx = max(idx,nbr_idx); 
                            params = spring_string(params, min_idx, max_idx, k_abs, rest_len, function_idx, output); 
                        end
                    end

                    if r_idx < N_r
                        rest_len = norm(points(:,theta_idx,r_idx,z_idx) - points(:,theta_idx,r_idx+1,z_idx)); 
                        k_abs = k_rel / rest_len;
                        nbr_idx = indices_global(theta_idx,r_idx+1,z_idx);
                        if ~isnan(nbr_idx)
                            min_idx = min(idx,nbr_idx); 
                            max_idx = max(idx,nbr_idx); 
                            params = spring_string(params, min_idx, max_idx, k_abs, rest_len, function_idx, output); 
                        end
                    end

                    if z_idx < N_z
                        rest_len = norm(points(:,theta_idx,r_idx,z_idx) - points(:,theta_idx,r_idx,z_idx+1)); 
                        k_abs = k_rel / rest_len;
                        nbr_idx = indices_global(theta_idx,r_idx,z_idx+1); 
                        if ~isnan(nbr_idx)
                            min_idx = min(idx,nbr_idx); 
                            max_idx = max(idx,nbr_idx); 
                            params = spring_string(params, min_idx, max_idx, k_abs, rest_len, function_idx, output);
                        end 
                    end
                
                end
            end            
        end 
    end 

end 




function params = write_inst_for_targets_as_bcs(params, leaflet)
    % if targets_for_bcs, then want to set the maximum vertices (which are the points connected to the exact ring location)
    % as the instrument points 

    j_max = leaflet.j_max; 
    k_max = leaflet.k_max; 
    
    % write the instrument file header here 
    fprintf(params.inst, '1   # num meters in file\n'); 
    fprintf(params.inst, 'meter_0   # name\n'); 
    fprintf(params.inst, '%d  # number of meter points\n', j_max); 
    

    k = k_max - 1; 
    instrument_idx = 0; 
    
    for j=1:j_max
    
        % pull the global index here 
        idx = leaflet.indices_global(j,k); 
        
        if ~leaflet.is_internal(j,k)
            error('placing non internal vertex in inst file')
        end 
        if ~leaflet.is_bc(j,k+1)
            error('placing an instrument point, but nbr up in k is not a bc')
        end 
        
        fprintf(params.inst, '%d \t0 \t %d\n', idx, instrument_idx); 
        instrument_idx = instrument_idx + 1; 
                    
    end 

end 

function params = place_rays(params, leaflet, ds, L, k_rel, k_target, ref_frac, eta, max_to_place)
    % 
    % Places rays of fibers emenating from the leaflet 
    % Angle of rays makes them (roughly) geodesics 
    % on the valve surface and plane surface  
    % 
    % Spacing is determined by reflection 
    %
    % Mesh is placed with a rectangular topology
    % Springs everywhere internal 
    % All points are made targets
    % 
    % Input
    %     r                Internal radius 
    %     h                Height of valve ring 
    %     L                Box width -- point must have max norm less than L to be included 
    %     N                Number of points placed on the initial circle (valve ring)
    %     spring           Files for writing, must be open already 
    %     vertex 
    %     target
    %     inst 
    %     params.global_idx       Running totals  
    %     total_vertices 
    %     total_springs 
    %     total_targets 
    % 

    X            = leaflet.X; 
    j_max        = leaflet.j_max; 
    k_max        = leaflet.k_max;
    is_bc        = leaflet.is_bc; 
    is_internal  = leaflet.is_internal; 
    
    function_idx = 0; 
    
    if isfield(leaflet, 'periodic_j')
        periodic_j = leaflet.periodic_j; 
    else
        periodic_j = zeros(k_max,1); 
    end 
    
    if ~isfield(leaflet, 'indices_global')
        error('Must place leaflets before placing rays'); 
    end 
    
    if ~exist('max_to_place', 'var')
        max_to_place = inf; 
    end 
    
    % output flag information 
    copy = params.copy; 
    if params.output.mesh(copy)
        output        = true; 
        
        % use leaflet stride here so fibers that continue as rays are plotted as such
        output_stride = params.output.stride_leaflet; 
    else 
        output        = false; 
    end 
    
    if strcmp(params.type, 'aortic')
        k_range = 1; 
    else 
        k_range = 1:k_max; 
    end 
    
    for j = 1:j_max
        for k = k_range
            if is_bc(j,k)
            
                % take mod one be 
                if output && ((mod(j,output_stride) == 1) || (output_stride == 1))
                    output_tmp = true;
                else 
                    output_tmp = false;
                end 
                
                pt_ring = X(:,j,k); 

                % only get a fiber if the previous point is included in the leaflet  
                neighbors = []; 
                
                % possible to have both directions of j nbr
                for j_nbr = [j-1,j+1]
                    k_nbr = k; 
                    if (j_nbr > 0) &&  (k_nbr > 0) && (j_nbr <= j_max) && (k_nbr <= k_max) && is_internal(j_nbr,k_nbr)
                        neighbors = [X(:,j_nbr,k_nbr), neighbors] ; 
                    end
                end 
                
                
                j_nbr = j; 
                if strcmp(params.type, 'aortic')
                    % k_nbr always up 
                    k_nbr = k+1; 
                else 
                    % k_nbr always down 
                    k_nbr = k-1; 
                end 

                if (j_nbr > 0) &&  (k_nbr > 0) && (j_nbr <= j_max) && (k_nbr <= k_max) && is_internal(j_nbr,k_nbr)
                    neighbors = [X(:,j_nbr,k_nbr), neighbors] ; 
                end        

                for x = neighbors 

                    % find the initial reflected point 
                    
                    % adjacent ring points determine local normal and tangent 
                    j_plus__1 = get_j_nbr(j+1, k, periodic_j, j_max); 
                    j_minus_1 = get_j_nbr(j-1, k, periodic_j, j_max);
                    
                    ring_nbr_plus  = X(:, j_plus__1, k); 
                    ring_nbr_minus = X(:, j_minus_1, k); 
                    
                    val = get_geodesic_continued_point(x, pt_ring, ring_nbr_minus, ring_nbr_plus, params.ring_center); 

                    % each point moves by this much from the initial point 
                    increment = val - x; 

                    % just zero this component, they are both near 3 but maybe not exactly 
                    increment(3) = 0.0; 
                    
                    % set to proper mesh width regardless of other spacing 
                    increment = ds * increment / norm(increment); 


                    % point = val; 
                    point_prev = pt_ring; 
                    point = pt_ring + increment; 
                    nbr_idx = leaflet.indices_global(j,k); 

                    n_placed = 0; 
                    
                    % just keep adding until points leave the domain 
                    % while norm(point(1:2), inf) < L   
                    while (params.x_min <= point(1))      && ...   
                          (point(1)     <= params.x_max)  && ...   
                          (params.y_min <= point(2))      && ...   
                          (point(2)     <= params.y_max)  && ...
                          (n_placed     <  max_to_place)
                      
                        % grab the index 
                        idx = params.global_idx;

                        % place point 
                        params.vertices(:,params.global_idx + 1) = point; 
                        
                        % it's a target too 
                        if exist('eta', 'var')
                            params = target_string(params, idx, k_target, eta);     
                        else
                            params = target_string(params, idx, k_target);     
                        end 

                        rest_len = ref_frac * norm(point - point_prev); 
                        k_abs = k_rel / rest_len;

                        params = spring_string(params, nbr_idx, idx, k_abs, rest_len, function_idx, output_tmp);

                        point_prev = point; 
                        point      = point + increment; 
                        params.global_idx = params.global_idx + 1; 
                        nbr_idx    = idx;
                        
                        n_placed = n_placed + 1; 
                         
                    end 

                end 
            end 
        end 
    end 
end 


function [val] = get_geodesic_continued_point(x, pt_ring, ring_nbr_minus, ring_nbr_plus, ring_center)
    %
    % Takes a point inside the the valve ring
    % Returns a point which allows a geodesic continuation 
    % Taking a segment from the ring point to the wall creates a geodesic 
    % 
    % Input: 
    %     x           Coordinates 
    %     pt_ring     Point on the valve ring to reflect relative to 
    %     r           Radius of valve ring 
    %     h           height of valve ring 
    %

    tol = 1e5 * eps; 
    
    if abs(pt_ring(3) - ring_center(3)) > tol
        error('Initial ring point is not near z ring plane'); 
    end 
    
    if abs(ring_nbr_plus(3) - ring_center(3)) > tol
        error('Initial ring point is not near z ring plane'); 
    end 
    
    if abs(ring_nbr_minus(3) - ring_center(3)) > tol
        error('Initial ring point is not near z ring plane'); 
    end 
    
    % local tangent implies local normal 
    tangent = ring_nbr_plus(1:2) - ring_nbr_minus(1:2); 
    normal  = [tangent(2); -tangent(1)]; 
    normal  = normal / norm(normal); 
    
    % translate ring point to origin (implicitly)
    % and x somewhere near the origin 
    val = x - pt_ring; 

    % rotate system around z such that normal now points in x axis direction 
    theta = atan2(normal(2), normal(1)); 
    
    val = rotation_matrix_z(-theta) * val; 
    
    % rotate val into the z = 0 plane, inside the transformed ring 
    % this is not inverted as it gets the image of the geodesic point before reflection 
    phi = atan2(val(3), val(1)); 
    val = rotation_matrix_y( -(phi - pi)) * val;
    
    if abs(val(3)) > tol
        error('should have landed in z=0 plane here...');
    end 

    % reflect, this would be the geodesic point if the system was flat 
    val = -val; 
    
    % undo rigid rotations and translations 
    val = rotation_matrix_z(theta) * val; 
    val = val + pt_ring; 
    
end 






function params = place_cartesian_net(params, leaflet, r_extra, L, ds, k_rel, k_target, ref_frac, eta)
    % 
    % Places a cartesian coordinate mesh in a box 
    % This is to avoid issues with the polar mesh at the edge 
    % Mesh is placed on [-L + ds/2, L-ds/2]
    % 
    % Springs everywhere internal in 2d arrangement 
    % All points are made targets
    % 
    % Input
    %     r                Internal radius, it is advisible to make this larger than the valve ring 
    %     h                Height of valve ring 
    %     L                Box width -- point must have max norm less than L to be included 
    %     ds               Mesh is placed starting ds from boundary 
    %     params.global_idx       Running totals  
    %     total_vertices 
    %     total_springs 
    %     total_targets 
    % 

    
    N = ceil(2*L / ds);  
    points = zeros(3,N,N);

    % compute the valve ring, adding r_extra to each point 
    X     = leaflet.X; 
    j_max = leaflet.j_max; 
    k_max = leaflet.k_max; 
   
    function_idx = 0; 
    
    ring_expanded = zeros(2,j_max); 
    
    for j=1:j_max
        % valve ring points from leaflet 
%         ring_pt    = X(1:2,j,k_max); 
%         increment  = r_extra * ring_pt / norm(ring_pt); 
% 
%         ring_expanded(:,j) = ring_pt + increment; 

        % valve ring points from leaflet 
        if strcmp(params.type, 'aortic')
            % ring_pt = X(1:2,j,1);            
            r = leaflet.skeleton.r; 
            ring_pt = [r*cos(2*pi*(j-1)/j_max); r*sin(2*pi*(j-1)/j_max)];
        else 
            ring_pt = X(1:2,j,k_max);
        end 
        
        increment = ring_pt - params.ring_center(1:2); 
        increment = r_extra * increment / norm(increment); 

        % expand in the direction of a vector from venter of ring to
        % current point 
        ring_expanded(:,j) = ring_pt + increment; 

    end 
    
    % output flag information 
    copy = params.copy; 
    if params.output.cartesian_mesh(copy)
        output        = true; 
        output_stride = params.output.stride_mesh; 
    else 
        output        = false; 
    end 
    
    
    % just keep a list of valid indices, mark NAN if out of physical bounds 
    indices_global = zeros(N,N);     
    
    
    % This loop should be one indexed, 
    % because we want to start not at the edge but ds in 
    for k=1:N
        for j=1:N
                        
            % coords_horiz = [ (j-1)*ds - L + ds/2; (k-1)*ds - L + ds/2]; 
            coords_horiz = [ (j-1)*ds + params.x_min + ds/2; (k-1)*ds + params.y_min + ds/2]; 
            
            % if one norm is less than L, then the point is within the domain  
            % if (norm(coords_horiz, inf) < L) && point_out_of_polygon(ring_expanded, coords_horiz)
            if (params.x_min    <= coords_horiz(1)) && ...   
               (coords_horiz(1) <= params.x_max   ) && ...   
               (params.y_min    <= coords_horiz(2)) && ...   
               (coords_horiz(2) <= params.y_max   ) && ... 
               point_out_of_polygon(ring_expanded, coords_horiz)
                
                points(:,j,k) = [coords_horiz; params.ring_center(3)]; 
                indices_global(j,k) = params.global_idx; 
                params.vertices(:,params.global_idx + 1) = points(:,j,k);

                % every valid vertex is a target point here 
                if exist('eta', 'var')
                    params = target_string(params, params.global_idx, k_target, eta);     
                else
                    params = target_string(params, params.global_idx, k_target);     
                end 

                params.global_idx = params.global_idx + 1; 
                
                
            else 
                points(:,j,k) = NaN * ones(3,1);
                indices_global(j,k) = NaN; 
            end     
            
        end 
        
    end 

    
    % below the first possible point 
    idx = -1; 
    
    for k=1:N
        for j=1:N
            
            % just ignore the NaNs 
            if ~isnan(indices_global(j,k))
 
                if output && ((mod(j,output_stride) == 1) || (output_stride == 1))
                    output_tmp_k = true; 
                else 
                    output_tmp_k = false; 
                end 
                
                if output && ((mod(k,output_stride) == 1) || (output_stride == 1))
                    output_tmp_j = true; 
                else 
                    output_tmp_j = false; 
                end    
                
                last_idx = idx; 
                idx = indices_global(j,k); 
                
                if last_idx >= idx
                    error('should always be placing points in order, something wrong'); 
                end 

                % check up directions for springs 
                if (j+1) <= N
                    if ~isnan(indices_global(j+1,k))
                        rest_len = ref_frac * norm(points(:,j,k) - points(:,j+1,k)); 
                        k_abs = k_rel / rest_len;
                        nbr_idx = indices_global(j+1,k); 
                        params = spring_string(params, idx, nbr_idx, k_abs, rest_len, function_idx, output_tmp_j); 
                    end 
                end 
                                
                % check up directions for springs 
                if (k+1) <= N
                    if ~isnan(indices_global(j,k+1))
                        rest_len = ref_frac * norm(points(:,j,k) - points(:,j,k+1)); 
                        k_abs = k_rel / rest_len; 
                        nbr_idx = indices_global(j,k+1); 
                        params = spring_string(params, idx, nbr_idx, k_abs, rest_len, function_idx, output_tmp_k); 
                    end 
                end 
                
            end 
           
        end 
    end 

end 


function outside = point_out_of_polygon(vert, test)
%
% Recklessly translated from 
% 
% https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
%
% Original code: 
% 
% int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy)
% {
%   int i, j, c = 0;
%   for (i = 0, j = nvert-1; i < nvert; j = i++) {
%     if ( ((verty[i]>testy) != (verty[j]>testy)) &&
% 	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
%        c = !c;
%   }
%   return c;
% }

    inside = false; 
    
    if size(vert,1) ~= 2
        error('Must pass two d vector of vertices to test');         
    end 
    
    nvert = size(vert,2); 
    
    i = 1; 
    j = nvert; 
    
    while i <= nvert 
        
        check_1 = ((vert(2,i) > test(2)) ~= (vert(2,j) > test(2))); 
        tmp     = (vert(1,j) - vert(1,i)) * (test(2)-vert(2,i)) / (vert(2,j) - vert(2,i)) + vert(1,i); 
        check_2 = (test(1) < tmp); 
        
        if check_1 && check_2
            inside = ~inside; 
        end 
        
        j = i; 
        i = i+1; 
    end 
    
    outside = ~inside; 

end 


function [params, total_lagrangian_placed] = place_lagrangian_tracers(params, n_lagrangian_tracers, double_z)
    % Places a uniform cartesian mesh of lagrangian particle tracers 
    % Simple lopp implementation 
    %
    %     params.global_idx               Running totals  
    %     total_vertices 
    %     vertex                   vertex file for writing  
    %     n_lagrangian_tracers
    %     L                        mesh placed in L/2

    % junk hack exit 
    if n_lagrangian_tracers == 0
        total_lagrangian_placed = 0; 
        return; 
    end 
    
    % includes ring_center
    % do not need to update this for nonzero ring_center
    x_min = params.x_min; 
    x_max = params.x_max; 
    y_min = params.y_min; 
    y_max = params.y_max; 
    z_min = params.z_min; 
    z_max = params.z_max;
    
    n_z_dir = n_lagrangian_tracers;  
    if double_z
        n_z_dir = 2*n_z_dir; 
    end 
    
    dx = (x_max - x_min) / n_lagrangian_tracers; 
    dy = (y_max - y_min) / n_lagrangian_tracers; 
    dz = (z_max - z_min) / n_z_dir; 
    
    total_lagrangian_placed = 0; 
    
    for i = 1:n_lagrangian_tracers
        for j = 1:n_lagrangian_tracers
            for k = 1:n_z_dir 
                
                x = (i - .5)*dx + x_min; 
                y = (j - .5)*dy + y_min; 
                z = (k - .5)*dz + z_min; 
                
                params.vertices(:,params.global_idx + 1) = [x y z]; 
    
                total_lagrangian_placed = total_lagrangian_placed + 1; 
                params.global_idx = params.global_idx + 1; 
            end 
        end 
    end 

end 





