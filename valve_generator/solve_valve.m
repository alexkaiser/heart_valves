function [valve valve_with_reference pass_all] = solve_valve(valve, interactive, from_history, build_reference)
% 
% Refines valve data structure to equilibrium 
% Applies auto-continuation to pressure and updates both leaflets 
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

pass_all = true; 

if length(valve.leaflets) ~= 1
    error('running with single leaflet assumption for now'); 
end 

if ~exist('from_history', 'var')
    from_history = false; 
end 

if ~exist('interactive', 'var')
    interactive = false; 
end 

if ~exist('build_reference', 'var')
    build_reference = true; 
end 

tol_global            = valve.tol_global; 
max_it                = valve.max_it; 
max_continuations     = valve.max_continuations; 
max_consecutive_fails = valve.max_consecutive_fails; 
max_total_fails       = valve.max_total_fails; 


if from_history
    if ~isfield(valve, 'tension_coeff_history')
        error('Must provide history to solve from history'); 
    end 
    
    [leaflet_current valve_current] = set_tension_coeffs(valve.leaflets(1), valve, valve.tension_coeff_history(1)); 
    valve = valve_current; 
    valve.leaflets(1) = leaflet_current; 
end 


for i=1:length(valve.leaflets)
    
    leaflet = valve.leaflets(i); 
    
    p_easy = 0; 
    p_goal    = leaflet.p_0; 
    goal_first = false; 
    
    [valve.leaflets(i) pass err] = solve_valve_pressure_auto_continuation(leaflet, tol_global, max_it, max_continuations, p_easy, p_goal, max_consecutive_fails, max_total_fails, goal_first); 

    if pass
        fprintf('Global solve passed, err = %e\n\n', err); 
    else 
        fprintf('Global solve failed, err = %e\n\n', err); 
    end 
    
%     fig = figure; 
%     surf_plot(valve.leaflets(i), fig); 
%     pause(0.01);
    
    pass_all = pass_all && pass; 
    
end 


if from_history
    
    fig = figure; 
    valve_plot(valve, fig); 
    title('History mode valve')
    
    for history_step = 2:length(valve.tension_coeff_history)
        [leaflet_current valve_current] = set_tension_coeffs(valve.leaflets(1), valve, valve.tension_coeff_history(history_step)); 
        
        valve = valve_current; 
        valve.leaflets(1) = leaflet_current;

        try
            [valve.leaflets(1) pass err] = newton_solve_valve(valve.leaflets(1), tol_global, max_it, max_consecutive_fails, max_total_fails);  
        catch 
            fprintf('This iteration in history failed. Maybe by some miracle the next one will pass. Moving on...\n'); 
            continue; 
        end 
        
        if ~pass 
            fprintf('This iteration in history failed. Maybe by some miracle the next one will pass. Moving on...\n'); 
            continue; 
        end 

        fig = figure; 
        [az el] = view; 
        clf(fig); 
        valve_plot(valve, fig);
        title(sprintf('History mode valve, it %d\n', history_step)); 
        view(az,el);
        
    end 
end 



if interactive && pass_all
    fprintf('Solve passed, interactive mode enabled.\n'); 
    
    num_passed = 1; 
    if from_history
        num_passed = length(valve.tension_coeff_history) + 1; 
    end 
    valve.tension_coeff_history(num_passed) = valve.leaflets(1).tension_coeffs; 
    
    fig = figure; 
    valve_plot(valve, fig); 
    title('Valve in interactive mode'); 
    
    fig_dissection_plot = figure; 
    if isfield(valve, 'name') && strcmp(valve.name, 'aortic') 
        fig_dissection_plot = dissection_plot_rest_height_aortic(valve, fig_dissection_plot); 
    else
        fig_dissection_plot = dissection_plot_rest_height(valve, fig_dissection_plot); 
    end 
    
%     if isfield(valve, 'name') && strcmp(valve.name, 'aortic')        
%         fiber_output    = true; 
%         fiber_stride    = 4; 
%         stride_offset_j = 0; 
%         circ  = false; 
%         rad   = false; 
%         ratio = true; 
%         height_plot = true; 
%         fig_ratio = figure; 
%         set(0, 'CurrentFigure', fig_ratio)
%         total_tension_surf_plot_aortic(valve.leaflets(1), fiber_output, fiber_stride, stride_offset_j, circ, rad, ratio, height_plot, fig_ratio)
%         title('ratio circ/radial tension')
%     end 
    
    while true 
    
        fprintf('Current tension_coeffs struct, which includes all valid variables:\n')
        fprintf('tension_coeffs = \n\n')
        disp(valve.leaflets(1).tension_coeffs)

        if isfield(valve, 'name') && strcmp(valve.name, 'aortic')
            % nothing 
        else
            % default mitral 
            fprintf('Current array variables:\n\n')
            fprintf('k_0_1 = \n\n')
            disp(valve.leaflets(1).tension_coeffs.k_0_1)
            fprintf('k_root = \n\n')
            disp(valve.leaflets(1).tension_coeffs.k_root)
            fprintf('c_dec_chordae_leaf = \n\n')
            disp(valve.leaflets(1).tension_coeffs.c_dec_chordae_leaf)
            fprintf('c_dec_chordae_root = \n\n')
            disp(valve.leaflets(1).tension_coeffs.c_dec_chordae_root)
            fprintf('\n\n'); 
        end 
        
        try 
            var_name = input('Enter the name of variable to change as a string (no quotes or whitespace, must match exactly):\n', 's'); 

            if ~ischar(var_name)
                fprintf('Must input a valid string for variable name.\n'); 
                continue; 
            end 

            if isempty(var_name)
                var_name = input('No variable name entered. Enter empty string again to leave interactive mode.\n', 's'); 
                if isempty(var_name)
                    break; 
                end 
            end 

            if isfield(valve.leaflets(1).tension_coeffs, var_name)

                tension_coeffs_current = valve.leaflets(1).tension_coeffs; 
                
                if length(tension_coeffs_current.(var_name)) > 1
                    fprintf('Found valid variable %s. You have selected an array variable. Contents of array are:', var_name);
                    tension_coeffs_current.(var_name)
                    idx = input('Enter the index you would like to change, ranges okay:\n');  
                    value = input('Input new value:\n'); 
                    tension_coeffs_current.(var_name)(idx) = value;
                else                    
                    value_old = tension_coeffs_current.(var_name);
                    fprintf('Found valid variable %s with old value %f.\n', var_name, value_old); 
                    value = input('Input new value:\n'); 
                    tension_coeffs_current.(var_name) = value;
                end 
                
                [leaflet_current valve_current] = set_tension_coeffs(valve.leaflets(1), valve, tension_coeffs_current); 

                try
                    [leaflet_current pass err] = newton_solve_valve(leaflet_current, tol_global, max_it, max_consecutive_fails, max_total_fails);  

                    if pass 
                        
                        % copy data from new version 
                        valve             = valve_current; 
                        valve.leaflets(1) = leaflet_current; 
                        
                        % save history 
                        num_passed = num_passed + 1; 
                        valve.tension_coeff_history(num_passed) = valve.leaflets(1).tension_coeffs; 
                        times = clock; 
                        save_name = sprintf('%s_tension_history_%d_%d_%d_%d.%d.%d.mat', valve.base_name, times(1), times(2), times(3), times(4), times(5), round(times(6))); 
                        history_tmp = valve.tension_coeff_history; 
                        save(save_name, 'history_tmp'); 

                        % update plot 
                        % if we do not have a figure add one 
                        if ~ishandle(fig)
                            fig = figure; 
                        end
                        set(0, 'CurrentFigure', fig)
                        [az el] = view;
                        clf(fig); 
                        valve_plot(valve, fig);
                        view(az,el);
                        
                        if ~ishandle(fig_dissection_plot)
                            fig_dissection_plot = figure; 
                        end
                        set(0, 'CurrentFigure', fig_dissection_plot)
                        if isfield(valve, 'name') && strcmp(valve.name, 'aortic') 
                            fig_dissection_plot = dissection_plot_rest_height_aortic(valve, fig_dissection_plot); 
                        else
                            fig_dissection_plot = dissection_plot_rest_height(valve, fig_dissection_plot); 
                        end 
                        % update aortic tension plots 
                        if isfield(valve, 'name') && strcmp(valve.name, 'aortic') 
                            
                            if ~exist('fig_ratio', 'var')
                                fig_ratio = figure; 
                            end
                            fiber_output    = true; 
                            fiber_stride    = 8; 
                            stride_offset_j = 4; 
                            circ = false; 
                            rad = false; 
                            ratio = true; 
                            height_plot = false; 
                            set(0, 'CurrentFigure', fig_ratio)
                            [az el] = view;
                            clf(fig_ratio); 
                            total_tension_surf_plot_aortic(valve.leaflets(1), fiber_output, fiber_stride, stride_offset_j, circ, rad, ratio, height_plot, fig_ratio);
                            title('ratio circ/radial tension'); 
                            view(az,el);
                        end 
                        
                    else 
                        fprintf('New parameters failed. Keeping old tension structure. No pass flag.\n'); 
                    end
                catch e %e is an MException struct                    
                    fprintf(1,'The identifier was: %s\n',e.identifier);
                    fprintf(1,'There was an error! The message was: %s\n',e.message);
                    fprintf('New parameters failed. Keeping old tension structure. Catch block, error called in Newton solve.\n'); 
                end 

            else 
                fprintf('Variable not found.\n'); 
            end 
            
         catch 
             fprintf('Error in interactive loop of some kind, restart loop.\n'); 
         end 
    end 
end 

    
    
if build_reference

    % constitutive law version 
    valve_with_reference = valve; 

    % kill off the old leaflet structure, new one has different fields, 
    % which makes matlab complain about assigning it to a structure array 
    valve_with_reference = rmfield(valve_with_reference, 'leaflets'); 

    for i=1:length(valve.leaflets)

        if isfield(valve, 'name') && strcmp(valve.name, 'aortic')
            valve_with_reference.leaflets(i) = set_rest_lengths_and_constants_aortic(valve.leaflets(i), valve); 
            
            plots = false; 
            fig = figure; 
            [sigma_circ, sigma_rad, sigma_circ_mean, sigma_rad_mean, ~, ...
            stress_circ, stress_rad, stress_circ_mean, stress_rad_mean, k_u_stress_mean, k_v_stress_mean] ...
                = estimate_tangent_modulus_aortic_with_reference(valve_with_reference.leaflets(i), valve.normal_thickness, fig);
            valve_with_reference.leaflets(i).sigma_circ = sigma_circ; 
            valve_with_reference.leaflets(i).sigma_rad = sigma_rad; 
            valve_with_reference.leaflets(i).sigma_circ_mean = sigma_circ_mean;  
            valve_with_reference.leaflets(i).sigma_rad_mean = sigma_rad_mean;                
            valve_with_reference.leaflets(i).stress_circ = stress_circ; 
            valve_with_reference.leaflets(i).stress_rad = stress_rad; 
            valve_with_reference.leaflets(i).stress_circ_mean = stress_circ_mean;  
            valve_with_reference.leaflets(i).stress_rad_mean = stress_rad_mean;                         
            valve_with_reference.leaflets(i).k_u_stress_mean = k_u_stress_mean;  
            valve_with_reference.leaflets(i).k_v_stress_mean = k_v_stress_mean; 
            sigma_circ_mean  
            sigma_rad_mean
            fprintf('\n')
            
            if ~plots
                close(fig); 
            end 
            
            if isfield(valve_with_reference, 'dilate_graft') && valve_with_reference.dilate_graft
                if isfield(valve_with_reference, 'dilation_dist')
                    if valve_with_reference.dilation_dist ~= 0
                        valve_with_reference = aortic_dilate_annulus(valve_with_reference); 
                    end 
                else 
                    error('must provide dilation_dist if valve.dilate_graft is true'); 
                end 
            end 
            
            if isfield(valve_with_reference, 'comm_raise_normal_height') 
                if valve_with_reference.comm_raise_normal_height ~= valve_with_reference.skeleton.normal_height 
                    valve_with_reference = aortic_raise_comm(valve_with_reference); 
                end
            end 
            
            
            if isfield(valve, 'dirichlet_free_edge_with_ref_only') && valve.dirichlet_free_edge_with_ref_only && ... 
               isfield(valve, 'dirichlet_free_edge_comm_ref_only') && valve.dirichlet_free_edge_comm_ref_only
                error('incompatible options')
            end 
                
            
            if isfield(leaflet, 'fused_commissure') && leaflet.fused_commissure                     
                if ~isfield(leaflet, 'fused_comm_idx')
                    error('must supply fused_comm_idx if leaflet.fused_commissure is true')
                end 
                
                if isfield(valve, 'extra_stretch_radial_dirichlet_free_edge')
                    extra_stretch_radial = valve.extra_stretch_radial_dirichlet_free_edge; 
                else
                    extra_stretch_radial = valve.strain_rad + 1.0; 
                end 
                
                leaflets_copy = aortic_free_edge_fuse_commissure(valve_with_reference.leaflets(i), extra_stretch_radial, leaflet.fused_comm_idx); 
                valve_with_reference = rmfield(valve_with_reference, 'leaflets'); 
                valve_with_reference.leaflets(1) = leaflets_copy; 
            
            elseif isfield(leaflet, 'pinch_commissure') && leaflet.pinch_commissure                     
                if ~isfield(leaflet, 'N_to_pinch')
                    error('must supply N_to_pinch if leaflet.pinch_commissure is true')
                end 
                
                if isfield(valve, 'extra_stretch_radial_dirichlet_free_edge')
                    extra_stretch_radial = valve.extra_stretch_radial_dirichlet_free_edge; 
                else
                    extra_stretch_radial = valve.strain_rad + 1.0; 
                end
                
                valve_with_reference.leaflets(i) = aortic_free_edge_pinch_comms(valve_with_reference.leaflets(i), extra_stretch_radial, leaflet.N_to_pinch);
%                 leaflets_copy = aortic_free_edge_pinch_comms(valve_with_reference.leaflets(i), extra_stretch_radial, leaflet.N_to_pinch); 
%                 valve_with_reference = rmfield(valve_with_reference, 'leaflets'); 
%                 valve_with_reference.leaflets(1) = leaflets_copy;    
                
            elseif isfield(valve, 'dirichlet_free_edge_with_ref_only') && valve.dirichlet_free_edge_with_ref_only
                % multiplicative stretch for setting aortic valve initial condition 
                % set to prescribed stain + 1 (prescribed stretch) for maintaining the leaflet height as the loaded height 
                if isfield(valve, 'extra_stretch_radial_dirichlet_free_edge')
                    extra_stretch_radial = valve.extra_stretch_radial_dirichlet_free_edge; 
                else
                    extra_stretch_radial = valve.strain_rad + 1.0; 
                end 
                if isfield(valve, 'extra_stretch_circ_dirichlet_free_edge')
                    extra_stretch_circ = valve.extra_stretch_circ_dirichlet_free_edge; 
                else
                    extra_stretch_circ = 1.0; 
                end 
                valve_with_reference.leaflets(i) = aortic_free_edge_to_dirichlet_bc(valve_with_reference.leaflets(i), extra_stretch_radial, extra_stretch_circ); 

            end 
            
        else
            % mitral default 
            valve_with_reference.leaflets(i) = set_rest_lengths_and_constants(valve.leaflets(i), valve); 
        end 

        % leave here to compute annulus force at systolic pressure, valve.p_physical 
        % comment to log diastolic or rest pressure 
        if isfield(valve_with_reference, 'log_annulus_force') && valve_with_reference.log_annulus_force
            if isfield(valve_with_reference, 'force_log_name') 
                [annulus_positions, forces_annulus] = compute_annulus_force(valve_with_reference.leaflets(i), valve_with_reference.force_log_name); 
            else 
                [annulus_positions, forces_annulus] = compute_annulus_force(valve_with_reference.leaflets(i)); 
            end 
        end 


        leaflet = valve_with_reference.leaflets(i); 

        p_easy = leaflet.p_0/10; 

        if isfield(valve, 'p_final')
            p_goal = valve.p_final; 
        else 
            p_goal = 0; 
        end 

        max_continuations_relaxed = 6; 

        fprintf('\n\nRefernece configuration initial solve:\n')
        [valve_with_reference.leaflets(i) pass err any_passed] = solve_valve_pressure_auto_continuation(leaflet, tol_global, max_it, max_continuations_relaxed, p_easy, p_goal, max_consecutive_fails, max_total_fails); 


        if isfield(valve, 'dirichlet_free_edge_comm_ref_only') && valve.dirichlet_free_edge_comm_ref_only
            if ~isfield(valve, 'n_fixed_comm')
                error('must provide n_fixed_comm if dirichlet_free_edge_comm_ref_only is true'); 
            end
            
            if ~isfield(valve, 'p_final_fixed_comm') 
                error('must provide p_final_fixed_comm if dirichlet_free_edge_comm_ref_only is true'); 
            end 
            
            plots_temp = true; 
            if plots_temp 
                fig = figure; 
                fig = valve_plot(valve_with_reference, fig); 
                title('after initial neumann bc solve before comm bc set');
            end     
            
            p_easy = p_goal; 
            p_goal = valve.p_final_fixed_comm; 
                        
            leaflet = aortic_comm_to_dirichlet_bc(valve_with_reference.leaflets(i), valve.n_fixed_comm);

            fprintf("\n\nRefernece config solve with fixed free edge near comm:\n")
            [valve_with_reference.leaflets(i) pass err any_passed] = solve_valve_pressure_auto_continuation(leaflet, tol_global, max_it, max_continuations_relaxed, p_easy, p_goal, max_consecutive_fails, max_total_fails);             
        end 
        
        % pre_extrude leaflet before  
        % this should probably be in a function 
        if isfield(valve_with_reference.leaflets(i), 'pre_extrude') && valve_with_reference.leaflets(i).pre_extrude
            
            ds_extrude = valve_with_reference.normal_thickness / (valve_with_reference.num_copies-1); 
            
            for copy = 1:valve_with_reference.num_copies
            
            
                if isfield(valve_with_reference, 'center_extrusion') && valve_with_reference.center_extrusion

                    if valve_with_reference.num_copies ~= 3
                        error('center extrusion only implemted with 3 layers') 
                    end 

                    extrude_length = (copy - 2) * ds_extrude; 
                else 
                    extrude_length = (copy - 1) * ds_extrude; 
                end 

                leaflet_temp = valve_with_reference.leaflets(i); 
                
                leaflet_temp.X = normal_extrude_aortic(leaflet_temp, extrude_length); 
            

                if extrude_length ~= 0
                    
                    % relax 
                    [leaflet_temp pass err any_passed] = solve_valve_pressure_auto_continuation(leaflet_temp, tol_global, max_it, max_continuations_relaxed, p_easy, p_goal, max_consecutive_fails, max_total_fails);             
                    
                    if ~pass
                        error('Extruded leaflet failed to converge'); 
                    end 
                    
                end 
                 
                valve_with_reference.leaflets_pre_extruded(copy) = leaflet_temp; 
                
                fig = figure; 
                surf_plot(leaflet_temp, fig); 
                title(sprintf('leaflet copy %d', copy));
                
            end 
                
        end 
        
                
        % reset the free edges of the aortic to neumann bc for final output 
        if isfield(valve, 'name') && strcmp(valve.name, 'aortic')
            if isfield(valve, 'dirichlet_free_edge_with_ref_only') && valve.dirichlet_free_edge_with_ref_only
                valve_with_reference.leaflets(i) = aortic_free_edge_to_neumann_bc(valve_with_reference.leaflets(i)); 
                
                neumann_visual_check = false; 
                if neumann_visual_check
                    leaflets_neumann = solve_valve_pressure_auto_continuation(valve_with_reference.leaflets(i), tol_global, max_it, max_continuations_relaxed, p_easy, p_goal, max_consecutive_fails, max_total_fails); 
                    surf_plot(leaflets_neumann); 
                    title('neumann bc with reference, not used')
                end 
                
            end 
        end 

        % uncomment to compute annulus force at diastolic pressure, valve.p_final
%         if isfield(valve_with_reference, 'log_annulus_force') && valve_with_reference.log_annulus_force
%             if isfield(valve_with_reference, 'force_log_name') 
%                 [annulus_positions, forces_annulus] = compute_annulus_force(valve_with_reference.leaflets(i), valve_with_reference.force_log_name); 
%             else 
%                 [annulus_positions, forces_annulus] = compute_annulus_force(valve_with_reference.leaflets(i)); 
%             end 
%         end 
        
        if isfield(valve, 'rotate_identical_leaflets') && valve.rotate_identical_leaflets 
            valve_with_reference = add_rotated_leaflets_aortic(valve_with_reference);
        end 
        
        if pass
            fprintf('Global solve passed, err = %e\n\n', err); 
        else 
            if any_passed
                fprintf('Global solve passed but with pressure %e, err = %e\n\n', valve_with_reference.leaflets(i).p_0, err); 
            else 
                fprintf('Global solve failed, err = %e\n\n', err); 
            end 
        end 

    %     fig = figure; 
    %     surf_plot(valve.leaflets(i), fig); 
    %     pause(0.01);

        pass_all = pass_all && pass; 

    end 
    
else 
    valve_with_reference = []; 

end 

