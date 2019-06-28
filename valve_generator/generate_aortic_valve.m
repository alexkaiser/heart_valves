% Script to build valve 

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

% Size parameter
% Number of points on free edge of each leaflet 
% 
N_each_range = 2^5; % 2.^(5:11); 

for N_each = N_each_range

    clearvars -except N_each
    
    N_each
    N = 3*N_each; 
    
    % Show some output 
    plots = false; 

    % Initialize structures 
    % Many parameters are in this script 

    valve = initialize_valve_data_structures_aortic_generic(N); 
    
    iteration_movie = false; 
    if iteration_movie
        valve.leaflets(1).iteration_movie = true;
        % valve.leaflets(1).movie_name = 'const_tension_newton'; 

        valve.leaflets(1).movie_name = 'iteration_newtons';

    end 


    interactive = false; 

    from_history = false; 
    if from_history 
        history_name = ''; 
        load(history_name); 
        valve.tension_coeff_history = history_tmp; 
    end 

    % if radial && bead_slip && attached
    %     valve = newton_solve_valve_attached(valve, valve.tol_global, valve.max_it); 
    % else 
    % %     valve = solve_valve(valve, p_range, repulsive_coeff_range); 

    [valve valve_with_reference pass_all] = solve_valve(valve, interactive, from_history); 
    % end 

    fig = figure; 
    fig = valve_plot(valve, fig); 

    title('Pressurized configuration fibers'); 
    saveas(fig, strcat(valve.base_name, '_pressurized'), 'fig'); 

    if ~isempty(valve_with_reference)
        fig = figure; 
        fig = valve_plot(valve_with_reference, fig); 
        title('Relaxed configuration radial fibers, reference config based constitutive law'); 
        saveas(fig, strcat(valve.base_name, '_relaxed'), 'fig'); 
    end 

    if pass_all 
        fprintf('Final solve passed.\n'); 
    else 
        fprintf('Final solve failed.\n'); 
    end 


    % Save current data 
    save(strcat(valve.base_name, '_final_data')); 

    % Write to simulation files 
    if ~isempty(valve_with_reference)
        output_to_ibamr_format(valve_with_reference); 
    end 

end 
