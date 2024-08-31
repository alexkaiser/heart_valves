function [annulus_positions, forces_annulus] = compute_annulus_force(leaflet, file_name)

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

    X_current          = leaflet.X; 
    j_max              = leaflet.j_max; 
    k_max              = leaflet.k_max; 
    R_v                = leaflet.R_v;
    k_v                = leaflet.k_v;
        
    annulus_positions = squeeze(X_current(:,:,k_max)); 
    forces_annulus = zeros(size(annulus_positions)); 

    if ~exist('file_name', 'var')
        logging = false; 
    else 
        logging = true;
        f = fopen(file_name, 'w'); 
    end 
    
    
    % Internal leaflet part 
    for j=1:j_max
        k=k_max; % constant k max 

        X = X_current(:,j,k); 

        F_tmp = zeros(3,1);

        % v type fibers 
        for k_nbr_tmp = [k-1,k+1]

            j_nbr_tmp = j; 

            [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

            if valid
                X_nbr = X_current(:,j_nbr,k_nbr); 

                tension = tension_with_reference(X, X_nbr, R_v(j_spr,k_spr), k_v(j_spr,k_spr), leaflet); 
                F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 
            end 

            if (k_nbr_tmp == (k-1)) && (~valid)
                error('point should be valid but not labeled as such'); 
            end 
            
            if (k_nbr_tmp == (k+1)) && valid
                error('finding valid spring away from annulus'); 
            end 
            
        end 

        forces_annulus(:,j) = F_tmp;

        if logging
            fprintf(f, '%.10f %.10f %.10f %.10f %.10f %.10f\n', ...
                annulus_positions(1,j), annulus_positions(2,j), annulus_positions(3,j), ...
                   forces_annulus(1,j),    forces_annulus(2,j),    forces_annulus(3,j)); 
        end 
        
    end

%     annulus_positions
%     forces_annulus

    if logging
        fclose(f); 
    end 

end 








