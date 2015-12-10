function [X_nbr R_nbr idx_chordae left] = get_neighbor(params, filter_params, j_nbr, k_nbr)
% 
% Returns the neighbor to a point on 
% By convention, if j==0 or k==0 then points off the leaflet are used 
% 
% If params.chordae is empty or does not exists, 
% then the papillary coordinate is returned 
% 
% If params.chordae is nonempty, 
% then the appropriate location on the chordae tree is returned 
% 

idx_chordae = 0; 
left = false; 


if (j_nbr > 0) && (k_nbr > 0) 
    
    X_nbr = params.X(:,j_nbr,k_nbr); 
    R_nbr = params.R(:,j_nbr,k_nbr); 
    
elseif (j_nbr == 0) && (k_nbr == 0) 
    
    error('there is no zero zero direction neighbor'); 
    
elseif (~isfield(params, 'chordae')) || isempty(params.chordae)
    
    % attach to the papillary muscle 
    left_papillary  = [0; -filter_params.a; 0]; 
    right_papillary = [0;  filter_params.a; 0]; 
    
    if j_nbr == 0
        X_nbr = left_papillary;
        R_nbr = left_papillary;
    elseif k_nbr == 0
        X_nbr = right_papillary;
        R_nbr = right_papillary;
    else 
        error('should be impossible to get here'); 
    end 
    
else
    
    % chordae part 
    if (j_nbr == 0) && (k_nbr > 0)
        
        % use the left chordae 
        [m max_internal] = size(params.chordae.C_left); 
        idx_chordae = max_internal + k_nbr;
        X_nbr = params.chordae.C_left(:,idx_chordae); 
        R_nbr = params.chordae.Ref_l(:,idx_chordae); 

    elseif (k_nbr == 0) && (j_nbr > 0)
        
        % use the right chordae 
        [m max_internal] = size(params.chordae.C_right); 
        idx_chordae = max_internal + j_nbr;
        X_nbr = params.chordae.C_right(:,idx_chordae); 
        R_nbr = params.chordae.Ref_r(:,idx_chordae); 
        
    else 
        error('should be impossible to get here'); 
    end 
end 





