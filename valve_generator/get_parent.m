function [p p_ref] = get_parent(C, i, papillary, ref)
% returns parent coordinates in binary tree of chordae 
% if i==1 then returns papillary 

    if i==1
        p = papillary; 
    else 
        p = C(:,floor(i/2));  
    end 
    
    if nargin > 3
        if i==1
            p_ref = papillary; 
        else 
            p_ref = ref(:,floor(i/2));  
        end        
    else 
        p_ref = []; 
    end 
    
    
end 
