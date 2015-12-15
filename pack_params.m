function params = pack_params(X,alpha,beta,N,p_0,R,ref_frac,chordae)
%
% Builds a single struct with all of the current parameters 
% 
% 
 
params.X         = X; 
params.alpha     = alpha; 
params.beta      = beta; 
params.p_0       = p_0; 
params.N         = N; 
params.R         = R; 
params.ref_frac  = ref_frac; 

if nargin > 7 
    params.chordae   = chordae; 
end 
    
    
    