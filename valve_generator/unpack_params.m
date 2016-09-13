function [X,alpha,beta,N,p_0,R,ref_frac,chordae] = unpack_params(params)
%
% Unpacks the struct into individual variables 
% 
% 
 
X         = params.X; 
alpha     = params.alpha; 
beta      = params.beta; 
p_0       = params.p_0; 
N         = params.N; 
R         = params.R; 
ref_frac  = params.ref_frac; 

if isfield(params, 'chordae') == 1 
    chordae = params.chordae; 
else 
    chordae = []; 
end
