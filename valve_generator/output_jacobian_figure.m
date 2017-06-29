function [] = output_jacobian_figure(N_jacobian)

data_str = sprintf('mitral_tree_%d_final_data.mat', N_jacobian); 

if ~exist(data_str, 'file') 
    error('cannot find data file')
end 
    
load(data_str)

J = build_jacobian_bead_slip(valve.leaflets(1)); 

fig = figure; 
spy(J,'k'); 
set(gcf,'color',[1 1 1])
title(sprintf('Jacobian nonzero pattern with N = %d', N_jacobian))

printfig(fig, 'jacobian_nnz')




