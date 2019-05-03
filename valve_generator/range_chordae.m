function range = range_chordae(chordae, idx_chordae, tree_idx)
% 
% Returns the range of indices corresponding to 
% 
% Input
%     total_internal    Total number of scalar points in the leaflet
%                       Does not include b.c.
%     N_chordae         Number of (vector) points in the tree of chordae 
%     idx_chordae       Current index in the (vector index) direction in the chordae 
%     left              Whether to use the left side of the chordae 
% 
% Output
%     range             Length three array of global indices 
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

if isfield(chordae(tree_idx), 'targets_for_bcs') && chordae(tree_idx).targets_for_bcs

    % everything up by one, because root is now an internal variable 
    range = chordae(tree_idx).min_global_idx + 3*idx_chordae + (0:2);

else 
    
    % default case with root as bc
    if idx_chordae == 0
        range = [];
    else 
        range = chordae(tree_idx).min_global_idx + 3*(idx_chordae-1) + (0:2);
    end 

end 

