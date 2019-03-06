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

data = [0, 1.6;
0.001, 1.6;
0.002, 1.6;
0.003, 1.6;
0.004, 1.6;
0.005, 1.6;
0.006, 1.6;
0.007, 1.6;
0.008, 1.6;
0.009, 1.6;
0.01, 1.6;
0.011, 1.6;
0.012, 1.6;
0.013, 1.6;
0.014, 1.6;
0.015, 1.6;
0.016, 1.6;
0.017, 1.6;
0.018, 1.6;
0.019, 1.6;
0.02, 1.6;
0.021, 1.6;
0.022, 1.6;
0.023, 1.6;
0.024, 1.6;
0.025, 1.6;
0.026, 1.6;
0.027, 1.6;
0.028, 1.6;
0.029, 1.6;
0.03, 1.6;
0.031, 1.6;
0.032, 1.6;
0.033, 1.6;
0.034, 1.6;
0.035, 1.6;
0.036, 1.6;
0.037, 1.6;
0.038, 1.6;
0.039, 1.6;
0.04, 1.6;
0.041, 1.6;
0.042, 1.6;
0.043, 1.6;
0.044, 1.6;
0.045, 1.6;
0.046, 1.6;
0.047, 1.6;
0.048, 1.6;
0.049, 1.6;
0.05, 1.6;
0.051, 1.6;
0.052, 1.6;
0.053, 1.6;
0.054, 1.6;
0.055, 1.6;
0.056, 1.6;
0.057, 1.6;
0.058, 1.6;
0.059, 1.6;
0.06, 1.6;
0.061, 1.6;
0.062, 1.6;
0.063, 1.6;
0.064, 1.6;
0.065, 1.6;
0.066, 1.6;
0.067, 1.6;
0.068, 1.6;
0.069, 1.6;
0.07, 1.6;
0.071, 1.6;
0.072, 1.6;
0.073, 1.6;
0.074, 1.6;
0.075, 1.6;
0.076, 1.6;
0.077, 1.6;
0.078, 1.6;
0.079, 1.6;
0.08, 1.6;
0.081, 1.6;
0.082, 1.6;
0.083, 1.6;
0.084, 1.6;
0.085, 1.6;
0.086, 1.6;
0.087, 1.6;
0.088, 1.6;
0.089, 1.6;
0.09, 1.6;
0.091, 1.6;
0.092, 1.6;
0.093, 1.6;
0.094, 1.6;
0.095, 1.6;
0.096, 1.6;
0.097, 1.6;
0.098, 1.6;
0.099, 1.6;
0.1, 1.6;
0.101, 1.6;
0.102, 1.6;
0.103, 1.6;
0.104, 1.6;
0.105, 1.6;
0.106, 1.6;
0.107, 1.6;
0.108, 1.6;
0.109, 1.6;
0.11, 1.6;
0.111, 1.6;
0.112, 1.6;
0.113, 1.6;
0.114, 1.6;
0.115, 1.6;
0.116, 1.6;
0.117, 1.6;
0.118, 1.6;
0.119, 1.6;
0.12, 1.6;
0.121, 1.6;
0.122, 1.6;
0.123, 1.6;
0.124, 1.6;
0.125, 1.6;
0.126, 1.6;
0.127, 1.6;
0.128, 1.6;
0.129, 1.6;
0.13, 1.6;
0.131, 1.6;
0.132, 1.6;
0.133, 1.6;
0.134, 1.6;
0.135, 1.6;
0.136, 1.6;
0.137, 1.6;
0.138, 1.6;
0.139, 1.6;
0.14, 1.6;
0.141, 1.6;
0.142, 1.6;
0.143, 1.6;
0.144, 1.6;
0.145, 1.6;
0.146, 1.6;
0.147, 1.6;
0.148, 1.6;
0.149, 1.6;
0.15, 1.6;
0.151, 1.6;
0.152, 1.6;
0.153, 1.6;
0.154, 1.6;
0.155, 1.6;
0.156, 1.6;
0.157, 1.6;
0.158, 1.6;
0.159, 1.6;
0.16, 1.6;
0.161, 1.6;
0.162, 1.6;
0.163, 1.6;
0.164, 1.6;
0.165, 1.6;
0.166, 1.6;
0.167, 1.6;
0.168, 1.6;
0.169, 1.6;
0.17, 1.6;
0.171, 1.6;
0.172, 1.6;
0.173, 1.6;
0.174, 1.6;
0.175, 1.6;
0.176, 1.6;
0.177, 1.6;
0.178, 1.6;
0.179, 1.6;
0.18, 1.6;
0.181, 1.6;
0.182, 1.6;
0.183, 1.6;
0.184, 1.6;
0.185, 1.6;
0.186, 1.6;
0.187, 1.6;
0.188, 1.6;
0.189, 1.6;
0.19, 1.6;
0.191, 1.6;
0.192, 1.6;
0.193, 1.6;
0.194, 1.6;
0.195, 1.6;
0.196, 1.6;
0.197, 1.6;
0.198, 1.6;
0.199, 1.6;
0.2, 1.6;
0.201, 1.6;
0.202, 1.6;
0.203, 1.6;
0.204, 1.6;
0.205, 1.6;
0.206, 1.6;
0.207, 1.6;
0.208, 1.6;
0.209, 1.6;
0.21, 1.6;
0.211, 1.6;
0.212, 1.6;
0.213, 1.6;
0.214, 1.6;
0.215, 1.6;
0.216, 1.6;
0.217, 1.6;
0.218, 1.6;
0.219, 1.6;
0.22, 1.6;
0.221, 1.6;
0.222, 1.6;
0.223, 1.6;
0.224, 1.6;
0.225, 1.6;
0.226, 1.6;
0.227, 1.6;
0.228, 1.6;
0.229, 1.6;
0.23, 1.6;
0.231, 1.6;
0.232, 1.6;
0.233, 1.6;
0.234, 1.6;
0.235, 1.6;
0.236, 1.6;
0.237, 1.6;
0.238, 1.6;
0.239, 1.6;
0.24, 1.6;
0.241, 1.6;
0.242, 1.6;
0.243, 1.6;
0.244, 1.6;
0.245, 1.6;
0.246, 1.6;
0.247, 1.6;
0.248, 1.6;
0.249, 1.6;
0.25, 1.6;
0.251, 1.6;
0.252, 1.6;
0.253, 1.6;
0.254, 1.6;
0.255, 1.6;
0.256, 1.6;
0.257, 1.6;
0.258, 1.6;
0.259, 1.6;
0.26, 1.6;
0.261, 1.6;
0.262, 1.6;
0.263, 1.6;
0.264, 1.6;
0.265, 1.6;
0.266, 1.6;
0.267, 1.6;
0.268, 1.6;
0.269, 1.6;
0.27, 1.6;
0.271, 1.6;
0.272, 1.6;
0.273, 1.6;
0.274, 1.6;
0.275, 1.6;
0.276, 1.6;
0.277, 1.6;
0.278, 1.6;
0.279, 1.6;
0.28, 1.6;
0.281, 1.6;
0.282, 1.6;
0.283, 1.6;
0.284, 1.6;
0.285, 1.6;
0.286, 1.6;
0.287, 1.6;
0.288, 1.6;
0.289, 1.6;
0.29, 1.6;
0.291, 1.6;
0.292, 1.6;
0.293, 1.6;
0.294, 1.6;
0.295, 1.6;
0.296, 1.6;
0.297, 1.6;
0.298, 1.6;
0.299, 1.6;
0.3, 1.6;
0.301, 1.6;
0.302, 1.6;
0.303, 1.6;
0.304, 1.6;
0.305, 1.6;
0.306, 1.6;
0.307, 1.6;
0.308, 1.6;
0.309, 1.6;
0.31, 1.6;
0.311, 1.6;
0.312, 1.6;
0.313, 1.6;
0.314, 1.6;
0.315, 1.6;
0.316, 1.6;
0.317, 1.6;
0.318, 1.6;
0.319, 1.6;
0.32, 1.6;
0.321, 1.6;
0.322, 1.6;
0.323, 1.6;
0.324, 1.6;
0.325, 1.6;
0.326, 1.6;
0.327, 1.6;
0.328, 1.6;
0.329, 1.6;
0.33, 1.6;
0.331, 1.6;
0.332, 1.6;
0.333, 1.6;
0.334, 1.6;
0.335, 1.6;
0.336, 1.6;
0.337, 1.6;
0.338, 1.6;
0.339, 1.6;
0.34, 1.6;
0.341, 1.6;
0.342, 1.6;
0.343, 1.6;
0.344, 1.6;
0.345, 1.6;
0.346, 1.6;
0.347, 1.6;
0.348, 1.6;
0.349, 1.6;
0.35, 1.6;
0.351, 1.6;
0.352, 1.6;
0.353, 1.6;
0.354, 1.6;
0.355, 1.6;
0.356, 1.6;
0.357, 1.6;
0.358, 1.6;
0.359, 1.6;
0.36, 1.6;
0.361, 1.6;
0.362, 1.6;
0.363, 1.6;
0.364, 1.6;
0.365, 1.6;
0.366, 1.6;
0.367, 1.6;
0.368, 1.6;
0.369, 1.6;
0.37, 1.6;
0.371, 1.6;
0.372, 1.6;
0.373, 1.6;
0.374, 1.6;
0.375, 1.6;
0.376, 1.6;
0.377, 1.6;
0.378, 1.6;
0.379, 1.6;
0.38, 1.6;
0.381, 1.6;
0.382, 1.6;
0.383, 1.6;
0.384, 1.6;
0.385, 1.6;
0.386, 1.6;
0.387, 1.6;
0.388, 1.6;
0.389, 1.6;
0.39, 1.6;
0.391, 1.6;
0.392, 1.6;
0.393, 1.6;
0.394, 1.6;
0.395, 1.6;
0.396, 1.6;
0.397, 1.6;
0.398, 1.6;
0.399, 1.6;
0.4, 1.6;
0.401, 1.59919;
0.402, 1.59678;
0.403, 1.59276;
0.404, 1.58714;
0.405, 1.57994;
0.406, 1.57117;
0.407, 1.56085;
0.408, 1.54899;
0.409, 1.53562;
0.41, 1.52078;
0.411, 1.50448;
0.412, 1.48676;
0.413, 1.46766;
0.414, 1.44721;
0.415, 1.42547;
0.416, 1.40246;
0.417, 1.37824;
0.418, 1.35285;
0.419, 1.32635;
0.42, 1.29879;
0.421, 1.27023;
0.422, 1.24072;
0.423, 1.21032;
0.424, 1.17909;
0.425, 1.14711;
0.426, 1.11442;
0.427, 1.0811;
0.428, 1.04721;
0.429, 1.01283;
0.43, 0.978017;
0.431, 0.942846;
0.432, 0.907387;
0.433, 0.871711;
0.434, 0.835892;
0.435, 0.8;
0.436, 0.764108;
0.437, 0.728289;
0.438, 0.692613;
0.439, 0.657154;
0.44, 0.621983;
0.441, 0.587171;
0.442, 0.552786;
0.443, 0.5189;
0.444, 0.48558;
0.445, 0.452893;
0.446, 0.420905;
0.447, 0.389681;
0.448, 0.359282;
0.449, 0.329772;
0.45, 0.301208;
0.451, 0.273649;
0.452, 0.24715;
0.453, 0.221764;
0.454, 0.197543;
0.455, 0.174535;
0.456, 0.152786;
0.457, 0.132341;
0.458, 0.113241;
0.459, 0.0955236;
0.46, 0.0792249;
0.461, 0.0643778;
0.462, 0.0510121;
0.463, 0.0391548;
0.464, 0.0288297;
0.465, 0.0200577;
0.466, 0.0128563;
0.467, 0.00724019;
0.468, 0.00322056;
0.469, 0.000805547;
0.47, 0;
0.471, 0.000805547;
0.472, 0.00322056;
0.473, 0.00724019;
0.474, 0.0128563;
0.475, 0.0200577;
0.476, 0.0288297;
0.477, 0.0391548;
0.478, 0.0510121;
0.479, 0.0643778;
0.48, 0.0792249;
0.481, 0.0955236;
0.482, 0.113241;
0.483, 0.132341;
0.484, 0.152786;
0.485, 0.174535;
0.486, 0.197543;
0.487, 0.221764;
0.488, 0.24715;
0.489, 0.273649;
0.49, 0.301208;
0.491, 0.329772;
0.492, 0.359282;
0.493, 0.389681;
0.494, 0.420905;
0.495, 0.452893;
0.496, 0.48558;
0.497, 0.5189;
0.498, 0.552786;
0.499, 0.587171;
0.5, 0.621983;
0.501, 0.657154;
0.502, 0.692613;
0.503, 0.728289;
0.504, 0.764108;
0.505, 0.8;
0.506, 0.835892;
0.507, 0.871711;
0.508, 0.907387;
0.509, 0.942846;
0.51, 0.978017;
0.511, 1.01283;
0.512, 1.04721;
0.513, 1.0811;
0.514, 1.11442;
0.515, 1.14711;
0.516, 1.17909;
0.517, 1.21032;
0.518, 1.24072;
0.519, 1.27023;
0.52, 1.29879;
0.521, 1.32635;
0.522, 1.35285;
0.523, 1.37824;
0.524, 1.40246;
0.525, 1.42547;
0.526, 1.44721;
0.527, 1.46766;
0.528, 1.48676;
0.529, 1.50448;
0.53, 1.52078;
0.531, 1.53562;
0.532, 1.54899;
0.533, 1.56085;
0.534, 1.57117;
0.535, 1.57994;
0.536, 1.58714;
0.537, 1.59276;
0.538, 1.59678;
0.539, 1.59919;
0.54, 1.6;
0.541, 1.6;
0.542, 1.6;
0.543, 1.6;
0.544, 1.6;
0.545, 1.6;
0.546, 1.6;
0.547, 1.6;
0.548, 1.6;
0.549, 1.6;
0.55, 1.6;
0.551, 1.6;
0.552, 1.6;
0.553, 1.6;
0.554, 1.6;
0.555, 1.6;
0.556, 1.6;
0.557, 1.6;
0.558, 1.6;
0.559, 1.6;
0.56, 1.6;
0.561, 1.6;
0.562, 1.6;
0.563, 1.6;
0.564, 1.6;
0.565, 1.6;
0.566, 1.6;
0.567, 1.6;
0.568, 1.6;
0.569, 1.6;
0.57, 1.6;
0.571, 1.6;
0.572, 1.6;
0.573, 1.6;
0.574, 1.6;
0.575, 1.6;
0.576, 1.6;
0.577, 1.6;
0.578, 1.6;
0.579, 1.6;
0.58, 1.6;
0.581, 1.6;
0.582, 1.6;
0.583, 1.6;
0.584, 1.6;
0.585, 1.6;
0.586, 1.6;
0.587, 1.6;
0.588, 1.6;
0.589, 1.6;
0.59, 1.6;
0.591, 1.6;
0.592, 1.6;
0.593, 1.6;
0.594, 1.6;
0.595, 1.6;
0.596, 1.6;
0.597, 1.6;
0.598, 1.6;
0.599, 1.6;
0.6, 1.6;
0.601, 1.6;
0.602, 1.6;
0.603, 1.6;
0.604, 1.6;
0.605, 1.6;
0.606, 1.6;
0.607, 1.6;
0.608, 1.6;
0.609, 1.6;
0.61, 1.6;
0.611, 1.6;
0.612, 1.6;
0.613, 1.6;
0.614, 1.6;
0.615, 1.6;
0.616, 1.6;
0.617, 1.6;
0.618, 1.6;
0.619, 1.6;
0.62, 1.6;
0.621, 1.6;
0.622, 1.6;
0.623, 1.6;
0.624, 1.6;
0.625, 1.6;
0.626, 1.6;
0.627, 1.6;
0.628, 1.6;
0.629, 1.6;
0.63, 1.6;
0.631, 1.6;
0.632, 1.6;
0.633, 1.6;
0.634, 1.6;
0.635, 1.6;
0.636, 1.6;
0.637, 1.6;
0.638, 1.6;
0.639, 1.6;
0.64, 1.6;
0.641, 1.6;
0.642, 1.6;
0.643, 1.6;
0.644, 1.6;
0.645, 1.6;
0.646, 1.6;
0.647, 1.6;
0.648, 1.6;
0.649, 1.6;
0.65, 1.6;
0.651, 1.6;
0.652, 1.6;
0.653, 1.6;
0.654, 1.6;
0.655, 1.6;
0.656, 1.6;
0.657, 1.6;
0.658, 1.6;
0.659, 1.6;
0.66, 1.6;
0.661, 1.6;
0.662, 1.6;
0.663, 1.6;
0.664, 1.6;
0.665, 1.6;
0.666, 1.6;
0.667, 1.6;
0.668, 1.6;
0.669, 1.6;
0.67, 1.6;
0.671, 1.6;
0.672, 1.6;
0.673, 1.6;
0.674, 1.6;
0.675, 1.6;
0.676, 1.6;
0.677, 1.6;
0.678, 1.6;
0.679, 1.6;
0.68, 1.6;
0.681, 1.6;
0.682, 1.6;
0.683, 1.6;
0.684, 1.6;
0.685, 1.6;
0.686, 1.6;
0.687, 1.6;
0.688, 1.6;
0.689, 1.6;
0.69, 1.6;
0.691, 1.6;
0.692, 1.6;
0.693, 1.6;
0.694, 1.6;
0.695, 1.6;
0.696, 1.6;
0.697, 1.6;
0.698, 1.6;
0.699, 1.6;
0.7, 1.6;
0.701, 1.6;
0.702, 1.6;
0.703, 1.6;
0.704, 1.6;
0.705, 1.6;
0.706, 1.6;
0.707, 1.6;
0.708, 1.6;
0.709, 1.6;
0.71, 1.6;
0.711, 1.6;
0.712, 1.6;
0.713, 1.6;
0.714, 1.6;
0.715, 1.6;
0.716, 1.6;
0.717, 1.6;
0.718, 1.6;
0.719, 1.6;
0.72, 1.6;
0.721, 1.6;
0.722, 1.6;
0.723, 1.6;
0.724, 1.6;
0.725, 1.6;
0.726, 1.6;
0.727, 1.6;
0.728, 1.6;
0.729, 1.6;
0.73, 1.6;
0.731, 1.6;
0.732, 1.6;
0.733, 1.6;
0.734, 1.6;
0.735, 1.6;
0.736, 1.6;
0.737, 1.6;
0.738, 1.6;
0.739, 1.6;
0.74, 1.6;
0.741, 1.6;
0.742, 1.6;
0.743, 1.6;
0.744, 1.6;
0.745, 1.6;
0.746, 1.6;
0.747, 1.6;
0.748, 1.6;
0.749, 1.6;
0.75, 1.6;
0.751, 1.6;
0.752, 1.6;
0.753, 1.6;
0.754, 1.6;
0.755, 1.6;
0.756, 1.6;
0.757, 1.6;
0.758, 1.6;
0.759, 1.6;
0.76, 1.6;
0.761, 1.6;
0.762, 1.6;
0.763, 1.6;
0.764, 1.6;
0.765, 1.6;
0.766, 1.6;
0.767, 1.6;
0.768, 1.6;
0.769, 1.6;
0.77, 1.6;
0.771, 1.6;
0.772, 1.6;
0.773, 1.6;
0.774, 1.6;
0.775, 1.6;
0.776, 1.6;
0.777, 1.6;
0.778, 1.6;
0.779, 1.6;
0.78, 1.6;
0.781, 1.6;
0.782, 1.6;
0.783, 1.6;
0.784, 1.6;
0.785, 1.6;
0.786, 1.6;
0.787, 1.6;
0.788, 1.6;
0.789, 1.6;
0.79, 1.6;
0.791, 1.6;
0.792, 1.6;
0.793, 1.6;
0.794, 1.6;
0.795, 1.6;
0.796, 1.6;
0.797, 1.6;
0.798, 1.6;
0.799, 1.6;
0.8, 1.6;
0.801, 1.6;
0.802, 1.6;
0.803, 1.6;
0.804, 1.6;
0.805, 1.6;
0.806, 1.6;
0.807, 1.6;
0.808, 1.6;
0.809, 1.6;
0.81, 1.6;
0.811, 1.6;
0.812, 1.6;
0.813, 1.6;
0.814, 1.6;
0.815, 1.6;
0.816, 1.6;
0.817, 1.6;
0.818, 1.6;
0.819, 1.6;
0.82, 1.6;
0.821, 1.6;
0.822, 1.6;
0.823, 1.6;
0.824, 1.6;
0.825, 1.6;
0.826, 1.6;
0.827, 1.6;
0.828, 1.6;
0.829, 1.6;
0.83, 1.6;
0.831, 1.6;
0.832, 1.6;
0.833, 1.6;
0.834, 1.6;
0.835, 1.6;
0.836, 1.6;
0.837, 1.6;
0.838, 1.6;
0.839, 1.6;
0.84, 1.6;
0.841, 1.6;
0.842, 1.6;
0.843, 1.6;
0.844, 1.6;
0.845, 1.6;
0.846, 1.6;
0.847, 1.6;
0.848, 1.6;
0.849, 1.6;
0.85, 1.6;
0.851, 1.6;
0.852, 1.6;
0.853, 1.6;
0.854, 1.6;
0.855, 1.6;
0.856, 1.6;
0.857, 1.6;
0.858, 1.6;
0.859, 1.6;
0.86, 1.6;
0.861, 1.6;
0.862, 1.6;
0.863, 1.6;
0.864, 1.6;
0.865, 1.6;
0.866, 1.6;
0.867, 1.6;
0.868, 1.6;
0.869, 1.6;
0.87, 1.6;
0.871, 1.6;
0.872, 1.6;
0.873, 1.6;
0.874, 1.6;
0.875, 1.6;
0.876, 1.6;
0.877, 1.6;
0.878, 1.6;
0.879, 1.6;
0.88, 1.6;
0.881, 1.6;
0.882, 1.6;
0.883, 1.6;
0.884, 1.6;
0.885, 1.6;
0.886, 1.6;
0.887, 1.6;
0.888, 1.6;
0.889, 1.6;
0.89, 1.6;
0.891, 1.6;
0.892, 1.6;
0.893, 1.6;
0.894, 1.6;
0.895, 1.6;
0.896, 1.6;
0.897, 1.6;
0.898, 1.6;
0.899, 1.6;
0.9, 1.6;
0.901, 1.6;
0.902, 1.6;
0.903, 1.6;
0.904, 1.6;
0.905, 1.6;
0.906, 1.6;
0.907, 1.6;
0.908, 1.6;
0.909, 1.6;
0.91, 1.6;
0.911, 1.6;
0.912, 1.6;
0.913, 1.6;
0.914, 1.6;
0.915, 1.6;
0.916, 1.6;
0.917, 1.6;
0.918, 1.6;
0.919, 1.6;
0.92, 1.6;
0.921, 1.6;
0.922, 1.6;
0.923, 1.6;
0.924, 1.6;
0.925, 1.6;
0.926, 1.6;
0.927, 1.6;
0.928, 1.6;
0.929, 1.6;
0.93, 1.6;
0.931, 1.6;
0.932, 1.6;
0.933, 1.6;
0.934, 1.6;
0.935, 1.6;
0.936, 1.6;
0.937, 1.6;
0.938, 1.6;
0.939, 1.6;
0.94, 1.6;
0.941, 1.6;
0.942, 1.6;
0.943, 1.6;
0.944, 1.6;
0.945, 1.6;
0.946, 1.6;
0.947, 1.6;
0.948, 1.6;
0.949, 1.6;
0.95, 1.6;
0.951, 1.6;
0.952, 1.6;
0.953, 1.6;
0.954, 1.6;
0.955, 1.6;
0.956, 1.6;
0.957, 1.6;
0.958, 1.6;
0.959, 1.6;
0.96, 1.6;
0.961, 1.6;
0.962, 1.6;
0.963, 1.6;
0.964, 1.6;
0.965, 1.6;
0.966, 1.6;
0.967, 1.6;
0.968, 1.6;
0.969, 1.6;
0.97, 1.6;
0.971, 1.6;
0.972, 1.6;
0.973, 1.6;
0.974, 1.6;
0.975, 1.6;
0.976, 1.6;
0.977, 1.6;
0.978, 1.6;
0.979, 1.6;
0.98, 1.6;
0.981, 1.6;
0.982, 1.6;
0.983, 1.6;
0.984, 1.6;
0.985, 1.6;
0.986, 1.6;
0.987, 1.6;
0.988, 1.6;
0.989, 1.6;
0.99, 1.6;
0.991, 1.6;
0.992, 1.6;
0.993, 1.6;
0.994, 1.6;
0.995, 1.6;
0.996, 1.6;
0.997, 1.6;
0.998, 1.6;
0.999, 1.6;
1, 1.6;
1.001, 1.6;
1.002, 1.6;
1.003, 1.6;
1.004, 1.6;
1.005, 1.6;
1.006, 1.6;
1.007, 1.6;
1.008, 1.6;
1.009, 1.6;
1.01, 1.6;
1.011, 1.6;
1.012, 1.6;
1.013, 1.6;
1.014, 1.6;
1.015, 1.6;
1.016, 1.6;
1.017, 1.6;
1.018, 1.6;
1.019, 1.6;
1.02, 1.6;
1.021, 1.6;
1.022, 1.6;
1.023, 1.6;
1.024, 1.6;
1.025, 1.6;
1.026, 1.6;
1.027, 1.6;
1.028, 1.6;
1.029, 1.6;
1.03, 1.6;
1.031, 1.6;
1.032, 1.6;
1.033, 1.6;
1.034, 1.6;
1.035, 1.6;
1.036, 1.6;
1.037, 1.6;
1.038, 1.6;
1.039, 1.6;
1.04, 1.6;
1.041, 1.6;
1.042, 1.6;
1.043, 1.6;
1.044, 1.6;
1.045, 1.6;
1.046, 1.6;
1.047, 1.6;
1.048, 1.6;
1.049, 1.6;
1.05, 1.6;
1.051, 1.6;
1.052, 1.6;
1.053, 1.6;
1.054, 1.6;
1.055, 1.6;
1.056, 1.6;
1.057, 1.6;
1.058, 1.6;
1.059, 1.6;
1.06, 1.6;
1.061, 1.6;
1.062, 1.6;
1.063, 1.6;
1.064, 1.6;
1.065, 1.6;
1.066, 1.6;
1.067, 1.6;
1.068, 1.6;
1.069, 1.6;
1.07, 1.6;
1.071, 1.6;
1.072, 1.6;
1.073, 1.6;
1.074, 1.6;
1.075, 1.6;
1.076, 1.6;
1.077, 1.6;
1.078, 1.6;
1.079, 1.6;
1.08, 1.6;
1.081, 1.6;
1.082, 1.6;
1.083, 1.6;
1.084, 1.6;
1.085, 1.6;
1.086, 1.6;
1.087, 1.6;
1.088, 1.6;
1.089, 1.6;
1.09, 1.6;
1.091, 1.6;
1.092, 1.6;
1.093, 1.6;
1.094, 1.6;
1.095, 1.6;
1.096, 1.6;
1.097, 1.6;
1.098, 1.6;
1.099, 1.6;
1.1, 1.6;
1.101, 1.6;
1.102, 1.6;
1.103, 1.6;
1.104, 1.6;
1.105, 1.6;
1.106, 1.6;
1.107, 1.6;
1.108, 1.6;
1.109, 1.6;
1.11, 1.6;
1.111, 1.6;
1.112, 1.6;
1.113, 1.6;
1.114, 1.6;
1.115, 1.6;
1.116, 1.6;
1.117, 1.6;
1.118, 1.6;
1.119, 1.6;
1.12, 1.6;
1.121, 1.6;
1.122, 1.6;
1.123, 1.6;
1.124, 1.6;
1.125, 1.6;
1.126, 1.6;
1.127, 1.6;
1.128, 1.6;
1.129, 1.6;
1.13, 1.6;
1.131, 1.6;
1.132, 1.6;
1.133, 1.6;
1.134, 1.6;
1.135, 1.6;
1.136, 1.6;
1.137, 1.6;
1.138, 1.6;
1.139, 1.6;
1.14, 1.6;
1.141, 1.6;
1.142, 1.6;
1.143, 1.6;
1.144, 1.6;
1.145, 1.6;
1.146, 1.6;
1.147, 1.6;
1.148, 1.6;
1.149, 1.6;
1.15, 1.6;
1.151, 1.6;
1.152, 1.6;
1.153, 1.6;
1.154, 1.6;
1.155, 1.6;
1.156, 1.6;
1.157, 1.6;
1.158, 1.6;
1.159, 1.6;
1.16, 1.6;
1.161, 1.6;
1.162, 1.6;
1.163, 1.6;
1.164, 1.6;
1.165, 1.6;
1.166, 1.6;
1.167, 1.6;
1.168, 1.6;
1.169, 1.6;
1.17, 1.6;
1.171, 1.6;
1.172, 1.6;
1.173, 1.6;
1.174, 1.6;
1.175, 1.6;
1.176, 1.6;
1.177, 1.6;
1.178, 1.6;
1.179, 1.6;
1.18, 1.6;
1.181, 1.6;
1.182, 1.6;
1.183, 1.6;
1.184, 1.6;
1.185, 1.6;
1.186, 1.6;
1.187, 1.6;
1.188, 1.6;
1.189, 1.6;
1.19, 1.6;
1.191, 1.6;
1.192, 1.6;
1.193, 1.6;
1.194, 1.6;
1.195, 1.6;
1.196, 1.6;
1.197, 1.6;
1.198, 1.6;
1.199, 1.6;
1.2, 1.6;
1.201, 1.59919;
1.202, 1.59678;
1.203, 1.59276;
1.204, 1.58714;
1.205, 1.57994;
1.206, 1.57117;
1.207, 1.56085;
1.208, 1.54899;
1.209, 1.53562;
1.21, 1.52078;
1.211, 1.50448;
1.212, 1.48676;
1.213, 1.46766;
1.214, 1.44721;
1.215, 1.42547;
1.216, 1.40246;
1.217, 1.37824;
1.218, 1.35285;
1.219, 1.32635;
1.22, 1.29879;
1.221, 1.27023;
1.222, 1.24072;
1.223, 1.21032;
1.224, 1.17909;
1.225, 1.14711;
1.226, 1.11442;
1.227, 1.0811;
1.228, 1.04721;
1.229, 1.01283;
1.23, 0.978017;
1.231, 0.942846;
1.232, 0.907387;
1.233, 0.871711;
1.234, 0.835892;
1.235, 0.8;
1.236, 0.764108;
1.237, 0.728289;
1.238, 0.692613;
1.239, 0.657154;
1.24, 0.621983;
1.241, 0.587171;
1.242, 0.552786;
1.243, 0.5189;
1.244, 0.48558;
1.245, 0.452893;
1.246, 0.420905;
1.247, 0.389681;
1.248, 0.359282;
1.249, 0.329772;
1.25, 0.301208;
1.251, 0.273649;
1.252, 0.24715;
1.253, 0.221764;
1.254, 0.197543;
1.255, 0.174535;
1.256, 0.152786;
1.257, 0.132341;
1.258, 0.113241;
1.259, 0.0955236;
1.26, 0.0792249;
1.261, 0.0643778;
1.262, 0.0510121;
1.263, 0.0391548;
1.264, 0.0288297;
1.265, 0.0200577;
1.266, 0.0128563;
1.267, 0.00724019;
1.268, 0.00322056;
1.269, 0.000805547;
1.27, 0;
1.271, 0.000805547;
1.272, 0.00322056;
1.273, 0.00724019;
1.274, 0.0128563;
1.275, 0.0200577;
1.276, 0.0288297;
1.277, 0.0391548;
1.278, 0.0510121;
1.279, 0.0643778;
1.28, 0.0792249;
1.281, 0.0955236;
1.282, 0.113241;
1.283, 0.132341;
1.284, 0.152786;
1.285, 0.174535;
1.286, 0.197543;
1.287, 0.221764;
1.288, 0.24715;
1.289, 0.273649;
1.29, 0.301208;
1.291, 0.329772;
1.292, 0.359282;
1.293, 0.389681;
1.294, 0.420905;
1.295, 0.452893;
1.296, 0.48558;
1.297, 0.5189;
1.298, 0.552786;
1.299, 0.587171;
1.3, 0.621983;
1.301, 0.657154;
1.302, 0.692613;
1.303, 0.728289;
1.304, 0.764108;
1.305, 0.8;
1.306, 0.835892;
1.307, 0.871711;
1.308, 0.907387;
1.309, 0.942846;
1.31, 0.978017;
1.311, 1.01283;
1.312, 1.04721;
1.313, 1.0811;
1.314, 1.11442;
1.315, 1.14711;
1.316, 1.17909;
1.317, 1.21032;
1.318, 1.24072;
1.319, 1.27023;
1.32, 1.29879;
1.321, 1.32635;
1.322, 1.35285;
1.323, 1.37824;
1.324, 1.40246;
1.325, 1.42547;
1.326, 1.44721;
1.327, 1.46766;
1.328, 1.48676;
1.329, 1.50448;
1.33, 1.52078;
1.331, 1.53562;
1.332, 1.54899;
1.333, 1.56085;
1.334, 1.57117;
1.335, 1.57994;
1.336, 1.58714;
1.337, 1.59276;
1.338, 1.59678;
1.339, 1.59919;
1.34, 1.6;
1.341, 1.6;
1.342, 1.6;
1.343, 1.6;
1.344, 1.6;
1.345, 1.6;
1.346, 1.6;
1.347, 1.6;
1.348, 1.6;
1.349, 1.6;
1.35, 1.6;
1.351, 1.6;
1.352, 1.6;
1.353, 1.6;
1.354, 1.6;
1.355, 1.6;
1.356, 1.6;
1.357, 1.6;
1.358, 1.6;
1.359, 1.6;
1.36, 1.6;
1.361, 1.6;
1.362, 1.6;
1.363, 1.6;
1.364, 1.6;
1.365, 1.6;
1.366, 1.6;
1.367, 1.6;
1.368, 1.6;
1.369, 1.6;
1.37, 1.6;
1.371, 1.6;
1.372, 1.6;
1.373, 1.6;
1.374, 1.6;
1.375, 1.6;
1.376, 1.6;
1.377, 1.6;
1.378, 1.6;
1.379, 1.6;
1.38, 1.6;
1.381, 1.6;
1.382, 1.6;
1.383, 1.6;
1.384, 1.6;
1.385, 1.6;
1.386, 1.6;
1.387, 1.6;
1.388, 1.6;
1.389, 1.6;
1.39, 1.6;
1.391, 1.6;
1.392, 1.6;
1.393, 1.6;
1.394, 1.6;
1.395, 1.6;
1.396, 1.6;
1.397, 1.6;
1.398, 1.6;
1.399, 1.6;
1.4, 1.6;
1.401, 1.6;
1.402, 1.6;
1.403, 1.6;
1.404, 1.6;
1.405, 1.6;
1.406, 1.6;
1.407, 1.6;
1.408, 1.6;
1.409, 1.6;
1.41, 1.6;
1.411, 1.6;
1.412, 1.6;
1.413, 1.6;
1.414, 1.6;
1.415, 1.6;
1.416, 1.6;
1.417, 1.6;
1.418, 1.6;
1.419, 1.6;
1.42, 1.6;
1.421, 1.6;
1.422, 1.6;
1.423, 1.6;
1.424, 1.6;
1.425, 1.6;
1.426, 1.6;
1.427, 1.6;
1.428, 1.6;
1.429, 1.6;
1.43, 1.6;
1.431, 1.6;
1.432, 1.6;
1.433, 1.6;
1.434, 1.6;
1.435, 1.6;
1.436, 1.6;
1.437, 1.6;
1.438, 1.6;
1.439, 1.6;
1.44, 1.6;
1.441, 1.6;
1.442, 1.6;
1.443, 1.6;
1.444, 1.6;
1.445, 1.6;
1.446, 1.6;
1.447, 1.6;
1.448, 1.6;
1.449, 1.6;
1.45, 1.6;
1.451, 1.6;
1.452, 1.6;
1.453, 1.6;
1.454, 1.6;
1.455, 1.6;
1.456, 1.6;
1.457, 1.6;
1.458, 1.6;
1.459, 1.6;
1.46, 1.6;
1.461, 1.6;
1.462, 1.6;
1.463, 1.6;
1.464, 1.6;
1.465, 1.6;
1.466, 1.6;
1.467, 1.6;
1.468, 1.6;
1.469, 1.6;
1.47, 1.6;
1.471, 1.6;
1.472, 1.6;
1.473, 1.6;
1.474, 1.6;
1.475, 1.6;
1.476, 1.6;
1.477, 1.6;
1.478, 1.6;
1.479, 1.6;
1.48, 1.6;
1.481, 1.6;
1.482, 1.6;
1.483, 1.6;
1.484, 1.6;
1.485, 1.6;
1.486, 1.6;
1.487, 1.6;
1.488, 1.6;
1.489, 1.6;
1.49, 1.6;
1.491, 1.6;
1.492, 1.6;
1.493, 1.6;
1.494, 1.6;
1.495, 1.6;
1.496, 1.6;
1.497, 1.6;
1.498, 1.6;
1.499, 1.6;
1.5, 1.6;
1.501, 1.6;
1.502, 1.6;
1.503, 1.6;
1.504, 1.6;
1.505, 1.6;
1.506, 1.6;
1.507, 1.6;
1.508, 1.6;
1.509, 1.6;
1.51, 1.6;
1.511, 1.6;
1.512, 1.6;
1.513, 1.6;
1.514, 1.6;
1.515, 1.6;
1.516, 1.6;
1.517, 1.6;
1.518, 1.6;
1.519, 1.6;
1.52, 1.6;
1.521, 1.6;
1.522, 1.6;
1.523, 1.6;
1.524, 1.6;
1.525, 1.6;
1.526, 1.6;
1.527, 1.6;
1.528, 1.6;
1.529, 1.6;
1.53, 1.6;
1.531, 1.6;
1.532, 1.6;
1.533, 1.6;
1.534, 1.6;
1.535, 1.6;
1.536, 1.6;
1.537, 1.6;
1.538, 1.6;
1.539, 1.6;
1.54, 1.6;
1.541, 1.6;
1.542, 1.6;
1.543, 1.6;
1.544, 1.6;
1.545, 1.6;
1.546, 1.6;
1.547, 1.6;
1.548, 1.6;
1.549, 1.6;
1.55, 1.6;
1.551, 1.6;
1.552, 1.6;
1.553, 1.6;
1.554, 1.6;
1.555, 1.6;
1.556, 1.6;
1.557, 1.6;
1.558, 1.6;
1.559, 1.6;
1.56, 1.6;
1.561, 1.6;
1.562, 1.6;
1.563, 1.6;
1.564, 1.6;
1.565, 1.6;
1.566, 1.6;
1.567, 1.6;
1.568, 1.6;
1.569, 1.6;
1.57, 1.6;
1.571, 1.6;
1.572, 1.6;
1.573, 1.6;
1.574, 1.6;
1.575, 1.6;
1.576, 1.6;
1.577, 1.6;
1.578, 1.6;
1.579, 1.6;
1.58, 1.6;
1.581, 1.6;
1.582, 1.6;
1.583, 1.6;
1.584, 1.6;
1.585, 1.6;
1.586, 1.6;
1.587, 1.6;
1.588, 1.6;
1.589, 1.6;
1.59, 1.6;
1.591, 1.6;
1.592, 1.6;
1.593, 1.6;
1.594, 1.6;
1.595, 1.6;
1.596, 1.6;
1.597, 1.6;
1.598, 1.6;
1.599, 1.6;
1.6, 1.6;
];
plot(data(:,1), data(:,2))
