function [k R] = get_rest_len_and_spring_constant_linear(X, X_nbr, tension, strain)
% 
% Given two points, a tension between them, and a desired specified strain 
% Compute rest length and spring constant that gives the current force 
% 
% Input: 
%     X          Point 
%     X_nbr      Other point 
%     tension    Force between the two (not that this must be an actual force, not a density)
%     strain     Desired strain 
% 
% Output: 
%     k          Spring constant k 
%                multiplies RELATIVE strain to get force 
%                This implies that k has units of force 
%     R          Rest length 
% 

k = tension / strain; 
R = norm(X - X_nbr) / (strain + 1); 

