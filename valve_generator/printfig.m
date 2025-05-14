function [] = printfig(fig, title)

%
% function [] = printfig(fig, title)
%
% Print figure to eps in current directory
% 
% Input: 
%     fig     Figure window. Set with
%                 fig = figure 
%     title   String for title 
% 
set(fig, 'Renderer', 'Painters');
print(fig, '-depsc', title);
