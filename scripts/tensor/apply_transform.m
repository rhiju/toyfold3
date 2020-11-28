function [next_t,next_R] = apply_transform( start_t, start_R, t, R );
%
% next_t = start_t + start_R * t;
% next_R = start_R * R;\
% 
% Inputs
%   start_t = (3x1) starting translation
%   start_R = (3x3) starting rotation 
%   t       = (3x1) translation to apply
%   R       = (3x3) rotation  to apply
%
% Outputs
%   next_t = (3x1) translation after applying transform
%   next_R = (3x3) rotation after applying transform
%
%
% (C) R. Das, Stanford, 2020

next_t = start_t + start_R * t;
next_R = start_R * R;