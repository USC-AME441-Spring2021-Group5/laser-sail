function [FDrag] = dragforce(rho,v,R)

%DRAGFORCE 
%   

%Check for too few inputs
if nargin < 1
    error('Too few inputs')
end


A = pi*R^2; % cross-sectional area of the sail
C_D = 2.667; % drag coefficient for diffuse reflection of a sphere from SMAD

FDrag = -0.5*rho*v.^2*C_D*A;
