function F_SRP = SRPforce(R,theta,q)
%SRPFORCE(R,THETA,Q) calculates the solar radiation pressure force depending on 
%the area of ilumination where R is the radius of the spherical sail,
%angle of incident solar radiation (theta [deg]) and the reflectace factor (q [unitless]),
%from 0 to 1 where 1 is perfect reflectance


c = 3e8;        %speed of light [m/s]
A_s = 2*pi*R^2;   %sunlit hemisphere area of the sail [m^2]
phi = 1366;     %solar constant adjusted for actual distance fro the sun [W/m^2]

F_SRP = (phi/c)*A_s*(1+q)*cosd(theta); %[N]