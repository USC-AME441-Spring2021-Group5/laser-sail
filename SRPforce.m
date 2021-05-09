function F_SRP = SRPforce(R,N,phi_prime,q)
%SRPFORCE(R,N,theta,Q) calculates the solar radiation pressure force depending on 
%the area of ilumination where R is the radius of the spherical sail, N the
%number of rays, theta [deg]  angle of the suns rays wrt the laser force
%vector [deg], (q [unitless]), from 0 to 1 where 1 is perfect reflectance 


c = 3e8;         %speed of light [m/s]
A_s = pi*R^2;    %cross sectional area of the sunlit hemisphere of the sail [m^2]
psi = 1366;      %solar constant adjusted for actual distance from the sun [W/m^2]

% finding the area per ray
A_dy = A_s/N;

SRPvec = zeros(2,N);
alpha = linspace(-90,90,N+1);

i = 1;
% For loop used to calculate the force due to the solar radiation
for j = alpha(1,:)
    % Calculating the force due to the ray of interest
    SRPvec(1,i) = (psi/c)*A_dy*(1+q).*cosd(j).*cosd(j); 
    SRPvec(2,i) = (psi/c)*A_dy*(1+q).*cosd(j).*sind(j); 
    i = i + 1;
end

% resultant SRP force vector wrt to the sun
SRPvecTOT = [sum(SRPvec(1,:)) sum(SRPvec(2,:))];

% resultnat SRP force vector wrt to sail refernce frame
F_SRP = [norm(SRPvecTOT)*cosd(phi_prime) norm(SRPvecTOT)*sind(phi_prime)];

end
