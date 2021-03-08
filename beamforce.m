function F = beamforce(R,P,PROFILE,X,Y)
%BEAMFORCE(R,P,PROFILE,X,Y) Calculates force on a spherical due to the
%photon momentum from an incoming laser beam. R is the radius of the
%spherical sail, P is the total power output of the laser, PROFILE is a
%string that idnetifies the laser beam profile


c = 3e8; % speed of light [m/s]

N = length(X);
% dy = 2*R/N;
std = R/sqrt(2*log(2));
switch PROFILE
    case 'uniform'
        profile = ones(1,N)./N;%./(2*R); % creating a uniform distribution for the laser profile
    case 'gaussian'
        profile = normpdf(Y,0,std);
        % Maybe ask on piazza about how to determine half beam width from std?
    case 'multi-mode gaussian'
        profile = normpdf(Y,.5*R,std) + normpdf(Y,-.5*R,std);
end

PVec = P*profile; % Assigning the power for each individual ray
% figure(2); plot(PVec,Y,P*normpdf(Y,.5*R,R),Y,P*normpdf(Y,-.5*R,R),Y)

bHat = [1 0]; % b^ vector. Ie, direction of beam propagation

% Finding the normal vector for every point on the sail
% dxVec = gradient(X); dyVec = gradient(Y);
   
Fvec = zeros(2,N);
nHatVec = Fvec;
% For loop used to calculate the force due to the beam
for j = 1:length(Y)-1
    % Calling normal vector for the point of interest
%     nHat = [dyVec(j) -dxVec(j)]/norm([dyVec(j) -dxVec(j)]);
    nHat = [(Y(j)-Y(j+1)) (X(j+1)-X(j))]/...
        norm([(Y(j)-Y(j+1)) (X(j+1)-X(j))]);
    nHatVec(1,j) = nHat(1); nHatVec(2,j) = nHat(2);
    % Calculating force due to the ray of interest
    Fvec(1,j) = 2*PVec(j)*dot(bHat,nHat)*nHat(1)/c; %*abs(dyVec(j))/c;
    Fvec(2,j) = 2*PVec(j)*dot(bHat,nHat)*nHat(2)/c;
    % Adding this ray's force to the total force
%     F = F + Fnow;
    
    %{
    May need to split the force calculatiosn up to the different cases as
    the integrals may be different... (see newport optics refs)
    %}
    
    %{
    THE FORCE EQUATION COULD BE WRONG...
    IT DOESN'T HAVE CONSISTENT UNITS!!!
    ...Yes it does... P is Power FLUX, not just Power.
    CORRECT UNITS WOULD NEED DIVISION BY SOME AREA.
    --> The power needs divsion by an area
    FROM ONLINE RESOURCES, THE FORCE WOULD BE AS FOLLOWS:
    FOR A PERFECT REFLECTOR
    F = 2*u;
    u = INT(I/A);
    I (FOR A GUASSIAN BEAM) = 2*P*S/(c*PI*W^2)
    WHERE W IS THE BEAM WIDTH AND S IS THE SURFACE AREA.
    %}
        
    % Plotting every 100th ray and its respective normal and force
    % vectors
%     figure(1); 
    hold on
    if mod(j,100) == 50 
        poi = [(.5*(X(j)+X(j+1))) (.5*(Y(j)+Y(j+1)))];
        ray = [-2*R poi(1); poi(2) poi(2)];
        plot(ray(1,:), ray(2,:), 'r')
%         quiver(poi(1), poi(2), nHat(1), nHat(2), .25, 'm')
        quiver(poi(1), poi(2), Fvec(1,j), Fvec(2,j), 50e3, 'b')
   end
end
% F = F/pi/((2*R)^2);
F = [sum(Fvec(1,:)) sum(Fvec(2,:))];


end

