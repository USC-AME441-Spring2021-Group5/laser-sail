function F = beamforce(R,P,LAMBDA,PROFILE,TOL,X,Y,D,DIVERGENCE,PLT)
%BEAMFORCE(R,P,PROFILE,X,Y) Calculates force on a spherical due to the
%photon momentum from an incoming laser beam. R is the radius of the
%spherical sail, P is the total power output of the laser, PROFILE is a
%string that idnetifies the laser beam profile

c = 3e8; % speed of light [m/s]

N = length(X);

% Pointing Error = +/- TOL as a length
PE = D*tan(TOL)*(2*rand(1) - 1); 

% Calculating beam width due at distance d due to beam divergence
if DIVERGENCE == 1
    W0 = R;
    W = W0*(sqrt(1 + (LAMBDA*D/pi/(W0^2))^2));
else
    W = R;
end

switch PROFILE
    case 'uniform'
        % creating a uniform distribution for the laser profile
        profile = ones(1,N)./N;
    case 'gaussian'      
        tempY = Y((length(Y)/2)+1:end);
        % Eq From NewpsOps: Power contained within radius R
        temp1 = 1 - exp(-2.*(tempY.^2)./(W.^2));
        temp2 = temp1 - [0 temp1(1:end-1)];
        profile = [temp2(end:-1:1) temp2]./2;

    case 'multi-mode gaussian'      
        tempY = Y((length(Y)/2)+1:end);
        % Eq From NewpsOps: Power contained within radius R
        temp1 = 1 - exp(-2.*(tempY.^2)./(W.^2));
        temp2 = temp1 - [0 temp1(1:end-1)];
        profile = [temp2 temp2(end:-1:1)]./2;

%{
Newport Optics gives an equation for power enclosed within a given radius.
Using this and iterate over the Y vec by solving for the power within
each radius and subtracting the last radius. Like a method of cylinderic
shells for guassian and multi-mode guassian
%}

end

% if PLT == 1
%     % This is specific to guassian
%     hold on
%     plot(diff(temp2).*1e4, tempY(2:end), 'm', ...
%         diff(temp2).*1e4, -tempY(2:end), 'm')
% end

%{
Recall that the beam waist W is defined to be the radius at which the
intensity drops to 1.35% of its axial intensity. Further recall that the
original beam waist is assumed to be the same maginatude as the sail's
radius. Hence, the sail will only sit on 86.5% of the entire beam. Thus,
only 86.5% of the beam's power is incident on the sail. Therefore, PVec
will only sum to .865*P.
%}
PVec = P*profile; % Assigning the power for each individual ray

bHat = [1 0]; % b^ vector. Ie, direction of beam propagation
   
Fvec = zeros(2,N);
nHatVec = Fvec;
% For loop used to calculate the force due to the beam
for j = 1:length(Y)-1
    % Calling normal vector for the point of interest
    nHat = [(Y(j)-Y(j+1)) (X(j+1)-X(j))]/...
        norm([(Y(j)-Y(j+1)) (X(j+1)-X(j))]);
    nHatVec(1,j) = nHat(1); nHatVec(2,j) = nHat(2);
    % Calculating force due to the ray of interest
    Fvec(1,j) = 2*PVec(j)*dot(bHat,nHat)*nHat(1)/c;
    Fvec(2,j) = 2*PVec(j)*dot(bHat,nHat)*nHat(2)/c;
         
    if PLT == 1
    % Plotting every 100th ray and its respective normal and force
    % vectors 
    hold on
    if mod(j,100) == 50 
        poi = [(.5*(X(j)+X(j+1))) (.5*(Y(j)+Y(j+1)))];
        ray = [-2*R poi(1); poi(2) poi(2)];
        plot(ray(1,:), ray(2,:), 'r')
%         quiver(poi(1), poi(2), nHat(1), nHat(2), .25, 'm')
%         quiver(poi(1), poi(2), Fvec(1,j), Fvec(2,j), 50e3, 'b')
    end
   end
end
F = [sum(Fvec(1,:)) sum(Fvec(2,:))];
end
