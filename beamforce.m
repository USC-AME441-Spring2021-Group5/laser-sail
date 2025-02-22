function [F,tau] = beamforce(R,P,LAMBDA,PROFILE,TOL,X,Y,D,DIVERGENCE,PLT,COM,SR)
%BEAMFORCE(R,P,PROFILE,X,Y) Calculates force on a spherical sail due to the
%photon momentum from an incoming laser beam. R is the radius of the
%spherical sail, P is the total power output of the laser, PROFILE is a
%string that idnetifies the laser beam profile

c = 3e8; % speed of light [m/s]
N = length(X);

% Pointing Error = +/- TOL as a length... May have to reevaluate
PE = D*tan(TOL)*(2*rand(1) - 1); 
% Q = [cos(TOL) -sin(TOL); sin(TOL) cos(TOL)];

% Calculating beam width due at distance d due to beam divergence
if DIVERGENCE == 1
    W0 = R;
    W = W0*(sqrt(1 + (LAMBDA*D/pi/(W0^2))^2));
else
    W = R;
end

if SR == 1
    sigh = 2e-6;
    T = 1e-3;
    siga = (1/sqrt(2))*atan(2*sigh/T);
end

switch PROFILE
    case 'uniform'
        % creating a uniform distribution for the laser profile
        profile = ones(1,N)./N;
    case 'gaussian'      
        tempY = Y((length(Y)/2)+1:end) + PE;
        % Eq From NewpsOps: Power contained within radius R
        temp1 = 1 - exp(-2.*(tempY.^2)./(W.^2));
        temp2 = temp1 - [0 temp1(1:end-1)];
        profile = [temp2(end:-1:1) temp2]./2;

    case 'multi-mode gaussian'      
        tempY = Y((length(Y)/2)+1:end) + PE;
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
if ~isreal(PVec)
    fprintf('PVec is imaginary!\n')
end

bHat = [1 0]; % b^ vector. Ie, direction of beam propagation
   
Fvec = zeros(2,N-1);
tauVec = zeros(3,N-1);
nHatVec = Fvec;
% For loop used to calculate the force due to the beam
for j = 1:length(Y)-1
    % Calling normal vector for the point of interest
    nHat = [(Y(j)-Y(j+1)) (X(j+1)-X(j))]/...
        norm([(Y(j)-Y(j+1)) (X(j+1)-X(j))]);
    if SR == 1
%         syms('a')
%         assume(a,'real')
%         random = randi(1,1);
%         prob = T/(2*sigh*sqrt(pi)*(cos(a)^2))*...
%             exp((T*(tan(a)^2/(2*sigh)))^2);
%         alph = double(solve(prob==random,a))
        alph = siga*randn;
        Q = [cosd(alph) -sind(alph); sind(alph) cosd(alph)];
        nHat = nHat*Q;
    end
    nHatVec(1,j) = nHat(1); nHatVec(2,j) = nHat(2);

    % Calculating force due to the ray of interest
    Fvec(1,j) = 2*(.5*(PVec(j)+PVec(j+1)))*dot(bHat,nHat)*nHat(1)/c;
    Fvec(2,j) = 2*(.5*(PVec(j)+PVec(j+1)))*dot(bHat,nHat)*nHat(2)/c;
    
    XY = [(.5*(X(j)+X(j+1))) (.5*(Y(j)+Y(j+1))) 0];
    r = XY - [COM 0];
    tauVec(:,j) = 2*(.5*(PVec(j)+PVec(j+1)))*dot(bHat,nHat)*...
        cross(r,[nHat 0])/c;
         
    if PLT == 1
    % Plotting every 100th ray and its respective normal and force
    % vectors 
    hold on
    if mod(j,100) == 50 
        poi = [(.5*(X(j)+X(j+1))) (.5*(Y(j)+Y(j+1)))];
        ray = [-1e3*R poi(1); poi(2) poi(2)];
        plot(ray(1,:), ray(2,:), 'r')
        quiver(poi(1), poi(2), nHat(1), nHat(2), .1, 'm')
        quiver(poi(1), poi(2), Fvec(1,j), Fvec(2,j), 1e6, 'b')
    end
    end
%    fprintf('\nray %i:\n%.16f\t%.16f\n%.16f\t%.16f\n', j, nHat(1), nHat(2),...
%        Fvec(1,j), Fvec(2,j))
end
F = [sum(Fvec(1,:)) sum(Fvec(2,:))];
tau = sum(tauVec(3,:)); 
% In the 2D, the cross product will always result in z-dir
% hence the above tau
end
