%{
AME 441 Group 5
Passive Stability of a Laser-Driven Light Sail
Spring 2021
%}
clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create "Real-Time" Plots? 1 == Yes, 0 == No
plt = 0;
PLTs = 1;
%% Create "Real-Time" Print Statements? 1 == Yes, 0 == No
stmnt = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deciding what perturbations to include
% 1 for include, 0 for don't
Drag = 0;
PointingTol = 0;
BeamDivergence = 1;

%% Setting up inital paramters and known constants

m = 1;                  % Mass of the sail [kg]
R = .1;                 % Radius of sail [m]
center = [0 0];         % Initial position of sail center
v = [0 0];              % Initial velocity of sail rel to beam sat [m/s]
P = 30;                 % Power of laser beam [W]
lambda = 980e-9;        % Wavelength of laser
N = 1e6;               % Number of rays --> KEEP EVEN
if mod(N,2) ~= 0
    error('N must be even to create power distribution.')
end
profile = 'multi-mode gaussian';    % Type of beam profile 
dt = 3600;                 % time differential for force calculation [s]
tTotal = 2*7*24*3600;      % total time to run the experiment
rho = 5.12e-19;         % atmospheric density at GEO (35,786km altitude) [kg/m^3]
if PointingTol == 0     % Pointing accuracy of laser sat (assumed)
    tol = .2*pi/180;    % Taken from ITU-R
else
    tol = 0;
end

if plt == 1
    myPlt = figure(1);      % Creating a figure to plot sail anf beam
end
t = 0;                  % initiating time
counter = 1;            % keep track of number of iterations

if PLTs == 1
    velVec = zeros(1,tTotal/dt);
    disVec = velVec;
    tVec = disVec;
    Fvec = tVec;
    offsetVec = Fvec;
end

%% While loop used to update sail posisiton over time
while t <= tTotal % letting code run for a "year". May change to a few weeks

    % Calculating x and y points on the sail for plotting and force
    % calulations
        
    yVec = linspace(center(2)-R,center(2)+R,N);

    xVec = center(1) - sqrt( R^2 - (yVec - center(2)).^2 ); %center(1) - R.*cos(thetaVec);
    xVec2 = center(1) + sqrt( R^2 - (yVec - center(2)).^2 ); %center(1) + R.*cos(thetaVec);
    % xVec2 is a dummy vector used only for plotting since the force is
    % only on the left side of the sail
    if ~isreal(yVec)
        fprintf('yVec isimaginary!\n')
    end
    if ~isreal(xVec)
        fprintf('xVec isimaginary!\n')
    end
    
    % Plotting most recent sail position
    if plt == 1
        myPlt; hold off
        plot([xVec xVec2(end:-1:1)],[yVec yVec(end:-1:1)],'k')%, xVec2, yVec, 'k'),
        axis equal
        xlim([center(1)-4*R, center(1)+4*R])
        ylim([center(2)-2*R, center(2)+2*R]) 
        xlabel('Distance from Source [m]')
        ylabel('Distance from Center [m]')
    end

%% Calculating total force and plotting   
    
%{  
So we have a seperate scripts that calculate different forces (drag and
from the beam) that we can call here. That way we make sure things are
kept seperate and can run by themselves first. Just so nothing gets too
hairy.
%}  

    FBeam = beamforce(R,P,lambda,profile,tol,xVec,yVec,center(1), ...
        BeamDivergence, plt);
    if Drag == 1
        FDrag = dragforce(rho,norm(v),R);
    else
        FDrag = 0;
    end
    F = FBeam + FDrag;
    
    if ~isreal(F)
        fprintf('F isimaginary!\n')
    end
    
    if plt == 1
        % Plotting total force vector at the sails "center of mass"
        hold on, quiver(center(1), center(2), F(1), F(2), ...
            R/(norm(F)), 'b', 'LineWidth', 2)
        hold on, plot(center(1),center(2), 'b*', 'LineWidth', 5)
        drawnow
    end
    
    % Calculating the accelleration due to the force
    a = F./m;
    % Updating the velocity due to the acceleration over the time
    % differential
    v = v + a.*dt;
    % Updating the position due to the velocity over the time differential
    center = center + v.*dt;
    
    if PLTs == 1
        velVec(counter) = v(1);
        disVec(counter) = center(1);
        tVec(counter) = t;
        Fvec(counter) = norm(F);
        offsetVec(counter) = center(2);
    end
        
    % A bunch of print statements just to see what is going on
    fprintf('\nt = %i:%i:%.2f\n', floor(t/3600),...
        floor((t-(floor(t/3600)*3600))/60), ...
        floor( (t-(floor(t/3600)*3600)-(floor((t-(floor(t/3600)*3600))/60)*60))))
    
    if stmnt == 1
        fprintf('F = %e N\n', F)
        fprintf('a = %e m/s^2\n', a)
        fprintf('v = %f m/s\n', v)
        fprintf('x = %f m\n', center)
    end
    
    %% Checking if the sail falls off of the beam and breaks if so
    if center(2)<-R || center(2)>R
        fprintf('\nSail has fallen out of potential well.\n')
        fprintf('Unstable at %ihrs %imins and %.2fs\n', floor(t/3600),...
        floor((t-(floor(t/3600)*3600))/60), ...
        floor( (t-(floor(t/3600)*3600)-(floor((t-(floor(t/3600)*3600))/60)*60))))
        break
    end
    
%     pause()
    t = t + dt; % Updating time so the while loop ends at some point.
    counter = counter + 1; % Updating interation number
end

if PLTs == 1
    figure, plot(tVec,disVec)
    xlabel('Time [s]'), ylabel('Distance [m]')
    figure, plot(tVec,velVec)
    xlabel('Time [s]'), ylabel('Velocity [m/s]')
    figure, plot(tVec,Fvec)
    xlabel('Time [s]'), ylabel('Force [N]')
    figure, plot(disVec, Fvec)
    xlabel('Distance [m]'), ylabel('Force [N]')
    figure, plot(tVec, offsetVec)
    xlabel('Time [s]'), ylabel('Sail Offset from Beam Axis [m]')
end
