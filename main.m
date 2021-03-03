%{
AME 441 Group 5
Passive Stability of a Laser-Driven Light Sail
Spring 2021
%}

clear; clc; close all;

%% Setting up inital paramters and known constants

m = 1;                  % Mass of the sail [kg]
R = 1;                  % Radius of sail [m]
center = [0 0];         % Initial position of sail center
v = [0 0];              % Initial velocity of sail [m/s]
P = 30;                 % Power of laser beam [W]
N = 501;                % Number of rays
profile = 'uniform';   % Type of beam profile 
dt = 1;                 % time differential for force calculation [s]

myPlt = figure(1);      % Creating a figure to plot sail anf beam
t = 0;                  % initiating time

%% While loop used to update sail posisiton over time
while t <= 3.15e7 % letting code run for a "year"

    % Calculating x and y points on the sail for plotting and force
    % calulations
    th = linspace(.5*pi,1.5*pi,N);
    yVec = R*sin(th) + center(2); %linspace(-R,R,N); Ignore these, just ideas
    xVec = R*cos(th) + center(1); %center(1) - sqrt( R^2 - (yVec - center(2)).^2 );
    xVec2 = R*cos(th+pi) + center(1); %center(1) + sqrt( R^2 - (yVec - center(2)).^2 );
    % xVec2 is a dummy vector used only for plotting since the force is
    % only on the left side of the sail
    
    % Plotting most recent sail position
    myPlt; hold off
    plot(xVec,yVec,'k', xVec2, yVec, 'k'), axis equal,
    xlim([center(1)-4*R, center(1)+4*R])
    ylim([center(2)-2*R, center(2)+2*R])    

%% Calculating total force and plotting   
    
%{  
So we have a seperate scripts that calculate different forces (drag and
 from the beam) that we can call here. That way we make sure things are
kept seperate and can run by themselves first. Just so nothing gets too
hairy.
%}  
    FBeam = beamforce(R,P,profile,xVec,yVec);
    FDrag = 0; % dragforce()
    F = FBeam + FDrag;
  
    % Plotting total force vector at the sails "center of mass"
    hold on, quiver(center(1), center(2), F(1), F(2), inv(norm(F)),...
        'b', 'LineWidth', 2)
    hold on, plot(center(1),center(2), 'b*', 'LineWidth', 5)
    drawnow
    
    % Calculating the accelleration due to the force
    a = F./m;
    % Updating the velocity due to the acceleration over the time
    % differential
    v = v + a.*dt;
    % Updating the position due to the velocity over the time differential
    center = center + v.*dt;
    
    % A bunch of print statements just to see what is going on
    fprintf('\nt = %i:%i:%.2f\n', floor(t/3600),...
        floor((t-(floor(t/3600)*3600))/60), ...
        floor( (t-(floor(t/3600)*3600)-(floor((t-(floor(t/3600)*3600))/60)*60))))        
    fprintf('F = %e N\n', F)
    fprintf('a = %e m/s^2\n', a)
    fprintf('v = %f m/s\n', v)
    fprintf('x = %f m\n', center)
    
%     pause()
    t = t + dt; % Updating time so the while loop ends at some point.
end
