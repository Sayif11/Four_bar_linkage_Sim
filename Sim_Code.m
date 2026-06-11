%% Four-Bar Linkage Simulation
clear;
clc;
close all;

%% Parameters

p.A = 1.0;      % Ground link height
p.C = 1.0;      % Coupler length
p.r = 0.5;      % Input crank
p.R = 0.8;      % Output rocker

p.torque = 0;

p.Ir = 60;
p.IR = 10;

%% Initial Conditions

theta0 = 0.0;
omega0 = 0.2;

%% Time

tFinal = 100;
nSteps = 400;

tspan = linspace(0,tFinal,nSteps);

%% Solve Dynamics

y0 = [theta0; omega0];

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

[t,y] = ode45(@(t,y) dynamics(t,y,p), tspan, y0, opts);

theta = y(:,1);
omega = y(:,2);

%% Angular Velocity Plot

figure('Name','Angular Velocity');

plot(t,omega,'LineWidth',1.5)

xlabel('Time')
ylabel('Angular Velocity (rad/s)')
title('Input Crank Angular Velocity')
grid on

%% Animation Figure

figure('Name','Four Bar Linkage');

axis equal
axis([-1.5 1.5 -0.6 1.6])

xlabel('X')
ylabel('Y')
title('Four-Bar Linkage')

hold on
grid on

origin = [0 0];
Apoint = [0 p.A];

% Initial geometry
B = point_r(theta(1),p);
D = point_R(theta(1),p);

% Ground link
hA = plot([origin(1) Apoint(1)], ...
          [origin(2) Apoint(2)], ...
          'b','LineWidth',3);

% Input crank
hr = plot([origin(1) B(1)], ...
          [origin(2) B(2)], ...
          'r','LineWidth',3);

% Output rocker
hR = plot([Apoint(1) D(1)], ...
          [Apoint(2) D(2)], ...
          'Color',[0.4 0.4 0.4], ...
          'LineWidth',3);

% Coupler
hC = plot([B(1) D(1)], ...
          [B(2) D(2)], ...
          'm','LineWidth',3);

% Joint markers
joints = scatter( ...
    [origin(1) Apoint(1) B(1) D(1)], ...
    [origin(2) Apoint(2) B(2) D(2)], ...
    100,'filled','y');

%% Animate

for k = 1:length(theta)

    th = theta(k);

    if ~isValidConfiguration(th,p)
        continue
    end

    B = point_r(th,p);
    D = point_R(th,p);

    set(hr,...
        'XData',[origin(1) B(1)],...
        'YData',[origin(2) B(2)]);

    set(hR,...
        'XData',[Apoint(1) D(1)],...
        'YData',[Apoint(2) D(2)]);

    set(hC,...
        'XData',[B(1) D(1)],...
        'YData',[B(2) D(2)]);

    joints.XData = [origin(1) Apoint(1) B(1) D(1)];
    joints.YData = [origin(2) Apoint(2) B(2) D(2)];

    drawnow
end

%% ============================================================
%% Local Functions
%% ============================================================

function dydt = dynamics(~,y,p)

theta = y(1);
omega = y(2);

if ~isValidConfiguration(theta,p)
    dydt = [0;0];
    return
end

dphi  = outputDerivative(theta,p);
ddphi = outputSecondDerivative(theta,p);

ratio = mechanicalAdvantage(theta,p);

alpha = ...
    (p.torque ...
    - p.IR*omega^2*ddphi*ratio) ...
    /(p.Ir + p.IR*dphi*ratio);

dydt = [omega; alpha];

end

% ------------------------------------------------------------

function s = sineInput(theta,p)

h = sqrt( ...
    p.A^2 + p.r^2 ...
    - 2*p.A*p.r*sin(theta));

s = (p.R^2 - p.C^2 + h^2)/(2*p.R*h);

end

% ------------------------------------------------------------

function valid = isValidConfiguration(theta,p)

valid = abs(sineInput(theta,p)) <= 1;

end

% ------------------------------------------------------------

function phi = outputAngle(theta,p)

s = sineInput(theta,p);

if abs(s) > 1
    error('Mechanism locked or impossible configuration.')
end

alpha = asin(s);

beta = atan2( ...
    p.r*cos(theta), ...
    p.A - p.r*sin(theta));

phi = beta - alpha;

end

% ------------------------------------------------------------

function P = point_r(theta,p)

P = [ ...
    p.r*cos(theta), ...
    p.r*sin(theta)];

end

% ------------------------------------------------------------

function P = point_R(theta,p)

phi = outputAngle(theta,p);

P = [ ...
    p.R*cos(phi), ...
    p.A + p.R*sin(phi)];

end

% ------------------------------------------------------------

function ratio = mechanicalAdvantage(theta,p)

phi = outputAngle(theta,p);

f = p.r*p.R*sin(phi-theta);

num = p.A*p.r*cos(theta) + f;
den = p.A*p.R*cos(phi)   + f;

ratio = num/den;

end

% ------------------------------------------------------------

function d = outputDerivative(theta,p)

h = 1e-6;

d = ( ...
    outputAngle(theta+h,p) ...
   -outputAngle(theta-h,p)) ...
   /(2*h);

end

% ------------------------------------------------------------

function d2 = outputSecondDerivative(theta,p)

h = 1e-6;

d2 = ( ...
     outputAngle(theta+h,p) ...
    -2*outputAngle(theta,p) ...
    +outputAngle(theta-h,p)) ...
    /(h^2);

end
