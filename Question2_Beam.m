%% QUESTION 2: MINIMUM-VOLUME SIMPLY SUPPORTED RECTANGULAR BEAM
% Select width b and depth h to minimize beam material volume.
clear; clc; close all;

% Given data
P = 200;              % central point load (kN)
L = 8;                % span (m)
E = 200e6;            % Young's modulus (kN/m^2), 200 GPa
sigmaAllow = 150e3;   % allowable bending stress (kN/m^2), 150 MPa
deltaAllow = L/360;   % allowable deflection (m)
rho = 1e4;            % quadratic penalty coefficient

% Design variables x = [b,h]
lb = [0.03 0.25];
ub = [0.15 0.70];
x0 = [0.08 0.50];

rawfun = @(x) L*x(1)*x(2); % beam volume (m^3)
violation = @(x) beamG(x,P,L,E,sigmaAllow,deltaAllow);
fun = @(x) rawfun(x) + rho*sum(max(violation(x),0).^2);

problem.name = 'Simply supported beam';
problem.fun = fun;
problem.rawfun = rawfun;
problem.violation = violation;
problem.lb = lb; problem.ub = ub; problem.x0 = x0;
problem.labels = {'b (m)','h (m)'};

Results = runSixSolvers(problem);

% Feasible design-space contour plot
bv = linspace(lb(1),ub(1),160);
hv = linspace(lb(2),ub(2),160);
[BB,HH] = meshgrid(bv,hv);
Z = L.*BB.*HH;
M = P*L/4;
sigma = 6*M./(BB.*HH.^2);
I = BB.*HH.^3/12;
delta = P*L^3./(48*E.*I);
feasible = (sigma<=sigmaAllow) & (delta<=deltaAllow);

figure('Color','w');
contourf(BB,HH,Z,20,'LineColor','none'); colorbar; hold on;
contour(BB,HH,double(feasible),[0.5 0.5],'k','LineWidth',2);
xlabel('Beam width b (m)'); ylabel('Beam depth h (m)');
title('Beam: material-volume contours and feasible boundary');
legend('Objective volume (m^3)','Feasible-region boundary','Location','best');
grid on;
exportgraphics(gcf,'Beam_design_space.png','Resolution',300);

function g = beamG(x,P,L,E,sigmaAllow,deltaAllow)
b=x(1); h=x(2);
M=P*L/4;
I=b*h^3/12;
sigma=6*M/(b*h^2);
delta=P*L^3/(48*E*I);
g=[sigma/sigmaAllow-1; delta/deltaAllow-1]; % feasible when g <= 0
end
