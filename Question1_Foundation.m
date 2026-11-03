%% QUESTION 1: OPTIMIZATION OF AN ISOLATED RECTANGULAR FOOTING
% Find footing width B and thickness D that minimize concrete volume.
clear; clc; close all;

% Given data
P = 1800;       % column load (kN)
qa = 200;       % allowable bearing pressure (kN/m^2)
Es = 20000;     % soil elastic modulus (kN/m^2)
nu = 0.30;      % Poisson's ratio
If = 1.12;      % simplified settlement influence factor
sAllow = 0.025; % allowable settlement (m)
rho = 1e4;      % quadratic penalty coefficient

% Design variables x = [B,D]; footing length L = 1.2B
lb = [2.5 0.45];
ub = [5.0 0.80];
x0 = [3.5 0.55];

rawfun = @(x) 1.2*x(1)^2*x(2);  % concrete volume (m^3)
violation = @(x) footingG(x,P,qa,Es,nu,If,sAllow);
fun = @(x) rawfun(x) + rho*sum(max(violation(x),0).^2);

problem.name = 'Foundation footing';
problem.fun = fun;
problem.rawfun = rawfun;
problem.violation = violation;
problem.lb = lb; problem.ub = ub; problem.x0 = x0;
problem.labels = {'B (m)','D (m)'};

Results = runSixSolvers(problem);

% Feasible design-space contour plot (B-D plane)
Bv = linspace(lb(1),ub(1),160);
Dv = linspace(lb(2),ub(2),160);
[BB,DD] = meshgrid(Bv,Dv);
Z = 1.2.*BB.^2.*DD;
G1 = P./(1.2.*BB.^2)./qa - 1;
q = P./(1.2.*BB.^2);
settlement = q.*BB.*(1-nu^2).*If./Es;
G2 = settlement./sAllow - 1;
feasible = (G1<=0) & (G2<=0);

figure('Color','w');
contourf(BB,DD,Z,20,'LineColor','none'); colorbar; hold on;
contour(BB,DD,double(feasible),[0.5 0.5],'k','LineWidth',2);
xlabel('Footing width B (m)'); ylabel('Footing thickness D (m)');
title('Foundation: concrete-volume contours and feasible boundary');
legend('Objective volume (m^3)','Feasible-region boundary','Location','best');
grid on;
exportgraphics(gcf,'Foundation_design_space.png','Resolution',300);

function g = footingG(x,P,qa,Es,nu,If,sAllow)
B=x(1);
q=P/(1.2*B^2);
s=q*B*(1-nu^2)*If/Es;
g=[q/qa-1; s/sAllow-1]; % feasible when each g <= 0
end
