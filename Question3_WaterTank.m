%% QUESTION 3: MINIMUM-SURFACE-AREA CLOSED CYLINDRICAL WATER TANK
% Find tank radius r and height h for a fixed capacity of 100 m^3.
% Eliminate h using the volume equality, so all solvers optimize r only.
clear; clc; close all;

V = 100;        % required capacity (m^3)
rMin = 2.0; rMax = 4.0;
hMin = 2.0; hMax = 10.0;

% h = V/(pi*r^2); derive radius bounds that also satisfy height bounds
lb = max(rMin,sqrt(V/(pi*hMax)));
ub = min(rMax,sqrt(V/(pi*hMin)));
x0 = 2.8;
rawfun = @(x) tankArea(x,V);

problem.name = 'Closed cylindrical water tank';
problem.fun = rawfun;
problem.rawfun = rawfun;
problem.violation = @(x) [];
problem.lb = lb; problem.ub = ub; problem.x0 = x0;
problem.labels = {'r (m)'};

Results = runSixSolvers(problem);

% Recover dimensions for each algorithm and print them
fprintf('\nTank dimensions from each solver:\n');
for k=1:height(Results)
    r = Results.X{k}(1);
    h = V/(pi*r^2);
    fprintf('%s: r = %.5f m, h = %.5f m, area = %.5f m^2\n', ...
        Results.Algorithm(k),r,h,Results.ObjectiveValue(k));
end

% Plot reduced objective function versus radius
rv = linspace(lb,ub,300);
Av = arrayfun(@(r) tankArea(r,V),rv);
figure('Color','w');
plot(rv,Av,'LineWidth',2);
xlabel('Tank radius r (m)'); ylabel('Closed-tank surface area (m^2)');
title('Water tank: surface area versus tank radius');
grid on;
exportgraphics(gcf,'WaterTank_design_curve.png','Resolution',300);

function A = tankArea(x,V)
r=x(1);
h=V/(pi*r^2);
A=2*pi*r*h+2*pi*r^2; % curved side plus top and bottom
end
