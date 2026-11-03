function Results = runSixSolvers(problem)
% RUNSIXSOLVERS Compare six MATLAB optimization algorithms on one problem.
% problem fields:
%   name, fun (penalized objective), rawfun (engineering objective),
%   violation (returns g(x)<=0), lb, ub, x0, labels
%
% Required toolboxes:
%   Optimization Toolbox: fmincon
%   Global Optimization Toolbox: ga, patternsearch, surrogateopt,
%       simulannealbnd, particleswarm

solvers = {'ga','fmincon','patternsearch','surrogateopt', ...
           'simulannealbnd','particleswarm'};
nRepeat = 3;
nSolver = numel(solvers);
times = nan(nSolver,1);
fvals = nan(nSolver,1);
violations = nan(nSolver,1);
solutions = cell(nSolver,1);
exitflags = nan(nSolver,1);

for k = 1:nSolver
    runTimes = nan(nRepeat,1);
    runX = cell(nRepeat,1);
    runF = nan(nRepeat,1);
    runV = nan(nRepeat,1);
    runExit = nan(nRepeat,1);

    for rep = 1:nRepeat
        rng(100*k + rep,'twister');
        n = numel(problem.lb);
        lb = problem.lb(:)';
        ub = problem.ub(:)';
        x0 = problem.x0(:)';
        objective = problem.fun;

        t0 = tic;
        switch solvers{k}
            case 'ga'
                opts = optimoptions('ga','Display','off', ...
                    'PopulationSize',40,'MaxGenerations',80);
                [x,~,ef] = ga(objective,n,[],[],[],[],lb,ub,[],opts);

            case 'fmincon'
                opts = optimoptions('fmincon','Display','off', ...
                    'Algorithm','sqp','MaxIterations',500, ...
                    'MaxFunctionEvaluations',5000);
                [x,~,ef] = fmincon(objective,x0,[],[],[],[],lb,ub,[],opts);

            case 'patternsearch'
                opts = optimoptions('patternsearch','Display','off', ...
                    'MaxIterations',500,'MaxFunctionEvaluations',5000);
                [x,~,ef] = patternsearch(objective,x0,[],[],[],[],lb,ub,[],opts);

            case 'surrogateopt'
                opts = optimoptions('surrogateopt','Display','off', ...
                    'MaxFunctionEvaluations',250);
                [x,~,ef] = surrogateopt(objective,lb,ub,[],[],[],[],[],opts);

            case 'simulannealbnd'
                opts = optimoptions('simulannealbnd','Display','off', ...
                    'MaxIterations',500,'MaxFunctionEvaluations',5000);
                [x,~,ef] = simulannealbnd(objective,x0,lb,ub,opts);

            case 'particleswarm'
                opts = optimoptions('particleswarm','Display','off', ...
                    'SwarmSize',40,'MaxIterations',100);
                [x,~,ef] = particleswarm(objective,n,lb,ub,opts);
        end
        runTimes(rep) = toc(t0);
        x = x(:)';
        runX{rep} = x;
        runF(rep) = problem.rawfun(x);
        g = problem.violation(x);
        if isempty(g)
            runV(rep) = 0;
        else
            runV(rep) = max([0; g(:)]);
        end
        runExit(rep) = ef;
    end

    % Choose the lowest-objective feasible result. If no repeat is
    % feasible, report the least-violating candidate and flag it.
    feasible = find(runV <= 1e-5);
    if ~isempty(feasible)
        [~,j] = min(runF(feasible));
        chosen = feasible(j);
    else
        [~,chosen] = min(runV);
        warning('%s did not find a feasible result for %s.',solvers{k},problem.name);
    end

    times(k) = median(runTimes);
    fvals(k) = runF(chosen);
    violations(k) = runV(chosen);
    solutions{k} = runX{chosen};
    exitflags(k) = runExit(chosen);
end

Algorithm = string(solvers(:));
MedianTime_s = times;
ObjectiveValue = fvals;
MaxConstraintViolation = violations;
ExitFlag = exitflags;
X = solutions;
Results = table(Algorithm,MedianTime_s,ObjectiveValue, ...
                MaxConstraintViolation,ExitFlag,X);
disp("Results for " + problem.name);
disp(Results);

safeName = regexprep(problem.name,'[^A-Za-z0-9]','_');
writetable(Results,[safeName '_results.csv']);

figure('Color','w');
bar(categorical(Algorithm),MedianTime_s);
ylabel('Median elapsed time (seconds)');
xlabel('Optimization algorithm');
title([problem.name ': computation-time comparison']);
grid on;
exportgraphics(gcf,[safeName '_time_comparison.png'],'Resolution',300);

figure('Color','w');
bar(categorical(Algorithm),ObjectiveValue);
ylabel('Best feasible objective value');
xlabel('Optimization algorithm');
title([problem.name ': objective-value comparison']);
grid on;
exportgraphics(gcf,[safeName '_objective_comparison.png'],'Resolution',300);

save([safeName '_results.mat'],'Results','problem');
end
