% Script to run the barrier survey numerical analysis

% Parameters 
tmax = 0.1;
level = 9;
lambda = 0.01;
idtype = 1;
idpar = [0.4, 0.075, 20.0];
vtype = 1;
lnv = linspace(-2,5, 251);
v_barrier = exp(lnv);

Fe = [];
tic
% Iterate through
for v = v_barrier
    vpar = [0.6, 0.8, v];

    [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
    i = round(length(x)*8/10);

    % Normalizing
    prob = prob/mean(prob(:,end));

    x1 = x(i);
    x2 = x(end);
    p1 = mean(prob(:,i));
    p2 = mean(prob(:,end));
    Fe = [Fe, (p2 - p1)/(x2-x1)];
end

plot(lnv, log(Fe), 'LineWidth', 2);
title("Barrier Survery")
xlabel("ln(V_0)")
ylabel("ln(Fe(0.8, 1.0))")
toc
