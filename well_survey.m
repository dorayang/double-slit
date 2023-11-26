% Script to run the well survey numerical analysis

% Parameters 
tmax = 0.1;
level = 9;
lambda = 0.01;
idtype = 1;
idpar = [0.4, 0.075, 0];
vtype = 1;
lnv = linspace(2, 10, 251);
v_well = -exp(lnv);
Fe_well = [];
tic
% Iterate through
for v = v_well
    vpar = [0.6, 0.8, v];

    [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
    i1 = round(length(x)*6/10);
    i2 = round(length(x)*8/10);

    % Normalizing
    prob = prob/mean(prob(:,end));

    x1 = x(i1);
    x2 = x(i2);
    p1 = mean(prob(:,i1));
    p2 = mean(prob(:,i2));
    Fe_well = [Fe_well, (p2 - p1)/(x2-x1)];
end

plot(lnv, log(Fe_well), 'LineWidth', 1.5  );
title("Well Survery")
xlabel("ln(V_0)")
ylabel("ln(Fe(0.6, 0.8))")
toc
