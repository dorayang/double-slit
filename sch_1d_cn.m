% Function that computes the numerical solution for the 1D Schrodinger
% equation following the Crank-Nicholson scheme
% Inputs
%   tmax: Maximum integration time
%   level: Discretization level
%   lamda: dt/dx
%   idtype: Selects initial data type
%   idpar: Vector of initial data parameters
%   vtype: Selects potential type
%   vpar: Vector of potential parameters
% Outputs
%   x: Vector of x coordinates [nx]
%   t: Vector of time [nt]
%   psi: Array of computed psi values [nt x nx]
%   psire: Array of computed psi real values [nt x nx]
%   psim: Array of computed psi imaginary values [nt x nx]
%   psimod: Array of computed sqrt(psi psi*) values [nt x nx]
%   prob: Array of computed running integral values [nt x nx]
%   v: Array of potential values [nx]

function [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar)
    % Initialize t and x vectors
    nx = 2^level + 1;
    x = linspace(0, 1, nx);
    dx = x(2) - x(1);
    dt = lambda*dx;
    nt = tmax / dt + 1;
    t = (0 : nt-1) * dt;
    
    % Initialize solution and initial data
    psi = zeros(nt, nx);
    if idtype == 0
        m = idpar(1);
        psi(1, :) = sin(m*pi*x);
    elseif idtype == 1
        x0 = idpar(1);
        delta = idpar(2);
        p = idpar(3);
        psi(1, :) = exp(1i*p*x).*exp(-((x-x0)./delta).^2);           
    else
      fprintf('sch_1d_cn: Invalid idtype=%d\n', idtype);
      return
    end
    % Fix boundary conditions
        psi(1,1) = 0;
        psi(1,end) = 0;

    % Initialize potential 
    if vtype == 0 
        v = zeros(nx,1);
    elseif vtype == 1
        v = zeros(nx,1);
        xmin = vpar(1);
        xmax = vpar(2);
        vc = vpar(3);
        for i = 1:length(x)
            if x(i)>xmin && x(i)<xmax
               v(i) = vc;
            end
        end
    else
      fprintf('sch_1d_cn: Invalid idtype=%d\n', vartype);
      return
    end

    % Set up triadigonal vectors
    dl = ones(nx,1)*1/(2*dx^2);
    du = dl;
    d = (1i/dt - 1/dx^2 - 1/2*v);

    % Fix boundary conditions
    d(1) = 1.0;
    du(2) = 0.0;
    dl(nx-1) = 0.0;
    d(nx) = 1.0;
    LHS = spdiags([dl d du], -1:1, nx, nx);
    
    % Initializing RHS and running integral arrays
    prob = zeros(nt,nx);
    f = zeros(nx,1);

    % Iterate through to solve the system
    for n = 1 : nt-1
        % Define RHS of linear system 
        f(2:nx-1) = -1/(2*dx^2).*psi(n, 3:nx) + (1i/dt + 1/dx^2 + 1/2*v(2:nx-1)').*psi(n,2:nx-1) -1/(2*dx^2).*psi(n, 1:nx-2);
        f(1) = 0.0;
        f(nx) = 0.0;

        % Solve system = updating approximation to next time 
        psi(n+1, :) = LHS \ f;       
    end

    % Computing running integral
    for j = 2:nx
        prob(:,j) = prob(:, j-1) + 1/2*(psi(:,j).*conj(psi(:,j)) + psi(:,j-1).*conj(psi(:,j-1)))*(x(j) - x(j-1));       
    end

    % Define other outputs
    psire = real(psi);
    psiim = imag(psi);
    psimod = sqrt(psi.*conj(psi));
end