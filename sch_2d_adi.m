% Function that computes the numerical solution of the 2D Schrodinger
% equation following the ADI scheme
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
%   y: Vector of y coordinates [ny]
%   t: Vector of time [nt]
%   psi: Array of computed psi values [nt x nx x ny]
%   psire: Array of computed psi real values [nt x nx x ny]
%   psim: Array of computed psi imaginary values [nt x nx x ny]
%   psimod: Array of computed sqrt(psi psi*) values [nt x nx x ny]
%   v: Array of potential values [nx x ny]
function [x, y, t, psi, psire, psiim, psimod, v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar)
    % Initialize x,y and t vectors
    nx = 2^level + 1;
    ny = nx;
    x = linspace(0, 1, nx); 
    y = x;
    dx = x(2) - x(1);
    dy = dx;
    dt = lambda*dx;
    nt = round(tmax / dt) + 1;
    t = (0 : nt-1) * dt;

    % Initialize solution and initial data
    psi = zeros(nt, ny, nx);
    if idtype == 0
        mx = idpar(1);
        my = idpar(2);
        psi(1, :, :) = sin(mx*pi*x).'*sin(my*pi*y);
    elseif idtype == 1
        x0 = idpar(1);
        y0 = idpar(2);
        deltax = idpar(3);
        deltay = idpar(4);
        px = idpar(5);
        py = idpar(6);
        psi(1, :, :) = exp(1i*px*x).'*exp(1i*py*y).*exp(-(((x-x0).'/deltax).^2+((y-y0)/deltay).^2));
    else
      fprintf('sch_2d_cn: Invalid idtype=%d\n', idtype);
      return
    end
    % Fix boundary conditions
    psi(1, 1, :) = zeros(nx,1); 
    psi(1, :, 1) = zeros(1,ny);
    psi(1, end, :) = zeros(nx, 1); 
    psi(1, :, end) = zeros(1, ny);
    
    % Initialize potential 
    if vtype == 0 
        v = zeros(nx);
    elseif vtype == 1
        v = zeros(nx);
        xmin = vpar(1);
        xmax = vpar(2);
        ymin = vpar(3);
        ymax = vpar(4);
        vc = vpar(5);
        for i = 1:length(x)
            if x(i)>xmin && x(i)<xmax
                for j = 1:length(y)
                    if y(j)>ymin && y(j)<ymax
                        v(i, j) = vc;
                    end
                end
            end
        end
    elseif vtype == 2
        v = zeros(nx);
        x1 = vpar(1);
        x2 = vpar(2);
        x3 = vpar(3);
        x4 = vpar(4);
        vc = vpar(5);

        j_prime = round((ny - 1)/4) + 1;

        for i = 1:length(x)
            if (x1 >= x2) || (x3 >= x4)
                fprintf('sch_2d_cn: Invalid vtype=%d\n', vpar);
                return
            end
            if (x1 <= x(i) && x(i) <= x2) || (x3 <= x(i) && x(i) <= x4) 
                v(i, j_prime) = 0;
            else
                v(i, j_prime) = vc;
                v(i, j_prime + 1) = vc;
            end
        end
    else
      fprintf('sch_2d_cn: Invalid vdtype=%d\n', vartype);
      return
    end

    % % SOLVE ADI
    % Set up tridiagonal vectors for step 1
    dl = -ones(nx,1)*(1i*dt/(2*dx^2));
    du = dl;
    d = ones(nx,1)*(1 + 1i*dt/(dx^2));
    % Fix boundary conditions
    d(1) = 1.0;
    du(2) = 0.0;
    dl(nx-1) = 0.0;
    d(nx) = 1.0;
    LHS1 = spdiags([dl d du], -1:1, nx, nx);

    % Iterate through to solve the system
    for n = 1 : nt-1
        % Initialize halfstep placeholder      
        f1 = zeros(nx,nx);
        psi_half = zeros(nx);

        % Loop through y values
        psi_temp = zeros(nx, ny);
           
        psi_temp(2:nx-1, 2:ny-1) = reshape(psi(n, 2:nx-1, 2:ny-1),[nx-2, ny-2]).*(1 - 1i*dt/(dy^2) - 1i*dt*v(2:nx-1, 2:ny-1)/2)...
                               + (reshape(psi(n, 2:nx-1, 3:ny ), [nx-2, ny-2])+ reshape(psi(n, 2:nx-1, 1:ny-2), [nx-2, ny-2]))*(1i*dt/(2*dy^2)) ;
        f1(2:nx-1, 2:ny-1) = psi_temp(2:nx-1, 2:ny-1)* (1 - 1i*dt/(dx^2)) + (psi_temp(3:nx, 2:ny-1) + psi_temp(1:nx-2, 2:ny-1))*(1i*dt/(2*dx^2));

        psi_half = LHS1 \ f1;

        for k = 2:nx-1
            % Set up tridiagonal vectors for step 2
            dl = -ones(nx,1)*(1i*dt/(2*dy^2));
            du = dl;
            d = (1 + 1i*dt/(dy^2) + 1i*dt*v(k, :).'/2);
            % Fix boundary conditions
            d(1) = 1.0;
            du(2) = 0.0;
            dl(nx-1) = 0.0;
            d(nx) = 1.0;
            LHS2 = spdiags([dl d du], -1:1, nx, nx);
            
            psi(n+1,k, :)  = LHS2 \ psi_half(k,:).';
        end 
    end
    psire = real(psi);
    psiim = imag(psi); 
    psimod = sqrt(psi.*conj(psi));
end