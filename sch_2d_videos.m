
% Gaussian initial condition with double slit
tmax = 0.003;
level = 9;
lambda = 0.01;
idtype = 1;
idpar = [0.5, 0.15, 0.5, 0.01, 0, 100];
vtype = 2;
slit_distance = 0.15;
slit_size = 0.02;
x2 = 0.5 - slit_distance/2;
x1 = x2 - slit_size;
x3 = 0.5 + slit_distance/2;
x4 = x3 +slit_size;
vpar = [x1, x2, x3, x4, 10^8];
vid_length = 10*4; 

% [x_v8, y_v8, t_v8, psi_v8, psire_v8, psiim_v8, psimod_v8, v_v8] = sch_2d_adi(tmax, level, lambda, 1, idpar, vtype, vpar);
video_generator("Double_slit_diffraction_pattern.avi", x_v8, y_v8, t_v8, psimod_v8, round(length(t_v8)/vid_length), "\surd{(\psi \psi^*)}", "Boosted Gaussian Double Slit Diffraction Pattern", 2, v_v8);


% Gaussian initial condition with barrier potential
tmax = 0.003;
level = 9;
lambda = 0.01;
idtype = 1;
idpar = [0.3, 0.3, 0.08, 0.08, 100, 100];
vtype = 1;
vpar = [0.5, 0.7, 0.5, 0.7, 10^8];
vid_length = 10*4; 
length(t_v6)

[x_v6, y_v6, t_v6, psi_v6, psire_v6, psiim_v6, psimod_v6, v_v6] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);

video_generator("Gaussian_barrier.avi", x_v6, y_v6, t_v6, psimod_v6, round(length(t_v6)/vid_length), "\surd{(\psi \psi^*)}", "Boosted Gaussian, Rectangular Barrier Potential", 2, v_v6);

% % Gaussian initial condition with well potential
tmax = 0.003;
level = 9;
lambda = 0.01;
idtype = 1;
idpar = [0.5, 0.7, 0.05, 0.05, 0, 100];
vtype = 1;
vpar = [0.3, 0.7, 0.8, 0.99, -10^8]; %-10^5 = bound state?
vid_length = 10*2; 

[x_v7, y_v7, t_v7, psi_v7, psire_v7, psiim_v7, psimod_v7, v_v7] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
video_generator("Gaussian_well_medium.avi", x_v7, y_v7, t_v7, psimod_v7, round(length(t_v7)/vid_length), "\surd{(\psi \psi^*)}", "Boosted Gaussian, Rectangular Well Potential (V = -10^5)", 2, v_v7);


