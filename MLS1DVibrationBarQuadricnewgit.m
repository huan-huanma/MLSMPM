%--------------------------------------------------------------------------
% This script implements the Moving Least Square Material Point Method (MLS-MPM)
% for a 1D elastic bar problem.
% Author: [Huanhuan Ma]
% Date: [2024.04]
%--------------------------------------------------------------------------
clear all;  % Clear workspace and variables

%--------------------------------------------------------------------------
% Material and simulation parameters
%--------------------------------------------------------------------------
E   = 100;        % Young's modulus
L   = 1;          % Length of the bar
rho = 1;          % Density
v0  = 0.1;        % Initial velocity amplitude
n   = 1;          % Mode number

% Compute wave speed and frequency
c      = sqrt(E/rho);       
beta1  = (2*n-1)*0.5*(pi/L);
omega1 = beta1*c;

% MLS and MPM parameters
alpha         = 1;     % Blending parameter: 1 (FLIP), 0 (PIC)
doubleMapping = 0;     % 1: standard particle update, 0: modified update
beta          = 1;     % Integration scheme parameter

%--------------------------------------------------------------------------
% Loop over mesh refinement (single case here)
%--------------------------------------------------------------------------
for ii = 6:6
    % Number of elements
    nn      = [5 10 20 30 40 60 80 160 320 640 1280];
    ne      = nn(ii);
    fprintf('Number of elements: %d \n', ne);

    % Build 1D computational grid
    mesh    = buildGrid1D(L, ne, 0);
    node    = mesh.node;
    numnode = length(node);

    % MLS shape function parameters
    shape = 'circle';       % Domain of influence shape
    dmax  = 2.5;            % Domain radius scaling factor
    
    form  = 'quartic_spline' ;
    %form  = 'cubic_spline' ;
    % form  = 'exp_spline';% Quartic spline weight function
    
    % Mesh information
    nodeCount = mesh.nodeCount;
    nodes     = mesh.node;
    deltax    = mesh.deltax;

    %--------------------------------------------------------------------------
    % Initialize material points
    %--------------------------------------------------------------------------
    ppc       = 2; % Particles per cell
    particles = buildParticlesGeom(mesh, ppc, rho); 

    % Extract particle properties
    xp  = particles.xp;      % Position
    vp  = particles.vp;      % Velocity
    Vp  = particles.Vp;      % Volume
    Vp0 = particles.Vp0;     % Initial volume
    Fp  = particles.Fp;      % Deformation gradient
    s   = particles.s;       % Stress
    eps = particles.eps;     % Strain
    Mp  = particles.Mp;      % Mass

    %--------------------------------------------------------------------------
    % Node domain of influence
    %--------------------------------------------------------------------------
    delta = mesh.deltax;
    di    = ones(length(node), 1) * dmax * delta;  % Radius for each node

    % Initialize particle velocities with a sine profile
    for p = 1:particles.pCount
        vp(p) = v0 * sin(beta1 * xp(p));
    end
    vp0 = vp;

    %--------------------------------------------------------------------------
    % Time integration parameters
    %--------------------------------------------------------------------------
    dtime = 0.2 * deltax / c;     % Time step
    time  = (2*pi/omega1) * 5;    % Final simulation time
    t     = 0;                    % Initial time

    % Storage arrays for analysis
    ta = [];     % Time history
    va = [];     % Center of mass velocity history
    ka = [];     % Kinetic energy
    sa = [];     % Strain energy

    % Initialize nodal quantities
    nmass      = zeros(nodeCount,1);
    nmomentum0 = zeros(nodeCount,1);
    nmomentum  = zeros(nodeCount,1);
    niforce    = zeros(nodeCount,1);
    nvelo      = zeros(nodeCount,1);

    %--------------------------------------------------------------------------
    % Main time-stepping loop
    %--------------------------------------------------------------------------
    while ( t < time )
        % Reset grid data for each time step
        nmass(:)      = 0;
        nmomentum0(:) = 0;
        nvelo(:)      = 0;
        niforce(:)    = 0;

        %-----------------------------------------
        % Particle-to-grid (P2G) projection
        %-----------------------------------------
        for i = 1:length(xp)
            pt    = xp(i);  % Particle position
            index = defineSupport(node, pt, di);  % Find nodes in support domain

            % Compute MLS shape functions and derivatives
            [phi, dphidx] = mlsQuadricBasis1D(pt, index, node, di, form);

            % Project particle mass, momentum, and stress to nodes
            mp    = Mp(i);
            vol   = Vp(i);
            velo  = vp(i);
            sigma = s(i);

            for j = 1:length(index)
                en = index(j);
                nmass(en)      = nmass(en) + mp * phi(j);
                nmomentum0(en) = nmomentum0(en) + mp * velo * phi(j);
                niforce(en)    = niforce(en) - vol * sigma * dphidx(j);
            end
        end

        % Update nodal momenta
        nmomentum = nmomentum0 + niforce * dtime;

        % Apply Dirichlet boundary condition at left end
        nmomentum0(1) = 0;
        nmomentum(1)  = 0;

        %-----------------------------------------
        % Grid-to-particle (G2P) update
        %-----------------------------------------
        k = 0;  % Kinetic energy
        u = 0;  % Strain energy
        for i = 1:length(xp)
            pt    = xp(i);
            index = defineSupport(node, pt, di);
            [phi, dphidx] = mlsQuadricBasis1D(pt, index, node, di, form);

            % Compute nodal velocities
            v10 = nmomentum0(index) ./ nmass(index);
            v1  = nmomentum(index) ./ nmass(index);

            % Update particle velocity
            if doubleMapping == 0
                vp(i) = alpha * vp(i) + phi' * (v1 - alpha * v10);
                Lp    = dphidx' * v1; % Velocity gradient
            else
                Lp = dphidx' * nvelo(index);
            end

            % Update particle position
            xp(i) = xp(i) + (1 - beta) * dtime * vp(i) + beta * dtime * (phi' * v1);

            % Update deformation gradient and volume
            Fp(i) = (1 + Lp * dtime) * Fp(i);
            Vp(i) = Fp(i) * Vp0(i);

            % Update stress and strain
            dEps   = dtime * Lp;
            s(i)   = s(i) + E * dEps;
            eps(i) = eps(i) + dEps;

            % Compute kinetic and strain energy
            k = k + 0.5 * vp(i)^2 * Mp(i);
            u = u + 0.5 * s(i) * eps(i) * Vp(i);
        end

        % Advance time
        t = t + dtime;

        % Store center of mass velocity and energies
        cv = 1 / sum(Mp) * dot(Mp, vp);
        ta = [ta; t];
        va = [va; cv];
        ka = [ka; k];
        sa = [sa; u];
    end

    %--------------------------------------------------------------------------
    % Post-processing: error analysis and plots
    %--------------------------------------------------------------------------
    vExact = v0 / (beta1 * L) * cos(omega1 * ta); % Analytical center of mass velocity

    rel_errorL2 = norm(va - vExact) / norm(vExact);
    fprintf('Relative L2 error = %d \n', rel_errorL2);

    % Plot: numerical vs analytical velocity
    figure;
    hold on;
    plot(ta, va, 'LineWidth', 2.2);
    plot(ta, vExact, 'LineWidth', 2.2, 'LineStyle', '--');
    xlabel('Time', 'FontSize', 18);
    ylabel('Velocity', 'FontSize', 18);
    legend('MLSMPM', 'Exact', 'FontSize', 15);
    grid on;
    set(gca, 'FontSize', 15);

    % Plot: absolute error
    figure;
    hold on;
    plot(ta, abs(va - vExact), 'LineWidth', 1.6);
    xlabel('Time', 'FontSize', 18);
    ylabel('$|v - v_h|$', 'Interpreter', 'latex', 'FontSize', 18);
    legend('Error', 'FontSize', 15);
    grid on;
    set(gca, 'FontSize', 15);

    % Plot: kinetic, strain, and total energy
    figure;
    hold on;
    plot(ta, ka, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.6); 
    plot(ta, sa, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.6, 'LineStyle', '--');
    plot(ta, ka + sa, 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.6);
    xlabel('Time', 'FontSize', 18);
    ylabel('Energy', 'FontSize', 18);
    legend('Kinetic', 'Strain', 'Total', 'FontSize', 15);
    grid on;
    set(gca, 'FontSize', 15);
end
