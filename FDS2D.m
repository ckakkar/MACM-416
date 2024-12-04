% 2D Heat Equation
% Author: Cyrus Kakkar
% Date: December 4, 2024
% Course: MACM 416

clear; clc; close all;

%% Parameters
alpha = 0.01;          % Thermal diffusivity
Lx = 1; Ly = 1;        % Domain size
resolutions = [50, 100, 200]; % Grid resolutions for convergence study
T_end = 0.1;           % Simulation end time
errors = [];           % Store L2 errors
runtimes = [];         % Store runtimes

%% Loop over different grid resolutions (Test Problem)
for N = resolutions
    % Grid setup
    Nx = N; Ny = N;
    dx = Lx / (Nx - 1);
    dy = Ly / (Ny - 1);
    dt = min(dx^2, dy^2) / (4 * alpha); % Time step (CFL condition)
    x = linspace(0, Lx, Nx);
    y = linspace(0, Ly, Ny);
    [X, Y] = meshgrid(x, y);

    % Initial condition (sinusoidal profile)
    u = sin(pi * X) .* sin(pi * Y);

    % Preallocate for new temperature values
    u_new = zeros(size(u));
    t = 0; % Initialize time

    % Initialize video writer for test problem movie
    v_test = VideoWriter(sprintf('test_problem_%dx%d.avi', Nx, Ny));
    open(v_test);

    % Time-stepping loop
    tic; % Start timer for runtime measurement
    while t < T_end
        % Update interior points using explicit scheme
        for i = 2:Nx-1
            for j = 2:Ny-1
                u_new(i, j) = u(i, j) + alpha * dt * ...
                    ((u(i+1, j) - 2*u(i, j) + u(i-1, j)) / dx^2 + ...
                     (u(i, j+1) - 2*u(i, j) + u(i, j-1)) / dy^2);
            end
        end

        % Update boundary conditions (Dirichlet: u = 0)
        u_new(:, 1) = 0; u_new(:, end) = 0;
        u_new(1, :) = 0; u_new(end, :) = 0;

        % Update time and solution
        u = u_new;
        t = t + dt;

        % Add frame to movie
        surf(X, Y, u, 'EdgeColor', 'none');
        xlabel('x'); ylabel('y'); zlabel('u');
        title(sprintf('Test Problem at t = %.3f', t));
        colorbar;
        axis([0 Lx 0 Ly -0.1 1]);
        frame = getframe(gcf);
        writeVideo(v_test, frame);
    end
    runtime = toc; % End timer for runtime measurement
    close(v_test);

    % Exact solution for validation
    u_exact = sin(pi * X) .* sin(pi * Y) .* exp(-pi^2 * alpha * T_end);

    % Compute L2 error
    error = sqrt(sum((u(:) - u_exact(:)).^2) / numel(u));
    errors = [errors, error];
    runtimes = [runtimes, runtime];

    % Display results for this resolution
    fprintf('Test Problem | Grid: %dx%d | L2 Error: %.5e | Runtime: %.2f s\n', Nx, Ny, error, runtime);
end

% Convergence plot
figure;
loglog(resolutions, errors, '-o', 'LineWidth', 1.5);
xlabel('Grid Resolution (Nx = Ny)');
ylabel('L2 Error');
title('Convergence of Numerical Method (Test Problem)');
grid on;

% Runtime plot
figure;
plot(resolutions, runtimes, '-o', 'LineWidth', 1.5);
xlabel('Grid Resolution (Nx = Ny)');
ylabel('Runtime (s)');
title('Algorithm Runtime vs Resolution (Test Problem)');
grid on;

% Visualization of Numerical and Exact Solutions (for finest grid)
figure;
subplot(1, 2, 1);
surf(X, Y, u, 'EdgeColor', 'none');
xlabel('x'); ylabel('y'); zlabel('u');
title('Numerical Solution (Test Problem)');
colorbar;

subplot(1, 2, 2);
surf(X, Y, u_exact, 'EdgeColor', 'none');
xlabel('x'); ylabel('y'); zlabel('u');
title('Exact Solution (Test Problem)');
colorbar;

% Error plot
figure;
error_plot = abs(u - u_exact); % Compute absolute error
surf(X, Y, error_plot, 'EdgeColor', 'none'); % Plot error
xlabel('x'); ylabel('y'); zlabel('|Error|');
title('Error (Numerical vs Exact Solution)');
colorbar;

%% Hard Problem Simulation
% Modify the parameters for the hard problem
f = @(x, y, t) 0.1 * sin(2 * pi * t) .* sin(pi * x) .* sin(pi * y); % Source term
alpha_hard = @(x, y) 0.01 + 0.005 * sin(pi * x) .* sin(pi * y);    % Non-uniform diffusivity

% Hard problem grid (finest resolution)
Nx = 200; Ny = 200;
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
dt = min(dx^2, dy^2) / (4 * 0.015); % Adjusted for max diffusivity
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);
u = zeros(size(X)); % Initial condition (uniform temperature)

% Initialize video writer for hard problem movie
v_hard = VideoWriter('hard_problem.avi');
open(v_hard);

t = 0;
while t < T_end
    % Update interior points with source term and non-uniform diffusivity
    for i = 2:Nx-1
        for j = 2:Ny-1
            alpha_loc = alpha_hard(X(i, j), Y(i, j));
            u_new(i, j) = u(i, j) + dt * ( ...
                alpha_loc * ((u(i+1, j) - 2*u(i, j) + u(i-1, j)) / dx^2 + ...
                             (u(i, j+1) - 2*u(i, j) + u(i, j-1)) / dy^2) + ...
                f(X(i, j), Y(i, j), t));
        end
    end

    % Update boundary conditions (Dirichlet: u = 0)
    u_new(:, 1) = 0; u_new(:, end) = 0;
    u_new(1, :) = 0; u_new(end, :) = 0;

    % Update time and solution
    u = u_new;
    t = t + dt;

    % Add frame to movie
    surf(X, Y, u, 'EdgeColor', 'none');
    xlabel('x'); ylabel('y'); zlabel('u');
    title(sprintf('Hard Problem at t = %.3f', t));
    colorbar;
    axis([0 Lx 0 Ly -0.1 1]);
    pause(0.1); 
    frame = getframe(gcf);
    writeVideo(v_hard, frame);
end
close(v_hard);

fprintf('Hard Problem simulation completei.avi.\n');
