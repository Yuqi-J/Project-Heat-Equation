% Parameter settings
L = 1;      % Interval length
T = 0.1;    % Total time
N = 100;    % Number of spatial grids
M = 1000;   % Time steps
alpha = 0.01; % Thermal diffusivity

dx = L / N; % Space step
dt = T / M; % Time step

% Check stability conditions
if dt > dx^2 / (2 * alpha)
    fprintf('Current time step dt = %f does not satisfy the stability condition (dt <= dx^2 / (2 * alpha))\n', dt);
else
    fprintf('Current time step dt = %f satisfies stability condition (dt <= dx^2 / (2 * alpha))\n', dt);
end

% Initial and boundary conditions
x = linspace(0, L, N+1); % Creating a spatial grid of points
u = sin(pi * x / L); % Initial temperature distribution (sine function)
u(1) = 0; % Left boundary condition (boundary temperature is 0)
u(end) = 0; % Right boundary condition (boundary temperature is 0)

% Explicit Euler method for solving the heat conduction equation
tic; % Start timing
u_matrix = zeros(M+1, N+1); % Store the results for each time step
u_matrix(1, :) = u; % Storing initial conditions
for n = 1:M
    u_new = zeros(1, N+1); % Initialize the temperature distribution for a new time step
    for i = 2:N
        % Explicit Euler method formula：u_new(i) = u(i) + alpha * dt / dx^2 * (u(i+1) - 2*u(i) + u(i-1))
        u_new(i) = u(i) + alpha * dt / dx^2 * (u(i+1) - 2*u(i) + u(i-1));
    end
    u = u_new; % Update temperature distribution
    u_matrix(n+1, :) = u; % Store the results of the current time step
end
elapsedTime = toc; % End timer
fprintf('Explicit Euler method computation time: %f seconds\n', elapsedTime);
fprintf('Average time per iteration: %e seconds\n', elapsedTime / M);

% Creating a time vector
t = linspace(0, T, M+1);

% Visualize the results
figure;
surf(x, t, u_matrix, 'EdgeColor', 'none'); % Plotting a 3D surface without showing grid lines
colormap jet; % Using the 'jet' colormap
c = colorbar; % Add a color bar
c.Label.String = 'temperature (u)'; % Color Bar Labels
xlabel('x'); % x-axis label
ylabel('t'); % y-axis label
zlabel('u(x,t)'); % z-axis label
title('Heat equation solution'); % Figure title
view(3); % 3D View

% 2D Visualization Results
figure;
contourf(x, t, u_matrix, 50); % Draw a contour map showing 50 contour levels
colormap jet; % Using the 'jet' colormap
c = colorbar; % Add a color bar
c.Label.String = 'temperature (u)'; % Color Bar Labels
xlabel('x'); % x-axis label
ylabel('t'); % y-axis label
title('2D Visualization of the Solution to the Heat Equation'); % Figure title

% Added graph visualization of numerical solutions to the heat conduction equation
figure;
plot(x, u_matrix(1, :), 'r', 'LineWidth', 2); % Initial conditions
hold on;
plot(x, u_matrix(round(M/4), :), 'g', 'LineWidth', 2); % 1/4 time step
plot(x, u_matrix(round(M/2), :), 'b', 'LineWidth', 2); % 1/2 time step
plot(x, u_matrix(end, :), 'm', 'LineWidth', 2); % Final time step
hold off;
xlabel('x'); % x-axis label
ylabel('u(x,t)'); % y-axis label
title('Plot visualization of numerical solutions to the heat equation'); % Figure title
legend('t = 0', 't = T/4', 't = T/2', 't = T'); % Add a legend

% Convergence Proof
errors = zeros(1, 5); % Store the maximum error at different grid numbers
timeSteps = zeros(1, 5); % Store the average time per iteration for different numbers of grids
for k = 1:5
    N = 20 * k; % Increase the number of spatial grids
    dx = L / N; % Update space step
    dt = T / M; % The time step remains constant
    x = linspace(0, L, N+1); % Create new spatial grid points
    u = sin(pi * x / L); % Initial temperature distribution
    u(1) = 0; % Left boundary condition
    u(end) = 0; % Right boundary condition
    tic; % Start timing
    for n = 1:M
        u_new = zeros(1, N+1); % Initialize the temperature distribution for a new time step
        for i = 2:N
            u_new(i) = u(i) + alpha * dt / dx^2 * (u(i+1) - 2*u(i) + u(i-1)); % Explicit Euler method formula
        end
        u = u_new; % Update temperature distribution
    end
    elapsedTime = toc; % End timer
    timeSteps(k) = elapsedTime / M; % Calculate the average time per iteration
    u_exact = sin(pi * x / L) * exp(-alpha * (pi / L)^2 * T); % Compute analytical solution
    errors(k) = max(abs(u - u_exact)); % Calculate the maximum error
end

% Visualizing the convergence proof
figure;
loglog(20:20:100, errors, '-o'); % Plot the error as a function of the number of spatial grid cells
title('Convergence Proof'); % Figure title
xlabel('Number of spatial grids'); % x-axis label
ylabel('Maximum error'); % y-axis label

% Visual efficiency indicators
figure;
loglog(20:20:100, timeSteps, '-o'); % Plot the average time per iteration as a function of the number of spatial grids
title('Average time per iteration'); % Figure title
xlabel('Number of spatial grids'); % x-axis label
ylabel('Average time per iteration (seconds)'); % y-axis label