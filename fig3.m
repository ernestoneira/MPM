% Parameters
alpha = [-0.9, -0.7, -0.8];  % Values for alpha1, alpha2, alpha3
f = [10, 15, 25];            % Values for f1, f2, f3 (frequencies)
c = [10, 2, 1];              % Coefficients for c1, c2, c3

% Define the time range
t = linspace(0, 1, 1000);  % Time range from 0 to 1 second with 1000 points

% Calculate dx(t)/dt (the time derivative of x(t))
dx_dt = c(1) * (2 * pi * alpha(1) * exp(2 * pi * alpha(1) * t) .* cos(2 * pi * f(1) * t) - ...
                2 * pi * f(1) * exp(2 * pi * alpha(1) * t) .* sin(2 * pi * f(1) * t)) + ...
        c(2) * (2 * pi * alpha(2) * exp(2 * pi * alpha(2) * t) .* cos(2 * pi * f(2) * t) - ...
                2 * pi * f(2) * exp(2 * pi * alpha(2) * t) .* sin(2 * pi * f(2) * t)) + ...
        c(3) * (2 * pi * alpha(3) * exp(2 * pi * alpha(3) * t) .* cos(2 * pi * f(3) * t) - ...
                2 * pi * f(3) * exp(2 * pi * alpha(3) * t) .* sin(2 * pi * f(3) * t));

% Create the figure window
figure('Position', [10 10 460 260], 'Color', [1 1 1]);

% Plot the graph
plot(t, dx_dt, 'k', 'LineWidth', 2);  % Plot the derivative of x(t) in black
xlabel('Time (s)');                   % Label for the x-axis (time)
ylabel('dx(t)/dt');                   % Label for the y-axis (derivative of x(t))
title('Derivative of x(t)');          % Title of the graph
grid on;                              % Enable grid for better visualization

% Save the figure as an EPS file
saveas(gcf, 'fig2', 'epsc');
