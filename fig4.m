% Parameters
alpha = [-0.9, -0.7, -0.8];  % Values of alpha1, alpha2, alpha3 (decay rates)
f = [10, 15, 25];            % Values of f1, f2, f3 (frequencies)
c = [10, 2, 1];              % Values of c1, c2, c3 (amplitudes)

% Define the time range
t = linspace(0, 1, 1000);  % Time from 0 to 1 second with 1000 points

% Calculate the original signal x(t)
x_t = c(1) * exp(2 * pi * alpha(1) * t) .* cos(2 * pi * f(1) * t) + ...  % First component
      c(2) * exp(2 * pi * alpha(2) * t) .* cos(2 * pi * f(2) * t) + ...  % Second component
      c(3) * exp(2 * pi * alpha(3) * t) .* cos(2 * pi * f(3) * t);       % Third component

% Define the poles as a complex array
poles = 2 * pi * [-0.8 + 1i*25, -0.8 - 1i*25, -0.7 + 1i*15, -0.7 - 1i*15, -0.9 + 1i*10, -0.9 - 1i*10];

% Define the residues as a complex array
residues = pi * [-0.8 + 1i*25, -0.8 - 1i*25, -1.4 + 1i*30, -1.4 - 1i*30, -9 + 1i*100, -9 - 1i*100];

% Initialize the reconstructed signal x_t_r as zeros (same size as t)
x_t_r = zeros(size(t));

% Loop over each pole and residue to reconstruct the signal x(t)
for n = 1:length(poles)
    x_t_r = x_t_r + (residues(n) / poles(n)) * exp(poles(n) * t);  % Add each term of the series
end

% Plot the original signal and the reconstructed signal
figure('Position', [10 10 460 260], 'Color', [1 1 1]);  % Create the figure window
plot(t, x_t, 'k', 'LineWidth', 3);  % Plot the original signal in black with thick line
hold on
plot(t, x_t_r, 'r', 'LineWidth', 1);  % Plot the reconstructed signal in red with thinner line
xlabel('Time (s)');  % Label the x-axis
ylabel('x(t)');  % Label the y-axis
title('Plot of x(t)');  % Title of the plot
grid on;  % Add grid to the plot

% Save the figure as an EPS file
saveas(gcf, 'fig4', 'epsc');
