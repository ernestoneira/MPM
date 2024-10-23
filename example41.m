% Example section 4.1: Analysis of a Noisy Signal

clear all
close all

alpha = [-0.9, -0.7, -0.8];  % Decay rates for each component
f = [10, 15, 25];            % Frequencies for each component
c = [10, 2, 1];              % Amplitudes for each component
noise_level = max(c) * 0.4;   % Noise level set to 20% of the max signal amplitude

% Define the time range
t = linspace(0, 1, 1000);  % Time range from 0 to 1 second with 1000 points

% Calculate the clean signal x(t)
x_t = c(1) * exp(2 * pi * alpha(1) * t) .* cos(2 * pi * f(1) * t) + ...
      c(2) * exp(2 * pi * alpha(2) * t) .* cos(2 * pi * f(2) * t) + ...
      c(3) * exp(2 * pi * alpha(3) * t) .* cos(2 * pi * f(3) * t);

% Add noise to the clean signal
noise = noise_level * (rand(1, length(t)) - 0.5);  % Generate random noise
x = x_t + noise;  % Noisy signal

% Compute the derivative of the noisy signal
dx_dt = diff(x) ./ diff(t);  % Numerical differentiation
newt = t(1:end-1) + (t(2) - t(1)) / 2;  % Adjust time vector for the derivative

% Apply the Matrix Pencil Method to reconstruct the signal
[p, r, x_r, intx] = matrixPencilMethod(newt, dx_dt, 'Accuracy', 0.02, 'method', 'CNRB', 'DC', 'y');

% Plot the original, noisy, and reconstructed signals
figure('Position', [10 10 460 260], 'Color', [1 1 1]);  % Create the figure window
plot(t, x, 'g', 'LineWidth', 1);  % Plot the noisy signal (green)
hold on
plot(t, x_t, 'k', 'LineWidth', 4);  % Plot the clean original signal (black)
plot(newt, intx, 'r', 'LineWidth', 1);  % Plot the reconstructed signal (red)
xlabel('Time (s)');  % Label the x-axis
ylabel('x(t)');  % Label the y-axis
title('Plot of x(t)');  % Title of the plot
grid on;  % Add grid to the plot
legend('Noisy Signal','Original Signal','Recovered Signal');  % Legend

% Calculate R^2 for Noisy vs Original, Reconstructed vs Original
% R^2 formula: R^2 = 1 - sum((y_real - y_pred)^2) / sum((y_real - mean(y_real))^2)

% R^2 between noisy signal and original signal
SS_res_noise = sum((x_t - x).^2);  % Residual sum of squares for noisy signal
SS_tot = sum((x_t - mean(x_t)).^2);  % Total sum of squares for the original signal
R2_noise = 1 - SS_res_noise / SS_tot;

% R^2 between reconstructed signal and original signal
SS_res_recon = sum((x_t(1:end-1) - intx).^2);  % Residual sum of squares for reconstructed signal
R2_recon = 1 - SS_res_recon / SS_tot;

% Display R^2 values
fprintf('R^2 (Noisy vs Original): %.4f\n', R2_noise);
fprintf('R^2 (Reconstructed vs Original): %.4f\n', R2_recon);

% Save the figure as an EPS file
saveas(gcf, 'fig5', 'epsc');
