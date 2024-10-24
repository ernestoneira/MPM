% Example section 4.1: Analysis of a Noisy Signal

clear all
close all

%variables for the F() calculation of the Montena D-dot probe
kbal=8;
katt=10;
kopt=1;
Rs=50;
Aeq=2e-4;

%distance to the sensor during the test meters
d=4;

data=dlmread('Mesoband_signal_raw.txt');

t=data(:,1);
V=data(:,2);

% Apply the Matrix Pencil Method to reconstruct the signal
[p, r, x_r, intx] = matrixPencilMethod(t, V, 'Accuracy', 0.05, 'method', 'CNRU', 'DC', 'n');

%compensation
Kf=10^((kbal+katt+kopt)/20)/(Rs*Aeq*e0);
E_int=cumtrapz(t, V);

E=intx*Kf*d;       
E_int=E_int*Kf*d;
% Plot the original, noisy, and reconstructed signals
figure('Position', [10 10 460 260], 'Color', [1 1 1]);  % Create the figure window
plot(t*1e9, E/1000, 'k', 'LineWidth', 2);  % Plot the noisy signal (green)
hold on
plot(t*1e9, E_int/1000, 'r', 'LineWidth', 2);  % Plot the noisy signal (green)

xlabel('Time (ns)');  % Label the x-axis
ylabel('Electric Field at 1 meter[kV/m]');  % Label the y-axis
title('Electric Field radiated by a MesoBand Source');  % Title of the plot
grid on;  % Add grid to the plot
axis([0 50 -18 25])
legend('Proposed Method','Integrated using trapezoidal rule')

saveas(gcf, 'fig7', 'epsc');
