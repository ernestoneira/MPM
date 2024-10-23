% Par?metros
alpha = [-0.9, -0.7, -0.8];  % Valores de alpha1, alpha2, alpha3
f = [10, 15, 25];            % Valores de f1, f2, f3
c = [10, 2, 1];              % Valores de c1, c2, c3

% Definir el rango de tiempo
t = linspace(0, 1, 1000);  % Tiempo de 0 a 1 segundo con 1000 puntos

% Calcular x(t)
x_t = c(1) * exp(2 * pi * alpha(1) * t) .* cos(2 * pi * f(1) * t) + ...
      c(2) * exp(2 * pi * alpha(2) * t) .* cos(2 * pi * f(2) * t) + ...
      c(3) * exp(2 * pi * alpha(3) * t) .* cos(2 * pi * f(3) * t);

% Dibujar la gr?fica
figure('Position',[10 10 460 260],'Color',[1 1 1]);
plot(t, x_t, 'k','LineWidth', 2);
hold on
xlabel('Time (s)');
ylabel('x(t)');
title('Plot of x(t)');
grid on;

% Guardar la figura como EPS
saveas(gcf, 'fig2', 'epsc');
