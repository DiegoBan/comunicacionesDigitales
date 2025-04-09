% Análisis del Pulso Coseno Alzado
% Este script genera gráficas de la respuesta al impulso y la respuesta en
% frecuencia del pulso coseno alzado para diferentes factores de roll-off.
clear;
close all;
clc;

% Parámetros
f0 = 1; % Frecuencia central (normalizada)
rolloff = [0, 0.25, 0.75, 1]; % Factores de roll-off (a)

% Vectores para las gráficas
t = linspace(0.01, 10, 1000); % Vector de tiempo (evitando t=0 para evitar división por cero)
f = linspace(-1.1, 1.1, 1000); % Vector de frecuencia ligeramente ampliado para mejor visualización

% Configuración de la figura para respuesta al impulso
figure('Name', 'Respuesta al Impulso', 'Position', [100, 100, 800, 600]);
hold on;

% Configuración de la figura para respuesta en frecuencia
figure('Name', 'Respuesta en Frecuencia', 'Position', [100, 700, 800, 600]);
hold on;

% Ciclo para cada factor de roll-off
color_codes = {'b', 'r', 'g', 'k'}; % Colores para cada valor de roll-off

for i = 1:length(rolloff)
    a = rolloff(i);
    
    % Cálculo de la respuesta al impulso usando la ecuación para pulso coseno alzado
    he_t = zeros(size(t));
    for j = 1:length(t)
        if abs(t(j) - 1/(2*f0)) < 1e-10 && a == 0 % Caso especial cuando a=0 y t=1/(2*f0)
            he_t(j) = pi/(2*f0);
        elseif abs(1 - (4*a*f0*t(j))^2) < 1e-10 % Caso especial cuando denominador es cercano a cero
            he_t(j) = (pi/(4*f0)) * (sin(pi*f0*t(j))/(pi*f0*t(j)));
        else
            term1 = sin(pi*f0*t(j)) / (pi*f0*t(j));
            term2 = cos(a*pi*f0*t(j)) / (1 - (2*a*f0*t(j))^2);
            he_t(j) = term1 * term2;
        end
    end
    
    % Cálculo de la respuesta en frecuencia
    He_f = zeros(size(f));
    for j = 1:length(f)
        % Respuesta en frecuencia del pulso coseno alzado
        if abs(f(j)) <= (1-a)*f0
            He_f(j) = 1; % Región plana
        elseif abs(f(j)) > (1-a)*f0 && abs(f(j)) <= (1+a)*f0
            % Región de transición con coseno alzado
            He_f(j) = 0.5 * (1 + cos((pi/(a*f0)) * (abs(f(j)) - (1-a)*f0)));
        else
            He_f(j) = 0; % Fuera del ancho de banda
        end
    end
    
    % Graficar respuesta al impulso
    figure(1);
    plot(t, he_t, color_codes{i}, 'LineWidth', 2, 'DisplayName', sprintf('a = %.2f', a));
    
    % Graficar respuesta en frecuencia
    figure(2);
    plot(f, He_f, color_codes{i}, 'LineWidth', 2, 'DisplayName', sprintf('a = %.2f', a));
end

% Configurar gráfica de respuesta al impulso
figure(1);
grid on;
title('Respuesta al Impulso del Pulso Coseno Alzado', 'FontSize', 14);
xlabel('Tiempo (t)', 'FontSize', 12);
ylabel('h_e(t)', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 12);
axis tight;

% Configurar gráfica de respuesta en frecuencia
figure(2);
grid on;
title('Respuesta en Frecuencia del Pulso Coseno Alzado', 'FontSize', 14);
xlabel('Frecuencia (f)', 'FontSize', 12);
ylabel('H_e(f)', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 12);
ylim([0, 1.1]); % Establecer límites de eje y de 0 a 1.1
xlim([-1.1, 1.1]); % Establecer límites de eje x entre -1.1 y 1.1

% Añadir líneas verticales para mostrar límites
figure(2);