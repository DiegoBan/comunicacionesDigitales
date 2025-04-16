%% Parámetros configurables
A = 1;              % Amplitud de la señal
fc = 1000;          % Frecuencia de la señal (Hz)
Ts_signal = 10e-6;  % Muestreo de la señal original (10 µs)
duracion = 0.1;     % Duración en segundos
fs = 5000;          % Frecuencia de muestreo PAM (Hz) - ¡Configurable!
d = 0.3;            % Ciclo de trabajo (τ/Ts_pulse) - ¡Configurable!

%% Cálculos derivados
Ts_pulse = 1/fs;    % Periodo de muestreo PAM
tau = d * Ts_pulse; % Ancho del pulso

%% Generar señal original
t = 0:Ts_signal:duracion - Ts_signal;  % Vector tiempo
m = A * sin(2*pi*fc*t);                % Señal original

%% Generar señal PAM con muestreo instantáneo
% 1. Tiempos de muestreo PAM
t_samples = 0:Ts_pulse:t(end);  % Tiempos de muestreo

% 2. Muestras instantáneas (valor exacto en t_samples)
m_samples = A * sin(2*pi*fc*t_samples);  % Sin interpolación

% 3. Crear tren de pulsos
s = zeros(size(t));
for n = 0:length(t_samples)-1
    % Definir ventana temporal para cada pulso
    ventana = (t >= n*Ts_pulse) & (t < n*Ts_pulse + tau);
    % Asignar amplitud de la muestra al pulso
    s(ventana) = m_samples(n+1);
end

%% Visualización
figure;
subplot(2,1,1);
plot(t, m, 'b');
xlim([0 0.005]);  % Primeros 5 ms
title('Señal original m(t)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(2,1,2);
plot(t, s, 'r');
xlim([0 0.005]);  % Primeros 5 ms
hold on;
stem(t_samples, m_samples, 'k', 'MarkerSize', 5);  % Muestras instantáneas
title(['PAM Instantáneo: fs = ', num2str(fs), ' Hz, d = ', num2str(d)]);
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;