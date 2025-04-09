
% Parámetros configurables
A = 1;              % Amplitud de la señal
fc = 1000;          % Frecuencia de la señal (Hz)
fs = 10000;         % Frecuencia de muestreo del tren de pulsos (Hz)
d = 0.1;            % Ciclo de trabajo (d = τ / Ts_pulse)
fs_signal = 1e5;    % Frecuencia de muestreo de m(t) (Hz)
duracion = 1;       % Duración en segundos

% Cálculo de parámetros derivados
Ts_signal = 1 / fs_signal;      % Periodo de muestreo de m(t)
Ts_pulse = 1 / fs;              % Periodo del tren de pulsos
tau = d * Ts_pulse;             % Ancho del pulso (τ)
t = 0:Ts_signal:duracion - Ts_signal; % Vector de tiempo

% Generar señal sinusoidal m(t)
m = A * sin(2 * pi * fc * t);

% Generar tren de pulsos p(t) con muestreo natural
p = (mod(t, Ts_pulse) < tau);   % 1 durante τ, 0 en otro caso

% Aplicar modulación AM
s = m .* p;  % Señal modulada

% Visualización
figure;
subplot(3,1,1);
plot(t, m);
xlim([0, 2/fc]);  % Mostrar 2 periodos de m(t)
title('Señal Original m(t)');
grid on;

subplot(3,1,2);
plot(t, p);
xlim([0, 2/fc]);  % Mostrar 2 periodos del tren de pulsos
ylim([-0.1, 1.1]);
title('Tren de Pulsos p(t)');
grid on;

subplot(3,1,3);
plot(t, s);
xlim([0, 2/fc]);  % Mostrar 2 periodos de la señal modulada
title('Señal Modulada s(t) (AM con muestreo natural)');
grid on;