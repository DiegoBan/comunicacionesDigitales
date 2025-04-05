%% Configuración de parámetros
A = 1;              % Amplitud de la señal
fc = 1000;          % Frecuencia de la señal (Hz)
Ts = 1/100000;      % Periodo de muestreo (segundos), configurable
T_total = 0.01;     % Duración total de la señal (segundos)
t = 0:Ts:T_total;   % Vector de tiempo

%% Generación de la señal sinusoidal m(t)
m_t = A * sin(2*pi*fc*t);

%% Visualización de la señal
figure;
plot(t, m_t, 'LineWidth',1.5);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Sinusoidal m(t)');
grid on;
