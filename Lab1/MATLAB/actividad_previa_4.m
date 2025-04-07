%% Parámetros configurables
A = 1;              % Amplitud de la señal
fc = 1000;         % Frecuencia de la señal (Hz)
Ts_signal = 10e-6;  % Muestreo de la señal original (10 µs)
duracion = 0.1;     % Duración en segundos
fs = 10000;          % Frecuencia de muestreo PAM (Hz)
d = 0.3;            % Ciclo de trabajo (τ/Ts_pulse)

%% Cálculos derivados
Ts_pulse = 1/fs;    % Periodo de muestreo PAM
tau = d * Ts_pulse; % Ancho del pulso

t_samples = 0:Ts_pulse:duracion - Ts_pulse;  % Vector de tiempos de muestreo

%% Generar señal original
t = 0:Ts_signal:duracion - Ts_signal;  % Vector tiempo (fila: 1xN)
m = A * sin(2*pi*fc*t); 

p = (mod(t, Ts_pulse) < tau);   % 1 durante τ, 0 en otro caso

% Aplicar modulación AM
s = m .* p;  % Señal modulada

% 3. Muestreo Instantáneo y Natural
m_samples_natural = A * sin(2*pi*fc*t_samples);   % PAM Natural
m_samples_instant = A * sin(2*pi*fc*t_samples);   % PAM Instantáneo

% 4. Crear tren de pulsos
s_natural = zeros(size(t));
s_instant = zeros(size(t));

for n = 1:length(t_samples)
    ventana = (t >= (n-1)*Ts_pulse) & (t < (n-1)*Ts_pulse + tau);
    s_natural(ventana) = m_samples_natural(n);
    s_instant(ventana) = m_samples_instant(n);
end

%% Graficación separada en subplots
figure;

% 1. Señal original
subplot(3,1,1);
plot(t, m, 'b', 'LineWidth', 1.5);
title('Señal Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
xlim([0 0.005]);
grid on;

% 2. Señal PAM Natural
subplot(3,1,2);
plot(t, s_natural, 'g--', 'LineWidth', 1.2);
hold on;
stem(t_samples, m_samples_natural, '^g', 'MarkerSize', 4, 'HandleVisibility', 'off');
title('PAM Natural');
xlabel('Tiempo (s)');
ylabel('Amplitud');
xlim([0 0.005]);
grid on;

% 3. Señal PAM Instantáneo
subplot(3,1,3);
plot(t, s_instant, 'r:', 'LineWidth', 1.2);
hold on;
stem(t_samples, m_samples_instant, 'or', 'MarkerSize', 4, 'HandleVisibility', 'off');
title('PAM Instantáneo');
xlabel('Tiempo (s)');
ylabel('Amplitud');
xlim([0 0.005]);
grid on;