%% Parámetros configurables
A = 1;              % Amplitud de la señal
fc = 1000;          % Frecuencia de la señal (Hz)
Ts_signal = 10e-6;  % Muestreo de la señal original (10 µs)
duracion = 0.1;     % Duración en segundos
fs = 5000;          % Frecuencia de muestreo PAM (Hz)
d = 0.3;            % Ciclo de trabajo (τ/Ts_pulse)
N = 3;              % Número de bits por palabra PCM

%% Cálculos derivados
Ts_pulse = 1/fs;    % Periodo de muestreo PAM
tau = d * Ts_pulse; % Ancho del pulso

%% Generar señal original
t = 0:Ts_signal:duracion - Ts_signal;
m = A * sin(2*pi*fc*t);

%% Generar señal PAM con muestreo instantáneo
t_samples = 0:Ts_pulse:t(end);
m_samples = A * sin(2*pi*fc*t_samples);

% Crear tren de pulsos PAM original
s = zeros(size(t));
for n = 0:length(t_samples)-1
    ventana = (t >= n*Ts_pulse) & (t < n*Ts_pulse + tau);
    s(ventana) = m_samples(n+1);
end

%% Generar señal PAM cuantizada
% Cuantización
quantized_indices = round((m_samples + A) * (2^N - 1) / (2*A));
quantized_indices = max(min(quantized_indices, 2^N - 1), 0);

% Convertir índices a valores de voltaje
samples_quantized = (quantized_indices/(2^N - 1)) * 2*A - A;

% Crear tren de pulsos PAM cuantizado
s_quant = zeros(size(t));
for n = 0:length(t_samples)-1
    ventana = (t >= n*Ts_pulse) & (t < n*Ts_pulse + tau);
    s_quant(ventana) = samples_quantized(n+1);
end

%% Visualización en mismo gráfico
figure;
hold on;

% Señal original
plot(t, m, 'b', 'LineWidth', 1.5, 'DisplayName', 'Señal original');

% PAM muestreado instantáneo
plot(t, s, 'r--', 'LineWidth', 1.2, 'DisplayName', 'PAM instantáneo');

% PAM cuantizado
plot(t, s_quant, 'g:', 'LineWidth', 2, 'DisplayName', ['PAM cuantizado (N=', num2str(N), ' bits)']);

xlim([0 0.005]);
title('Comparación de señales');
xlabel('Tiempo (s)');
ylabel('Amplitud');
legend('show', 'Location', 'southeast');
grid on;
hold off;

%% Visualización adicional de niveles de cuantización (opcional)
figure;
stem(t_samples, m_samples, 'r', 'DisplayName', 'Muestras originales');
hold on;
stem(t_samples, samples_quantized, 'g', 'LineWidth', 2, 'DisplayName', 'Muestras cuantizadas');
title('Niveles de cuantización');
xlabel('Tiempo (s)');
ylabel('Amplitud');
xlim([0 0.005]);
legend('show');
grid on;