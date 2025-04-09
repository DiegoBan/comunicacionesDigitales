%% Parámetros configurables
A = 1;              % Amplitud de la señal
fc = 1000;          % Frecuencia de la señal (Hz)
Ts_signal = 10e-6;  % Muestreo de la señal original (10 µs)
duracion = 0.1;     % Duración en segundos
fs = 5000;          % Frecuencia de muestreo PAM (Hz)
d = 0.3;            % Ciclo de trabajo (τ/Ts_pulse)

%% Cálculos derivados
Ts_pulse = 1/fs;    % Periodo de muestreo PAM
tau = d * Ts_pulse; % Ancho del pulso

%% Generar señal original
t = 0:Ts_signal:duracion - Ts_signal;  % Vector tiempo (fila: 1xN)
m = A * sin(2*pi*fc*t); 

p = (mod(t, Ts_pulse) < tau);   % 1 durante τ, 0 en otro caso

% Aplicar modulación AM
s = m .* p;  % Señal modulada


% 3. Muestreo Instantáneo 
m_samples_instant = A * sin(2*pi*fc*t_samples);

% 4. Crear tren de pulsos
s_natural = zeros(size(t));
s_instant = zeros(size(t));

for n = 1:length(t_samples)
    ventana = (t >= (n-1)*Ts_pulse) & (t < (n-1)*Ts_pulse + tau);
    s_natural(ventana) = m_samples_natural(n);
    s_instant(ventana) = m_samples_instant(n);
end

%% Graficación (sin error)
figure;
hold on;

% Señal original
plot(t, m, 'b', 'LineWidth', 1.5, 'DisplayName', 'Señal original');

% Señal PAM natural
plot(t, s_natural, 'g--', 'LineWidth', 1.2, 'DisplayName', 'PAM Natural');

% Señal PAM instantáneo
plot(t, s_instant, 'r:', 'LineWidth', 1.2, 'DisplayName', 'PAM Instantáneo');

% Marcadores de muestras (mismas longitudes)
stem(t_samples, m_samples_natural, '^g', 'MarkerSize', 4, 'HandleVisibility', 'off');
stem(t_samples, m_samples_instant, 'or', 'MarkerSize', 4, 'HandleVisibility', 'off');

xlim([0 0.005]); 
xlabel('Tiempo (s)'); 
ylabel('Amplitud');
title(['Comparación PAM: fs = ', num2str(fs), ' Hz, d = ', num2str(d)]);
legend('Location', 'best');
grid on;