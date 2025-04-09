%actividad presencial
% Código MATLAB para generar diagramas de ojo para pulso coseno alzado
% con codificación NRZ-L y canal AWGN

clear;
close all;
clc;

% Parámetros
num_bits = 10000;               % Número de bits a transmitir
Tb = 1;                         % Tiempo de bit (normalizado)
samples_per_bit = 16;           % Muestras por bit
fs = samples_per_bit/Tb;        % Frecuencia de muestreo
T = Tb/samples_per_bit;         % Periodo de muestreo
f0 = 1/(2*Tb);                  % Frecuencia fundamental
B = 2*f0;                       % Ancho de banda absoluto
SNR_dB = 20;                    % Relación señal a ruido en dB

% Factores de roll-off a evaluar
rolloff = [0, 0.25, 0.75, 1];   
num_rolloff = length(rolloff);

% Generación de bits aleatorios (0 y 1)
bits = randi([0, 1], 1, num_bits);

% Codificación NRZ-L: 0 -> -1, 1 -> +1
symbols = 2*bits - 1;

% Tiempo para la visualización de diagrama de ojo (2 periodos de bit)
eye_span = 2*Tb;
eye_samples = eye_span/T;

% Creación del tren de pulsos (secuencia de símbolos sobremuestreada)
pulse_train = zeros(1, num_bits * samples_per_bit);
for i = 1:num_bits
    pulse_train((i-1)*samples_per_bit + 1:i*samples_per_bit) = symbols(i);
end

% Duración del filtro (en símbolos)
filter_span = 10;
t_filter = -filter_span*Tb:T:filter_span*Tb;

% Procesamiento para cada factor de roll-off
for i = 1:num_rolloff
    a = rolloff(i);  % Factor de roll-off actual
    f_delta = a*f0;  % Ancho de banda de roll-off
    
    % Respuesta al impulso del filtro coseno alzado
    h_rc = zeros(size(t_filter));
    for j = 1:length(t_filter)
        t = t_filter(j);
        if t == 0
            h_rc(j) = 2*f0;  % Para t=0, usar regla de L'Hôpital
        else
            sinc_term = sin(2*pi*f0*t)/(2*pi*f0*t);
            if a == 0
                cos_term = 1;  % Para roll-off = 0 (Nyquist)
            else
                cos_term = cos(2*pi*f_delta*t)/(1-(4*f_delta*t)^2);
            end
            h_rc(j) = 2*f0 * sinc_term * cos_term;
        end
    end
    
    % Normalizar el filtro para mantener la energía
    h_rc = h_rc / sqrt(sum(h_rc.^2));
    
    % Aplicar el filtro al tren de pulsos mediante convolución
    tx_signal = conv(pulse_train, h_rc, 'same');
    
    % Aplicar ruido gaussiano blanco aditivo (AWGN)
    signal_power = mean(abs(tx_signal).^2);
    SNR_linear = 10^(SNR_dB/10);
    noise_power = signal_power / SNR_linear;
    noise = sqrt(noise_power) * randn(size(tx_signal));
    rx_signal = tx_signal + noise;
    
    % Crear el diagrama de ojo
    figure('Name', ['Diagrama de Ojo (α = ', num2str(a), ')']);
    
    % Índices para dividir la señal en segmentos para el diagrama de ojo
    total_samples = length(rx_signal);
    eye_segments = floor(total_samples/samples_per_bit) - 1;
    
    hold on;
    for k = 1:eye_segments
        start_idx = (k-1)*samples_per_bit + 1;
        end_idx = start_idx + eye_samples - 1;
        
        if end_idx <= total_samples
            t_eye = (0:eye_samples-1)*T;
            plot(t_eye, rx_signal(start_idx:end_idx), 'b', 'LineWidth', 0.2);
        end
    end
    
    % Configuración del gráfico
    grid on;
    title(['Diagrama de Ojo - Pulso Coseno Alzado (α = ', num2str(a), ')']);
    xlabel('Tiempo (t/Tb)');
    ylabel('Amplitud');
    xlim([0, eye_span]);
    
    % Marcadores para el centro del ojo
    line([Tb, Tb], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
    
    % Mostrar información adicional en el diagrama
    text(0.1*Tb, max(ylim)*0.9, ['SNR = ', num2str(SNR_dB), ' dB'], 'FontSize', 8);
    text(0.1*Tb, max(ylim)*0.8, ['Muestras/bit = ', num2str(samples_per_bit)], 'FontSize', 8);
    text(0.1*Tb, max(ylim)*0.7, ['Total bits = ', num2str(num_bits)], 'FontSize', 8);
    
    hold off;
    
    % Analizar apertura del ojo
    mid_point = samples_per_bit/2;
    eye_openings = zeros(1, eye_segments);
    
    for k = 1:eye_segments
        start_idx = (k-1)*samples_per_bit + mid_point;
        eye_value = rx_signal(start_idx);
        if bits(k) == 1
            eye_openings(k) = eye_value;
        else
            eye_openings(k) = -eye_value;
        end
    end
    
    % Estadísticas de apertura del ojo
    avg_opening = mean(eye_openings);
    std_opening = std(eye_openings);
    
    fprintf('Factor de roll-off α = %0.2f:\n', a);
    fprintf('  Apertura media del ojo: %0.4f\n', avg_opening);
    fprintf('  Desviación estándar: %0.4f\n', std_opening);
    fprintf('  Relación apertura/desviación: %0.4f\n\n', avg_opening/std_opening);
end

% Opcional: Crear un diagrama de ojo adicional que muestre todos los factores juntos
figure('Name', 'Comparación de Diagramas de Ojo para Diferentes Factores de Roll-off');

% Para la comparación, generaremos una señal específica para cada factor y tramaremos 2 períodos
time_span = 2*Tb;
t_compare = 0:T:time_span-T;
sample_count = length(t_compare);

% Datos para la comparación: patrón alternante (peor caso)
test_bits = [1 0 1 0 1 0];
test_symbols = 2*test_bits - 1;

% Arrays para almacenar las señales para cada factor de roll-off
signals = zeros(num_rolloff, sample_count);

% Crear subplot para cada factor de roll-off
for i = 1:num_rolloff
    a = rolloff(i);
    f_delta = a*f0;
    
    % Crear filtro
    filter_span = 10;
    t_filter = -filter_span*Tb:T:filter_span*Tb;
    h_rc = zeros(size(t_filter));
    
    for j = 1:length(t_filter)
        t = t_filter(j);
        if t == 0
            h_rc(j) = 2*f0;
        else
            sinc_term = sin(2*pi*f0*t)/(2*pi*f0*t);
            if a == 0
                cos_term = 1;
            else
                cos_term = cos(2*pi*f_delta*t)/(1-(4*f_delta*t)^2);
            end
            h_rc(j) = 2*f0 * sinc_term * cos_term;
        end
    end
    
    % Normalizar
    h_rc = h_rc / max(h_rc);
    
    % Crear señal de prueba (con símbolos espaciados)
    pulse_length = length(t_filter);
    middle = (pulse_length + 1)/2;
    test_signal = zeros(1, 3 * pulse_length);
    
    for j = 1:length(test_symbols)
        position = middle + (j-1)*samples_per_bit;
        if position + samples_per_bit <= length(test_signal)
            test_signal(position:position+samples_per_bit-1) = test_symbols(j);
        end
    end
    
    % Aplicar filtro
    filtered_signal = conv(test_signal, h_rc, 'same');
    
    % Crear subplot
    subplot(2, 2, i);
    
    % Mostrar diagrama de ojo para este caso específico
    hold on;
    for k = 1:3
        start_idx = (k-1)*samples_per_bit + 1;
        if start_idx + sample_count - 1 <= length(filtered_signal)
            plot(t_compare, filtered_signal(start_idx:start_idx+sample_count-1), 'b', 'LineWidth', 1);
        end
    end
    
    grid on;
    title(['α = ', num2str(a)]);
    xlabel('Tiempo (t/Tb)');
    ylabel('Amplitud');
    xlim([0, time_span]);
    
    % Marcadores
    line([Tb, Tb], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
    hold off;
end

% Ajustar diseño de la figura
sgtitle('Comparación de Diagramas de Ojo para Diferentes Factores de Roll-off');