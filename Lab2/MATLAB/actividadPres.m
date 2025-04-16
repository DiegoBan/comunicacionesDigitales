% Parámetros principales
num_bits = 10^4;           % Número de bits
Rs = 1e3;                  % Tasa de símbolos (1 kHz)
sps = 8;                   % Muestras por símbolo
Fs = Rs * sps;             % Frecuencia de muestreo 
span = 10;                 % Span del filtro en símbolos
snr_dB = 20;               % Relación señal/ruido en dB

rolloff_values = [0, 0.25, 0.75, 1];  % Valores de roll-off α

figure;

for k = 1:length(rolloff_values)
    rolloff = rolloff_values(k);
    
    % Generación de bits aleatorios
    bits = randi([0, 1], 1, num_bits);
    
    % Codificación NRZ-L (0 -> -1, 1 -> +1)
    symbols = 2 * bits - 1;
    
    % Sobremuestreo (interpolación manual)
    symbolsUp = zeros(1, length(symbols) * sps);
    symbolsUp(1:sps:end) = symbols;
    
    % Crear y aplicar filtro Raised Cosine (coseno alzado)
    rrcFilter = rcosdesign(rolloff, span, sps, 'normal');
    filteredSignal = filter(rrcFilter, 1, symbolsUp);
    
    % Calcular potencia de la señal para ajustar el nivel de ruido
    filteredSignalPower = sum(abs(filteredSignal).^2) / length(filteredSignal);
    noisePower = filteredSignalPower / (10^(snr_dB / 10));
    noise = sqrt(noisePower) * randn(size(filteredSignal));
    
    % Agregar ruido gaussiano blanco (AWGN)
    receivedSignal = filteredSignal + noise;
    
    % Generar diagrama de ojo
    numTraces = floor(length(receivedSignal) / (2 * sps));
    eyeData = zeros(numTraces, 2 * sps);
    
    for i = 1:numTraces
        idx = (i - 1) * sps + 1;
        if (idx + 2 * sps - 1) <= length(receivedSignal)
            eyeData(i, :) = receivedSignal(idx : idx + 2 * sps - 1);
        end
    end
    
    % Subgráfico para cada α
    subplot(2, 2, k);
    plot(eyeData', 'b');
    title(['Diagrama de ojo - \alpha = ', num2str(rolloff)]);
    xlabel('Muestras');
    ylabel('Amplitud');
    grid on;
end

sgtitle('Diagramas de ojo para diferentes valores de \alpha (coseno alzado)');
