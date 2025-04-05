%-------------------------------
% Parámetros iniciales (ya definidos anteriormente)
fc = 1000;       % Frecuencia de la señal sinosoidal
fm = 100000;     % Frecuencia de muestreo 
tm = 1/fm;       % Periodo de muestreo
tiempo_max = 0.01; % Duración total (se puede ajustar)
ls = 200;        % Número de muestras

% Vector de tiempo
Tiempo = (0:ls-1)*tm;

% Señal sinosoidal (portadora)
y = sin(2*pi*fc*Tiempo);

%----- PAM: Muestreo Instantáneo ---------------------
fs = 5000;       % Frecuencia de muestreo para PAM
ts = 1/fs;       % Periodo de muestreo PAM
tau = 0.5*ts;    % Duración del pulso
r = floor(ts/tm); % Número de muestras entre dos instantes de muestreo
s = floor(tau/tm); % Número de muestras que dura el pulso activo

Vector_instantaneo_muestral = zeros(1,length(Tiempo));
for i = 1:length(y)
    if mod(i, r) == 0
       Vector_instantaneo_muestral(i:i+s) = y(i);
    end
end
Muestreo_Instantaneo = Vector_instantaneo_muestral(1:length(Tiempo));

%-------------------------------
% PCM: Cuantización de la señal instantánea
% Configurar el número de bits para cada palabra PCM (N es configurable)
N = 4;   % Ejemplo: 4 bits (se puede cambiar a cualquier valor entero positivo)

% Número de niveles de cuantización
L = 2^N;

% Determinar los valores mínimo y máximo de la señal muestreada
x = Muestreo_Instantaneo;
x_min = min(x);
x_max = max(x);

% Calcular el tamaño de paso (delta) de cuantización
delta = (x_max - x_min) / (L - 1);

% Cuantizar cada muestra: calcular el índice de cuantización y reconstruir el valor cuantizado
PCM_indices = round((x - x_min) / delta);
PCM_quantized = PCM_indices * delta + x_min;

% (Opcional) Convertir cada índice a una cadena binaria de N bits
PCM_bin = cell(size(PCM_indices));
for i = 1:length(PCM_indices)
   PCM_bin{i} = dec2bin(PCM_indices(i), N);
end

%-------------------------------
% Graficar resultados: Señal original vs. Señal cuantizada (PCM)
figure;
plot(Tiempo, x, 'b.-'); hold on;
stairs(Tiempo, PCM_quantized, 'r', 'LineWidth', 1.5);
title('Señal Instantánea y Cuantizada (PCM)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
legend('Señal Instantánea', 'PCM Cuantizado');
grid on;

% Mostrar algunos ejemplos de codificación PCM en consola
disp('Ejemplos de palabras PCM (en binario):');
for i = 1:10:length(PCM_bin)
    fprintf('Tiempo = %.6f s, Valor cuantizado = %.4f, Código PCM = %s\n', Tiempo(i), PCM_quantized(i), PCM_bin{i});
end
