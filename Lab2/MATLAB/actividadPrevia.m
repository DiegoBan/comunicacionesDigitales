% Parámetros
Frecuencia_Base = 1;            % Frecuencia base
Valores_Alpha = [0, 0.25, 0.75, 1];  % Factores de roll-off
t = linspace(0, 5/Frecuencia_Base, 1000);  % Vector de tiempo
f = linspace(-2*Frecuencia_Base, 2*Frecuencia_Base, 1000); % Vector de frecuencia

%% Graficar todas las respuestas al impulso en una misma figura con líneas más gruesas
figure;
hold on;
for Alpha = Valores_Alpha
    fd = Alpha * Frecuencia_Base;
    % Cálculo de la respuesta al impulso (he(t))
    he_t = 2*Frecuencia_Base * (sin(2*pi*Frecuencia_Base*t)./(2*pi*Frecuencia_Base*t)) ...
           .* (cos(2*pi*fd*t)./(1 - (4*fd*t).^2));
    % Corrección en t = 0 para evitar indeterminación
    he_t(t==0) = 2*Frecuencia_Base;
    
    % Graficar con línea de ancho 2 y asignar etiqueta para la leyenda
    plot(t, he_t, 'LineWidth', 2, 'DisplayName', ['\alpha = ' num2str(Alpha)]);
end
hold off;
grid on;
title('Respuesta al impulso para diferentes \alpha');
xlabel('Tiempo t');
ylabel('he(t)');
legend show;

%% Graficar todas las respuestas en frecuencia en una misma figura con líneas más gruesas
figure;
hold on;
for Alpha = Valores_Alpha
    H = zeros(size(f));
    for i = 1:length(f)
        Frecuencia_absoluta = abs(f(i));
        if Frecuencia_absoluta < Frecuencia_Base * (1 - Alpha)
            H(i) = 1;
        elseif Frecuencia_absoluta <= Frecuencia_Base * (1 + Alpha)
            H(i) = 0.5 * (1 + cos(pi/(2*Alpha*Frecuencia_Base) * (Frecuencia_absoluta - Frecuencia_Base * (1 - Alpha))));
        else
            H(i) = 0;
        end
    end
    
    % Graficar la respuesta en frecuencia con línea de ancho 2 y asignar etiqueta
    plot(f, H, 'LineWidth', 2, 'DisplayName', ['\alpha = ' num2str(Alpha)]);
end
hold off;
grid on;
title('Respuesta en frecuencia para diferentes \alpha');
xlabel('Frecuencia f');
ylabel('He(f)');
legend show;
