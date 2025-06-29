clear; clc; close all;

%% 1. Cargar señal de audio original (x(t))
[x, fs] = audioread('mu.wav');              % Leer archivo de audio
x = x(:,1);                                  % Asegurar que la señal sea mono (un solo canal)
x = x / max(abs(x));                         % Normalizar amplitud de la señal
N = length(x);                               % Número total de muestras
t = (0:N-1)/fs;                              % Vector de tiempo en segundos

figure;
plot(t, x, 'Color', [1 0.6 0.6], 'LineWidth', 1.2);
xlabel('Tiempo [s]');
ylabel('Amplitud');
title('Señal de Entrada x(t)');
grid on;
saveas(gcf, '1_entrada_xt.png');

%% 2. Generar ruido agudo y añadirlo a la señal original
ruido_blanco = 0.1 * randn(size(x));         % Ruido blanco de baja amplitud
fc_ruido = 6000;                             % Frecuencia de corte para el filtro pasa alto (Hz)
[bh, ah] = butter(4, fc_ruido/(fs/2), 'high'); % Filtro pasa alto de 4to orden (Butterworth)
ruido_agudo = filter(bh, ah, ruido_blanco);  % Aplicar filtro al ruido blanco

x_ruido = x + ruido_agudo;                   % Señal contaminada con ruido agudo
x_ruido = x_ruido / max(abs(x_ruido));       % Normalizar para mantener escala uniforme

figure;
plot(t, x_ruido, 'Color', [1 0.7 0.7], 'LineWidth', 1.1);
xlabel('Tiempo [s]');
ylabel('Amplitud');
title('Señal Contaminada con Ruido Agudo');
grid on;
saveas(gcf, '2_ruido_agudo_xt.png');

%% 3. Transformada de Fourier de la señal contaminada
X = fftshift(fft(x_ruido));                 % FFT centrada en cero
f = linspace(-fs/2, fs/2, N);               % Vector de frecuencias
mag_X = abs(X);                             % Magnitud del espectro

figure;
plot(f, mag_X / max(mag_X), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.3);
hold on;
fc = 6000;                                  % Frecuencia de corte del filtro ideal
xline(fc, '--k', ['f_c = ' num2str(fc) ' Hz'], 'LabelVerticalAlignment', 'bottom');
xlabel('Frecuencia [Hz]');
ylabel('Magnitud Normalizada');
title('Espectro de Magnitud de X(f)');
xlim([-8000 8000]); ylim([0 1.1]);
grid on;
saveas(gcf, '3_espectro_xf.png');

%% 4. Filtro pasa bajo ideal en frecuencia
H = double(abs(f) <= fc);                   % Filtro rectangular pasa bajo

figure;
plot(f, H, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.8);
xlabel('Frecuencia [Hz]');
ylabel('Magnitud');
title(['Respuesta en Frecuencia H(f) del Filtro Pasa Bajo Ideal (f_c = ', num2str(fc), ' Hz)']);
xlim([-8000 8000]); ylim([-0.1 1.1]);
grid on;
saveas(gcf, '4_respuesta_frecuencia_Hf.png');

%% 5. Comparación entre espectro y filtro
figure;
plot(f, mag_X / max(mag_X), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.3, 'DisplayName', 'Espectro X(f)');
hold on;
plot(f, H, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.6, 'LineStyle', '--', 'DisplayName', 'Filtro Ideal H(f)');
xlabel('Frecuencia [Hz]');
ylabel('Magnitud Normalizada');
title('Espectro X(f) y Filtro Pasa Bajo Ideal H(f)');
xlim([-8000 8000]); ylim([-0.1 1.2]);
legend('Location', 'northeast');
grid on;
saveas(gcf, '5_espectro_y_filtro.png');

%% 6. Respuesta al impulso h(t) del filtro ideal
t_h = linspace(-3/fc, 3/fc, 1000);           % Vector tiempo para h(t)
h_t = 2 * fc * sinc(2 * fc * t_h);           % h(t) ideal = sinc escalada

figure;
plot(t_h, h_t, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.7);
xlabel('Tiempo [s]');
ylabel('Amplitud');
title('Respuesta al Impulso h(t) del Filtro Pasa Bajo Ideal');
grid on;
saveas(gcf, '6_respuesta_impulso_ht.png');

%% 7. Aplicar filtro en frecuencia
Y = X .* H';                                % Aplicar filtro (multiplicación espectral)
mag_Y = abs(Y);                             % Magnitud del espectro filtrado

figure;
plot(f, mag_Y / max(mag_Y), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.3);
xlabel('Frecuencia [Hz]');
ylabel('Magnitud Normalizada');
title('Espectro de Y(f) después del Filtrado');
xlim([-8000 8000]); ylim([0 1.1]);
grid on;
saveas(gcf, '7_espectro_filtrado.png');

% Comparación espectral directa
figure;
plot(f, mag_X / max(mag_X), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Antes del filtrado');
hold on;
plot(f, mag_Y / max(mag_Y), 'b', 'LineWidth', 1.2, 'DisplayName', 'Después del filtrado');
xline(fc, '--k', ['f_c = ' num2str(fc) ' Hz'], 'LabelVerticalAlignment', 'bottom');
xlabel('Frecuencia [Hz]');
ylabel('Magnitud Normalizada');
title('Comparación de espectros antes y después del filtrado');
legend; grid on;
saveas(gcf, '7b_comparacion_espectros.png');

%% 8. Transformada inversa para recuperar y(t)
y = real(ifft(ifftshift(Y)));              % Transformada inversa centrada
y = y / max(abs(y));                        % Normalizar
t = (0:length(y)-1)/fs;

figure;
plot(t, y, 'Color', [0.6 1 0.6], 'LineWidth', 1.4);
xlabel('Tiempo [s]');
ylabel('Amplitud');
title('Señal Filtrada y(t)');
grid on;
saveas(gcf, '8_salida_yt.png');

%% 9. Filtrado por convolución (verificación)
h_t_discreta = 2*fc * sinc(2*fc * (-N/2:N/2)/fs);    % h(t) discreta en ventana larga
y_conv = conv(x_ruido, h_t_discreta, 'same');        % Convolución en tiempo
y_conv = y_conv / max(abs(y_conv));                  % Normalización

figure;
plot(t, y_conv, 'Color', [0.6 0.8 1], 'LineWidth', 1.3);
xlabel('Tiempo [s]');
ylabel('Amplitud');
title('Señal Filtrada mediante Convolución Temporal');
grid on;
saveas(gcf, '9_filtrado_convolucion.png');

%% 10. Guardar señal filtrada
audiowrite('salida_filtrada.wav', y, fs);

%% 11. Calcular relación SNR antes y después
ruido_inicial = x_ruido - x;
ruido_final = y - x;
snr_antes = 20*log10(norm(x)/norm(ruido_inicial));
snr_despues = 20*log10(norm(x)/norm(ruido_final));

fprintf('\\n--- Métricas de SNR ---\\n');
fprintf('SNR antes del filtrado: %.2f dB\\n', snr_antes);
fprintf('SNR después del filtrado: %.2f dB\\n', snr_despues);

%% 12. Reproducir señales
disp("Reproduciendo señal ORIGINAL...");
sound(x, fs); pause(length(x)/fs + 1);

disp("Reproduciendo señal CONTAMINADA...");
sound(x_ruido, fs); pause(length(x)/fs + 1);

disp("Reproduciendo señal FILTRADA (frecuencia)...");
sound(y, fs); pause(length(y)/fs + 1);

disp("Reproduciendo señal FILTRADA (convolución)...");
sound(y_conv, fs); pause(length(y_conv)/fs + 1);

%% 13. Comparar SNR con diferentes frecuencias de corte
fc_values = [3000, 4000, 5000, 6000, 7000, 8000, 9000];
snrs = zeros(size(fc_values));

for i = 1:length(fc_values)
    fc_test = fc_values(i);
    H_test = double(abs(f) <= fc_test);
    Y_test = X .* H_test';
    y_test = real(ifft(ifftshift(Y_test)));
    y_test = y_test / max(abs(y_test));
    ruido_final = y_test - x;
    snrs(i) = 20 * log10(norm(x)/norm(ruido_final));
end

fprintf('\n--- Resultados de SNR por frecuencia de corte ---\n');
fprintf('f_c (Hz)\tSNR (dB)\n');
for i = 1:length(fc_values)
    fprintf('%d\t\t%.2f\n', fc_values(i), snrs(i));
end

figure;
plot(fc_values, snrs, '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Frecuencia de Corte f_c [Hz]');
ylabel('SNR [dB]');
title('SNR vs. Frecuencia de Corte del Filtro Pasa Bajo');
grid on;
saveas(gcf, 'snr_vs_fc.png');
