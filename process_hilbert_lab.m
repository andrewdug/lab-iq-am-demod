%% ========================================================================
% Лабораторная работа:
% "Формирование низкочастотных квадратурных составляющих сигнала
%  с помощью цифрового фильтра Гильберта"
% Скрипт цифровой обработки сигналов (GNU Octave / MATLAB) для пп. 1-8
% Бригада M = 1: f0 = 430 МГц, Fm = 3 МГц, fd = 1720 МГц
%% ========================================================================
clear; clc; close all;

%% 1. Исходные параметры (Бригада M = 1)
M  = 1;
f0 = (400 + 30*M)*1e6;   % Несущая частота: 430 МГц
Fm = (3*M)*1e6;          % Частота модуляции: 3 МГц
m_mod = 1.0;             % Глубина АМ-модуляции (100%)
U0 = 1.0;                % Амплитуда несущей, В

% Выбор частоты дискретизации из условия работы на центре плоской вершины АЧХ:
% f0 = fd / 4  =>  fd = 4 * f0
fd = 4 * f0;             % 1720 МГц (1.72 ГГц)
dt = 1 / fd;             % Период дискретизации: ~581.395 пс

fprintf('====================================================================\n');
fprintf(' ИСХОДНЫЕ ДАННЫЕ (Бригада M = %d)\n', M);
fprintf('====================================================================\n');
fprintf(' Несущая частота f0          = %.2f МГц\n', f0/1e6);
fprintf(' Частота модуляции Fm        = %.2f МГц\n', Fm/1e6);
fprintf(' Частота дискретизации fd    = %.2f МГц (fd = 4*f0)\n', fd/1e6);
fprintf(' Период дискретизации dt     = %.6f нс (%.6e с)\n', dt*1e9, dt);
fprintf(' Элементарная задержка 2*dt  = %.6f нс\n', 2*dt*1e9);

%% 2. Расчет коэффициентов КИХ-фильтра Гильберта 10-го и 14-го порядков
% Фильтр 10-го порядка (N = 11 отсчетов, n0 = 6)
N10  = 11;
n0_10 = 6;
n10  = 1:N10;
k10  = n10 - n0_10;
h10_rect = zeros(1, N10);
nz10 = (k10 ~= 0);
h10_rect(nz10) = (1 - cos(pi * k10(nz10))) ./ (pi * k10(nz10));

% Окно Хэмминга для N = 11
w10_hamm = 0.54 - 0.46 * cos(2*pi*(n10 - 1)/(N10 - 1));
h10_hamm = h10_rect .* w10_hamm;
t0_10 = (n0_10 - 1) * dt; % 5 * dt = 2.906977 нс

fprintf('\n--------------------------------------------------------------------\n');
fprintf(' КОЭФФИЦИЕНТЫ ФИЛЬТРА 10-ГО ПОРЯДКА (N = 11, n0 = 6, t0 = %.6f нс)\n', t0_10*1e9);
fprintf('--------------------------------------------------------------------\n');
fprintf('  n | n-n0 |  h(n) (без окна) |   w(n) (Хэмминг) | hw(n) (с окном)\n');
fprintf('--------------------------------------------------------------------\n');
for i = 1:N10
    fprintf(' %2d |  %+2d  |   %+12.8f   |    %10.8f    |  %+12.8f\n', ...
        n10(i), k10(i), h10_rect(i), w10_hamm(i), h10_hamm(i));
end

% Фильтр 14-го порядка (N = 15 отсчетов, n0 = 8) с окном Хэмминга (п. 4)
N14  = 15;
n0_14 = 8;
n14  = 1:N14;
k14  = n14 - n0_14;
h14_rect = zeros(1, N14);
nz14 = (k14 ~= 0);
h14_rect(nz14) = (1 - cos(pi * k14(nz14))) ./ (pi * k14(nz14));
w14_hamm = 0.54 - 0.46 * cos(2*pi*(n14 - 1)/(N14 - 1));
h14_hamm = h14_rect .* w14_hamm;
t0_14 = (n0_14 - 1) * dt; % 7 * dt = 4.069767 нс

fprintf('\n--------------------------------------------------------------------\n');
fprintf(' КОЭФФИЦИЕНТЫ ФИЛЬТРА 14-ГО ПОРЯДКА (N = 15, n0 = 8, t0 = %.6f нс)\n', t0_14*1e9);
fprintf('--------------------------------------------------------------------\n');
for i = 1:2:N14
    fprintf(' %2d |  %+2d  |   %+12.8f   |    %10.8f    |  %+12.8f\n', ...
        n14(i), k14(i), h14_rect(i), w14_hamm(i), h14_hamm(i));
end

%% 3-4. Построение и сравнение АЧХ и ФЧХ (загрузка из Qucs-S .dat или расчет)
if exist('1_Hilbert_FIR10_Rect.dat', 'file') && exist('2_Hilbert_FIR10_Hamming.dat', 'file')
    ds1 = loadQucsDataSet('1_Hilbert_FIR10_Rect.dat');
    ds2 = loadQucsDataSet('2_Hilbert_FIR10_Hamming.dat');
    f_ac = getQucsVariable(ds1, 'acfrequency');
    H10_rect_ac = getQucsVariable(ds1, 's_hilb.v') ./ getQucsVariable(ds1, 's_in.v');
    H10_hamm_ac = getQucsVariable(ds2, 's_hilb.v') ./ getQucsVariable(ds2, 's_in.v');
    H10_del_ac  = getQucsVariable(ds2, 's_del.v')  ./ getQucsVariable(ds2, 's_in.v');
    if exist('3_Hilbert_FIR14_Hamming_IncreasedOrder.dat', 'file')
        ds3 = loadQucsDataSet('3_Hilbert_FIR14_Hamming_IncreasedOrder.dat');
        H14_hamm_ac = getQucsVariable(ds3, 's_hilb.v') ./ getQucsVariable(ds3, 's_in.v');
    else
        H14_hamm_ac = zeros(size(f_ac));
        for idx = 1:N14
            H14_hamm_ac = H14_hamm_ac + h14_hamm(idx) * exp(-1j*2*pi*f_ac*(idx-1)*dt);
        end
    end
else
    f_ac = linspace(1e6, fd/2, 860);
    H10_rect_ac = zeros(size(f_ac));
    H10_hamm_ac = zeros(size(f_ac));
    H14_hamm_ac = zeros(size(f_ac));
    for idx = 1:N10
        H10_rect_ac = H10_rect_ac + h10_rect(idx) * exp(-1j*2*pi*f_ac*(idx-1)*dt);
        H10_hamm_ac = H10_hamm_ac + h10_hamm(idx) * exp(-1j*2*pi*f_ac*(idx-1)*dt);
    end
    for idx = 1:N14
        H14_hamm_ac = H14_hamm_ac + h14_hamm(idx) * exp(-1j*2*pi*f_ac*(idx-1)*dt);
    end
    H10_del_ac = exp(-1j*2*pi*f_ac*t0_10);
end

figure('Name', 'Пункты 3-4: АЧХ и ФЧХ фильтров Гильберта', 'Color', 'w', 'Position', [80, 80, 1000, 680]);
subplot(2,1,1);
plot(f_ac/1e6, abs(H10_rect_ac), 'b--', 'LineWidth', 1.5); hold on;
plot(f_ac/1e6, abs(H10_hamm_ac), 'r-',  'LineWidth', 2.0);
plot(f_ac/1e6, abs(H14_hamm_ac), 'g-.', 'LineWidth', 1.8);
xline(f0/1e6, 'k:', 'LineWidth', 1.2);
grid on;
xlabel('Частота f, МГц');
ylabel('Модуль коэффициента передачи |H(f)|');
title('Сравнение АЧХ КИХ-фильтров Гильберта на одной координатной плоскости (п. 3, 4)');
legend('10-й порядок без окна (прямоуг.)', ...
       '10-й порядок с окном Хэмминга', ...
       '14-й порядок с окном Хэмминга', ...
       'Несущая частота f_0 = 430 МГц', 'Location', 'south');
ylim([0 1.25]);

subplot(2,1,2);
phi_rect_deg = unwrap(angle(H10_rect_ac))*180/pi;
phi_hamm_deg = unwrap(angle(H10_hamm_ac))*180/pi;
phi_del_deg  = unwrap(angle(H10_del_ac))*180/pi;
phi_diff_deg = angle(H10_hamm_ac ./ H10_del_ac)*180/pi;
plot(f_ac/1e6, phi_rect_deg, 'b--', 'LineWidth', 1.5); hold on;
plot(f_ac/1e6, phi_hamm_deg, 'r-',  'LineWidth', 1.8);
plot(f_ac/1e6, phi_del_deg,  'k:',  'LineWidth', 1.5);
plot(f_ac/1e6, phi_diff_deg, 'm-',  'LineWidth', 2.0);
grid on;
xlabel('Частота f, МГц');
ylabel('Фаза, градусы');
title('ФЧХ фильтра Гильберта, линии задержки t_0 и разность фаз \Delta\phi = -90^\circ');
legend('ФЧХ фильтра (без окна)', 'ФЧХ фильтра (Хэмминг)', ...
       'ФЧХ линии задержки s''(t)', 'Разность фаз \phi_{hilb} - \phi_{del} = -90^\circ', ...
       'Location', 'southwest');

%% 5. Временные диаграммы s(t), s'(t), s_hilb(t) и проверка задержки t0
% Формируем дискретные отсчеты с частотой дискретизации fd = 1720 МГц
% на интервале 3 периодов модуляции (T_sim = 3 / Fm = 1 мкс, N_samp = 1720 отсчетов)
N_samp = round(3 * fd / Fm); % 1720 отсчетов
n_idx  = 1:N_samp;
t_n    = (n_idx - 1) * dt;   % Дискретные моменты времени (n-1)*dt

% Входной АМ-сигнал (фаза несущей 30 град., как в схеме Qucs-S):
phi0 = 30 * pi / 180;
u_mod_n = sin(2*pi*Fm*t_n);                         % Модулирующее напряжение
A_env_n = U0 * (1 + m_mod * u_mod_n);               % Истинная огибающая АМ-сигнала
s_in_n  = A_env_n .* sin(2*pi*f0*t_n + phi0);       % Входной АМ-сигнал s(n)

% Задержанный на t0 = 5*dt сигнал s'(n) = s(n - 5)
delay_taps = n0_10 - 1; % 5 отсчетов
s_del_n  = [zeros(1, delay_taps), s_in_n(1:end-delay_taps)];

% Сигнал на выходе КИХ-фильтра Гильберта (с окном Хэмминга и без окна)
s_hilb_hamm_n = filter(h10_hamm, 1, s_in_n);
s_hilb_rect_n = filter(h10_rect, 1, s_in_n);

figure('Name', 'Пункт 5: Синхронизация s(t), s''(t) и s_hilb(t)', 'Color', 'w', 'Position', [100, 100, 1000, 680]);
subplot(2,1,1);
plot(t_n*1e9, s_in_n,        'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 1.2); hold on;
plot(t_n*1e9, s_del_n,       'b-', 'LineWidth', 1.8);
plot(t_n*1e9, s_hilb_hamm_n, 'r-', 'LineWidth', 1.8);
xline(t0_10*1e9, 'k--', 'LineWidth', 1.2);
xlim([0, 25]);
grid on;
xlabel('Время t, нс');
ylabel('Напряжение, В');
title(sprintf('Начальный участок осциллограмм: задержка сигнала t_0 = 5\\Delta t_д = %.3f нс и сдвиг фаз 90^\\circ', t0_10*1e9));
legend('Входной сигнал s(t)', 'Задержанный сигнал s''(t) = s(t - t_0)', ...
       'Сигнал на выходе фильтра \hat{s}(t)', 'Время групповой задержки t_0 = 2.907 нс');

subplot(2,1,2);
plot(t_n*1e6, s_del_n,       'b-', 'LineWidth', 0.8); hold on;
plot(t_n*1e6, s_hilb_hamm_n, 'r-', 'LineWidth', 0.8);
plot(t_n*1e6, A_env_n,       'k--', 'LineWidth', 1.5);
plot(t_n*1e6, -A_env_n,      'k--', 'LineWidth', 1.5);
grid on;
xlabel('Время t, мкс');
ylabel('Напряжение, В');
title('Проверка синхронности огибающих задержанного s''(t) и сопряженного по Гильберту \hat{s}(t) сигналов');
legend('s''(t) = s(t - t_0)', '\hat{s}(t)', 'Огибающая \pm A(t)');

%% 6. Формирование аналитического сигнала z(n) = s'(n) + j*s_hilb(n), годограф и спектры
% Отбрасываем переходный процесс первых 10 отсчетов фильтра
valid_idx = (N10 + 1):N_samp;
z_n = s_del_n + 1j * s_hilb_hamm_n;

% Годограф аналитического сигнала на комплексной плоскости
figure('Name', 'Пункт 6: Годограф и спектры s(t) и z(t)', 'Color', 'w', 'Position', [120, 120, 1050, 500]);
subplot(1,2,1);
plot(real(z_n(valid_idx)), imag(z_n(valid_idx)), 'b-', 'LineWidth', 1.0); hold on;
plot(real(z_n(valid_idx(1:20))), imag(z_n(valid_idx(1:20))), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
axis equal; grid on;
xlim([-2.2 2.2]); ylim([-2.2 2.2]);
xlabel('Действительная часть Re\{z(n)\} = s''(n), В');
ylabel('Мнимая часть Im\{z(n)\} = \hat{s}(n), В');
title('Годограф аналитического сигнала z(n) на комплексной плоскости');

% Двухсторонние спектры вещественного сигнала s(t) и аналитического сигнала z(t)
N_fft = N_samp;
f_axis = (-N_fft/2 : N_fft/2 - 1) * (fd / N_fft); % от -860 МГц до +860 МГц
S_in_spec = fftshift(abs(fft(s_in_n)) / N_fft);
Z_spec    = fftshift(abs(fft(z_n))    / N_fft);

subplot(1,2,2);
plot(f_axis/1e6, S_in_spec, 'b-',  'LineWidth', 1.5); hold on;
plot(f_axis/1e6, Z_spec,    'r--', 'LineWidth', 1.5);
grid on;
xlim([-600 600]);
xlabel('Частота f, МГц');
ylabel('Амплитуда спектральных составляющих, В');
title('Спектры входного вещественного s(t) и аналитического z(t) сигналов');
legend('Спектр |S(f)| (двухсторонний)', 'Спектр |Z(f)| (односторонний, f > 0)');

%% 7. Формирование низкочастотных квадратурных составляющих I(n) и Q(n)
% По формулам задания:
% I(n) =  s'(n)*cos(w0*(n-1)*dt) + s_hilb(n)*sin(w0*(n-1)*dt)
% Q(n) = -s'(n)*sin(w0*(n-1)*dt) + s_hilb(n)*cos(w0*(n-1)*dt)
w0 = 2 * pi * f0;
lo_cos_n = cos(w0 * t_n);
lo_sin_n = sin(w0 * t_n);

I_n =  s_del_n .* lo_cos_n + s_hilb_hamm_n .* lo_sin_n;
Q_n = -s_del_n .* lo_sin_n + s_hilb_hamm_n .* lo_cos_n;

% Амплитудные спектры квадратурных составляющих I(n) и Q(n)
I_spec = fftshift(abs(fft(I_n(valid_idx))) / length(valid_idx));
Q_spec = fftshift(abs(fft(Q_n(valid_idx))) / length(valid_idx));
f_iq   = (-length(valid_idx)/2 : length(valid_idx)/2 - 1) * (fd / length(valid_idx));

figure('Name', 'Пункт 7: Квадратурные составляющие I(t), Q(t) и их спектры', 'Color', 'w', 'Position', [140, 140, 1000, 680]);
subplot(2,1,1);
plot(t_n(valid_idx)*1e6, I_n(valid_idx), 'b-', 'LineWidth', 2.0); hold on;
plot(t_n(valid_idx)*1e6, Q_n(valid_idx), 'g-', 'LineWidth', 2.0);
grid on;
xlabel('Время t, мкс');
ylabel('Напряжение, В');
title('Осциллограммы низкочастотных квадратурных составляющих I(t) и Q(t)');
legend('Синфазная составляющая I(t)', 'Квадратурная составляющая Q(t)');

subplot(2,1,2);
stem(f_iq/1e6, I_spec, 'b', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
stem(f_iq/1e6, Q_spec, 'g--', 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([-15 15]);
grid on;
xlabel('Частота f, МГц');
ylabel('Амплитуда спектра, В');
title('Амплитудные спектры низкочастотных квадратурных составляющих |I(f)| и |Q(f)| (пики на 0 и \pm 3 МГц)');
legend('Спектр |I(f)|', 'Спектр |Q(f)|');

%% 8. Цифровой амплитудный демодулятор U(n) = sqrt(I^2(n) + Q^2(n)) и оценка искажений
U_demod_hamm10 = sqrt(I_n.^2 + Q_n.^2);

% Для сравнения: демодуляция с фильтром без окна (прямоугольное окно)
I_rect_n =  s_del_n .* lo_cos_n + s_hilb_rect_n .* lo_sin_n;
Q_rect_n = -s_del_n .* lo_sin_n + s_hilb_rect_n .* lo_cos_n;
U_demod_rect10 = sqrt(I_rect_n.^2 + Q_rect_n.^2);

% Истинная задержанная на t0 огибающая модулирующего сигнала:
A_env_del_n = [zeros(1, delay_taps), A_env_n(1:end-delay_taps)];
u_mod_del_n = [zeros(1, delay_taps), u_mod_n(1:end-delay_taps)];

figure('Name', 'Пункт 8: Демодуляция АМ-сигнала и оценка искажений', 'Color', 'w', 'Position', [160, 160, 1000, 680]);
subplot(2,1,1);
plot(t_n(valid_idx)*1e6, A_env_del_n(valid_idx),    'k--', 'LineWidth', 2.0); hold on;
plot(t_n(valid_idx)*1e6, U_demod_hamm10(valid_idx), 'r-',  'LineWidth', 1.8);
plot(t_n(valid_idx)*1e6, U_demod_rect10(valid_idx), 'b:',  'LineWidth', 1.5);
plot(t_n(valid_idx)*1e6, u_mod_del_n(valid_idx),    'Color', [0 0.6 0], 'LineStyle', '-.', 'LineWidth', 1.5);
grid on;
xlabel('Время t, мкс');
ylabel('Напряжение, В');
title('Сравнение восстановленной огибающей U(n) = \sqrt{I^2(n)+Q^2(n)} и модулирующего сигнала (п. 8)');
legend('Идеальная огибающая 1 + u_{мод}(t - t_0)', ...
       'Демодулированная огибающая U(n) (окно Хэмминга, 10 пор.)', ...
       'Демодулированная огибающая U(n) (без окна, 10 пор.)', ...
       'Модулирующее напряжение u_{мод}(t - t_0)');

subplot(2,1,2);
err_hamm = U_demod_hamm10(valid_idx) - A_env_del_n(valid_idx);
err_rect = U_demod_rect10(valid_idx) - A_env_del_n(valid_idx);
plot(t_n(valid_idx)*1e6, err_rect, 'b--', 'LineWidth', 1.2); hold on;
plot(t_n(valid_idx)*1e6, err_hamm, 'r-',  'LineWidth', 1.5);
grid on;
xlabel('Время t, мкс');
ylabel('Абсолютная погрешность \Delta U(t), В');
title('Погрешность восстановления огибающей \Delta U(t) = U(t) - [1 + u_{мод}(t - t_0)]');
legend(sprintf('Без окна (макс. откл. %.1f%% из-за выброса Гиббса |H(f_0)|=1.103)', max(abs(err_rect))/2*100), ...
       sprintf('С окном Хэмминга (макс. откл. %.2f%%, |H_w(f_0)|=1.013)', max(abs(err_hamm))/2*100));

fprintf('\n====================================================================\n');
fprintf(' ОЦЕНКА ТОЧНОСТИ ДЕМОДУЛЯЦИИ (Пункт 8):\n');
fprintf('  - Без весового окна (10 пор.): |H(f0)| = %.5f (погрешность огибающей +%.2f%%)\n', ...
    abs(sum(h10_rect .* exp(-1j*2*pi*f0*(0:N10-1)*dt))), (abs(sum(h10_rect .* exp(-1j*2*pi*f0*(0:N10-1)*dt)))-1)*100);
fprintf('  - С окном Хэмминга  (10 пор.): |H(f0)| = %.5f (погрешность огибающей +%.2f%%)\n', ...
    abs(sum(h10_hamm .* exp(-1j*2*pi*f0*(0:N10-1)*dt))), (abs(sum(h10_hamm .* exp(-1j*2*pi*f0*(0:N10-1)*dt)))-1)*100);
fprintf('  - С окном Хэмминга  (14 пор.): |H(f0)| = %.5f (погрешность огибающей %.2f%%)\n', ...
    abs(sum(h14_hamm .* exp(-1j*2*pi*f0*(0:N14-1)*dt))), (abs(sum(h14_hamm .* exp(-1j*2*pi*f0*(0:N14-1)*dt)))-1)*100);
fprintf('====================================================================\n');
