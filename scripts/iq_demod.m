%% Цифровая амплитудная демодуляция по I/Q-составляющим (п. 2 и п. 4 задания)
% Вход: два CSV, экспортированные из графиков Qucs-S (I_out и Q_out).
% Огибающая: U[n] = sqrt(I[n]^2 + Q[n]^2).
% Работает в MATLAB (R2016b и новее) и в Octave.

clear; close all;

%% Исходные данные
fileI   = 'I_out.csv';   % CSV с графика I (время; напряжение)
fileQ   = 'Q_out.csv';   % CSV с графика Q (время; напряжение)
N       = 1;             % номер бригады
Fm      = 3 * N * 1e6;   % частота модуляции, Гц
m       = 1;             % коэффициент модуляции по заданию
Nsp     = 20;            % отсчётов на период модуляции (по заданию >= 10)
tSkip   = 0.3e-6;        % отбросить переходный процесс фильтров, с
outName = 'demod';       % префикс выходных файлов (PNG и CSV)

%% Чтение данных Qucs
[tI, I] = readQucsCsv(fileI);
[tQ, Q] = readQucsCsv(fileQ);

%% Дискретизация I[n], Q[n]
Ts    = 1 / (Nsp * Fm);
tFrom = max(tI(1), tQ(1));
tTo   = min(tI(end), tQ(end));
tn    = (ceil(tFrom / Ts) * Ts : Ts : tTo).';
In    = interp1(tI, I, tn, 'linear');
Qn    = interp1(tQ, Q, tn, 'linear');

%% Демодуляция
Un = sqrt(In.^2 + Qn.^2);

%% Сравнение с модулирующим напряжением
% На установившемся участке U[n] ~ U0 + a*cos(w*t) + b*sin(w*t);
% U0 - коэффициент передачи тракта, фаза гармоники - задержка фильтров.
w      = 2 * pi * Fm;
steady = tn >= tSkip;
if nnz(steady) < Nsp        % нужен хотя бы один период модуляции
    error('Мало отсчётов после tSkip = %.3g с: увеличьте длительность моделирования.', tSkip);
end
A     = [ones(nnz(steady), 1), cos(w * tn(steady)), sin(w * tn(steady))];
coef  = A \ Un(steady);
U0    = coef(1);
theta = atan2(coef(3), coef(2));
tau   = mod(theta, 2 * pi) / w;                % задержка огибающей, с
mMeas = hypot(coef(2), coef(3)) / U0;          % измеренный коэффициент модуляции

Uref  = U0 * (1 + m * cos(w * (tn - tau)));    % модулирующее напряжение в масштабе U[n]
err   = Un(steady) - Uref(steady);
delta = sqrt(mean(err.^2)) / sqrt(mean(Uref(steady).^2)) * 100;

fprintf('Отсчётов I[n], Q[n]: %d, шаг Ts = %.2f нс (%d на период Fм)\n', numel(tn), Ts * 1e9, Nsp);
fprintf('Коэффициент передачи U0 = %.4f В\n', U0);
fprintf('Коэффициент модуляции: задан %.2f, измерен %.3f\n', m, mMeas);
fprintf('Задержка огибающей tau = %.1f нс\n', tau * 1e9);
fprintf('СКО отклонения U[n] от модулирующего напряжения: %.2f %%\n', delta);

%% Графики
tus = tn * 1e6;
fig = figure('Name', 'IQ-демодуляция', 'Position', [100 100 900 650]);
if exist('theme', 'file')           % MATLAB R2025a+: светлый фон для отчёта
    theme(fig, 'light');
end

subplot(2, 1, 1);
stem(tus, In, 'r', 'filled', 'MarkerSize', 3); hold on;
stem(tus, Qn, 'b', 'filled', 'MarkerSize', 3);
grid on; xlabel('t, мкс'); ylabel('U, В');
title(sprintf('Отсчёты I[n], Q[n], Ts = %.1f нс', Ts * 1e9));
legend('I[n]', 'Q[n]', 'Location', 'best');

subplot(2, 1, 2);
plot(tus, Uref, 'k--', 'LineWidth', 1.2); hold on;
plot(tus, Un, 'mo-', 'MarkerSize', 4, 'LineWidth', 1);
xline(tSkip * 1e6, ':', 'установление');
grid on; xlabel('t, мкс'); ylabel('U, В');
title(sprintf('Огибающая U[n] = sqrt(I^2 + Q^2) и модулирующее напряжение (СКО %.1f %%)', delta));
legend('U_0(1 + m cos 2\piF_м(t - \tau))', 'U[n]', 'Location', 'best');

saveas(fig, [outName '.png']);

%% Сохранение отсчётов
fid = fopen([outName '.csv'], 'w');
if fid < 0
    error('Не удалось создать файл %s.csv', outName);
end
fprintf(fid, 't;I;Q;U;Uref\n');
fprintf(fid, '%.6e;%.6e;%.6e;%.6e;%.6e\n', [tn, In, Qn, Un, Uref].');
fclose(fid);
fprintf('Сохранено: %s.png, %s.csv\n', outName, outName);

%% Чтение CSV, экспортированного из Qucs-S
function [t, y] = readQucsCsv(fname)
    if ~isfile(fname)
        error('Файл не найден: %s', fname);
    end
    lines = regexp(fileread(fname), '\r?\n', 'split');
    data  = zeros(numel(lines), 2);
    n = 0;
    for i = 1:numel(lines)
        s = strtrim(lines{i});
        if contains(s, ';')
            parts = strsplit(strrep(s, ',', '.'), ';');   % ";" и десятичная запятая
        else
            parts = regexp(s, '[,\t ]+', 'split');
        end
        v = str2double(parts);
        if numel(v) >= 2 && ~any(isnan(v(1:2)))
            n = n + 1;
            data(n, :) = v(1:2);
        end
    end
    if n < 2
        error('В файле %s не найдено числовых данных (время; значение).', fname);
    end
    [t, idx] = unique(data(1:n, 1));
    y = data(idx, 2);
end
