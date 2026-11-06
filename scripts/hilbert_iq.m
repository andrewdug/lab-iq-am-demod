%% Фильтр Гильберта: задержка, аналитический сигнал, НЧ квадратуры I/Q, АМ-демодуляция (п. 5-8)
% Вход: LR2k5_N1_hilbert_TR.dat - результат моделирования схемы LR2k5_N1_hilbert_TR.sch в Qucs-S
% (файл появляется рядом со схемой после F2; экспорт графиков в CSV не нужен).
%   x0   - входной АМ-сигнал s(t);
%   xd   - s'(t) = s(t - 5Dt), xd14 - s(t - 7Dt), xd18 - s(t - 9Dt)  (линии задержки);
%   y, yw, y14, y18 - сопряжённый по Гильберту сигнал: порядок 10 без окна, порядки 10, 14 и 18 с окном Хэмминга.
% MATLAB R2018b и новее (Signal Processing Toolbox не нужен).

clear; close all;

%% Исходные данные
fileDat = 'LR2k5_N1_hilbert_TR.dat';
N       = 1;                        % номер бригады
f0      = (400 + 30 * N) * 1e6;     % несущая, Гц
Fm      = 3 * N * 1e6;              % частота модуляции, Гц
U0      = 1;                        % амплитуда несущей, В
m       = 1;                        % коэффициент модуляции
fs      = 4 * f0;                   % частота дискретизации (f0 = fs/4)
Dt      = 1 / fs;
tAn     = [1e-6, 2e-6];             % окно спектрального анализа: 3 периода Fм = 1720 отсчётов (шаг ДПФ 1 МГц)
tSkip   = 0.05e-6;                  % отбросить заполнение линии задержки, с
outName = 'hilbert_iq';             % префикс PNG и CSV

variants = struct( ...
    'name',  {'порядок 10, без окна', 'порядок 10, окно Хэмминга', 'порядок 14, окно Хэмминга', 'порядок 18, окно Хэмминга'}, ...
    'ref',   {'xd.Vt', 'xd.Vt', 'xd14.Vt', 'xd18.Vt'}, ...
    'out',   {'y.Vt', 'yw.Vt', 'y14.Vt', 'y18.Vt'}, ...
    'order', {10, 10, 14, 18});
iMain = 3;                          % вариант для графиков (итоговый фильтр после выбора порядка)

%% Чтение результатов Qucs-S: текущая папка, иначе ../qucs рядом со скриптом
if ~isfile(fileDat)
    alt = fullfile(fileparts(mfilename('fullpath')), '..', 'qucs', fileDat);
    if isfile(alt)
        fileDat = alt;
    end
end
fprintf('Данные: %s\n', fileDat);
ds = readQucsDat(fileDat);
t  = ds('time');
s  = ds('x0.Vt');
if t(end) < tAn(2) - Dt / 2
    error('Длительность моделирования %.3g с меньше конца окна анализа %.3g с.', t(end), tAn(2));
end

%% Отсчёты с шагом Dt: n = 1, 2, ...;  t_n = (n - 1) Dt
tn   = (0 : floor(t(end) / Dt + 1e-9)).' * Dt;
iAn  = find(tn >= tAn(1) - Dt / 2 & tn < tAn(2) - Dt / 2);
fAn  = (-numel(iAn) / 2 : numel(iAn) / 2 - 1).' * fs / numel(iAn);
sn   = interp1(t, s, tn, 'spline');
Hs   = idealHilbert(t, s);                % эталон H{s}(t) для измерения задержки

Uall = zeros(numel(tn), numel(variants));
res = repmat(struct('name', '', 't0', 0, 't0Meas', 0, 'K', 0, 'dPhi', 0, ...
                    'imgZ', 0, 'imgIQ', 0, 'errU', 0), 1, numel(variants));
for v = 1 : numel(variants)
    V  = variants(v);
    t0 = V.order / 2 * Dt;                % теория: (n0 - 1) Dt
    steady = t >= tSkip;
    [t0Meas, K] = measureDelay(ds(V.out), Hs, steady, 20 * Dt);

    sd = interp1(t, ds(V.ref), tn, 'spline');   % s'(n)
    sh = interp1(t, ds(V.out), tn, 'spline');   % s^(n)
    z  = sd + 1i * sh;                          % аналитический сигнал z(n)

    w0t = 2 * pi * f0 * tn;                     % w0 (n - 1) Dt
    I   = sd .* cos(w0t) + sh .* sin(w0t);
    Q   = -sd .* sin(w0t) + sh .* cos(w0t);
    U   = sqrt(I.^2 + Q.^2);
    Uenv = U0 * (1 + m * sin(2 * pi * Fm * (tn - t0)));

    Sz  = spectrum(z(iAn));
    Siq = spectrum(I(iAn) + 1i * Q(iAn));
    ok  = tn >= tSkip;
    res(v) = struct('name', V.name, 't0', t0, 't0Meas', t0Meas, 'K', K, ...
        'dPhi', phaseAt(sh(iAn), sd(iAn), tn(iAn), f0), ...
        'imgZ', db(at(Sz, fAn, -f0) / at(Sz, fAn, f0)), ...
        'imgIQ', db(at(Siq, fAn, -fs / 2) / at(Siq, fAn, 0)), ...
        'errU', rms(U(ok) - Uenv(ok)) / rms(Uenv(ok)) * 100);

    if v == iMain
        M = struct('V', V, 't0', t0, 'sd', sd, 'sh', sh, 'z', z, 'I', I, 'Q', Q, 'U', U, 'Uenv', Uenv);
    end
    Uall(:, v) = U - Uenv;
end

%% Результаты
fprintf('fs = %.0f МГц, Dt = %.4f нс, w0*Dt = %.1f град\n', fs / 1e6, Dt * 1e9, 360 * f0 * Dt);
fprintf('%-28s %8s %8s %8s %9s %9s %9s %8s\n', 'Фильтр', 't0 теор', 't0 изм', '|K|', 'фаза,°', 'z(-f0),дБ', 'I/Q обр.', 'СКО U');
for v = 1 : numel(res)
    r = res(v);
    fprintf('%-28s %6.3f нс %6.3f нс %8.4f %9.2f %9.1f %9.1f %6.2f %%\n', r.name, r.t0 * 1e9, ...
        r.t0Meas * 1e9, r.K, r.dPhi, r.imgZ, r.imgIQ, r.errU);
end

%% Графики для итогового фильтра
plotSignals(t, ds, M, tn, Dt, outName);
plotHodograph(t, ds, M, tn, tAn, outName);
plotSpectra(sn(iAn), M, iAn, fAn, f0, outName);
plotIQ(M, tn, iAn, fAn, f0, Fm, outName);
plotDemod(M, tn, Uall, {variants.name}, tSkip, outName);

fid = fopen([outName '.csv'], 'w');
if fid < 0
    error('Не удалось создать файл %s.csv', outName);
end
fprintf(fid, 't;s;s_d;s_h;I;Q;U;Uenv\n');
fprintf(fid, '%.6e;%.6e;%.6e;%.6e;%.6e;%.6e;%.6e;%.6e\n', [tn, sn, M.sd, M.sh, M.I, M.Q, M.U, M.Uenv].');
fclose(fid);
fprintf('Сохранено: %s_*.png, %s.csv (%s)\n', outName, outName, M.V.name);

%% ---------- обработка ----------
function Hs = idealHilbert(t, s)
    % Эталонное преобразование Гильберта по БПФ. На интервале моделирования 2 мкс
    % укладывается целое число периодов f0 и Fм, поэтому сигнал периодичен и БПФ точно.
    n  = numel(t) - 1;                       % последняя точка = первая следующего периода
    X  = fft(s(1 : n));
    k  = [0 : ceil(n / 2) - 1, -floor(n / 2) : -1].';
    Hs = struct('X', X .* (-1i * sign(k)), 'f', k / (t(end) - t(1)));
end

function [tau, K] = measureDelay(y, Hs, steady, tauMax)
    % Задержка tau и коэффициент передачи K: y(t) ~ K * H{s}(t - tau) по МНК;
    % сначала перебор по сетке (шаг Dt/12), затем уточнение fminbnd
    y    = y(steady);
    cost = @(tau) fitCost(tau, y, Hs, steady);
    grid = linspace(0, tauMax, 241);
    [~, i] = min(arrayfun(cost, grid));
    tau  = fminbnd(cost, grid(max(i - 1, 1)), grid(min(i + 1, end)), optimset('TolX', 1e-15));
    [~, K] = fitCost(tau, y, Hs, steady);
end

function [c, K] = fitCost(tau, y, Hs, steady)
    h = real(ifft(Hs.X .* exp(-2i * pi * Hs.f * tau)));
    h = [h; h(1)];
    h = h(steady);
    K = (h.' * y) / (h.' * h);
    c = sum((y - K * h).^2);
end

function S = spectrum(x)
    S = fftshift(abs(fft(x))) / numel(x);
end

function a = at(S, f, f1)
    [~, i] = min(abs(f - f1));
    a = S(i);
end

function d = db(x)
    d = 20 * log10(max(x, 1e-12));
end

function r = rms(x)
    r = sqrt(mean(x.^2));
end

function dphi = phaseAt(a, b, tn, f)
    % разность фаз a и b на частоте f, град
    e = exp(-2i * pi * f * tn);
    dphi = angle(sum(a .* e) / sum(b .* e)) * 180 / pi;
end

%% ---------- графики ----------
function fig = newFig(name)
    fig = figure('Name', name, 'Position', [100 100 950 650], 'Color', 'w');
end

function plotSignals(t, ds, M, tn, Dt, outName)
    fig = newFig('Осциллограммы');
    win = [1.000e-6, 1.010e-6];
    k   = t >= win(1) & t <= win(2);
    kn  = tn >= win(1) & tn <= win(2);
    s = ds('x0.Vt'); sd = ds(M.V.ref); sh = ds(M.V.out);
    subplot(2, 1, 1);
    plot(t(k) * 1e9, s(k), 'Color', [0 .6 0]); hold on;
    plot(t(k) * 1e9, sd(k), 'k', 'LineWidth', 1.2);
    plot(t(k) * 1e9, sh(k), 'b--', 'LineWidth', 1.2);
    plot(tn(kn) * 1e9, M.sd(kn), 'ko', tn(kn) * 1e9, M.sh(kn), 'bs', 'MarkerSize', 4);
    grid on; xlabel('t, нс'); ylabel('U, В');
    title(sprintf('Фильтр Гильберта (%s): t_0 = %.3f нс = %d\\Deltat', M.V.name, M.t0 * 1e9, round(M.t0 / Dt)));
    legend('s(t)', 's''(t) = s(t - t_0)', 'ŝ(t)', 's''(n)', 'ŝ(n)', 'Location', 'eastoutside');
    subplot(2, 1, 2);
    plot(t * 1e6, s, 'Color', [0 .6 0]);
    grid on; xlabel('t, мкс'); ylabel('U, В'); title('Входной АМ-сигнал s(t): f_0 = 430 МГц, F_м = 3 МГц, m = 1');
    saveas(fig, [outName '_signals.png']);
end

function plotHodograph(t, ds, M, tn, tAn, outName)
    fig = newFig('Годограф z(n)');
    k  = t >= tAn(1) & t <= tAn(1) + (tAn(2) - tAn(1)) / 3;     % один период Fм
    kn = tn >= tAn(1) & tn < tAn(2);
    sd = ds(M.V.ref); sh = ds(M.V.out);
    plot(sd(k), sh(k), 'Color', [.75 .75 .75]); hold on;
    plot(real(M.z(kn)), imag(M.z(kn)), 'b.', 'MarkerSize', 8);
    a = linspace(0, 2 * pi, 361);
    plot(2 * cos(a), 2 * sin(a), 'k:');
    axis equal; grid on; xlabel('s''(n), В'); ylabel('ŝ(n), В');
    title({'Годограф аналитического сигнала z(n) = s''(n) + j ŝ(n)', ...
           'серый - непрерывный z(t) за период Fм; f_д = 4f_0: фаза растёт на 90° за отсчёт'});
    legend('z(t)', 'z(n)', '|z| = U_0(1 + m)', 'Location', 'eastoutside');
    saveas(fig, [outName '_hodograph.png']);
end

function plotSpectra(sAn, M, iAn, fAn, f0, outName)
    fig = newFig('Спектры s и z');
    Ss = spectrum(sAn);
    Sz = spectrum(M.z(iAn));
    ref = max(Sz);
    subplot(2, 1, 1);
    stem(fAn / 1e6, db(Ss / ref), 'BaseValue', -100, 'Marker', 'none', 'Color', [0 .6 0]);
    grid on; ylim([-100 5]); xlim([-1.1 1.1] * f0 / 1e6 * 2);
    xlabel('f, МГц'); ylabel('дБ'); title('Амплитудный спектр исходного сигнала s(n): составляющие на \pmf_0, \pm(f_0 \pm F_м)');
    subplot(2, 1, 2);
    stem(fAn / 1e6, db(Sz / ref), 'BaseValue', -100, 'Marker', 'none', 'Color', 'b');
    grid on; ylim([-100 5]); xlim([-1.1 1.1] * f0 / 1e6 * 2);
    xlabel('f, МГц'); ylabel('дБ');
    title(sprintf('Спектр аналитического сигнала z(n): отрицательные частоты подавлены на %.1f дБ', ...
        -db(at(Sz, fAn, -f0) / at(Sz, fAn, f0))));
    saveas(fig, [outName '_spectra.png']);
end

function plotIQ(M, tn, iAn, fAn, f0, Fm, outName)
    fig = newFig('НЧ квадратуры I, Q');
    B = -1i * exp(-2i * pi * f0 * M.t0) * M.Uenv;       % теория: z e^{-jw0 t} = -j U(t - t0) e^{-jw0 t0}
    subplot(2, 1, 1);
    plot(tn * 1e6, M.I, 'r', tn * 1e6, M.Q, 'b'); hold on;
    plot(tn * 1e6, real(B), 'k:', tn * 1e6, imag(B), 'k--');
    grid on; xlabel('t, мкс'); ylabel('U, В'); xlim([0 tn(end) * 1e6]);
    title('Осциллограммы НЧ квадратурных составляющих I(n), Q(n)');
    legend('I(n)', 'Q(n)', 'I теор.', 'Q теор.', 'Location', 'eastoutside');
    subplot(2, 1, 2);
    SI = spectrum(M.I(iAn)); SQ = spectrum(M.Q(iAn)); ref = max([SI; SQ]);
    stem(fAn / 1e6, db(SI / ref), 'BaseValue', -100, 'Marker', 'none', 'Color', 'r'); hold on;
    stem(fAn / 1e6, db(SQ / ref), 'BaseValue', -100, 'Marker', 'none', 'Color', 'b');
    grid on; ylim([-100 5]); xlabel('f, МГц'); ylabel('дБ');
    title(sprintf('Спектры I и Q: составляющие 0, \\pm%g МГц; max|Q(n)| = %.1g В; на \\pmf_д/2 - след |K(f_0)| \\neq 1', ...
        Fm / 1e6, max(abs(M.Q(iAn)))));
    legend('|I(f)|', '|Q(f)|', 'Location', 'eastoutside');
    saveas(fig, [outName '_iq.png']);
end

function plotDemod(M, tn, Uall, names, tSkip, outName)
    fig = newFig('АМ-демодуляция');
    subplot(2, 1, 1);
    plot(tn * 1e6, M.Uenv, 'k--', 'LineWidth', 1.2); hold on;
    plot(tn * 1e6, M.U, 'm');
    xline(tSkip * 1e6, ':', 'заполнение');
    grid on; xlabel('t, мкс'); ylabel('U, В'); xlim([0 tn(end) * 1e6]);
    title(sprintf('Огибающая U(n) = sqrt(I^2 + Q^2) и модулирующее U_0(1 + m sin 2\\piF_м(t - t_0)), %s', M.V.name));
    legend('модулирующее', 'U(n)', 'Location', 'eastoutside');
    subplot(2, 1, 2);
    plot(tn * 1e6, Uall);
    grid on; xlabel('t, мкс'); ylabel('\DeltaU, В'); xlim([tSkip, tn(end)] * 1e6);
    title('Ошибка демодуляции U(n) - модулирующее (чередование знака от отсчёта к отсчёту: |K(f_0)| \neq 1)');
    legend(names, 'Location', 'eastoutside');
    saveas(fig, [outName '_demod.png']);
end

%% ---------- чтение .dat Qucs ----------
function ds = readQucsDat(fname)
    % Набор данных Qucs: блоки <indep имя N> ... </indep> и <dep имя зависимость> ... </dep>
    if ~isfile(fname)
        error('Файл не найден: %s (запустите моделирование LR2k5_N1_hilbert_TR.sch в Qucs-S)', fname);
    end
    txt = fileread(fname);
    tok = regexp(txt, '<(?:indep|dep) (\S+)[^>]*>(.*?)</(?:indep|dep)>', 'tokens');
    if isempty(tok)
        error('В %s нет блоков данных Qucs.', fname);
    end
    ds = containers.Map();
    for i = 1 : numel(tok)
        body = tok{i}{2};
        if any(body == 'j')
            continue;                        % комплексные величины (AC) здесь не нужны
        end
        ds(tok{i}{1}) = sscanf(body, '%f');
    end
    need = {'time', 'x0.Vt', 'xd.Vt', 'xd14.Vt', 'xd18.Vt', 'y.Vt', 'yw.Vt', 'y14.Vt', 'y18.Vt'};
    miss = need(~isKey(ds, need));
    if ~isempty(miss)
        error('В %s нет переменных: %s', fname, strjoin(miss, ', '));
    end
end
