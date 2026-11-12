%% Фильтр Гильберта: задержка, аналитический сигнал, НЧ квадратуры I/Q, АМ-демодуляция (п. 5-8)
% Скрипт сам формирует z(n) = s'(n) + j s^(n), I(n), Q(n) и U(n) из сигналов, полученных в Qucs-S.
% Два варианта входных данных (переменная inputMode ниже):
%  'dat' - результаты отдельных схем, по одной на фильтр (файл .dat появляется рядом со схемой после F2):
%            LR2k5_N1_hilbert_R10_TR.dat - порядок 10, без окна;
%            LR2k5_N1_hilbert_W10_TR.dat - порядок 10, окно Хэмминга;
%            LR2k5_N1_hilbert_W14_TR.dat - порядок 14, окно Хэмминга;
%            LR2k5_N1_hilbert_W18_TR.dat - порядок 18, окно Хэмминга.
%          В каждом: x0 - входной АМ-сигнал s(t), xd - s'(t) = s(t - t0), y - s^(t) этого фильтра.
%          Обрабатываются те файлы, что найдены (хотя бы один).
%  'csv' - три CSV, экспортированные из графиков Qucs-S (правый клик по графику -> Export to CSV),
%          каждый "время; напряжение": s(t), s'(t) и s^(t) одного фильтра (любой схемы, в т.ч.
%          LR2k5_N1_hilbert_W18_TR.sch: x0.Vt, xd.Vt, y.Vt). Порядок фильтра задаётся в csvOrder.
% MATLAB R2018b и новее (Signal Processing Toolbox не нужен).

clear; close all;

%% Исходные данные
datFiles = {'LR2k5_N1_hilbert_R10_TR.dat', 'LR2k5_N1_hilbert_W10_TR.dat', 'LR2k5_N1_hilbert_W14_TR.dat', ...
            'LR2k5_N1_hilbert_W18_TR.dat'};
datNames = {'порядок 10, без окна', 'порядок 10, окно Хэмминга', 'порядок 14, окно Хэмминга', ...
            'порядок 18, окно Хэмминга'};
datOrder = [10 10 14 18];
N       = 1;                        % номер бригады
f0      = (400 + 30 * N) * 1e6;     % несущая, Гц
Fm      = 3 * N * 1e6;              % частота модуляции, Гц
U0      = 1;                        % амплитуда несущей, В
m       = 1;                        % коэффициент модуляции
usDiv   = 1;                        % субдискретизация: fs = 4 f0 / usDiv (1 - без неё; 43 - схемы
                                    % LR2k5_N1_hilbert_{R10,W10,W18}_D43_TR.sch, fs = 40 МГц)
fs      = 4 * f0 / usDiv;           % частота дискретизации (f0 = fs/4 или её образ)
Dt      = 1 / fs;
fa      = f0 - round(f0 / fs) * fs; % где несущая оказывается после дискретизации: +fs/4 или -fs/4
tAn     = [1e-6, 2e-6];             % окно спектрального анализа: 3 периода Fм (шаг ДПФ 1 МГц)
tSkip   = max(0.05e-6, 20 * Dt);    % отбросить заполнение линии задержки, с
outName = 'hilbert_iq';             % префикс PNG и CSV
ordList = 6 : 4 : 30;              % порядки для сравнения (цифровая модель фильтра, без Qucs);
                                    % L/2 нечётное, как у 10 и 18: при чётном L/2 крайние отводы нулевые
afcOrd  = [6 10 18 30];             % порядки, для которых рисуются АЧХ
if usDiv ~= 1
    datFiles = strrep(datFiles, '_TR.dat', sprintf('_D%d_TR.dat', usDiv));
    outName  = sprintf('hilbert_iq_D%d', usDiv);
end

%% Источник данных: 'dat' или 'csv' (см. шапку)
inputMode = 'dat';
csvS      = 's.csv';                % s(t)  - входной АМ-сигнал (x0.Vt)
csvSd     = 'sd.csv';               % s'(t) - выход линии задержки (xd.Vt / xd18.Vt)
csvSh     = 'sh.csv';               % s^(t) - выход фильтра Гильберта (y.Vt / yw.Vt / y18.Vt)
csvOrder  = 18;                     % порядок фильтра, с которого сняты s' и s^ (10 или 18)
csvName   = 'порядок 18, окно Хэмминга';   % подпись на графиках

%% Чтение результатов Qucs-S (файлы ищутся в текущей папке, иначе в ../qucs рядом со скриптом)
switch inputMode
    case 'dat'
        [ds, variants] = readFilterSet(datFiles, datNames, datOrder);
    case 'csv'
        variants = struct('name', csvName, 'ref', 'sd', 'out', 'sh', 'order', csvOrder);
        ds = readCsvSet(findInput(csvS), findInput(csvSd), findInput(csvSh));
    otherwise
        error('inputMode должен быть ''dat'' или ''csv''.');
end
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

R = struct([]);                           % результаты и сигналы всех фильтров
for v = 1 : numel(variants)
    V  = variants(v);
    t0 = V.order / 2 * Dt;                % теория: (n0 - 1) Dt
    steady = t >= tSkip;
    [t0Meas, K] = measureDelay(ds(V.out), Hs, steady, 20 * Dt, f0);

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
    r = struct('V', V, 'name', V.name, 't0', t0, 't0Meas', t0Meas, 'K', K, ...
        'dPhi', phaseAt(sh(iAn), sd(iAn), tn(iAn), f0), ...
        'imgZ', db(at(Sz, fAn, -fa) / at(Sz, fAn, fa)), ...
        'imgIQ', db(at(Siq, fAn, -fs / 2) / at(Siq, fAn, 0)), ...
        'errU', rms(U(ok) - Uenv(ok)) / rms(Uenv(ok)) * 100, ...
        'sd', sd, 'sh', sh, 'z', z, 'I', I, 'Q', Q, 'U', U, 'Uenv', Uenv);
    if isempty(R)
        R = r;
    else
        R(v) = r;
    end
end
M = R(end);                                  % итоговый фильтр (обычно порядок 18)

%% Результаты
fprintf('fs = %.4g МГц, Dt = %.4f нс, w0*Dt = %.1f град, несущая после дискретизации %+.4g МГц\n', ...
    fs / 1e6, Dt * 1e9, mod(360 * f0 * Dt, 360), fa / 1e6);
fprintf('%-28s %8s %8s %8s %9s %9s %9s %8s\n', 'Фильтр', 't0 теор', 't0 изм', '|K|', 'фаза,°', 'z(-f0),дБ', 'I/Q обр.', 'СКО U');
for v = 1 : numel(R)
    r = R(v);
    fprintf('%-28s %6.3f нс %6.3f нс %8.4f %9.2f %9.1f %9.1f %6.2f %%\n', r.name, r.t0 * 1e9, ...
        r.t0Meas * 1e9, r.K, r.dPhi, r.imgZ, r.imgIQ, r.errU);
end

%% Графики: на каждом - все фильтры для сравнения
% plotSignals(t, ds, M, tn, Dt, outName);
plotHodograph(R, tn, tAn, sprintf('f_д = %.4g МГц: фаза меняется на %+d° за отсчёт', fs / 1e6, ...
    round(mod(360 * f0 * Dt + 180, 360) - 180)), outName);
plotSpectra(sn(iAn), R, iAn, fAn, fa, fs, outName);
plotIQ(R, tn, iAn, fAn, Fm, outName);
plotDemod(R, tn, tSkip, outName);

%% Сравнение порядков фильтра: цифровая модель с теми же коэффициентами, что в схемах Qucs
% h(k) = 2/(pi k) при нечётном k = n - n0 (k = -L/2 ... L/2), с окном Хэмминга
% w = 0,54 - 0,46 cos(2 pi (k + L/2) / L); при спектре, инвертированном субдискретизацией, знак меняется.
% На вход подаются те же отсчёты s(n), что и в Qucs, поэтому для порядков 10 и 18 модель
% должна совпасть с результатами схем (они отмечены на графиках).
sgnH = sign(fa);
O = struct('L', {}, 'win', {}, 'K', {}, 'imgZ', {}, 'errU', {});
for win = [0 1]
    for L = ordList
        [k, h] = hilbertCoef(L, win, sgnH);
        sh = filter(h, 1, sn);                          % s^(n), задержан на L/2 отсчётов
        sd = [zeros(L / 2, 1); sn(1 : end - L / 2)];    % s'(n) = s(n - L/2)
        z  = sd + 1i * sh;
        w0t = 2 * pi * f0 * tn;
        I  = sd .* cos(w0t) + sh .* sin(w0t);
        Q  = -sd .* sin(w0t) + sh .* cos(w0t);
        U  = sqrt(I.^2 + Q.^2);
        Uenv = U0 * (1 + m * sin(2 * pi * Fm * (tn - L / 2 * Dt)));
        ok = tn >= max(tSkip, (L + 1) * Dt);
        Sz = spectrum(z(iAn));
        O(end + 1) = struct('L', L, 'win', win, ...
            'K', abs(sum(h .* exp(-2i * pi * fa * k * Dt))), ...   % |K| на частоте несущей
            'imgZ', db(at(Sz, fAn, -fa) / at(Sz, fAn, fa)), ...
            'errU', rms(U(ok) - Uenv(ok)) / rms(Uenv(ok)) * 100); %#ok<SAGROW>
    end
end
fprintf('\nСравнение порядков (цифровая модель):\n%-8s %-14s %9s %11s %8s\n', 'порядок', 'окно', '|K(f0)|', 'z(-f0), дБ', 'СКО U');
for j = 1 : numel(O)
    o = O(j);
    fprintf('%-8d %-14s %9.4f %11.1f %6.2f %%\n', o.L, winName(o.win), o.K, o.imgZ, o.errU);
end
plotOrderAFC(afcOrd, sgnH, fs, fa, Fm, outName);
plotOrderMetrics(O, R, outName);

fid = fopen([outName '.csv'], 'w');
if fid < 0
    error('Не удалось создать файл %s.csv', outName);
end
% столбцы: t; s; затем для каждого фильтра k = 1, 2, ...: s_d, s_h, I, Q, U, Uenv
hdr = {'t', 's'};
dat = [tn, sn];
for v = 1 : numel(R)
    hdr = [hdr, strcat({'s_d', 's_h', 'I', 'Q', 'U', 'Uenv'}, sprintf('_%d', v))]; %#ok<AGROW>
    dat = [dat, R(v).sd, R(v).sh, R(v).I, R(v).Q, R(v).U, R(v).Uenv];         %#ok<AGROW>
    fprintf('  столбцы *_%d в CSV - %s\n', v, R(v).name);
end
fprintf(fid, '%s\n', strjoin(hdr, ';'));
fprintf(fid, [strjoin(repmat({'%.6e'}, 1, size(dat, 2)), ';') '\n'], dat.');
fclose(fid);
fprintf('Сохранено: %s_{hodograph,spectra,iq,demod,orders,orders_afc}.png, %s.csv\n', outName, outName);

%% ---------- обработка ----------
function Hs = idealHilbert(t, s)
    % Эталонное преобразование Гильберта по БПФ. На интервале моделирования 2 мкс
    % укладывается целое число периодов f0 и Fм, поэтому сигнал периодичен и БПФ точно.
    n  = numel(t) - 1;                       % последняя точка = первая следующего периода
    X  = fft(s(1 : n));
    k  = [0 : ceil(n / 2) - 1, -floor(n / 2) : -1].';
    Hs = struct('X', X .* (-1i * sign(k)), 'f', k / (t(end) - t(1)));
end

function [tau, K] = measureDelay(y, Hs, steady, tauMax, f0)
    % Задержка tau и коэффициент передачи K: y(t) ~ K * H{s}(t - tau) по МНК;
    % сначала перебор по сетке (шаг не больше 1/20 периода несущей - иначе при большом Dt
    % минимум можно спутать с соседним через полпериода), затем уточнение fminbnd
    y    = y(steady);
    cost = @(tau) fitCost(tau, y, Hs, steady);
    grid = linspace(0, tauMax, max(241, ceil(tauMax * f0 * 20) + 1));
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

function [k, h] = hilbertCoef(L, win, sgn)
    % коэффициенты КИХ-фильтра Гильберта порядка L (L + 1 отвод), k = -L/2 ... L/2
    k = (-L / 2 : L / 2).';
    h = zeros(size(k));
    odd = mod(k, 2) ~= 0;
    h(odd) = sgn * 2 ./ (pi * k(odd));
    if win
        h = h .* (0.54 - 0.46 * cos(2 * pi * (k + L / 2) / L));
    end
end

function s = winName(win)
    if win
        s = 'Хэмминг';
    else
        s = 'без окна';
    end
end

%% ---------- графики ----------
function fig = newFig(name)
    fig = figure('Name', name, 'Position', [100 100 950 650], 'Color', 'w');
end

% function plotSignals(t, ds, M, tn, Dt, outName)
%     fig = newFig('Осциллограммы');
%     win = [1.000e-6, 1.010e-6];
%     k   = t >= win(1) & t <= win(2);
%     kn  = tn >= win(1) & tn <= win(2);
%     s = ds('x0.Vt'); sd = ds(M.V.ref); sh = ds(M.V.out);
%     subplot(2, 1, 1);
%     plot(t(k) * 1e9, s(k), 'Color', [0 .6 0]); hold on;
%     plot(t(k) * 1e9, sd(k), 'k', 'LineWidth', 1.2);
%     plot(t(k) * 1e9, sh(k), 'b--', 'LineWidth', 1.2);
%     plot(tn(kn) * 1e9, M.sd(kn), 'ko', tn(kn) * 1e9, M.sh(kn), 'bs', 'MarkerSize', 4);
%     grid on; xlabel('t, нс'); ylabel('U, В');
%     title(sprintf('Фильтр Гильберта (%s): t_0 = %.3f нс = %d\\Deltat', M.V.name, M.t0 * 1e9, round(M.t0 / Dt)));
%     legend('s(t)', 's''(t) = s(t - t_0)', 'ŝ(t)', 's''(n)', 'ŝ(n)', 'Location', 'eastoutside');
%     subplot(2, 1, 2);
%     plot(t * 1e6, s, 'Color', [0 .6 0]);
%     grid on; xlabel('t, мкс'); ylabel('U, В'); title('Входной АМ-сигнал s(t): f_0 = 430 МГц, F_м = 3 МГц, m = 1');
%     saveas(fig, [outName '_signals.png']);
% end

function c = filtColors(n)
    c = [0 .45 .74; .85 .33 .10; .47 .67 .19; .49 .18 .56];
    c = c(1 + mod(0 : n - 1, size(c, 1)), :);
end

function plotHodograph(R, tn, tAn, subtitle, outName)
    % годографы z(n) всех фильтров рядом, в одном масштабе
    nV  = numel(R);
    fig = figure('Name', 'Годограф z(n)', 'Position', [100 100 420 * nV + 100 520], 'Color', 'w');
    kn  = tn >= tAn(1) & tn < tAn(2);
    a   = linspace(0, 2 * pi, 361);
    c   = filtColors(nV);
    for v = 1 : nV
        r = R(v);
        subplot(1, nV, v);
        plot(real(r.z(kn)), imag(r.z(kn)), '.', 'Color', c(v, :), 'MarkerSize', 8); hold on;
        plot(2 * cos(a), 2 * sin(a), 'k:');
        axis equal; grid on; xlim([-2.4 2.4]); ylim([-2.4 2.4]);
        xlabel('s''(n), В'); ylabel('ŝ(n), В');
        title({r.name, sprintf('|K(f_0)| = %.4f, \\Delta\\phi = %.1f°', r.K, r.dPhi)});
    end
    addSuperTitle({'Годограф аналитического сигнала z(n) = s''(n) + j ŝ(n); пунктир - |z| = U_0(1 + m)', subtitle});
    saveas(fig, [outName '_hodograph.png']);
end

function plotSpectra(sAn, R, iAn, fAn, fa, fs, outName)
    % вверху спектр s(n), ниже спектры z(n) всех фильтров
    % fa - частота несущей после дискретизации: f0 = fs/4 без субдискретизации, иначе её образ +-fs/4
    nV  = numel(R);
    fig = figure('Name', 'Спектры s и z', 'Position', [100 50 950 230 * (nV + 1)], 'Color', 'w');
    c   = filtColors(nV);
    Ss  = spectrum(sAn);
    subplot(nV + 1, 1, 1);
    stem(fAn / 1e6, db(Ss / max(Ss)), 'BaseValue', -100, 'Marker', 'none', 'Color', [0 .6 0]);
    grid on; ylim([-100 5]); xlim([-1.1 1.1] * fs / 2e6); ylabel('дБ');
    title(sprintf('Амплитудный спектр s(n): несущая на \\pm%.4g МГц и боковые \\pmF_м', abs(fa) / 1e6));
    for v = 1 : nV
        Sz = spectrum(R(v).z(iAn));
        subplot(nV + 1, 1, v + 1);
        stem(fAn / 1e6, db(Sz / max(Sz)), 'BaseValue', -100, 'Marker', 'none', 'Color', c(v, :));
        grid on; ylim([-100 5]); xlim([-1.1 1.1] * fs / 2e6); ylabel('дБ');
        title(sprintf('Спектр z(n), %s: несущая на %+.4g МГц, зеркальная подавлена на %.1f дБ', ...
            R(v).name, fa / 1e6, -R(v).imgZ));
    end
    xlabel('f, МГц');
    saveas(fig, [outName '_spectra.png']);
end

function plotIQ(R, tn, iAn, fAn, Fm, outName)
    % строка на фильтр: слева осциллограммы I(n), Q(n), справа их спектры
    nV  = numel(R);
    fig = figure('Name', 'НЧ квадратуры I, Q', 'Position', [100 50 1300 260 * nV], 'Color', 'w');
    lim = [0, 0];                            % общий масштаб по вертикали
    for v = 1 : nV
        lim = [min(lim(1), min([R(v).I; R(v).Q])), max(lim(2), max([R(v).I; R(v).Q]))];
    end
    for v = 1 : nV
        r = R(v);
        subplot(nV, 2, 2 * v - 1);
        plot(tn * 1e6, r.I, 'r', tn * 1e6, r.Q, 'b');
        grid on; ylabel('U, В'); xlim([0 tn(end) * 1e6]); ylim(lim + [-0.1 0.1]);
        title(sprintf('I(n), Q(n), %s', r.name));
        legend('I(n)', 'Q(n)', 'Location', 'eastoutside');
        subplot(nV, 2, 2 * v);
        SI = spectrum(r.I(iAn)); SQ = spectrum(r.Q(iAn)); ref = max([SI; SQ]);
        stem(fAn / 1e6, db(SI / ref), 'BaseValue', -100, 'Marker', 'none', 'Color', 'r'); hold on;
        stem(fAn / 1e6, db(SQ / ref), 'BaseValue', -100, 'Marker', 'none', 'Color', 'b');
        grid on; ylim([-100 5]); ylabel('дБ');
        title(sprintf('Спектры I, Q: 0, \\pm%g МГц; max|Q| = %.1g В; на \\pmf_д/2: %.1f дБ', ...
            Fm / 1e6, max(abs(r.Q(iAn))), r.imgIQ));
        legend('|I(f)|', '|Q(f)|', 'Location', 'eastoutside');
    end
    subplot(nV, 2, 2 * nV - 1); xlabel('t, мкс');
    subplot(nV, 2, 2 * nV);     xlabel('f, МГц');
    saveas(fig, [outName '_iq.png']);
end

function plotDemod(R, tn, tSkip, outName)
    % строка на фильтр: U(n) и модулирующее; внизу - ошибка демодуляции всех фильтров вместе
    nV  = numel(R);
    fig = figure('Name', 'АМ-демодуляция', 'Position', [100 50 1100 220 * (nV + 1)], 'Color', 'w');
    c   = filtColors(nV);
    for v = 1 : nV
        r = R(v);
        subplot(nV + 1, 1, v);
        plot(tn * 1e6, r.Uenv, 'k--', 'LineWidth', 1.2); hold on;
        plot(tn * 1e6, r.U, 'Color', c(v, :));
        grid on; ylabel('U, В'); xlim([0 tn(end) * 1e6]); ylim([0 2.4]);
        title(sprintf('U(n) = sqrt(I^2 + Q^2), %s: СКО от модулирующего %.2f %%', r.name, r.errU));
        legend('U_0(1 + m sin 2\piF_м(t - t_0))', 'U(n)', 'Location', 'eastoutside');
    end
    subplot(nV + 1, 1, nV + 1);
    ok = tn >= tSkip;
    eMax = 0;
    names = cell(1, nV);
    for v = 1 : nV
        e = R(v).U - R(v).Uenv;
        plot(tn * 1e6, e, 'Color', c(v, :)); hold on;
        eMax = max(eMax, max(abs(e(ok))));
        names{v} = R(v).name;
    end
    yl = 1.15 * eMax * [-1 1];               % масштаб по установившемуся режиму
    plot(tSkip * 1e6 * [1 1], yl, 'k:');
    grid on; xlim([0 tn(end) * 1e6]); ylim(yl);
    xlabel('t, мкс'); ylabel('\DeltaU, В');
    title('Ошибка демодуляции U(n) - U_0(1 + m sin 2\piF_м(t - t_0)) всех фильтров');
    legend([names, {'заполнение'}], 'Location', 'eastoutside');
    saveas(fig, [outName '_demod.png']);
end

function plotOrderAFC(ords, sgn, fs, fa, Fm, outName)
    % АЧХ фильтра Гильберта |K(f)| для разных порядков: без окна и с окном Хэмминга
    fig = figure('Name', 'АЧХ для разных порядков', 'Position', [100 80 1000 700], 'Color', 'w');
    f   = linspace(0, fs / 2, 2001).';
    c   = filtColors(numel(ords));
    for win = [0 1]
        subplot(2, 1, win + 1);
        names = cell(1, numel(ords));
        for i = 1 : numel(ords)
            [k, h] = hilbertCoef(ords(i), win, sgn);
            K = abs(exp(-2i * pi * f * k.' / fs) * h);
            plot(f / 1e6, K, 'Color', c(i, :), 'LineWidth', 1.2); hold on;
            names{i} = sprintf('порядок %d', ords(i));
        end
        yl = [0 1.25];
        plot(abs(fa) / 1e6 * [1 1], yl, 'k--');
        plot((abs(fa) + Fm * [-1 1]) / 1e6, [1 1], 'k', 'LineWidth', 3);
        grid on; ylim(yl); xlim([0 fs / 2e6]);
        xlabel('f, МГц'); ylabel('|K(f)|');
        title(sprintf('АЧХ фильтра Гильберта, %s (f_д = %.4g МГц)', winName(win), fs / 1e6));
        legend([names, {'несущая', 'полоса АМ-сигнала'}], 'Location', 'eastoutside');
    end
    saveas(fig, [outName '_orders_afc.png']);
end

function plotOrderMetrics(O, R, outName)
    % показатели качества в зависимости от порядка; квадраты - результаты схем Qucs-S
    fig = figure('Name', 'Сравнение порядков', 'Position', [100 80 1000 800], 'Color', 'w');
    mk  = {'o-', 's-'};
    cl  = [0 .45 .74; .85 .33 .10];
    qW  = arrayfun(@(r) double(~isempty(strfind(r.name, 'Хэмминг'))), R);
    qL  = arrayfun(@(r) r.V.order, R);
    fld = {'imgZ', 'errU', 'K'};
    ttl = {'Подавление зеркальной составляющей z(n)', 'СКО огибающей U(n) от модулирующего', ...
           'Отклонение |K(f_0)| от 1'};
    yl  = {'подавление, дБ', 'СКО, %', '||K(f_0)| - 1|, %'};
    for j = 1 : 3
        subplot(3, 1, j);
        for win = [0 1]
            o = O([O.win] == win);
            y = [o.(fld{j})];
            q = [R(qW == win).(fld{j})];
            switch fld{j}
                case 'imgZ', y = -y;              q = -q;
                case 'K',    y = abs(y - 1) * 100; q = abs(q - 1) * 100;
            end
            if j == 1
                plot([o.L], y, mk{win + 1}, 'Color', cl(win + 1, :), 'LineWidth', 1.2); hold on;
            else
                semilogy([o.L], y, mk{win + 1}, 'Color', cl(win + 1, :), 'LineWidth', 1.2); hold on;
            end
            if ~isempty(q)
                plot(qL(qW == win), q, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', cl(win + 1, :));
            end
        end
        grid on; xlim([min([O.L]) - 2, max([O.L]) + 2]); xticks(unique([O.L]));
        ylabel(yl{j}); title(ttl{j});
        if j == 1
            legend('без окна (модель)', 'без окна (Qucs-S)', 'Хэмминг (модель)', 'Хэмминг (Qucs-S)', ...
                   'Location', 'eastoutside');
        else
            legend('без окна', 'Qucs-S', 'Хэмминг', 'Qucs-S', 'Location', 'eastoutside');
        end
    end
    xlabel('порядок фильтра L');
    saveas(fig, [outName '_orders.png']);
end

function addSuperTitle(txt)
    if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
        sgtitle(txt);
    end
end

%% ---------- чтение .dat Qucs ----------
function [ds, variants] = readFilterSet(files, names, orders)
    % Результаты отдельных схем (по одному фильтру): из каждой берутся xd (s') и y (s^);
    % s(t) и время - из первой найденной, остальные при необходимости пересчитываются на её сетку.
    ds = containers.Map();
    variants = struct('name', {}, 'ref', {}, 'out', {}, 'order', {});
    for k = 1:numel(files)
        f = findInput(files{k});
        if ~isfile(f)
            fprintf('Нет файла %s - фильтр "%s" пропущен.\n', files{k}, names{k});
            continue;
        end
        fprintf('Данные: %s\n', f);
        d = readQucsDat(f);
        if ~isKey(ds, 'time')
            ds('time')  = d('time');
            ds('x0.Vt') = d('x0.Vt');
        end
        tag = sprintf('f%d', k);
        ds([tag '.xd']) = onGrid(d('time'), d('xd.Vt'), ds('time'));
        ds([tag '.y'])  = onGrid(d('time'), d('y.Vt'),  ds('time'));
        variants(end + 1) = struct('name', names{k}, 'ref', [tag '.xd'], 'out', [tag '.y'], 'order', orders(k)); %#ok<AGROW>
    end
    if isempty(variants)
        error(['Не найден ни один файл: %s.\nОткройте схемы LR2k5_N1_hilbert_{R10,W10,W18}_TR.sch ' ...
               'в Qucs-S и запустите моделирование (F2).'], strjoin(files, ', '));
    end
end

function v = onGrid(t, x, tRef)
    if numel(t) == numel(tRef) && max(abs(t - tRef)) < 1e-15
        v = x;
    else
        v = interp1(t, x, tRef, 'spline');
    end
end

function ds = readQucsDat(fname)
    % Набор данных Qucs: блоки <indep имя N> ... </indep> и <dep имя зависимость> ... </dep>
    if ~isfile(fname)
        error('Файл не найден: %s (запустите моделирование схемы в Qucs-S)', fname);
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
    need = {'time', 'x0.Vt', 'xd.Vt', 'y.Vt'};
    miss = need(~isKey(ds, need));
    if ~isempty(miss)
        error('В %s нет переменных: %s', fname, strjoin(miss, ', '));
    end
end

%% ---------- чтение CSV, экспортированных из Qucs-S ----------
function f = findInput(f)
    % файл в текущей папке, иначе в ../qucs рядом со скриптом
    if ~isfile(f)
        alt = fullfile(fileparts(mfilename('fullpath')), '..', 'qucs', f);
        if isfile(alt)
            f = alt;
        end
    end
end

function ds = readCsvSet(fS, fSd, fSh)
    % Три CSV "время; напряжение" -> общий равномерный отсчёт времени (как у TR в Qucs).
    [t, s]   = readQucsCsv(fS);
    [t1, sd] = readQucsCsv(fSd);
    [t2, sh] = readQucsCsv(fSh);
    fprintf('Данные: %s, %s, %s\n', fS, fSd, fSh);
    tu = linspace(t(1), t(end), numel(t)).';          % БПФ-эталон H{s} требует равномерной сетки
    if max(abs(t - tu)) > 1e-3 * (t(2) - t(1))
        warning('Шаг времени в %s неравномерный: сигналы пересчитаны на равномерную сетку.', fS);
    end
    if t1(1) > tu(1) || t1(end) < tu(end) || t2(1) > tu(1) || t2(end) < tu(end)
        error('Интервалы времени в CSV не совпадают: экспортируйте графики целиком (0...2 мкс).');
    end
    ds = containers.Map();
    ds('time')  = tu;
    ds('x0.Vt') = interp1(t, s, tu, 'spline');
    ds('sd')    = interp1(t1, sd, tu, 'spline');
    ds('sh')    = interp1(t2, sh, tu, 'spline');
end

function [t, y] = readQucsCsv(fname)
    % Qucs-S: заголовок с именами, затем "время;значение" (разделитель ; , или табуляция,
    % десятичная точка или запятая при разделителе ;). Берутся первые два числовых столбца.
    if ~isfile(fname)
        error('Файл не найден: %s (экспортируйте график из Qucs-S: правый клик -> Export to CSV)', fname);
    end
    lines = regexp(fileread(fname), '\r?\n', 'split');
    data  = zeros(numel(lines), 2);
    n = 0;
    for i = 1:numel(lines)
        s = strtrim(lines{i});
        if contains(s, ';')
            parts = strsplit(strrep(s, ',', '.'), ';');
        else
            parts = regexp(s, '[,\t ]+', 'split');
        end
        v = str2double(strrep(parts, '"', ''));
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
