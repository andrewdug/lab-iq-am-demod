%% Фильтр Гильберта: обработка результатов Qucs-S и сравнение вариантов фильтра (п. 5-8)
% Скрипт читает результаты моделирования схем LR2k5_N1_hilbert_<вариант>[_D43]_{TR,AC}.sch
% (файлы .dat появляются рядом со схемой после F2) и сравнивает все найденные варианты:
%   R10 - порядок 10, без окна;      W10 - порядок 10, окно Хэмминга;
%   W14 - порядок 14, окно Хэмминга; W18 - порядок 18, окно Хэмминга.
% Две группы: fд = 4 f0 = 1720 МГц (без суффикса) и субдискретизация fд = 4 f0 / 43 = 40 МГц (_D43).
% Из TR берутся x0 - входной АМ-сигнал s(t), xd - s'(t) = s(t - t0), y - выход фильтра s^(t);
% из AC - АЧХ K, разность фаз dF с элементом задержки и значения на f0 (K_f0, t0_ns).
% Варианты без .dat пропускаются. Файлы ищутся в текущей папке, иначе в ../qucs.
% MATLAB R2018b и новее (Signal Processing Toolbox не нужен).

clear; close all;

%% Исходные данные
N       = 1;                        % номер бригады
f0      = (400 + 30 * N) * 1e6;     % несущая, Гц
Fm      = 3 * N * 1e6;              % частота модуляции, Гц
U0      = 1;                        % амплитуда несущей, В
m       = 1;                        % коэффициент модуляции
divList = [1 43];                   % fд = 4 f0 / div: 1 - 1720 МГц, 43 - 40 МГц (субдискретизация)
tAn     = [1e-6, 2e-6];             % окно спектрального анализа: 3 периода Fм (шаг ДПФ 1 МГц)
filt    = struct('tag',   {'R10', 'W10', 'W14', 'W18'}, ...
                 'name',  {'порядок 10, без окна', 'порядок 10, Хэмминг', ...
                           'порядок 14, Хэмминг', 'порядок 18, Хэмминг'}, ...
                 'order', {10, 10, 14, 18});

%% Обработка обеих групп
S = struct([]);                                  % сводка по всем вариантам
for div = divList
    P = struct('f0', f0, 'Fm', Fm, 'U0', U0, 'm', m, 'fs', 4 * f0 / div, 'tAn', tAn);
    P.Dt    = 1 / P.fs;
    P.fa    = f0 - round(f0 / P.fs) * P.fs;      % несущая после дискретизации: +fs/4 или -fs/4
    P.tSkip = max(0.05e-6, 20 * P.Dt);           % отбросить заполнение линии задержки, с
    sfx = '';
    if div ~= 1
        sfx = sprintf('_D%d', div);
    end
    outName = ['hilbert_iq' sfx];

    [ds, V] = readFilterSet(filt, sfx);
    if isempty(V)
        fprintf('\nfд = %.4g МГц: результатов TR нет (схемы LR2k5_N1_hilbert_*%s_TR.sch не промоделированы).\n', ...
            P.fs / 1e6, sfx);
        continue;
    end
    R = processFilters(ds, V, P);
    A = readAfcSet(V, sfx);

    fprintf('\nfд = %.4g МГц, Dt = %.4f нс, w0*Dt = %.1f град, несущая после дискретизации %+.4g МГц\n', ...
        P.fs / 1e6, P.Dt * 1e9, mod(360 * f0 * P.Dt, 360), P.fa / 1e6);
    fprintf('%-24s %9s %9s %9s %8s %8s %10s %9s %8s\n', 'Фильтр', 't0 теор', 't0 TR', 't0 AC', ...
        '|K| TR', '|K| AC', 'фаза TR,°', 'z(-f0),дБ', 'СКО U');
    for v = 1 : numel(R)
        r = R(v); a = A(v);
        fprintf('%-24s %6.3f нс %6.3f нс %6.3f нс %8.4f %8.4f %10.2f %9.1f %6.2f %%\n', r.name, ...
            r.t0 * 1e9, r.t0Meas * 1e9, a.t0 * 1e9, r.K, a.K0, r.dPhi, r.imgZ, r.errU);
        S = [S, struct('fs', P.fs, 'name', r.name, 't0', r.t0, 't0TR', r.t0Meas, 't0AC', a.t0, ...
            'KTR', r.K, 'KAC', a.K0, 'dPhi', r.dPhi, 'imgZ', r.imgZ, 'errU', r.errU)]; %#ok<AGROW>
    end

    plotAFC(A, R, P, outName);
    plotHodograph(R, P, outName);
    plotSpectra(ds, R, P, outName);
    plotIQ(R, P, outName);
    plotDemod(R, outName);
    fprintf('Сохранено: %s_{afc,hodograph,spectra,iq,demod}.png\n', outName);
end
if isempty(S)
    error(['Не найдено ни одного результата TR. Откройте схемы LR2k5_N1_hilbert_{R10,W10,W14,W18}_TR.sch ' ...
           'в Qucs-S и запустите моделирование (F2).']);
end
writeSummary(S, 'hilbert_iq_summary.csv');

%% ---------- обработка ----------
function R = processFilters(ds, V, P)
    % для каждого фильтра: задержка, |K|, z(n), I(n), Q(n), U(n) и показатели качества
    t   = ds('time');
    tn  = (0 : floor(t(end) / P.Dt + 1e-9)).' * P.Dt;   % отсчёты с шагом Dt
    iAn = find(tn >= P.tAn(1) - P.Dt / 2 & tn < P.tAn(2) - P.Dt / 2);
    fAn = (-numel(iAn) / 2 : numel(iAn) / 2 - 1).' * P.fs / numel(iAn);
    if t(end) < P.tAn(2) - P.Dt / 2
        error('Длительность моделирования %.3g с меньше конца окна анализа %.3g с.', t(end), P.tAn(2));
    end
    Hs = idealHilbert(t, ds('x0.Vt'));                  % эталон H{s}(t) для измерения задержки
    R  = struct([]);
    for v = 1 : numel(V)
        t0 = V(v).order / 2 * P.Dt;                     % теория: (n0 - 1) Dt
        % поиск до L*Dt (теория L/2*Dt): больший диапазон при fд = 40 МГц превышает период
        % огибающей 1/Fм, и подгонка может уйти на соседний период
        [t0Meas, K] = measureDelay(ds(V(v).out), Hs, t >= P.tSkip, V(v).order * P.Dt, P.f0);
        sd  = interp1(t, ds(V(v).ref), tn, 'spline');   % s'(n)
        sh  = interp1(t, ds(V(v).out), tn, 'spline');   % s^(n)
        z   = sd + 1i * sh;                             % аналитический сигнал z(n)
        w0t = 2 * pi * P.f0 * tn;
        I   = sd .* cos(w0t) + sh .* sin(w0t);
        Q   = -sd .* sin(w0t) + sh .* cos(w0t);
        U   = sqrt(I.^2 + Q.^2);
        Uenv = P.U0 * (1 + P.m * sin(2 * pi * P.Fm * (tn - t0)));
        Sz  = spectrum(z(iAn));
        Siq = spectrum(I(iAn) + 1i * Q(iAn));
        ok  = tn >= P.tSkip;
        r = struct('name', V(v).name, 'tag', V(v).tag, 't0', t0, 't0Meas', t0Meas, 'K', K, ...
            'dPhi', phaseAt(sh(iAn), sd(iAn), tn(iAn), P.f0), ...
            'imgZ', db(at(Sz, fAn, -P.fa) / at(Sz, fAn, P.fa)), ...
            'imgIQ', db(at(Siq, fAn, -P.fs / 2) / at(Siq, fAn, 0)), ...
            'errU', rms(U(ok) - Uenv(ok)) / rms(Uenv(ok)) * 100, ...
            'tn', tn, 'iAn', iAn, 'fAn', fAn, 'sd', sd, 'sh', sh, 'z', z, 'I', I, 'Q', Q, 'U', U, 'Uenv', Uenv);
        R = [R, r]; %#ok<AGROW>
    end
end

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
    if K <= 0
        c = Inf;                             % фильтр Гильберта не инвертирует сигнал
    end
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

function writeSummary(S, fname)
    fid = fopen(fname, 'w');
    if fid < 0
        error('Не удалось создать файл %s', fname);
    end
    fprintf(fid, 'fs_MHz;filter;t0_theory_ns;t0_TR_ns;t0_AC_ns;K_TR;K_AC;dPhi_deg;image_dB;errU_pct\n');
    for s = S
        fprintf(fid, '%.4g;%s;%.4f;%.4f;%.4f;%.5f;%.5f;%.2f;%.1f;%.3f\n', s.fs / 1e6, s.name, ...
            s.t0 * 1e9, s.t0TR * 1e9, s.t0AC * 1e9, s.KTR, s.KAC, s.dPhi, s.imgZ, s.errU);
    end
    fclose(fid);
    fprintf('\nСводка по всем вариантам: %s\n', fname);
end

%% ---------- графики: на каждом все фильтры группы ----------
function fig = newFig(varargin)
    fig = figure(varargin{:});
    if exist('theme', 'file')               % MATLAB R2025a+: светлый фон для отчёта
        theme(fig, 'light');
    end
end

function c = filtColors(n)
    c = [0 .45 .74; .85 .33 .10; .47 .67 .19; .49 .18 .56];
    c = c(1 + mod(0 : n - 1, size(c, 1)), :);
end

function plotAFC(A, R, P, outName)
    % АЧХ и разность фаз с элементом задержки по результатам AC Qucs-S
    has = arrayfun(@(a) ~isempty(a.f), A);
    if ~any(has)
        return;
    end
    fig   = newFig('Name', 'АЧХ и ФЧХ (Qucs-S, AC)', 'Position', [100 50 1100 850], 'Color', 'w');
    c     = filtColors(numel(A));
    names = {R(has).name};
    band  = (P.f0 + P.Fm * [-1 1]) / 1e6;
    for j = 1 : 3
        subplot(3, 1, j);
        for v = find(has)
            if j < 3
                y = A(v).K;
            else
                y = A(v).dF;
                y(A(v).K < 1e-3) = NaN;          % при |K| = 0 (f = 0, fд/2) фаза не определена
            end
            plot(A(v).f / 1e6, y, 'Color', c(v, :), 'LineWidth', 1.2); hold on;
        end
        grid on;
        switch j
            case 1
                yl = [0 1.25]; ylabel('|K(f)|');
                title(sprintf('АЧХ фильтров Гильберта (f_д = %.4g МГц); пунктир - несущая f_0', P.fs / 1e6));
            case 2
                xlim(P.f0 / 1e6 + [-60 60]); yl = ylim; ylabel('|K(f)|');
                plot(band, [1 1], 'k', 'LineWidth', 3);
                title('АЧХ вблизи f_0; жирная линия - полоса АМ-сигнала f_0 \pm F_м');
            case 3
                yl = [-180 180]; ylabel('\Delta\phi, град'); xlabel('f, МГц');
                title('Разность фаз выхода фильтра и элемента задержки t_0');
        end
        plot(P.f0 / 1e6 * [1 1], yl, 'k--'); ylim(yl);
        legend(names, 'Location', 'eastoutside');
    end
    saveas(fig, [outName '_afc.png']);
end

function plotHodograph(R, P, outName)
    % годографы z(n) всех фильтров рядом, в одном масштабе
    nV  = numel(R);
    fig = newFig('Name', 'Годограф z(n)', 'Position', [100 100 420 * nV + 100 520], 'Color', 'w');
    a   = linspace(0, 2 * pi, 361);
    c   = filtColors(nV);
    for v = 1 : nV
        r  = R(v);
        kn = r.tn >= P.tAn(1) & r.tn < P.tAn(2);
        subplot(1, nV, v);
        plot(real(r.z(kn)), imag(r.z(kn)), '.', 'Color', c(v, :), 'MarkerSize', 8); hold on;
        plot(2 * cos(a), 2 * sin(a), 'k:');
        axis equal; grid on; xlim([-2.4 2.4]); ylim([-2.4 2.4]);
        xlabel('s''(n), В'); ylabel('ŝ(n), В');
        title({r.name, sprintf('|K(f_0)| = %.4f, \\Delta\\phi = %.1f°', r.K, r.dPhi)});
    end
    addSuperTitle({'Годограф аналитического сигнала z(n) = s''(n) + j ŝ(n); пунктир - |z| = U_0(1 + m)', ...
        sprintf('f_д = %.4g МГц: фаза меняется на %+d° за отсчёт', P.fs / 1e6, ...
        round(mod(360 * P.f0 * P.Dt + 180, 360) - 180))});
    saveas(fig, [outName '_hodograph.png']);
end

function plotSpectra(ds, R, P, outName)
    % вверху спектр s(n), ниже спектры z(n) всех фильтров
    nV  = numel(R);
    fig = newFig('Name', 'Спектры s и z', 'Position', [100 50 950 230 * (nV + 1)], 'Color', 'w');
    c   = filtColors(nV);
    r1  = R(1);
    sn  = interp1(ds('time'), ds('x0.Vt'), r1.tn, 'spline');
    Ss  = spectrum(sn(r1.iAn));
    subplot(nV + 1, 1, 1);
    stem(r1.fAn / 1e6, db(Ss / max(Ss)), 'BaseValue', -100, 'Marker', 'none', 'Color', [0 .6 0]);
    grid on; ylim([-100 5]); xlim([-1.1 1.1] * P.fs / 2e6); ylabel('дБ');
    title(sprintf('Амплитудный спектр s(n): несущая на \\pm%.4g МГц и боковые \\pmF_м', abs(P.fa) / 1e6));
    for v = 1 : nV
        Sz = spectrum(R(v).z(R(v).iAn));
        subplot(nV + 1, 1, v + 1);
        stem(R(v).fAn / 1e6, db(Sz / max(Sz)), 'BaseValue', -100, 'Marker', 'none', 'Color', c(v, :));
        grid on; ylim([-100 5]); xlim([-1.1 1.1] * P.fs / 2e6); ylabel('дБ');
        title(sprintf('Спектр z(n), %s: несущая на %+.4g МГц, зеркальная подавлена на %.1f дБ', ...
            R(v).name, P.fa / 1e6, -R(v).imgZ));
    end
    xlabel('f, МГц');
    saveas(fig, [outName '_spectra.png']);
end

function plotIQ(R, P, outName)
    % строка на фильтр: слева осциллограммы I(n), Q(n), справа их спектры
    nV  = numel(R);
    fig = newFig('Name', 'НЧ квадратуры I, Q', 'Position', [100 50 1300 260 * nV], 'Color', 'w');
    lim = [0, 0];                            % общий масштаб по вертикали
    for v = 1 : nV
        lim = [min(lim(1), min([R(v).I; R(v).Q])), max(lim(2), max([R(v).I; R(v).Q]))];
    end
    for v = 1 : nV
        r = R(v);
        subplot(nV, 2, 2 * v - 1);
        plot(r.tn * 1e6, r.I, 'r', r.tn * 1e6, r.Q, 'b');
        grid on; ylabel('U, В'); xlim([0 r.tn(end) * 1e6]); ylim(lim + [-0.1 0.1]);
        title(sprintf('I(n), Q(n), %s', r.name));
        legend('I(n)', 'Q(n)', 'Location', 'eastoutside');
        subplot(nV, 2, 2 * v);
        SI = spectrum(r.I(r.iAn)); SQ = spectrum(r.Q(r.iAn)); ref = max([SI; SQ]);
        stem(r.fAn / 1e6, db(SI / ref), 'BaseValue', -100, 'Marker', 'none', 'Color', 'r'); hold on;
        stem(r.fAn / 1e6, db(SQ / ref), 'BaseValue', -100, 'Marker', 'none', 'Color', 'b');
        grid on; ylim([-100 5]); ylabel('дБ');
        title(sprintf('Спектры I, Q: 0, \\pm%g МГц; max|Q| = %.1g В; на \\pmf_д/2: %.1f дБ', ...
            P.Fm / 1e6, max(abs(r.Q(r.iAn))), r.imgIQ));
        legend('|I(f)|', '|Q(f)|', 'Location', 'eastoutside');
    end
    subplot(nV, 2, 2 * nV - 1); xlabel('t, мкс');
    subplot(nV, 2, 2 * nV);     xlabel('f, МГц');
    saveas(fig, [outName '_iq.png']);
end

function plotDemod(R, outName)
    % строка на фильтр: U(n) и модулирующее
    nV  = numel(R);
    fig = newFig('Name', 'АМ-демодуляция', 'Position', [100 50 1100 220 * nV], 'Color', 'w');
    c   = filtColors(nV);
    for v = 1 : nV
        r = R(v);
        subplot(nV, 1, v);
        plot(r.tn * 1e6, r.Uenv, 'k--', 'LineWidth', 1.2); hold on;
        plot(r.tn * 1e6, r.U, 'Color', c(v, :));
        grid on; ylabel('U, В'); xlim([0 r.tn(end) * 1e6]); ylim([0 2.4]);
        title(sprintf('U(n) = sqrt(I^2 + Q^2), %s: СКО от модулирующего %.2f %%', r.name, r.errU));
        legend('U_0(1 + m sin 2\piF_м(t - t_0))', 'U(n)', 'Location', 'eastoutside');
    end
    xlabel('t, мкс');
    saveas(fig, [outName '_demod.png']);
end

function addSuperTitle(txt)
    if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
        sgtitle(txt);
    end
end

%% ---------- чтение .dat Qucs ----------
function [ds, V] = readFilterSet(filt, sfx)
    % Результаты TR отдельных схем (по одному фильтру): из каждой берутся xd (s') и y (s^);
    % s(t) и время - из первой найденной, остальные при необходимости пересчитываются на её сетку.
    ds = containers.Map();
    V  = struct('tag', {}, 'name', {}, 'order', {}, 'ref', {}, 'out', {});
    for k = 1 : numel(filt)
        f = findInput(sprintf('LR2k5_N1_hilbert_%s%s_TR.dat', filt(k).tag, sfx));
        if ~isfile(f)
            continue;
        end
        fprintf('Данные TR: %s\n', f);
        d = readQucsDat(f, {'time', 'x0.Vt', 'xd.Vt', 'y.Vt'});
        if ~isKey(ds, 'time')
            ds('time')  = d('time');
            ds('x0.Vt') = d('x0.Vt');
        end
        tag = filt(k).tag;
        ds([tag '.xd']) = onGrid(d('time'), d('xd.Vt'), ds('time'));
        ds([tag '.y'])  = onGrid(d('time'), d('y.Vt'),  ds('time'));
        V(end + 1) = struct('tag', tag, 'name', filt(k).name, 'order', filt(k).order, ...
            'ref', [tag '.xd'], 'out', [tag '.y']); %#ok<AGROW>
    end
end

function A = readAfcSet(V, sfx)
    % Результаты AC тех же вариантов: АЧХ K(f), разность фаз dF(f) и значения на f0
    A = struct('f', {}, 'K', {}, 'dF', {}, 'K0', {}, 't0', {});
    for v = 1 : numel(V)
        f = findInput(sprintf('LR2k5_N1_hilbert_%s%s_AC.dat', V(v).tag, sfx));
        a = struct('f', [], 'K', [], 'dF', [], 'K0', NaN, 't0', NaN);
        if isfile(f)
            fprintf('Данные AC: %s\n', f);
            d = readQucsDat(f, {'acfrequency', 'K', 'dF', 'K_f0', 't0_ns'});
            a = struct('f', d('acfrequency'), 'K', d('K'), 'dF', d('dF'), 'K0', d('K_f0'), ...
                't0', d('t0_ns') * 1e-9);
        end
        A(v) = a;
    end
end

function f = findInput(f)
    % файл в текущей папке, иначе в ../qucs рядом со скриптом
    if ~isfile(f)
        alt = fullfile(fileparts(mfilename('fullpath')), '..', 'qucs', f);
        if isfile(alt)
            f = alt;
        end
    end
end

function v = onGrid(t, x, tRef)
    if numel(t) == numel(tRef) && max(abs(t - tRef)) < 1e-15
        v = x;
    else
        v = interp1(t, x, tRef, 'spline');
    end
end

function ds = readQucsDat(fname, need)
    % Набор данных Qucs: блоки <indep имя N> ... </indep> и <dep имя зависимость> ... </dep>;
    % комплексные величины пропускаются (нужные - вещественные: напряжения TR, mag/angle AC)
    txt = fileread(fname);
    tok = regexp(txt, '<(?:indep|dep) (\S+)[^>]*>(.*?)</(?:indep|dep)>', 'tokens');
    ds  = containers.Map();
    for i = 1 : numel(tok)
        body = tok{i}{2};
        if ~any(body == 'j')
            ds(tok{i}{1}) = sscanf(body, '%f');
        end
    end
    miss = need(~isKey(ds, need));
    if ~isempty(miss)
        error('В %s нет переменных: %s', fname, strjoin(miss, ', '));
    end
end
