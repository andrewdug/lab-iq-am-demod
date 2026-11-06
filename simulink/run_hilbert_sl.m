function run_hilbert_sl()
%RUN_HILBERT_SL Запускает hilbert_iq_sl для всех вариантов фильтра, строит графики п. 6-8 и ГВЗ.
%  Рисунки sl_*.png и sl_results.json сохраняются в report/figures.

here = fileparts(mfilename('fullpath'));
addpath(here);
outDir = fullfile(here, '..', 'report', 'figures');
d = load(fullfile(here, 'qucs_data.mat'));
mdl = 'hilbert_iq_sl';
load_system(mdl);

f0 = 430e6; Fm = 3e6; fs = 4 * f0; Dt = 1 / fs;
iMain = 3;
tSkip = 0.05e-6;
ian = d.tn >= 1e-6 - Dt / 2 & d.tn < 2e-6 - Dt / 2;       % 1720 отсчётов, шаг ДПФ 1 МГц
na  = nnz(ian);
fan = (-na / 2 : na / 2 - 1).' * fs / na;

results = struct([]);
for v = 1 : numel(d.names)
    assignin('base', 'iVar', v);
    out = sim(mdl);
    tt = out.get('tOut'); I = out.get('IOut'); Q = out.get('QOut');
    U  = out.get('UOut'); Ur = out.get('UrefOut'); z = squeeze(out.get('zOut'));
    if norm(tt - d.tn) > 1e-12
        error('Сетка времени модели не совпала с отсчётами Qucs-S');
    end
    ok = tt >= tSkip;
    Sz = spec(z(ian)); Si = spec(I(ian) + 1i * Q(ian));
    [tg, fg] = groupDelay(d.f, d.Hf(:, v));
    tgF0 = interp1(fg, tg, f0);
    r = struct('name', d.names{v}, 'order', d.orders(v), ...
        'image_z_db', dB(at(Sz, fan, -f0) / at(Sz, fan, f0)), ...
        'image_iq_db', dB(at(Si, fan, -fs / 2) / at(Si, fan, 0)), ...
        'u_err_rms_percent', rmsv(U(ok) - Ur(ok)) / rmsv(Ur(ok)) * 100, ...
        'q_max', max(abs(Q(ian))), ...
        'tg_f0_ns', tgF0 * 1e9, 'tg_ripple_ps', (max(tg) - min(tg)) * 1e12, ...
        'tg_theory_ns', d.orders(v) / 2 * Dt * 1e9);
    results = [results, r]; %#ok<AGROW>
    if v == iMain
        M = struct('t', tt, 'z', z, 'I', I, 'Q', Q, 'U', U, 'Ur', Ur, 'r', r);
    end
end
close_system(mdl, 0);

plotHodograph(d, M, ian, outDir);
plotSpectra(d, M, ian, fan, f0, outDir);
plotIQ(M, ian, fan, outDir);
plotDemod(M, outDir);
plotGroupDelay(d, Dt, f0, outDir);

fid = fopen(fullfile(outDir, 'sl_results.json'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s', jsonencode(results, 'PrettyPrint', true));
fclose(fid);
for r = results
    fprintf('%-24s z(-f0) %6.1f дБ | I/Q(fs/2) %6.1f дБ | СКО U %5.2f %% | max|Q| %.1e | ГВЗ %.3f нс (теор. %.3f, разброс %.2e пс)\n', ...
        r.name, r.image_z_db, r.image_iq_db, r.u_err_rms_percent, r.q_max, r.tg_f0_ns, r.tg_theory_ns, r.tg_ripple_ps);
end
end

function [tg, fg] = groupDelay(f, H)
% ГВЗ: t_гр = -dФ/dw по ФЧХ Qucs-S; края, где |K| близок к нулю, отбрасываются
keep = f >= 5e6 & f <= 855e6;
ph = unwrap(angle(H(keep)));
fg = f(keep);
tg = -gradient(ph, 2 * pi * fg);
fg = fg(2 : end - 1);
tg = tg(2 : end - 1);
end

function S = spec(x)
S = fftshift(abs(fft(x))) / numel(x);
end

function a = at(S, f, f1)
[~, i] = min(abs(f - f1));
a = S(i);
end

function y = dB(x)
y = 20 * log10(max(x, 1e-12));
end

function r = rmsv(x)
r = sqrt(mean(x .^ 2));
end

function saveFig(fig, outDir, name)
exportgraphics(fig, fullfile(outDir, name), 'Resolution', 170);
close(fig);
end

function fig = newFig(w, h)
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 w h], 'Color', 'w');
end

function plotHodograph(d, M, ian, outDir)
fig = newFig(21, 10);
tiledlayout(fig, 1, 2, 'Padding', 'compact');
nexttile;
plot(real(M.z(ian)), imag(M.z(ian)), '.', 'Color', [1 .4 0], 'MarkerSize', 5); hold on;
a = linspace(0, 2 * pi, 361); plot(2 * cos(a), 2 * sin(a), 'k:');
axis equal; grid on; xlim([-2.2 2.2]); ylim([-2.2 2.2]);
xlabel('s''(n), В'); ylabel('ŝ(n), В'); title('Годограф z(n) (Simulink)');
nexttile;
plot(M.I(ian), M.Q(ian), '.', 'Color', [1 .4 0], 'MarkerSize', 5);
axis equal; grid on; xlim([-2.2 2.2]); ylim([-2.2 2.2]);
xlabel('I(n), В'); ylabel('Q(n), В'); title('Годограф I-Q (Simulink)');
saveFig(fig, outDir, 'sl_hodograph.png');
end

function plotSpectra(d, M, ian, fan, f0, outDir)
fig = newFig(21, 13);
tiledlayout(fig, 2, 1, 'Padding', 'compact');
Ss = spec(d.s(ian)); Sz = spec(M.z(ian)); ref = max(Sz);
nexttile;
stem(fan / 1e6, max(dB(Ss / ref), -100), 'BaseValue', -100, 'Marker', 'none', 'Color', [.17 .63 .17]);
grid on; ylim([-100 5]); xlim([-860 860]); xlabel('f, МГц'); ylabel('дБ'); title('Спектр s(n)');
nexttile;
stem(fan / 1e6, max(dB(Sz / ref), -100), 'BaseValue', -100, 'Marker', 'none', 'Color', [1 .4 0]);
grid on; ylim([-100 5]); xlim([-860 860]); xlabel('f, МГц'); ylabel('дБ');
title(sprintf('Спектр z(n): подавление -f_0 на %.1f дБ', -dB(at(Sz, fan, -f0) / at(Sz, fan, f0))));
saveFig(fig, outDir, 'sl_spectra.png');
end

function plotIQ(M, ian, fan, outDir)
fig = newFig(21, 13);
tiledlayout(fig, 2, 1, 'Padding', 'compact');
nexttile;
plot(M.t * 1e6, M.I, 'r', M.t * 1e6, M.Q, 'b'); grid on;
xlabel('t, мкс'); ylabel('U, В'); legend('I(n)', 'Q(n)'); title('НЧ квадратуры (Simulink)');
nexttile;
Si = spec(M.I(ian)); ref = max(Si);
stem(fan / 1e6, max(dB(Si / ref), -100), 'BaseValue', -100, 'Marker', 'none', 'Color', 'r');
grid on; xlim([-15 15]); ylim([-100 5]); xlabel('f, МГц'); ylabel('дБ'); title('Спектр I(n) вблизи нуля');
saveFig(fig, outDir, 'sl_iq.png');
end

function plotDemod(M, outDir)
fig = newFig(21, 9);
plot(M.t * 1e6, M.Ur, 'k--', 'LineWidth', 1.3); hold on;
plot(M.t * 1e6, M.U, 'Color', [1 .4 0]);
grid on; xlabel('t, мкс'); ylabel('U, В'); legend('модулирующее', 'U(n)');
title('Выход АМ-демодулятора (Simulink)');
saveFig(fig, outDir, 'sl_demod.png');
end

function plotGroupDelay(d, Dt, f0, outDir)
colors = [.84 .15 .16; .12 .47 .71; 1 .4 0; .17 .63 .17];
fig = newFig(21, 10);
hold on;
for v = 1 : numel(d.names)
    [tg, fg] = groupDelay(d.f, d.Hf(:, v));
    plot(fg / 1e6, tg * 1e9, 'Color', colors(v, :), 'LineWidth', 1.4, 'DisplayName', d.names{v});
end
xline(f0 / 1e6, ':', 'f_0', 'HandleVisibility', 'off');
grid on; xlim([0 860]); ylim([0 7]);
xlabel('f, МГц'); ylabel('t_{гр}, нс');
title('ГВЗ фильтров Гильберта: t_{гр} = -dФ/d\omega по ФЧХ Qucs-S');
legend('Location', 'southeast');
saveFig(fig, outDir, 'sl_gvz.png');
end

