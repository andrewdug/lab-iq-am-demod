%HILBERT_SL_INIT Параметры модели hilbert_iq_sl (запускается автоматически при открытии и старте).
%  Вариант фильтра задаётся переменной iVar: 1 - порядок 10 без окна, 2 - порядок 10 с окном
%  Хэмминга, 3 - порядок 14 с окном Хэмминга (итоговый), 4 - порядок 18 с окном Хэмминга.

M   = 1;                      % номер бригады
f0  = (400 + 30 * M) * 1e6;   % несущая, Гц
Fm  = 3 * M * 1e6;            % частота модуляции, Гц
U0  = 1;                      % амплитуда несущей, В
m   = 1;                      % коэффициент модуляции
fs  = 4 * f0;                 % частота дискретизации, Гц
Dt  = 1 / fs;                 % интервал дискретизации, с
w0  = 2 * pi * f0;            % несущая, рад/с

if ~exist('iVar', 'var')
    iVar = 3;
end

hilbertData = load(fullfile(fileparts(get_param(bdroot, 'FileName')), 'qucs_data.mat'));
t0     = hilbertData.orders(iVar) / 2 * Dt;      % задержка фильтра (n0 - 1) Δt
Tend   = hilbertData.tn(end);
sdIn   = [hilbertData.tn, hilbertData.sd(:, iVar)];   % s'(n)
shIn   = [hilbertData.tn, hilbertData.sh(:, iVar)];   % ŝ(n)
sIn    = [hilbertData.tn, hilbertData.s];             % s(n)
