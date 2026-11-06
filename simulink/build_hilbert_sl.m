function build_hilbert_sl()
%BUILD_HILBERT_SL Собирает модель hilbert_iq_sl.slx: цифровая обработка п. 6-8.
%  Вход - отсчёты s(n), s'(n), ŝ(n) из Qucs-S (qucs_data.mat).

here = fileparts(mfilename('fullpath'));
addpath(here);
mdl = 'hilbert_iq_sl';
if bdIsLoaded(mdl)
    close_system(mdl, 0);
end
new_system(mdl);
% Запуск из папки модели: работает и при открытии .slx двойным щелчком из любой текущей папки
initCmd = 'run(fullfile(fileparts(get_param(bdroot, ''FileName'')), ''hilbert_sl_init.m''))';
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'FixedStepDiscrete', 'FixedStep', 'Dt', ...
    'StopTime', 'Tend', 'PostLoadFcn', initCmd, 'InitFcn', initCmd, ...
    'SaveOutput', 'off', 'SaveTime', 'off', 'SaveState', 'off', 'SignalLogging', 'off');

% --- источники: результаты Qucs-S ---
fromWs(mdl, 's(n) из Qucs-S', 'sIn', [40 40]);
fromWs(mdl, 'sd(n) из Qucs-S', 'sdIn', [40 140]);
fromWs(mdl, 'sh(n) из Qucs-S', 'shIn', [40 240]);
blk(mdl, 'simulink/Sources/Digital Clock', 'Время t(n)', [40 340 80 370], 'SampleTime', 'Dt');

% --- п. 6: аналитический сигнал ---
blk(mdl, 'simulink/Math Operations/Real-Imag to Complex', 'z(n) = sd(n) + j sh(n)', [260 170 320 230]);
add_line(mdl, 'sd(n) из Qucs-S/1', 'z(n) = sd(n) + j sh(n)/1', 'autorouting', 'smart');
add_line(mdl, 'sh(n) из Qucs-S/1', 'z(n) = sd(n) + j sh(n)/2', 'autorouting', 'smart');
sink(mdl, 'zOut', [400 185 460 215], 'z(n) = sd(n) + j sh(n)/1');

% --- п. 7: НЧ квадратуры ---
iq = [mdl '/НЧ квадратуры I, Q'];
add_block('built-in/Subsystem', iq, 'Position', [260 300 380 420]);
buildIQ(iq);
add_line(mdl, 'sd(n) из Qucs-S/1', 'НЧ квадратуры I, Q/1', 'autorouting', 'smart');
add_line(mdl, 'sh(n) из Qucs-S/1', 'НЧ квадратуры I, Q/2', 'autorouting', 'smart');
add_line(mdl, 'Время t(n)/1', 'НЧ квадратуры I, Q/3', 'autorouting', 'smart');

% --- п. 8: демодулятор ---
dm = [mdl '/АМ-демодулятор'];
add_block('built-in/Subsystem', dm, 'Position', [520 300 640 420]);
buildDemod(dm);
add_line(mdl, 'НЧ квадратуры I, Q/1', 'АМ-демодулятор/1', 'autorouting', 'smart');
add_line(mdl, 'НЧ квадратуры I, Q/2', 'АМ-демодулятор/2', 'autorouting', 'smart');
sink(mdl, 'IOut', [520 440 580 470], 'НЧ квадратуры I, Q/1');
sink(mdl, 'QOut', [520 490 580 520], 'НЧ квадратуры I, Q/2');
sink(mdl, 'UOut', [860 300 920 330], 'АМ-демодулятор/1');

% --- модулирующее напряжение для сравнения ---
blk(mdl, 'simulink/User-Defined Functions/Fcn', 'Модулирующее U0(1 + m sin 2pi Fm (t - t0))', [520 540 640 590], ...
    'Expr', 'U0*(1+m*sin(2*pi*Fm*(u(1)-t0)))');
add_line(mdl, 'Время t(n)/1', 'Модулирующее U0(1 + m sin 2pi Fm (t - t0))/1', 'autorouting', 'smart');
sink(mdl, 'UrefOut', [860 550 920 580], 'Модулирующее U0(1 + m sin 2pi Fm (t - t0))/1');
sink(mdl, 'tOut', [40 430 100 460], 'Время t(n)/1');

blk(mdl, 'simulink/Math Operations/Sum', 'Ошибка U - Uмод', [760 420 790 450], 'Inputs', '+-', 'ShowName', 'off');
add_line(mdl, 'АМ-демодулятор/1', 'Ошибка U - Uмод/1', 'autorouting', 'smart');
add_line(mdl, 'Модулирующее U0(1 + m sin 2pi Fm (t - t0))/1', 'Ошибка U - Uмод/2', 'autorouting', 'smart');
sink(mdl, 'dUOut', [860 420 920 450], 'Ошибка U - Uмод/1');

% --- годографы (XY Graph) ---
blk(mdl, 'simulink/Sinks/XY Graph', 'Годограф z(n)', [400 60 480 120]);
add_line(mdl, 'sd(n) из Qucs-S/1', 'Годограф z(n)/1', 'autorouting', 'smart');
add_line(mdl, 'sh(n) из Qucs-S/1', 'Годограф z(n)/2', 'autorouting', 'smart');
blk(mdl, 'simulink/Sinks/XY Graph', 'Годограф I-Q', [760 180 840 240]);
add_line(mdl, 'НЧ квадратуры I, Q/1', 'Годограф I-Q/1', 'autorouting', 'smart');
add_line(mdl, 'НЧ квадратуры I, Q/2', 'Годограф I-Q/2', 'autorouting', 'smart');

% --- осциллографы ---
blk(mdl, 'simulink/Sinks/Scope', 'Осциллограф I, Q', [760 490 800 530], 'NumInputPorts', '2');
add_line(mdl, 'НЧ квадратуры I, Q/1', 'Осциллограф I, Q/1', 'autorouting', 'smart');
add_line(mdl, 'НЧ квадратуры I, Q/2', 'Осциллограф I, Q/2', 'autorouting', 'smart');
blk(mdl, 'simulink/Sinks/Scope', 'Осциллограф U, Uмод', [960 400 1000 440], 'NumInputPorts', '2');
add_line(mdl, 'АМ-демодулятор/1', 'Осциллограф U, Uмод/1', 'autorouting', 'smart');
add_line(mdl, 'Модулирующее U0(1 + m sin 2pi Fm (t - t0))/1', 'Осциллограф U, Uмод/2', 'autorouting', 'smart');
blk(mdl, 'simulink/Sinks/Scope', 'Осциллограф s, sd, sh', [260 40 300 80], 'NumInputPorts', '3');
add_line(mdl, 's(n) из Qucs-S/1', 'Осциллограф s, sd, sh/1', 'autorouting', 'smart');
add_line(mdl, 'sd(n) из Qucs-S/1', 'Осциллограф s, sd, sh/2', 'autorouting', 'smart');
add_line(mdl, 'sh(n) из Qucs-S/1', 'Осциллограф s, sd, sh/3', 'autorouting', 'smart');

% --- спектр z(n) ---
load_system('dspsnks4');
blk(mdl, 'dspsnks4/Spectrum Analyzer', 'Спектр z(n)', [400 228 460 262]);
add_line(mdl, 'z(n) = sd(n) + j sh(n)/1', 'Спектр z(n)/1', 'autorouting', 'smart');

save_system(mdl, fullfile(here, [mdl '.slx']));
close_system(mdl, 0);
fprintf('Сохранено: %s.slx\n', mdl);
end

function fromWs(mdl, name, var, xy)
blk(mdl, 'simulink/Sources/From Workspace', name, [xy, xy + [120 30]], 'VariableName', var, ...
    'SampleTime', 'Dt', 'Interpolate', 'on', 'OutputAfterFinalValue', 'Holding final value');
end

function h = blk(sys, lib, name, pos, varargin)
h = add_block(lib, [sys '/' name], 'Position', pos, varargin{:});
end

function sink(mdl, var, pos, src)
blk(mdl, 'simulink/Sinks/To Workspace', var, pos, 'VariableName', var, 'SaveFormat', 'Array', 'SampleTime', 'Dt');
add_line(mdl, src, [var '/1'], 'autorouting', 'smart');
end

function buildIQ(s)
% I(n) =  s'(n) cos(w0 (n-1) Dt) + ŝ(n) sin(w0 (n-1) Dt)
% Q(n) = -s'(n) sin(w0 (n-1) Dt) + ŝ(n) cos(w0 (n-1) Dt)
clearDefault(s);
blk(s, 'simulink/Sources/In1', 'sd', [30 60 60 74]);
blk(s, 'simulink/Sources/In1', 'sh', [30 220 60 234], 'Port', '2');
blk(s, 'simulink/Sources/In1', 't', [30 380 60 394], 'Port', '3');
blk(s, 'simulink/Math Operations/Gain', 'w0', [110 370 170 404], 'Gain', 'w0');
add_line(s, 't/1', 'w0/1', 'autorouting', 'smart');
blk(s, 'simulink/Math Operations/Trigonometric Function', 'cos', [210 330 260 360], 'Operator', 'cos');
blk(s, 'simulink/Math Operations/Trigonometric Function', 'sin', [210 410 260 440], 'Operator', 'sin');
add_line(s, 'w0/1', 'cos/1', 'autorouting', 'smart');
add_line(s, 'w0/1', 'sin/1', 'autorouting', 'smart');
blk(s, 'simulink/Math Operations/Product', 'sd cos', [330 40 370 80]);
blk(s, 'simulink/Math Operations/Product', 'sh sin', [330 110 370 150]);
blk(s, 'simulink/Math Operations/Product', 'sd sin', [330 200 370 240]);
blk(s, 'simulink/Math Operations/Product', 'sh cos', [330 270 370 310]);
pairs = {'sd/1', 'sd cos/1'; 'cos/1', 'sd cos/2'; 'sh/1', 'sh sin/1'; 'sin/1', 'sh sin/2'; ...
         'sd/1', 'sd sin/1'; 'sin/1', 'sd sin/2'; 'sh/1', 'sh cos/1'; 'cos/1', 'sh cos/2'};
for k = 1 : size(pairs, 1)
    add_line(s, pairs{k, 1}, pairs{k, 2}, 'autorouting', 'smart');
end
blk(s, 'simulink/Math Operations/Sum', 'Сумма I', [440 80 470 110], 'Inputs', '++', 'ShowName', 'off');
blk(s, 'simulink/Math Operations/Sum', 'Сумма Q', [440 240 470 270], 'Inputs', '-+', 'ShowName', 'off');
add_line(s, 'sd cos/1', 'Сумма I/1', 'autorouting', 'smart');
add_line(s, 'sh sin/1', 'Сумма I/2', 'autorouting', 'smart');
add_line(s, 'sd sin/1', 'Сумма Q/1', 'autorouting', 'smart');
add_line(s, 'sh cos/1', 'Сумма Q/2', 'autorouting', 'smart');
blk(s, 'simulink/Sinks/Out1', 'I', [530 88 560 102]);
blk(s, 'simulink/Sinks/Out1', 'Q', [530 248 560 262], 'Port', '2');
add_line(s, 'Сумма I/1', 'I/1', 'autorouting', 'smart');
add_line(s, 'Сумма Q/1', 'Q/1', 'autorouting', 'smart');
end

function buildDemod(s)
% U(n) = sqrt(I^2(n) + Q^2(n))
clearDefault(s);
blk(s, 'simulink/Sources/In1', 'I', [30 40 60 54]);
blk(s, 'simulink/Sources/In1', 'Q', [30 120 60 134], 'Port', '2');
blk(s, 'simulink/Math Operations/Math Function', 'I^2', [110 30 170 64], 'Operator', 'square');
blk(s, 'simulink/Math Operations/Math Function', 'Q^2', [110 110 170 144], 'Operator', 'square');
add_line(s, 'I/1', 'I^2/1', 'autorouting', 'smart');
add_line(s, 'Q/1', 'Q^2/1', 'autorouting', 'smart');
blk(s, 'simulink/Math Operations/Sum', 'Сумма', [220 70 250 100], 'Inputs', '++', 'ShowName', 'off');
add_line(s, 'I^2/1', 'Сумма/1', 'autorouting', 'smart');
add_line(s, 'Q^2/1', 'Сумма/2', 'autorouting', 'smart');
blk(s, 'simulink/Math Operations/Sqrt', 'Корень', [300 70 340 100]);
add_line(s, 'Сумма/1', 'Корень/1', 'autorouting', 'smart');
blk(s, 'simulink/Sinks/Out1', 'U', [400 78 430 92]);
add_line(s, 'Корень/1', 'U/1', 'autorouting', 'smart');
end

function clearDefault(s)
% built-in/Subsystem создаётся пустым; функция оставлена на случай шаблона по умолчанию
lines = find_system(s, 'FindAll', 'on', 'SearchDepth', 1, 'Type', 'line');
if ~isempty(lines)
    delete_line(lines);
end
blocks = find_system(s, 'SearchDepth', 1, 'Type', 'Block');
for k = 1 : numel(blocks)
    if ~strcmp(blocks{k}, s)
        delete_block(blocks{k});
    end
end
end
