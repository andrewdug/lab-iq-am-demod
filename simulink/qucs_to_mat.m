function qucs_to_mat()
%QUCS_TO_MAT Передаёт результаты Qucs-S в цифровой обработчик (п. 6 задания).
%  Читает qucs/LR2k5_N1_hilbert_TR.dat и _AC.dat, берёт отсчёты с шагом Δt
%  и сохраняет simulink/qucs_data.mat. Файл .mat хранится в репозитории,
%  поэтому модель Simulink запускается без повторного моделирования в Qucs-S.

here = fileparts(mfilename('fullpath'));
qucs = fullfile(here, '..', 'qucs');

M  = 1;                    % номер бригады
f0 = (400 + 30 * M) * 1e6;
fs = 4 * f0;
Dt = 1 / fs;

tr = readDat(fullfile(qucs, 'LR2k5_N1_hilbert_TR.dat'), ...
    {'time', 'x0.Vt', 'xd.Vt', 'xd14.Vt', 'xd18.Vt', 'y.Vt', 'yw.Vt', 'y14.Vt', 'y18.Vt'});
tt = tr('time');
tn = (0 : floor(tt(end) / Dt + 1e-9)).' * Dt;      % t_n = (n - 1) Δt
s  = sample(tr, 'x0.Vt', tn);
sd = [sample(tr, 'xd.Vt', tn), sample(tr, 'xd.Vt', tn), sample(tr, 'xd14.Vt', tn), sample(tr, 'xd18.Vt', tn)];
sh = [sample(tr, 'y.Vt', tn), sample(tr, 'yw.Vt', tn), sample(tr, 'y14.Vt', tn), sample(tr, 'y18.Vt', tn)];

ac = readDat(fullfile(qucs, 'LR2k5_N1_hilbert_AC.dat'), ...
    {'acfrequency', 'xd.v', 'xd14.v', 'xd18.v', 'y.v', 'yw.v', 'y14.v', 'y18.v'});
f  = ac('acfrequency');
Hd = [ac('xd.v'), ac('xd.v'), ac('xd14.v'), ac('xd18.v')];     % элементы задержки
Hf = [ac('y.v'), ac('yw.v'), ac('y14.v'), ac('y18.v')];        % фильтры

names  = {'порядок 10, без окна', 'порядок 10, Хэмминг', 'порядок 14, Хэмминг', 'порядок 18, Хэмминг'};
orders = [10 10 14 18];
save(fullfile(here, 'qucs_data.mat'), 'tn', 's', 'sd', 'sh', 'f', 'Hd', 'Hf', 'names', 'orders', '-v7');
fprintf('qucs_data.mat: %d отсчётов, %d частот\n', numel(tn), numel(f));
end

function x = sample(ds, name, tn)
x = interp1(ds('time'), ds(name), tn, 'linear');
end

function ds = readDat(fname, want)
% Набор данных Qucs: блоки <indep имя N> ... </indep> и <dep имя зависимость> ... </dep>
txt = fileread(fname);
tok = regexp(txt, '<(?:indep|dep) (\S+)[^>]*>(.*?)</(?:indep|dep)>', 'tokens');
ds  = containers.Map();
for i = 1 : numel(tok)
    name = tok{i}{1};
    if ~any(strcmp(name, want))
        continue
    end
    body = tok{i}{2};
    if any(body == 'j')
        v = sscanf(regexprep(body, '([+-])j', ' $1'), '%f');   % пары Re Im
        ds(name) = complex(v(1 : 2 : end), v(2 : 2 : end));
    else
        ds(name) = sscanf(body, '%f');
    end
end
miss = want(~isKey(ds, want));
if ~isempty(miss)
    error('В %s нет переменных: %s', fname, strjoin(miss, ', '));
end
end
