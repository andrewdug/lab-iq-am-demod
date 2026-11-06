function save_model_images()
%SAVE_MODEL_IMAGES Сохраняет изображения модели для отчёта. Запускать при наличии дисплея (не в -batch).
here = fileparts(mfilename('fullpath'));
outDir = fullfile(here, '..', 'report', 'figures');
mdl = 'hilbert_iq_sl';
open_system(fullfile(here, [mdl '.slx']));
print(['-s' mdl], '-dpng', '-r170', fullfile(outDir, 'sl_model.png'));
print(['-s' mdl '/НЧ квадратуры I, Q'], '-dpng', '-r170', fullfile(outDir, 'sl_iq_block.png'));
print(['-s' mdl '/АМ-демодулятор'], '-dpng', '-r170', fullfile(outDir, 'sl_demod_block.png'));
close_system(mdl, 0);
end
