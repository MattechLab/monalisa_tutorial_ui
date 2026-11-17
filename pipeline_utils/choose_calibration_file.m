function candidate = choose_calibration_file(folder)
%CHOOSE_CALIBRATION_FILE Include timestamped folders; never pick first silently.
candidate = fullfile(folder,'C','C.mat');
if isfile(candidate), return; end
d = dir(fullfile(folder,'*_recon_C*','C','C.mat'));
assert(~isempty(d),'No calibration C/C.mat found under %s',folder);
paths = arrayfun(@(entry) fullfile(entry.folder,entry.name),d,'UniformOutput',false);
if numel(paths)==1, candidate=paths{1}; return; end
[index,ok] = listdlg('PromptString','Choose the calibration for this task:', ...
    'SelectionMode','single','ListString',paths,'ListSize',[750 250]);
assert(ok,'Calibration selection cancelled.');
candidate = paths{index};
end
