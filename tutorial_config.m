function config = tutorial_config()
%TUTORIAL_CONFIG Visible, editable default paths for main_script_3TPrismaXA60.
% Edit these assignments, then rerun main_script_3TPrismaXA60 cell 1.
% Setup does not read/write saved settings or silently override this file.
% For a different task within this dataset, edit tutorial_task_config_3T.m:
% main/prescan inputs, calibration method, task name, ROVir, and recon settings.
repoRoot = fileparts(mfilename('fullpath'));
workspaceRoot = fileparts(repoRoot);

% Dependencies: default to sibling checkouts. Change to your installed paths.
config.monalisaRoot = fullfile(workspaceRoot,'monalisa');
config.pulseqRoot = fullfile(workspaceRoot,'pulseq');

% No assumed dataset. Empty opens a folder picker at each setup.
% Set an explicit dataset path here only if you want to reuse that choice.
config.baseFolder = '';

% Results: relative to the selected raw dataset folder, not this code repo.
% Default: <baseFolder>/recon_results. An absolute path overrides this.
config.outputRoot = 'recon_results';
% Generated MEX binaries stay with the code, separate from scan results.
config.mexCacheRoot = fullfile(repoRoot,'runtime_cache','mex');
config.sourceFile = [mfilename('fullpath') '.m']; % Printed by setup for clarity.
% macOS: setup recreates missing cache files from your compiled MONALISA MEXs,
% linking copies to MATLAB's libomp. See runtime_cache/README.md for sources
% and rebuild instructions. Windows/Linux use MONALISA's original binaries.
end
