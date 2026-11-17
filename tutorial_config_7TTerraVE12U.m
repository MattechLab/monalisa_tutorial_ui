function config = tutorial_config_7TTerraVE12U()
%TUTORIAL_CONFIG_7TTERRAVE12U Independent, visible default paths for 7T.
% Edit these paths, then run main_script_7TTerraVE12U cell 1.
repoRoot = fileparts(mfilename('fullpath'));
workspaceRoot = fileparts(repoRoot);
config.monalisaRoot = fullfile(workspaceRoot,'monalisa');
config.pulseqRoot = fullfile(workspaceRoot,'pulseq');
% No assumed dataset. Empty opens a folder picker at each setup.
% Set an explicit dataset path here only if you want to reuse that choice.
config.baseFolder = '';
% Relative to the selected raw dataset: <baseFolder>/recon_results.
% Set an absolute path here only if you want a different destination.
config.outputRoot = 'recon_results';
config.mexCacheRoot = fullfile(repoRoot,'runtime_cache','mex');
config.sourceFile = [mfilename('fullpath') '.m']; % Printed by setup for clarity.
end
