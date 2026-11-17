function config = setup_tutorial(config)
%SETUP_TUTORIAL Configure local dependencies and the MONALISA MEX runtime.
% Run from the repository root, or let main_script_3TPrismaXA60 call this function.
% Defaults are visible in tutorial_config.m; selections apply to this run only.
% Generated MEX files go to config.mexCacheRoot. No settings are saved.
repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'monalisa_scripts'), ...
    fullfile(repoRoot,'pipeline_utils'), fullfile(repoRoot,'seq_utils'));
defaults = tutorial_config();
if nargin == 0
    config = defaults;
end
assert(isstruct(config) && isscalar(config),'Configuration must be a scalar struct.');
fields = fieldnames(defaults);
for k = 1:numel(fields)
    if ~isfield(config,fields{k}), config.(fields{k}) = defaults.(fields{k}); end
end
parent = fileparts(repoRoot);
config.monalisaRoot = dependency_root(config,'monalisaRoot', ...
    fullfile(parent,'monalisa'),fullfile('src','rawDataReader','createRawDataReader.m'),'MONALISA');
config.pulseqRoot = dependency_root(config,'pulseqRoot', ...
    fullfile(parent,'pulseq'),fullfile('matlab','+mr','@Sequence','Sequence.m'),'Pulseq');
if ~isfield(config,'baseFolder') || ~isfolder(config.baseFolder)
    selection_reminder('Setup: raw-data folder', ...
        ['Select the subject/acquisition folder containing the main raw-data file (.dat) ' ...
         'and its matching Pulseq sequence (.seq).']);
    config.baseFolder = uigetdir(repoRoot,'Select dataset folder (edit the entry script''s config file to keep this choice)');
    assert(~isequal(config.baseFolder,0),'Dataset selection cancelled.');
end
if ~isfield(config,'outputRoot') || isempty(config.outputRoot)
    config.outputRoot = 'recon_results';
end
% Resolve AFTER dataset selection so a new dataset gets its own results.
config.outputRoot = resolve_output_root(config.baseFolder,config.outputRoot);
addpath(genpath(fullfile(config.monalisaRoot,'src')));
addpath(fullfile(config.monalisaRoot,'third_part','twix_for_monalisa'));
addpath(fullfile(config.pulseqRoot,'matlab'));
required = {'mapVBVD_JH_for_monalisa','bmTraj2SparseMat','bmMathilda', ...
    'fitgmdist','graythresh','regionprops'};
for k = 1:numel(required)
    assert(~isempty(which(required{k})),'Missing dependency: %s',required{k});
end
requiredMex = {'bmImCrossMean3_mex','bmImLaplacian3_mex', ...
    'bmGriddingMatrix_secret_length3_mex','bmSparseMat_r_nJump2index_mex', ...
    'bmGridder_n2u_leight3_mex'};
for k = 1:numel(requiredMex)
    assert(exist(requiredMex{k},'file')==3, ...
        'Build MONALISA MEX files for this MATLAB/platform first. Missing: %s',requiredMex{k});
end
runtime = setup_monalisa_runtime(config.monalisaRoot,config.mexCacheRoot);
fprintf('Default paths: %s (edit this file to change future runs)\n',config.sourceFile);
disp(config);
if ismac
    fprintf('MEX cache: %s\n',fullfile(config.mexCacheRoot,computer('arch'),version('-release')));
    fprintf('Missing cache files are recreated by this setup from MONALISA source MEXs.\n');
end
fprintf('MEX runtime: %s\n',runtime);
end

function root = dependency_root(config,field,candidate,marker,label)
if isfield(config,field), root = config.(field); else, root = candidate; end
if ~isfile(fullfile(root,marker))
    selection_reminder(['Setup: ' label ' installation folder'], ...
        sprintf('Select the %s repository folder.\nIt must contain:\n%s',label,marker));
    root = uigetdir(fileparts(candidate),['Select the ' label ' repository']);
    assert(~isequal(root,0),'%s selection cancelled.',label);
end
assert(isfile(fullfile(root,marker)),'Invalid %s repository: %s',label,root);
end
