function [mDir,reconstructionCFile] = mitosius_nobin_ui(defaultDatasetRoot, defaultReconRoot, rovirOptions, preparationOptions)
%RECON_MITOSIUS_UI  User-selection wrapper for recon_mitosius().
%
% Lets user pick:
%   - RECON_FOLDER_DATE (dataset subfolder)
%   - measurement file (twix/raw)
%   - sequence file (.seq)
%   - C source (either pick MID folder OR directly pick C.mat)
% Results default to <selected dataset folder>/recon_results.
% The configured absolute output override is honored when supplied.
% Optional rovirOptions: weightsFile (MAT containing W), nVirtualCoils.
% Omit/leave empty to retain physical coils. Outputs pair Mitosius with its C.
% Optional preparationOptions: matrixSize/nShotOff prefill the settings dialog;
% sequenceParams reuses a resolved header only for the same sequence file.
% files(mainRaw/mainSeq), datasetDir and CfilePath bypass repeated pickers.
% showSettings=false uses the supplied numerical settings without a dialog.
% requireProvenance=true requires matching task C/ROVir provenance.

clc;
mDir = ''; reconstructionCFile = '';
if nargin<3 || isempty(rovirOptions), rovirOptions = struct(); end
if nargin<4 || isempty(preparationOptions), preparationOptions = struct(); end
assert(isstruct(preparationOptions) && isscalar(preparationOptions), ...
    'preparationOptions must be a scalar struct.');
assert(isstruct(rovirOptions) && isscalar(rovirOptions),'rovirOptions must be a scalar struct.');
useROVir = isfield(rovirOptions,'weightsFile') && ~isempty(rovirOptions.weightsFile);
if useROVir
    assert(isfile(rovirOptions.weightsFile),'ROVir weights file not found: %s',rovirOptions.weightsFile);
    assert(isfield(rovirOptions,'nVirtualCoils'),'Specify rovirOptions.nVirtualCoils.');
    weights = load(rovirOptions.weightsFile,'W');
    assert(isfield(weights,'W'),'ROVir weights file must contain W.');
    validateattributes(rovirOptions.nVirtualCoils,{'numeric'}, ...
        {'scalar','integer','positive','<=',min(size(weights.W))});
end
% ---- Defaults (edit to your environment) ----
if nargin<2 || isempty(defaultDatasetRoot) || isempty(defaultReconRoot)
    config = setup_tutorial();
    if nargin<1 || isempty(defaultDatasetRoot), defaultDatasetRoot = config.baseFolder; end
    if nargin<2 || isempty(defaultReconRoot)
        if strcmp(defaultDatasetRoot,config.baseFolder)
            defaultReconRoot = config.outputRoot;
        else
            defaultReconRoot = 'recon_results';
        end
    end
end
% Keep the default tied to the dataset chosen below, including if it changes.
if strcmp(defaultReconRoot,resolve_output_root(defaultDatasetRoot,'recon_results'))
    defaultReconRoot = 'recon_results';
end

% Explicit task paths bypass pickers; legacy calls keep the interactive UI.
if isfield(preparationOptions,'datasetDir') && ~isempty(preparationOptions.datasetDir)
    datasetDir = preparationOptions.datasetDir;
else
disp( 'Select dataset folder containing raw data')
selection_reminder('Mitosius 1/4: raw-data folder', ...
    ['Select the RAW-DATA folder containing the main .dat file and matching .seq file ' ...
     '(for example, 0p25mm). This step selects the acquisition folder, not a calibration or Mitosius result folder.']);
datasetDatePath = uigetdir(defaultDatasetRoot, 'Select RECON_FOLDER_DATE folder (inside datasets)');
if isequal(datasetDatePath, 0); disp('Cancelled.'); return; end



datasetDir = datasetDatePath;

end
if isfield(preparationOptions,'files') && ~isempty(preparationOptions.files.mainRaw)
    measureFile = preparationOptions.files.mainRaw;
else
% ---- Pick measurement file ----
disp('Select measurement file (meas_name)');
selection_reminder('Mitosius 2/4: main raw-data file', ...
    ['Select the MAIN measurement .dat file to reconstruct. ' ...
     'Use the same acquisition used for the C maps and ROVir weights.']);
[meas_name, meas_path] = uigetfile( ...
    {'*.*','All files (*.*)'}, ...
    'Select measurement file (meas_name)', ...
    datasetDir);
if isequal(meas_name, 0); disp('Cancelled.'); return; end
measureFile = fullfile(meas_path, meas_name);

end
if isfield(preparationOptions,'files') && ~isempty(preparationOptions.files.mainSeq)
    seqFile = preparationOptions.files.mainSeq;
else
% ---- Pick sequence file ----
disp('Select sequence file (seqName)');
selection_reminder('Mitosius 3/4: matching Pulseq sequence', ...
    sprintf('Select the Pulseq sequence (.seq) matching this MAIN raw-data file:\n%s',measureFile));
[seqName, seq_path] = uigetfile( ...
    {'*.seq','Pulseq (*.seq)'; '*.*','All files (*.*)'}, ...
    'Select sequence file (seqName)', ...
    datasetDir);
if isequal(seqName, 0); disp('Cancelled.'); return; end
seqFile = fullfile(seq_path, seqName);

end
C_measMID = '';
if isfield(preparationOptions,'CfilePath') && ~isempty(preparationOptions.CfilePath)
    CfilePath = preparationOptions.CfilePath;
else
% ---- Choose C source: pick C.mat directly OR pick a C_measMID folder ----
choice = questdlg('How do you want to provide sensitivity map C?', ...
    'C selection', ...
    'Pick C.mat directly', 'Pick C MID folder', 'Pick C.mat directly');

if isempty(choice); disp('Cancelled.'); return; end

C_measMID  = '';

switch choice
    case 'Pick C.mat directly'
        disp('Select C.mat');
        selection_reminder('Mitosius 4/4: original C.mat', ...
            ['Select the ORIGINAL physical-coil C.mat saved by the C estimation cell. ' ...
             'It is inside the calibration result folder, under C/C.mat. ' ...
             'ROVir compression, if enabled, will be applied to this C automatically.']);
        [cfile, cpath] = uigetfile({'C.mat','C.mat'; '*.mat','MAT-files (*.mat)'}, ...
            'Select C.mat', datasetDir);
        if isequal(cfile, 0); disp('Cancelled.'); return; end
        CfilePath = fullfile(cpath, cfile);

    case 'Pick C MID folder'
        disp('Select <C_measMID>_recon_C folder (or its parent)');
        selection_reminder('Mitosius 4/4: C calibration folder', ...
            ['Select the calibration RESULT folder that contains C/C.mat, ' ...
             'for example MID00019_recon_C_main. Choose the original physical-coil calibration.']);
        CmidPath = uigetdir(datasetDir, 'Select <C_measMID>_recon_C folder (or its parent)');
        if isequal(CmidPath, 0); disp('Cancelled.'); return; end

        CfilePath = choose_calibration_file(CmidPath);

    otherwise
        disp('Cancelled.');
        return;
end

end
% ---- Resolve results beside the selected raw dataset ----
reconRoot = resolve_output_root(datasetDir,defaultReconRoot);
fprintf('Reconstruction output root: %s\n',reconRoot);

% ---- Ask for any “small knobs” you typically tweak ----
prompt  = {'Normalization matrix_size (60 or 120 recommended):', ...
           'nShot_off (initial shots to discard):'};
defans  = {'120','50'};
if isfield(preparationOptions,'matrixSize'), defans{1} = num2str(preparationOptions.matrixSize); end
if isfield(preparationOptions,'nShotOff'), defans{2} = num2str(preparationOptions.nShotOff); end
if isfield(preparationOptions,'showSettings') && ~preparationOptions.showSettings
    answer = defans;
else
    answer = inputdlg(prompt, 'Reconstruction settings', [1 60], defans);
end
if isempty(answer); disp('Cancelled.'); return; end

matrix_size = str2double(answer{1});
nShot_off   = str2double(answer{2});
validateattributes(matrix_size,{'numeric'},{'scalar','integer','positive','even','>=',4});
validateattributes(nShot_off,{'numeric'},{'scalar','integer','nonnegative'});

% ---- Call a modified core that accepts full paths + roots ----
opts = struct();
opts.datasetDir   = datasetDir;    % full date folder
opts.reconRoot    = reconRoot;
opts.measureFile  = measureFile;
opts.seqFile      = seqFile;
opts.CfilePath    = CfilePath;
opts.C_measMID    = C_measMID;
opts.matrix_size  = matrix_size;
opts.nShot_off    = nShot_off;
opts.rovirOptions = rovirOptions;
if isfield(preparationOptions,'sequenceParams'), opts.sequenceParams = preparationOptions.sequenceParams; end
opts.requireProvenance = isfield(preparationOptions,'requireProvenance') && preparationOptions.requireProvenance;

[mDir,reconstructionCFile] = recon_mitosius_from_paths(opts);

end
