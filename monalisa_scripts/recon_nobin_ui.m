function [x, x0, outDir] = recon_nobin_ui(defaultBaseFolder, defaultDatasetFolder, flagCS, preparedCFile, preparedMitosiusDir, reconstructionOptions)
%RECON_NOBIN_UI  User-selection wrapper for recon_nobin().
%
% Outputs:
%   x, x0  (same as recon_nobin)
%
% Requirements:
%   - The recon directory should already exist and contain:
%       <reconDir>/T1_LIBRE_woBinning/mitosius/  (y,t,ve saved by bmMitosius_create)
%   - C.mat should be accessible (selectable)
% Optional preparedCFile/preparedMitosiusDir bypass pickers with matched outputs
% returned by mitosius_nobin_ui (including ROVir-compressed C).
% Optional reconstructionOptions: Matrix_size, nIter, delta, rho, outputRoot,
% showSettings (false uses these settings directly). outDir is the new run.

clc; outDir='';
if nargin<6, reconstructionOptions=struct(); end

% ---- Defaults (edit to your machine) ----
if nargin<2 || isempty(defaultBaseFolder) || isempty(defaultDatasetFolder)
    config = setup_tutorial();
    if nargin<1 || isempty(defaultBaseFolder), defaultBaseFolder = config.outputRoot; end
    if nargin<2 || isempty(defaultDatasetFolder), defaultDatasetFolder = config.baseFolder; end
end

if nargin<3 || isempty(flagCS)
    flagCS = 1;
end
if flagCS
    disp('Will perform compressed sensing...')
else
     disp('Will NOT perform compressed sensing...')
end
% ---- Choose reconDir directly (most robust) ----
disp('Select <MIDxxxxx_recon> folder (reconDir)');
if nargin>=5 && ~isempty(preparedMitosiusDir)
    assert(isfile(fullfile(preparedMitosiusDir,'mitosius_size.mat')), ...
        'Prepared Mitosius not found: %s',preparedMitosiusDir);
    reconDir = fileparts(fileparts(preparedMitosiusDir));
else
    preparedMitosiusDir = '';
    selection_reminder('Reconstruction: prepared result folder', ...
        ['Select the MIDxxxxx_recon folder created by Mitosius preparation. ' ...
         'It contains T1_LIBRE_woBinning/mitosius. This step selects the prepared-data result folder.']);
    reconDir = uigetdir(defaultBaseFolder, 'Select <MIDxxxxx_recon> folder (reconDir)');
end
if isequal(reconDir, 0); disp('Cancelled.'); x=[]; x0=[]; return; end

% infer MID from folder name if possible
[~, reconFolderName] = fileparts(reconDir);
token = regexp(reconFolderName, 'MID\d+', 'match');
if isempty(token)
    measMID = '';
else
    measMID = token{1};
end

% infer RECON_FOLDER_DATE from parent folder (base/date/MID_recon)
parentDir = fileparts(reconDir);
[~, RECON_FOLDER_DATE_guess] = fileparts(parentDir);

% ---- Optional: choose measurement file (only used to extract MID if not found) ----
if isempty(measMID)
    pickMeas = questdlg('MID not detected from reconDir name. Pick measurement file to extract MID?', ...
        'MID missing', 'Yes', 'No', 'Yes');
    if strcmp(pickMeas, 'Yes')
        disp('Select measurement file (to parse MID)');
        selection_reminder('Reconstruction: identify the acquisition', ...
            'Select the original MAIN .dat file belonging to this reconstruction. Its name supplies the measurement MID.');
        [meas_name, ~] = uigetfile({'*.*','All files'}, 'Select measurement file (to parse MID)');
        if isequal(meas_name, 0); disp('Cancelled.'); x=[]; x0=[]; return; end
        token = regexp(meas_name, 'MID\d+', 'match');
        assert(~isempty(token), 'Could not parse MID\d+ from measurement filename.');
        measMID = token{1};

        % If reconDir name didn’t include MID, make sure user selected correct folder
        if ~contains(reconFolderName, measMID)
            warning('Selected reconDir name does not contain %s. Continuing anyway.', measMID);
        end
    end
end

% ---- Select C source ----
datasetFolder = defaultDatasetFolder;
baseFolder    = defaultBaseFolder;
RECON_FOLDER_DATE = RECON_FOLDER_DATE_guess;

% allow user to override datasetFolder (helpful if running elsewhere)
if ~isfolder(datasetFolder)
    disp('Select datasetFolder (datasets root)');
    selection_reminder('Reconstruction: raw-data folder', ...
        'Select the RAW-DATA folder containing the main .dat file and matching .seq for this reconstruction.');
    datasetFolder = uigetdir(pwd, 'Select datasetFolder (datasets root)');
    if isequal(datasetFolder, 0); disp('Cancelled.'); x=[]; x0=[]; return; end
end

if nargin>=4 && ~isempty(preparedCFile)
    choice = 'Prepared C';
else
    choice = questdlg('How do you want to provide sensitivity map C?', ...
    'C selection', ...
    'Pick C.mat directly', 'Pick C MID folder', 'Pick C.mat directly');
end

if isempty(choice); disp('Cancelled.'); x=[]; x0=[]; return; end

C_measMID = '';

switch choice
    case 'Prepared C'
        CfilePath = preparedCFile;
    case 'Pick C.mat directly'
        disp('Select C.mat');
        selection_reminder('Reconstruction: C matching the prepared data', ...
            ['Select the C file matching the exported Mitosius coil channels. ' ...
             'For ROVir-compressed data, choose C_compressed_N.mat inside the mitosius folder. ' ...
             'For data exported without compression, choose the original C.mat.']);
        [cfile, cpath] = uigetfile({'C.mat','C.mat'; '*.mat','MAT-files (*.mat)'}, ...
            'Select C.mat', datasetFolder);
        if isequal(cfile, 0); disp('Cancelled.'); x=[]; x0=[]; return; end
        CfilePath = fullfile(cpath, cfile);

        token = regexp(CfilePath, 'MID\d+', 'match');
        if ~isempty(token)
            C_measMID = token{1};
        end

    case 'Pick C MID folder'
        disp('Select <C_measMID>_recon_C folder (or parent containing it)');
        selection_reminder('Reconstruction: original C calibration folder', ...
            ['Select the calibration RESULT folder containing the original C/C.mat. ' ...
             'This folder option is for uncompressed data. For compressed data, ' ...
             'cancel and use Pick C.mat directly to select C_compressed_N.mat.']);
        CmidPath = uigetdir(datasetFolder, 'Select <C_measMID>_recon_C folder (or parent containing it)');
        if isequal(CmidPath, 0); disp('Cancelled.'); x=[]; x0=[]; return; end

        CfilePath = choose_calibration_file(CmidPath);

    otherwise
        disp('Cancelled.'); x=[]; x0=[]; return;
end

assert(isfile(CfilePath), 'C.mat not found: %s', CfilePath);

% ---- User-settable numeric knobs ----
prompt = { ...
    'Matrix_size:', ...
    'nIter:', ...
    'delta:', ...
    'rho (leave empty to use rho = 10*delta):' ...
    };
defans = { '120', '20', '1', '' };

keys={'Matrix_size','nIter','delta','rho'};
for k=1:numel(keys)
    if isfield(reconstructionOptions,keys{k}), defans{k}=num2str(reconstructionOptions.(keys{k}),17); end
end
if isfield(reconstructionOptions,'showSettings') && ~reconstructionOptions.showSettings
    answ=defans;
else
    answ = inputdlg(prompt, 'Recon parameters', [1 60], defans);
end
if isempty(answ); disp('Cancelled.'); x=[]; x0=[]; return; end

Matrix_size = str2double(answ{1});
nIter       = str2double(answ{2});
delta       = str2double(answ{3});
rho_in      = str2double(answ{4});

assert(isfinite(Matrix_size) && Matrix_size>0, 'Matrix_size invalid.');
assert(isfinite(nIter) && nIter>0, 'nIter invalid.');
assert(isfinite(delta) && delta>0, 'delta invalid.');

if ~isfinite(rho_in)
    rho = 10*delta;
else
    rho = rho_in;
end

% ---- Pack opts and run ----
opts = struct();
opts.reconDir      = reconDir;
opts.RECON_FOLDER_DATE = RECON_FOLDER_DATE;
opts.baseFolder    = baseFolder;
opts.datasetFolder = datasetFolder;
opts.measMID       = measMID;
opts.C_measMID     = C_measMID;
opts.CfilePath     = CfilePath;
opts.mitosiusDir   = preparedMitosiusDir;

opts.Matrix_size   = Matrix_size;
opts.nIter         = nIter;
opts.delta         = delta;
opts.rho           = rho;
opts.flagCS = flagCS;
if isfield(reconstructionOptions,'outputRoot'), opts.outputRoot=reconstructionOptions.outputRoot; end
[x, x0, outDir] = recon_nobin_from_paths(opts);

end
