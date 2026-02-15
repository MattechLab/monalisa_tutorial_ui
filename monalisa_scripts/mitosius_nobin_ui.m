function mitosius_nobin_ui(defaultDatasetRoot, defaultReconRoot)
%RECON_MITOSIUS_UI  User-selection wrapper for recon_mitosius().
%
% Lets user pick:
%   - RECON_FOLDER_DATE (dataset subfolder)
%   - measurement file (twix/raw)
%   - sequence file (.seq)
%   - C source (either pick MID folder OR directly pick C.mat)
%   - output recon root (recon_folder base)

clc;

% ---- Defaults (edit to your environment) ----
if nargin<1 || isempty(defaultDatasetRoot)
defaultDatasetRoot = '/Users/cag/Documents/Dataset/datasets';
end
if nargin<2 || isempty(defaultReconRoot)
defaultReconRoot   = '/Users/cag/Documents/Dataset/recon_results';
end

% ---- Select dataset date folder (subfolder under datasets) ----
disp( 'Select dataset folder containing raw data')
datasetDatePath = uigetdir(defaultDatasetRoot, 'Select RECON_FOLDER_DATE folder (inside datasets)');
if isequal(datasetDatePath, 0); disp('Cancelled.'); return; end



datasetDir = datasetDatePath;  % full path to date folder

% ---- Pick measurement file ----
disp('Select measurement file (meas_name)');
[meas_name, meas_path] = uigetfile( ...
    {'*.*','All files (*.*)'}, ...
    'Select measurement file (meas_name)', ...
    datasetDir);
if isequal(meas_name, 0); disp('Cancelled.'); return; end
measureFile = fullfile(meas_path, meas_name);

% ---- Pick sequence file ----
disp('Select sequence file (seqName)');
[seqName, seq_path] = uigetfile( ...
    {'*.seq','Pulseq (*.seq)'; '*.*','All files (*.*)'}, ...
    'Select sequence file (seqName)', ...
    datasetDir);
if isequal(seqName, 0); disp('Cancelled.'); return; end
seqFile = fullfile(seq_path, seqName);

% ---- Choose C source: pick C.mat directly OR pick a C_measMID folder ----
choice = questdlg('How do you want to provide sensitivity map C?', ...
    'C selection', ...
    'Pick C.mat directly', 'Pick C MID folder', 'Pick C.mat directly');

if isempty(choice); disp('Cancelled.'); return; end

CfilePath = '';
C_measMID  = '';

switch choice
    case 'Pick C.mat directly'
        disp('Select C.mat');
        [cfile, cpath] = uigetfile({'C.mat','C.mat'; '*.mat','MAT-files (*.mat)'}, ...
            'Select C.mat', datasetDir);
        if isequal(cfile, 0); disp('Cancelled.'); return; end
        CfilePath = fullfile(cpath, cfile);

    case 'Pick C MID folder'
        disp('Select <C_measMID>_recon_C folder (or its parent)');
        CmidPath = uigetdir(datasetDir, 'Select <C_measMID>_recon_C folder (or its parent)');
        if isequal(CmidPath, 0); disp('Cancelled.'); return; end

        % Try to infer C_measMID from folder name
        % Accept either ".../<MID>_recon_C" or a parent folder containing it
        [~, folderName] = fileparts(CmidPath);

        if contains(folderName, '_recon_C')
            token = regexp(folderName, 'MID\d+', 'match');
            if ~isempty(token); C_measMID = token{1}; end
            candidate = fullfile(CmidPath, 'C', 'C.mat');
        else
            % look for "*_recon_C/C/C.mat" under chosen dir
            d = dir(fullfile(CmidPath, '*_recon_C', 'C', 'C.mat'));
            if isempty(d)
                error('Could not find *_recon_C/C/C.mat under: %s', CmidPath);
            end
            candidate = fullfile(d(1).folder, d(1).name);

            token = regexp(d(1).folder, 'MID\d+', 'match');
            if ~isempty(token); C_measMID = token{1}; end
        end

        CfilePath = candidate;

    otherwise
        disp('Cancelled.');
        return;
end

% ---- Select recon output root (optional) ----
disp('Select recon output root (recon_folder)');
reconRoot = uigetdir(defaultReconRoot, 'Select recon output root (recon_folder)');
if isequal(reconRoot, 0)
    reconRoot = defaultReconRoot; % user cancelled -> use default
end

% ---- Ask for any “small knobs” you typically tweak ----
prompt  = {'matrix_size (<=240 recommended unless you know memory is OK):', ...
           'nShot_off (Bern manual):'};
defans  = {'480','14'};
answer  = inputdlg(prompt, 'Reconstruction settings', [1 60], defans);
if isempty(answer); disp('Cancelled.'); return; end

matrix_size = str2double(answer{1});
nShot_off   = str2double(answer{2});
if isnan(matrix_size) || isnan(nShot_off)
    error('matrix_size and nShot_off must be numeric.');
end

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

recon_mitosius_from_paths(opts);

end

function recon_mitosius_from_paths(opts)
%RECON_MITOSIUS_FROM_PATHS  Core recon using explicit paths from UI/batch.
%
% Required fields:
%   opts.datasetDir   (full path to .../datasets/<RECON_FOLDER_DATE>)
%   opts.reconRoot    (full path to recon root)
%   opts.measureFile  (full path to measurement)
%   opts.seqFile      (full path to .seq)
%   opts.CfilePath    (full path to C.mat)
%
% Optional:
%   opts.matrix_size
%   opts.nShot_off

arguments
    opts struct
end

disp(opts.measureFile);
disp(opts.seqFile);

% ---- Extract MID from meas filename ----
token = regexp(opts.measureFile, 'MID\d+', 'match');
if isempty(token)
    % try from the basename
    [~, base] = fileparts(opts.measureFile);
    token = regexp(base, 'MID\d+', 'match');
end
assert(~isempty(token), 'Could not extract MID\d+ from measurement file path/name.');
measMID = token{1};
disp(['measMID = ', measMID]);

% ---- Output recon directory ----
reconDir = fullfile(opts.reconRoot, [measMID '_recon']);
disp(['reconDir = ', reconDir]);

if ~isfolder(reconDir)
    mkdir(reconDir);
    disp(['Directory created: ', reconDir]);
else
    disp(['Directory already exists: ', reconDir]);
end
cd(reconDir);

otherDir = fullfile(reconDir, 'T1_LIBRE_woBinning', 'other');
if ~isfolder(otherDir)
    mkdir(otherDir);
end

% ---- Load and Configure Data ----
saveflag = 1; %#ok<NASGU>
seqParams = extract_seq_params(opts.seqFile);

flagSS = 1; %#ok<NASGU>
flagExcludeSI = 1; %#ok<NASGU>

reader = createRawDataReader(opts.measureFile, 1);

% Manual part exposed as opts.nShot_off
if isfield(opts, 'nShot_off') && ~isempty(opts.nShot_off)
    reader.acquisitionParams.nShot_off = opts.nShot_off;
else
    reader.acquisitionParams.nShot_off = 14; % fallback
end

reader.acquisitionParams.traj_type = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name = strcat(opts.seqFile);

% Hash check
isMatch = check_hash(opts.measureFile, reader.acquisitionParams.pulseqTrajFile_name); %#ok<NASGU>

if isfield(seqParams, 'nshot')
    reader.acquisitionParams.nShot = seqParams.nshot;
end
if isfield(seqParams, 'nseg')
    reader.acquisitionParams.nSeg = seqParams.nseg;
end

y_tot = reader.readRawData(true, true);               % Filter nShotOff and SI  [nCh, nSample, nLine]
t_tot = bmTraj(reader.acquisitionParams);             % Compute trajectory
ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3'); %#ok<NASGU>

load(opts.CfilePath, 'C');  % Load sensitivity maps
disp(['C loaded from: ', opts.CfilePath]);

% ---- Reconstruction parameters ----
if isfield(opts, 'matrix_size') && ~isempty(opts.matrix_size)
    matrix_size = opts.matrix_size;
else
    matrix_size = 480;
end

N_u = [matrix_size, matrix_size, matrix_size];
dK_u  = 1 ./ (seqParams.fov*2e3);

C = bmImResize(C, [48, 48, 48], N_u);
x_tot = bmMathilda(y_tot, t_tot, bmVolumeElement(t_tot,'voronoi_full_radial3'), C, N_u, N_u, dK_u); %#ok<NASGU>

temp_im = x_tot(round(matrix_size/4):round(matrix_size/4*3), ...
                round(matrix_size/4):round(matrix_size/4*3), ...
                round(matrix_size/2));
normalize_val = mean(temp_im(:));
disp(['normalize_val = ', num2str(normalize_val)]);

% NOTE: your original "if real(y_tot)<1" is ambiguous (y_tot array).
% I keep the intent but make it deterministic:
testval = real(y_tot(1));  % choose one sample
if testval < 1
    y_tot = y_tot / normalize_val;
end

% ---- eMask ----
eMask = true(1, seqParams.nshot*seqParams.nseg);
eMask(1:reader.acquisitionParams.nShot_off*seqParams.nseg) = false;

eMaskFilePath = fullfile(otherDir, 'eMask_woBin.mat');
save(eMaskFilePath, 'eMask');
disp(['eMask saved: ', eMaskFilePath]);

% ---- Load eye mask (you used a different filename without .mat) ----
% If you truly want to load the .mat you just saved:
eyeMaskStruct = load(eMaskFilePath, 'eMask');
eyeMask = eyeMaskStruct.eMask;

size_Mask = size(eyeMask);
nbins = size_Mask(1); %#ok<NASGU>

% WARNING: In your snippet, you reshape eMask into [nbins, nseg, nshot].
% That only makes sense if eMask is actually nbins x ... (not a 1xN vector).
% I’m leaving your operations in place, but you should confirm eMask dims.
eyeMask = reshape(eyeMask, [numel(eyeMask)/ (seqParams.nseg*seqParams.nshot), seqParams.nseg, seqParams.nshot]); 
eyeMask(:, 1, :) = [];
eyeMask(:, :, 1:reader.acquisitionParams.nShot_off ) = [];
eyeMask = bmPointReshape(eyeMask);

% ---- mitosius saving ----
mDir = fullfile(reconDir, 'T1_LIBRE_woBinning', 'mitosius');
if ~isfolder(mDir); mkdir(mDir); end
disp(['mDir = ', mDir]);

[y, t] = bmMitosis(y_tot, t_tot, eyeMask);
y = bmPermuteToCol(y);
ve  = bmVolumeElement(t, 'voronoi_full_radial3');

bmMitosius_create(mDir, y, t, ve);
disp(['Mitosius files saved: ', mDir]);

end