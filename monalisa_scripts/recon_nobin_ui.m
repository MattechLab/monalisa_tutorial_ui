function [x, x0] = recon_nobin_ui(defaultBaseFolder, defaultDatasetFolder, flagCS)
%RECON_NOBIN_UI  User-selection wrapper for recon_nobin().
%
% Outputs:
%   x, x0  (same as recon_nobin)
%
% Requirements:
%   - The recon directory should already exist and contain:
%       <reconDir>/T1_LIBRE_woBinning/mitosius/  (y,t,ve saved by bmMitosius_create)
%   - C.mat should be accessible (selectable)

clc;

% ---- Defaults (edit to your machine) ----
if nargin<1 || isempty(defaultBaseFolder)
defaultBaseFolder    = '/Users/cag/Documents/Dataset/recon_results/';
end
if nargin<2 || isempty(defaultDatasetFolder)
defaultDatasetFolder = '/Users/cag/Documents/Dataset/datasets/';
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
reconDir = uigetdir(defaultBaseFolder, 'Select <MIDxxxxx_recon> folder (reconDir)');
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
[baseFolderGuess, RECON_FOLDER_DATE_guess] = fileparts(parentDir); %#ok<ASGLU>

% ---- Optional: choose measurement file (only used to extract MID if not found) ----
if isempty(measMID)
    pickMeas = questdlg('MID not detected from reconDir name. Pick measurement file to extract MID?', ...
        'MID missing', 'Yes', 'No', 'Yes');
    if strcmp(pickMeas, 'Yes')
        disp('Select measurement file (to parse MID)');
        [meas_name, meas_path] = uigetfile({'*.*','All files'}, 'Select measurement file (to parse MID)');
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
    datasetFolder = uigetdir(pwd, 'Select datasetFolder (datasets root)');
    if isequal(datasetFolder, 0); disp('Cancelled.'); x=[]; x0=[]; return; end
end

choice = questdlg('How do you want to provide sensitivity map C?', ...
    'C selection', ...
    'Pick C.mat directly', 'Pick C MID folder', 'Pick C.mat directly');

if isempty(choice); disp('Cancelled.'); x=[]; x0=[]; return; end

CfilePath = '';
C_measMID = '';

switch choice
    case 'Pick C.mat directly'
        disp('Select C.mat');
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
        CmidPath = uigetdir(datasetFolder, 'Select <C_measMID>_recon_C folder (or parent containing it)');
        if isequal(CmidPath, 0); disp('Cancelled.'); x=[]; x0=[]; return; end

        [~, folderName] = fileparts(CmidPath);
        if contains(folderName, '_recon_C')
            token = regexp(folderName, 'MID\d+', 'match');
            if ~isempty(token); C_measMID = token{1}; end
            candidate = fullfile(CmidPath, 'C', 'C.mat');
        else
            d = dir(fullfile(CmidPath, '*_recon_C', 'C', 'C.mat'));
            assert(~isempty(d), 'Could not find *_recon_C/C/C.mat under: %s', CmidPath);
            candidate = fullfile(d(1).folder, d(1).name);
            token = regexp(candidate, 'MID\d+', 'match');
            if ~isempty(token); C_measMID = token{1}; end
        end
        CfilePath = candidate;

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

answ = inputdlg(prompt, 'Recon parameters', [1 60], defans);
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

opts.Matrix_size   = Matrix_size;
opts.nIter         = nIter;
opts.delta         = delta;
opts.rho           = rho;
opts.flagCS = flagCS;
[x, x0] = recon_nobin_from_ui(opts);

end

function [x, x0] = recon_nobin_from_ui(opts)
%RECON_NOBIN_FROM_UI  Core recon_nobin with explicit selectable paths/params.

arguments
    opts struct
end

reconDir = opts.reconDir;
disp(['reconDir: ', reconDir]);
cd(reconDir);

mDir = fullfile(reconDir, 'T1_LIBRE_woBinning', 'mitosius');
disp(['mDir: ', mDir]);

assert(isfolder(mDir), 'Mitosius folder not found: %s', mDir);

CfilePath = opts.CfilePath;
disp(['CfilePath: ', CfilePath]);
assert(isfile(CfilePath), 'C.mat not found: %s', CfilePath);

Matrix_size = opts.Matrix_size;
nIter       = opts.nIter;
delta       = opts.delta;
rho         = opts.rho;

% ---- Load mitosius ----
y   = bmMitosius_load(mDir, 'y');
t   = bmMitosius_load(mDir, 't');
ve  = bmMitosius_load(mDir, 've');

disp('Mitosius has been loaded!');
disp(mDir);

% ---- Recon grid params (same as your code) ----
ReconFov = 240; % mm
N_u   = [Matrix_size, Matrix_size, Matrix_size];
n_u   = [Matrix_size, Matrix_size, Matrix_size];
dK_u  = [1, 1, 1] ./ ReconFov;

nFr = size(y, 1);

% ---- Load / resize C ----
S = load(CfilePath);
assert(isfield(S, 'C'), 'C.mat does not contain variable C.');
C = S.C;

C = bmImResize(C, [48, 48, 48], N_u);

% ---- Compute x0 (same logic) ----
serial = false;
if serial
    nCh = 44;
    x0 = bmZero([N_u, nCh], 'complex_single');
    for i = 1:nCh
        x0(:, :, :, i) = bmMathilda(y{1}(:, i), t{1}, ve{1}, [], N_u, n_u, dK_u, [], [], [], []);
    end
    x0 = bmCoilSense_pinv(C, x0, N_u);
else
    x0 = cell(nFr, 1);
    for i = 1:nFr
        x0{i} = bmMathilda(y{i}, t{i}, ve{i}, C, N_u, n_u, dK_u, [], [], [], []);
    end
end

% ---- Save x0 ----
outDir = fullfile(reconDir, 'T1_LIBRE_woBinning', 'output');
if ~isfolder(outDir); mkdir(outDir); end

x0Path = fullfile(outDir, ['x0_' num2str(Matrix_size) '.mat']);
save(x0Path, 'x0', '-v7.3');
disp(['x0 saved: ', x0Path]);

% ---- Sparse matrices ----
    if opts.flagCS
        [Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);
        disp('Gridding is done!')
        % ---- Steva/Teva params ----
        witness_ind = [1, 12, 15];
        nCGD    = 4;
        ve_max  = 10 * prod(dK_u(:));
        
        % ---- Run recon ----
        if nFr <= 1
            disp('Preparing Steva...')
            witness_info = bmWitnessInfo('stevaMorphosia_custom', witness_ind);
            witness_info.save_witnessIm_flag = 0;
        
            x = bmSteva( ...
                x0{1}, [], [], y{1}, ve{1}, C, Gu{1}, Gut{1}, n_u, ...
                delta, rho, nCGD, ve_max, ...
                nIter, ...
                witness_info);
        
        else
            witness_info = bmWitnessInfo('tevaMorphosia_custom', witness_ind);
            witness_info.save_witnessIm_flag = true;
        
            x = bmTevaMorphosia_chain( ...
                x0, [], [], ...
                y, ve, C, ...
                Gu, Gut, n_u, ...
                [], [], ...
                delta, rho, 'normal', ...
                nCGD, ve_max, ...
                nIter, ...
                witness_info);
        end
        
        % ---- Save x ----
        xPath = fullfile(outDir, sprintf('x_nIter%d_delta_%.3f_rho_%.3f_%d.mat', nIter, delta, rho, Matrix_size));
        save(xPath, 'x', '-v7.3');
        disp(['x saved: ', xPath]);
    else 
        x=[];
     end
end