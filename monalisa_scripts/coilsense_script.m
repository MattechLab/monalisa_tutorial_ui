%% =========================================================================
%  COILSENSE_SCRIPT
%
%  Purpose:
%  --------
%  Estimate non-Cartesian coil sensitivity maps using:
%      - Body coil prescan (reference)
%      - Head coil prescan
%      - Main measurement (for orientation validation)
%
%  Workflow:
%  ---------
%  1. Select prescan sequence + raw data (body + head coils)
%  2. Configure raw data readers
%  3. Validate orientation consistency with main measurement
%  4. Compute non-Cartesian gridding operators
%  5. Generate automatic mask
%  6. Estimate coil sensitivity maps (two-stage method)
%
%  Inputs:
%  -------
%  baseFolder : starting directory for file selection dialogs
%
%  Outputs:
%  --------
%  C1                  : Final coil sensitivity maps
%  reconDir            : Reconstruction directory path
%  prescan_seqParams   : Extracted prescan sequence parameters
%  xsos_all            : Sum-of-squares prescan image (orientation check)
%  ncc_all             : Normalized cross-correlation values
%
% =========================================================================
function [C1, reconDir, prescan_seqParams, xsos_all, ncc_all] = coilsense_script(baseFolder)

clc;

%% =========================================================================
%% === Dependency Reminder ===
% The following toolboxes must be added to MATLAB path:
%   - monalisa
%   - pulseq
%   - mapVBVD
% =========================================================================
warning('Reminder: please make sure monalisa, pulseq and mapVBVD are added to the path!');

%% =========================================================================
%% === 1. User Input: Prescan Files ===
% User selects:
%   - Prescan sequence file (*.seq)
%   - Body coil raw data (*.dat)
%   - Head coil raw data (*.dat)
% =========================================================================

disp('Select prescan sequence file');
[seqName, seqFolder] = uigetfile('*.seq', ...
    'Select prescan sequence file', baseFolder);

if seqName == 0
    error('Sequence file selection was cancelled');
end

disp('Select bc data file');
[bodyCoilFilename, rawDir] = uigetfile('*.dat', ...
    'Select bc data file', baseFolder);

if bodyCoilFilename == 0
    error('bodyCoil data file selection was cancelled');
end

disp('Select hc data file');
[arrayCoilFilename, rawDir] = uigetfile('*.dat', ...
    'Select hc data file', baseFolder);

if arrayCoilFilename == 0
    error('arrayCoil data file selection was cancelled');
end


%% =========================================================================
%% === 2. Initialize Reconstruction Directory ===
% Construct reconstruction path:
%   datasets  →  recon_results
% Extract MID from filename
% =========================================================================

saveflag = 1;  % (Currently unused but kept for compatibility)

% Full paths to selected raw data files
bodyCoilFile  = fullfile(rawDir, bodyCoilFilename);
arrayCoilFile = fullfile(rawDir, arrayCoilFilename);

% Replace "datasets" folder with "recon_results"
[parentDir, ~] = fileparts(rawDir);
parts = strsplit(parentDir, filesep);
parts(strcmp(parts,'datasets')) = {'recon_results'};
parentDir = fullfile(filesep, parts{:});

% Extract MID number (e.g. MID00012)
token = regexp(arrayCoilFile, 'MID\d+', 'match');
MID   = token{1};

% Final reconstruction directory
reconDir = fullfile(parentDir, [MID '_recon_C']);


%% =========================================================================
%% === 3. Load Prescan Sequence Parameters ===
% Extract pulseq sequence parameters
% =========================================================================

seqFile = fullfile(seqFolder, seqName);
seqParams = extract_seq_params(seqFile);
prescan_seqParams = seqParams;   % Return to caller


%% =========================================================================
%% === 4. Configure Body Coil Reader ===
% createRawDataReader unifies Siemens / ISMRMRD formats
% =========================================================================

bodyCoilreader = createRawDataReader(bodyCoilFile, true);

% Optional sequence parameters (if available)
if isfield(seqParams, 'nshot')
    bodyCoilreader.acquisitionParams.nShot = seqParams.nshot;
end

if isfield(seqParams, 'nseg')
    bodyCoilreader.acquisitionParams.nSeg = seqParams.nseg;
end

% Prescan configuration
bodyCoilreader.acquisitionParams.nShot_off = 14;
bodyCoilreader.acquisitionParams.traj_type = 'pulseq';
bodyCoilreader.acquisitionParams.pulseqTrajFile_name = ...
    strcat(seqFolder, seqName);

% Validate pulseq hash consistency
isMatch = check_hash(bodyCoilFile, ...
    bodyCoilreader.acquisitionParams.pulseqTrajFile_name);


%% =========================================================================
%% === 5. Configure Head (Array) Coil Reader ===
% Same configuration as body coil
% =========================================================================

arrayCoilReader = createRawDataReader(arrayCoilFile, true);

if isfield(seqParams, 'nshot')
    arrayCoilReader.acquisitionParams.nShot = seqParams.nshot;
end

if isfield(seqParams, 'nseg')
    arrayCoilReader.acquisitionParams.nSeg = seqParams.nseg;
end

arrayCoilReader.acquisitionParams.nShot_off = 14;
arrayCoilReader.acquisitionParams.traj_type = 'pulseq';
arrayCoilReader.acquisitionParams.pulseqTrajFile_name = ...
    strcat(seqFolder, seqName);

isMatch = check_hash(arrayCoilFile, ...
    arrayCoilReader.acquisitionParams.pulseqTrajFile_name);

% Store for consistency reference
nShotOff = arrayCoilReader.acquisitionParams.nShot_off;


%% =========================================================================
%% === 6. Orientation Consistency Check ===
% Ensure prescan orientation matches main measurement
% =========================================================================

disp('Select main meas sequence file:');
[measSeqName, seqFolder] = uigetfile('*.seq', ...
    'Select main meas sequence file', baseFolder);

if measSeqName == 0
    error('Sequence file selection was cancelled');
end

seqFile = fullfile(seqFolder, measSeqName);
measSeqParams = extract_seq_params(seqFile);

disp('Select main raw data file:');
[measFilename, rawDir] = uigetfile('*.dat', ...
    'Select main raw data file', baseFolder);

if measFilename == 0
    error('measurement data file selection was cancelled');
end

measureFile = fullfile(rawDir, measFilename);
reader = createRawDataReader(measureFile, 1);

if isfield(measSeqParams, 'nshot')
    reader.acquisitionParams.nShot = measSeqParams.nshot;
end

if isfield(measSeqParams, 'nseg')
    reader.acquisitionParams.nSeg = measSeqParams.nseg;
end

reader.acquisitionParams.nShot_off = 20;
reader.acquisitionParams.traj_type = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name = seqFile;

isMatch = check_hash(measureFile, ...
    reader.acquisitionParams.pulseqTrajFile_name);

% Compare orientations
[orientation_ok, xsos_all, ncc_all] = ...
    checkPrescanOrientation(bodyCoilreader, ...
                            arrayCoilReader, ...
                            reader, 64);

bmImage(xsos_all);

if orientation_ok
    disp('Orientation check passed!');
else
    error('Orientation check failed!')
end


%% =========================================================================
%% === 7. Define Reconstruction Grid ===
% Compute Cartesian grid parameters from FoV
% =========================================================================

dK_u = [1, 1, 1] ./ arrayCoilReader.acquisitionParams.FoV;
N_u  = [48, 48, 48];   % Reconstruction grid size


%% =========================================================================
%% === 8. Compute Non-Cartesian Data & Trajectory ===
% Extract k-space data and trajectory information
% =========================================================================

[y_body, t, ve] = bmCoilSense_nonCart_data(bodyCoilreader, N_u);
y_surface       = bmCoilSense_nonCart_data(arrayCoilReader, N_u);


%% =========================================================================
%% === 9. Construct Gridding Operators ===
% Gn  : Uniform → Non-uniform
% Gu  : Non-uniform → Uniform
% Gut : Transpose of Gu
% =========================================================================

[Gn, Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);


%% =========================================================================
%% === 10. Generate Object Mask ===
% Automatically generate mask from body coil image
% (Note: reported to fail on macOS in some cases)
% =========================================================================

mask = bmCoilSense_nonCart_mask_automatic(y_body, Gn, false);


%% =========================================================================
%% === 11. Coil Sensitivity Estimation ===
%
% Two-stage estimation:
%
%   Stage 1: Body coil reference
%   Stage 2: Head coil primary estimate
%   Stage 3: Iterative refinement
%
% =========================================================================

close all;

% --- Reference from body coil ---
[y_ref, C_ref] = bmCoilSense_nonCart_ref(y_body, Gn, mask, []);

% --- Primary estimate (array coil vs body reference) ---
C_array_prime = bmCoilSense_nonCart_primary( ...
                    y_surface, y_ref, C_ref, Gn, ve, mask);

% --- Secondary refinement ---
nIter = 5;
[C1, x] = bmCoilSense_nonCart_secondary( ...
                    y_surface, C_array_prime, ...
                    y_ref, C_ref, ...
                    Gn, Gu, Gut, ve, ...
                    nIter, false);

% Display final sensitivity maps
bmImage(C1);

%% =========================================================================
end