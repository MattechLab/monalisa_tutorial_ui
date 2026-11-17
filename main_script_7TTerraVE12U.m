%% 1. Independent 7T setup (edit tutorial_config_7TTerraVE12U.m)
% Run cell by cell. This script estimates relative C from the MAIN scan.
% It does not use the Prisma BC/HC prescan estimator.
% Default scan outputs: <selected raw dataset folder>/recon_results/.
close all; clc; clear;
repoRoot = fileparts(which('main_script_7TTerraVE12U'));
assert(~isempty(repoRoot),'Open this repository as MATLAB''s Current Folder first.');
addpath(repoRoot);
config = tutorial_config_7TTerraVE12U();
config = setup_tutorial(config);
baseFolder = config.baseFolder;

%% 2. Select inputs and inspect scanner software (header only)
rawFile = fullfile(baseFolder,'meas_MID00019_FID05422_0_25_7T_LIBRE.dat');
seqFile = fullfile(baseFolder, ...
    'aa_seq701_t1_libre_main_TR8.5ms_TE3.8ms_swap1_FA6_Traj1_nSeg44_nShot1055_Fid0_7T_0.25_smaller_FOV_LIBRE.seq');
if ~isfile(rawFile)
    selection_reminder('Terra inputs 1/2: main raw data', ...
        'Select the Terra MAIN raw-data file (.dat) for this acquisition.');
    [name,folder] = uigetfile('*.dat','Select Terra main raw data',baseFolder);
    assert(~isequal(name,0),'Raw-data selection cancelled.');
    rawFile = fullfile(folder,name);
end
if ~isfile(seqFile)
    selection_reminder('Terra inputs 2/2: matching sequence', ...
        sprintf('Select the Pulseq sequence (.seq) used for this raw-data file:\n%s',rawFile));
    [name,folder] = uigetfile('*.seq','Select the matching Terra sequence',baseFolder);
    assert(~isequal(name,0),'Sequence selection cancelled.');
    seqFile = fullfile(folder,name);
end
systemInfo = read_twix_system_info(rawFile);
disp(systemInfo);
% Canonical field names/units match Yiwei's sequence style. [] means absent.
sequenceHeader = read_sequence_header(seqFile);
disp(sequenceHeader);
% Supplied scan: Terra, syngo MR E12, N4_VE12U_LATEST_20181126.
assert(strcmpi(strtrim(systemInfo.scannerModel),'Terra') && ...
    systemInfo.fieldStrengthT>6 && systemInfo.fieldStrengthT<8, ...
    'This entry script expects a 7T Terra acquisition.');

%% 3. Verify sequence and reconstruct low-resolution main-scan coil images
% Verifies the signed payload; added nSeg/nShot definitions are recognized
% only when removing them restores the exact stored/TWIX signature.
assert(check_hash(rawFile,seqFile,systemInfo.headerSignature),'Sequence/header mismatch.');
gridSize = [48 48 48]; % Sensitivity grid, compatible with tutorial C consumers.
nShotOff = 50; % Explicit initial choice: steady-state shots to be discarded.
% Missing nSam/flagSelfNav/etc. are requested in an input dialog by the loader.
% flagSelfNav=1 removes the first readout of each shot; 0 keeps every readout.
[y,t,ve,acq,rawInfo] = load_siemens_lowres(rawFile,seqFile,gridSize,nShotOff);
sequenceHeader = rawInfo.sequenceParams.header;
disp(sequenceHeader); % Includes required values confirmed in the input dialog.
x0 = complex(zeros([gridSize acq.nCh],'single'));
for channel = 1:acq.nCh
    x0(:,:,:,channel) = bmMathilda(y(:,channel),t,ve,[], ...
        gridSize,gridSize,[1 1 1]./acq.FoV);
end
clear y t ve;
xrss = sqrt(sum(abs(x0).^2,4));
xrms = xrss/sqrt(acq.nCh);
bmImage(xrms);

%% 4. Review the signal-support mask (adjust threshold and rerun this cell)
supportThreshold = 0.1; % Fraction of the maximum RSS; excludes weak background.
mask = xrss>supportThreshold*max(xrss(:));
bmImage({xrms,single(mask)});
%%
answer = questdlg('Use this signal-support mask for main-scan C estimation?', ...
    'Review Terra calibration mask','Use mask','Cancel','Cancel');
assert(strcmp(answer,'Use mask'),'Mask not approved; adjust threshold and rerun this cell.');

%% 5. Estimate relative C from the main scan
% Direct calibration: smooth complex images / RSS, then Laplace extension.
% These maps include a common object phase and have no BC intensity reference.
smoothSigma = 1; % Gaussian standard deviation in calibration-grid voxels.
laplaceIterations = 1000;
laplaceTolerance = 1e-4;
[C,calibrationDetails] = estimate_relative_sensitivities( ...
    x0,mask,smoothSigma,laplaceIterations,laplaceTolerance);
bmImage(C);

%% 6. Save C and its main-scan provenance
[~,stem] = fileparts(rawFile);
token = regexp(stem,'MID\d+','match','once');
if ~isempty(token), stem = token; end
runName = [stem '_recon_C_main'];
reconDir = fullfile(config.outputRoot,runName);
assert(~isfolder(reconDir),'Output folder already exists.');
saveCDir = fullfile(reconDir,'C'); mkdir(saveCDir);
CfilePath = fullfile(saveCDir,'C.mat');
save(CfilePath,'C','-v7.3');
save(fullfile(reconDir,'mask.mat'),'mask');
save(fullfile(reconDir,'calibration.mat'),'systemInfo','rawInfo','rawFile', ...
    'seqFile','gridSize','nShotOff','calibrationDetails');
header = rawInfo.sequenceParams.header;
provenance = rawInfo.sequenceParams.headerprovenance;
save(fullfile(reconDir,'main_sequence_header.mat'),'header','provenance');
fprintf('Main-scan relative sensitivities saved: %s\n',CfilePath);
return;

%% Optional: Compute eye-ROI ROVir weights from the same main-scan images
% Run this cell separately after C estimation; draw eye ROI and select z range.
coilComDir = fullfile(reconDir,'coilCom'); %#ok<UNRCH> Optional standalone cell.
[W,lambda,masksOut_rovir] = computeROVirFromEyeROI( ...
    x0,[],[],xrms,'IRMode','complement','DilateROIPixels',1, ...
    'SelectZRange',true,'IRGapPixels',1,'CheckMasks',true);
% Archive the previous results only after the new estimation succeeds.
if isfolder(coilComDir)
    archiveDir = [coilComDir '_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))];
    assert(~isfolder(archiveDir) && ~isfile(archiveDir),'Archive destination already exists.');
    [moved,message] = movefile(coilComDir,archiveDir);
    assert(moved,'Could not archive previous ROVir results: %s',message);
    fprintf('Previous ROVir results archived: %s\n',archiveDir);
end
mkdir(coilComDir);
save(fullfile(coilComDir,'rovir_W.mat'),'W','-v7.3');
save(fullfile(coilComDir,'lambda.mat'),'lambda');
save(fullfile(coilComDir,'masksOut_rovir.mat'),'masksOut_rovir');
save(fullfile(coilComDir,'xrms.mat'),'xrms');
save(fullfile(coilComDir,'x0.mat'),'x0','-v7.3');
fprintf('\nROVir results saved successfully:\n  %s\n',coilComDir);
fprintf('Files: rovir_W.mat, lambda.mat, masksOut_rovir.mat, xrms.mat, x0.mat\n');

%% Optional: Prepare main data for reconstruction (Mitosius)
% Run after saving C. If useROVir=true, first compute/save the ROVir weights.
% Select this same main raw file, sequence, and the physical-coil C.mat.
useROVir = true;
nVirtualCoils = 15; % EDIT HERE: number of retained ROVir channels (1..32).
rovirOptions = struct();
if useROVir
    rovirOptions.weightsFile = fullfile(reconDir,'coilCom','rovir_W.mat');
    rovirOptions.nVirtualCoils = nVirtualCoils;
end
preparationOptions = struct();
preparationOptions.matrixSize = 120; % Normalization grid; final grid is set below.
preparationOptions.nShotOff = nShotOff; % Reuse the main-scan choice from cell 3.
preparationOptions.sequenceParams = rawInfo.sequenceParams;
[mitosiusDir,reconCfilePath] = mitosius_nobin_ui( ...
    baseFolder,config.outputRoot,rovirOptions,preparationOptions);

%% Optional: Reconstruct the prepared main data
% flagCS=0: gridded reconstruction; flagCS=1: compressed sensing.
% The dialog sets reconstruction matrix, nIter, delta, and rho.
% For FOV=200 mm, matrix 120 gives a 1.67-mm preview grid; 400 gives 0.5 mm.
mitosiusDir = '/Volumes/Extreme SSD/Arina/260806/0p25mm_actual0p5mm/recon_results/MID00019_recon/T1_LIBRE_woBinning/mitosius_nch15';
reconCfilePath ="/Volumes/Extreme SSD/Arina/260806/0p25mm_actual0p5mm/recon_results/MID00019_recon/T1_LIBRE_woBinning/mitosius_nch15/C_compressed_15.mat";
flagCS = 1;
assert(exist('mitosiusDir','var')==1 && isfolder(mitosiusDir) && ...
    exist('reconCfilePath','var')==1 && isfile(reconCfilePath), ...
    'Run the Mitosius cell first to obtain its matching data and C.');
[x_recon,x0_recon] = recon_nobin_ui( ...
    config.outputRoot,baseFolder,flagCS,reconCfilePath,mitosiusDir);

%% Anisotropic scanner xRMS: seq712-714, independent of calibration/setup cells
% Run this cell alone. The separate scanner pipeline needs no phantom/calibration.
anisotropicRepoRoot = fileparts(which('main_script_7TTerraVE12U'));
assert(~isempty(anisotropicRepoRoot),'Open this repository as MATLAB''s Current Folder first.');
addpath(anisotropicRepoRoot);
anisotropicDataset = '/Volumes/Extreme SSD/261110';
anisotropicShotOff = 20;
anisotropicKernelSigmas = [.61 .32 .20 .10]; % Retained baseline/comparison settings.
anisotropicShowImages = true;
anisotropicPlan = run_anisotropic_scanner(anisotropicDataset, ...
    anisotropicKernelSigmas,anisotropicShotOff,[],anisotropicShowImages);
anisotropicRunDir = anisotropicPlan.outDir;

%% Isotropic scanner xRMS: seq701-705, independent of calibration/setup cells
% Run this cell alone: 240^3 matrix, 240-mm FOV, 1-mm isotropic output voxels.
isotropicRepoRoot = fileparts(which('main_script_7TTerraVE12U'));
assert(~isempty(isotropicRepoRoot),'Open this repository as MATLAB''s Current Folder first.');
addpath(isotropicRepoRoot);
isotropicDataset = '/Volumes/Extreme SSD/261110';
isotropicShotOff = 50; % Same initial-shot exclusion as the 7T main-scan cell.
isotropicPlan = run_isotropic_scanner(isotropicDataset,[],isotropicShotOff);
isotropicRunDir = isotropicPlan.outDir;
