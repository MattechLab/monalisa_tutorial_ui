%% 1. Setup a NEW task: edit tutorial_task_config_3T.m before running
% Run cell by cell. Re-running setup clears old variables and makes a fresh
% task folder; different tasks in one dataset cannot overwrite each other.
% Machine/dependency/dataset paths: tutorial_config.m.
% Task inputs, calibration, compression, and grids: tutorial_task_config_3T.m.
% For a quick raw-data preview only, run the LAST cell independently.
close all; clc; clear;
repoRoot = fileparts(which('main_script_3TPrismaXA60'));
assert(~isempty(repoRoot),'Open this repository as MATLAB''s Current Folder first.');
addpath(repoRoot);
task = tutorial_task_config_3T();
calibrationMethods = {'prescan','main','reuse'};
calibrationLabels = {'BC/HC prescans','Main scan only','Reuse saved C'};
defaultMethod = find(strcmp(task.calibrationMethod,calibrationMethods),1);
assert(~isempty(defaultMethod),'Unknown calibration method in tutorial_task_config_3T.m.');
calibrationChoice = questdlg('How do you want to obtain coil sensitivities for this task?', ...
    'Use BC/HC prescans?',calibrationLabels{:},calibrationLabels{defaultMethod});
assert(~isempty(calibrationChoice),'Calibration selection cancelled.');
task.calibrationMethod = calibrationMethods{strcmp(calibrationChoice,calibrationLabels)};
config = setup_tutorial(tutorial_config());
baseFolder = config.baseFolder;
task = resolve_task_inputs(task,baseFolder);
taskRoot = fullfile(config.outputRoot, ...
    [task.name '_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
assert(~isfolder(taskRoot),'Task output already exists.');
mkdir(taskRoot);
save(fullfile(taskRoot,'task_config.mat'),'task','config');
fprintf('Task output: %s\n',taskRoot);
disp(task);

%% 2. Estimate C (or load existing physical-coil C)
switch task.calibrationMethod
    case 'prescan'
        coilOptions = struct('outputRoot',taskRoot,'files',task.files, ...
            'maskFile',task.maskFile,'showFigures',true, ...
            'prescanShotOff',task.prescanShotOff,'mainShotOff',task.mainShotOff, ...
            'nIter',task.prescanIterations,'orientation',task.orientation);
        [C,reconDir,prescan_seqParams,xsos_all,ncc_all,mainInfo] = ...
            coilsense_script(baseFolder,coilOptions);
    case 'main'
        % Same direct smoothed-RSS/Laplace estimator as 7T; not ESPIRiT.
        mainOptions = task.mainCalibration;
        mainOptions.files = task.files;
        mainOptions.maskFile = task.maskFile;
        [C,reconDir,mainRms,coilImages,mainInfo] = ...
            coilsense_main_script(baseFolder,taskRoot,task.mainShotOff,mainOptions);
    case 'reuse'
        S = load(task.existingCFile,'C','mainIdentity');
        assert(isfield(S,'C'),'Existing C file must contain C.');
        C = S.C;
        reconDir = fullfile(taskRoot,'reused_C'); mkdir(reconDir);
        [mainRms,coilImages,~,mainInfo] = check_orient_xrms( ...
            baseFolder,task.previewMatrix,task.mainShotOff,task.files);
        if isfield(S,'mainIdentity')
            assert(isequal(S.mainIdentity,main_scan_identity(mainInfo)), ...
                'Existing C belongs to another main acquisition.');
        else
            bmImage(C);
            answer=questdlg('This C has no acquisition provenance. Confirm its physical coil order and spatial correspondence to this main scan.', ...
                'Review reused calibration','Accept','Cancel','Cancel');
            assert(strcmp(answer,'Accept'),'Reused C not approved.');
        end
end
mainIdentity = main_scan_identity(mainInfo);
assert(isequal([size(C,1),size(C,2),size(C,3)],[48 48 48]) && ...
    size(C,4)==mainInfo.nCh && all(isfinite(C(:))) && any(abs(C(:))>0), ...
    'C must be finite, nonzero, physical-coil 48^3 maps matching the main scan.');

%% 3. Save physical-coil C with task/acquisition provenance
saveCDir = fullfile(reconDir,'C');
if ~isfolder(saveCDir), mkdir(saveCDir); end
CfilePath = fullfile(saveCDir,'C.mat');
assert(~isfile(CfilePath),'C.mat already exists; start a new task.');
save(CfilePath,'C','mainIdentity','task','-v7.3');
fprintf('Coil sensitivity saved: %s\n',CfilePath);
% Run subsequent cells individually. Whole-script execution stops here.
return;

%% Optional: Main-data inspection / physical-coil images for ROVir
% Run after setup; C estimation is not required for this preview.
% For inspection without task setup, use the LAST cell.
assert(exist('task','var')==1 && exist('baseFolder','var')==1, ...
    'Run setup cell 1 first, or use the last cell for a standalone preview.');
previewParams = [];
if exist('mainInfo','var')==1
    previewParams = mainInfo.sequenceParams;
end
[mainRms,coilImages,~,previewInfo] = check_orient_xrms( ...
    baseFolder,task.previewMatrix,task.mainShotOff,task.files,previewParams);
if exist('mainIdentity','var')==1
    assert(isequal(main_scan_identity(previewInfo),mainIdentity),'Main scan changed.');
end
bmImage(mainRms); bmImage(coilImages);

%% Optional: Estimate ROVir weights from this task's physical-coil images
% Prescan mode: run the inspection cell first. Main/reuse mode returns images.
assert(exist('coilImages','var')==1 && exist('mainRms','var')==1, ...
    'Run the main-data inspection cell first.');
coilComDir = fullfile(reconDir,'coilCom');
[W,lambda,masksOut_rovir] = computeROVirFromEyeROI( ...
    coilImages,[],[],mainRms,'IRMode','complement','DilateROIPixels',1, ...
    'SelectZRange',true,'ObjectThresh',0.01,'IRGapPixels',1,'CheckMasks',true);
% Preserve an earlier review only after the new estimation succeeds.
if isfolder(coilComDir)
    archiveDir = [coilComDir '_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))];
    assert(~isfolder(archiveDir),'Archive already exists.');
    [ok,message] = movefile(coilComDir,archiveDir); assert(ok,'%s',message);
end
mkdir(coilComDir);
rovirProvenance = struct('mainIdentity',mainIdentity,'nShotOff',task.mainShotOff, ...
    'sourceCFile',CfilePath,'sourceCMD5',file_md5(CfilePath),'taskName',task.name);
save(fullfile(coilComDir,'rovir_W.mat'),'W','rovirProvenance','-v7.3');
save(fullfile(coilComDir,'masksOut_rovir.mat'),'masksOut_rovir');
save(fullfile(coilComDir,'lambda.mat'),'lambda');
% Preserve established saved variable names for consumers of these files.
xrms = mainRms; x0 = coilImages;
save(fullfile(coilComDir,'xrms.mat'),'xrms');
save(fullfile(coilComDir,'x0.mat'),'x0','-v7.3');

%% Optional: Prepare this task's exact main data and physical-coil C
rovirOptions = struct();
if task.useROVir
    rovirOptions.weightsFile = fullfile(reconDir,'coilCom','rovir_W.mat');
    rovirOptions.nVirtualCoils = task.nVirtualCoils;
end
preparationOptions = struct('matrixSize',task.normalizationMatrix, ...
    'nShotOff',task.mainShotOff,'sequenceParams',mainInfo.sequenceParams, ...
    'files',task.files,'CfilePath',CfilePath,'datasetDir',baseFolder, ...
    'showSettings',false,'requireProvenance',true);
[mitosiusDir,reconCfilePath] = mitosius_nobin_ui( ...
    baseFolder,taskRoot,rovirOptions,preparationOptions);

%% Optional: Reconstruct prepared data with explicit task settings
if ~isempty(task.preparedMitosiusDir)
    mitosiusDir = task.preparedMitosiusDir;
    reconCfilePath = task.preparedCFile;
end
assert(exist('mitosiusDir','var')==1 && isfolder(mitosiusDir) && ...
    exist('reconCfilePath','var')==1 && isfile(reconCfilePath), ...
    'Prepare data first, or configure both prepared paths for an existing export.');
reconstructionOptions = task.reconstruction;
reconstructionOptions.showSettings = false;
reconstructionOptions.outputRoot = fullfile(taskRoot,'reconstructions');
[x_recon,x0_recon,reconstructionRunDir] = recon_nobin_ui( ...
    taskRoot,baseFolder,task.reconstruction.flagCS,reconCfilePath, ...
    mitosiusDir,reconstructionOptions);

%% Quick raw-data preview ONLY: run this cell independently, even after clear
% No task configuration, BC/HC prescans, or coil sensitivities needed.
% Choose the dataset folder, main .seq, and main .dat when prompted.
% Missing sequence metadata may still require confirmation.
% These variables are separate from an existing reconstruction task.
quickPreviewRepoRoot = fileparts(which('main_script_3TPrismaXA60'));
assert(~isempty(quickPreviewRepoRoot),'Open this repository as MATLAB''s Current Folder first.');
addpath(quickPreviewRepoRoot);
quickPreviewConfig = setup_tutorial(tutorial_config());
quickPreviewMatrix = 120; % Same preview setting as the task default
quickPreviewShotOff = 20; % Initial shots to discard; edit for this scan
[quickPreviewRms,quickPreviewCoilImages,~,quickPreviewInfo] = check_orient_xrms( ...
    quickPreviewConfig.baseFolder,quickPreviewMatrix,quickPreviewShotOff);
bmImage(quickPreviewRms); bmImage(quickPreviewCoilImages);
