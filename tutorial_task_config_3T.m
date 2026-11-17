function task = tutorial_task_config_3T()
%TUTORIAL_TASK_CONFIG_3T EDIT HERE to start a different task in one dataset.
% Machine/dependency and dataset-root paths remain in tutorial_config.m.
% Paths below are absolute or relative to config.baseFolder. Empty input
% paths open a picker ONCE in setup; resolved paths are reused by all stages.
task.name = 'prisma_task'; % e.g. static_fullcoil, motion_rovir25, mainC_test
task.files.mainRaw = ''; % Main Siemens .dat
task.files.mainSeq = ''; % Matching Pulseq .seq
task.files.prescanSeq = ''; % BC/HC .seq (only for prescan calibration)
task.files.bodyRaw = ''; % BC .dat (only for prescan calibration)
task.files.headRaw = ''; % HC .dat (only for prescan calibration)

task.calibrationMethod = 'prescan'; % Default dialog choice: 'prescan', 'main', or 'reuse'
task.existingCFile = ''; % Physical-coil C.mat; empty opens a picker for 'reuse'
task.maskFile = ''; % Optional approved logical 48^3 mask MAT
task.mainShotOff = 50;
task.prescanShotOff = 14;
task.previewMatrix = 60; % Main coil images for inspection/ROVir
task.mainCalibration.supportThreshold = 0.1;
task.mainCalibration.smoothSigma = 1; % Calibration voxels
task.mainCalibration.laplaceIterations = 1000; % Batch, not total limit
task.mainCalibration.laplaceTolerance = 1e-4;
task.prescanIterations = 5;
% Review metrics are array-axis diagnostics, not head-motion registration.
% These optional gates need acquisition-specific validation; [] disables.
task.orientation.minimumCorrelation = [];
task.orientation.minimumMargin = [];

task.useROVir = true;
task.nVirtualCoils = 25;
task.normalizationMatrix = 120; % 60 or 120; existing normalization retained
task.reconstruction.Matrix_size = 120; % Final grid, independent of normalization
task.reconstruction.nIter = 20;
task.reconstruction.delta = 1;
task.reconstruction.rho = 10;
task.reconstruction.flagCS = 1; % 0: gridded; 1: CS
% Optional: reconstruct an existing matched export instead of preparing anew.
task.preparedMitosiusDir = '';
task.preparedCFile = ''; % Exact C paired with that export (compressed if used)
end
