%% =========================================================================
%  FULL RECONSTRUCTION PIPELINE – EDUCATIONAL ENTRY SCRIPT
%
%  This script demonstrates the complete workflow:
%
%   1) (Optional) Quick data quality check (RMS reconstruction)
%   2) Coil sensitivity estimation
%   3) Save coil sensitivity maps
%   4) Prepare data using Mitosius wrapper
%   5) Perform gridded or CS reconstruction
%
%  Recommended folder structure:
%  --------------------------------
%  Place all related files under ONE dataset folder:
%       - Prescan .seq
%       - Raw data .dat
%       - Main measurement .seq
%       - Main measurement .dat
%
%  Place all intermediate files and outputs under ONE recon_result folder:
%       - mitosius
%       - outputs/x.mat
%       - etc
%
% =========================================================================


%% === 1. Define Dataset Base Folder ======================================

% Option A (Recommended):
% Set your dataset directory here.
baseFolder = '/path/to/your/datasets/';

% Option B:
% Leave empty to browse manually each time.
% baseFolder = [];


%% === 1.5 Optional: Quick Data Inspection ================================
% Before estimating coil sensitivities, you may want to quickly inspect
% the raw data quality and orientation consistency.
%
% This reconstructs:
%   xrms : root-mean-square image across coils
%   x0   : initial reconstruction (non-optimized)
%
% This step is useful for:
%   - Checking FOV coverage
%   - Verifying orientation
%   - Detecting severe artifacts
% -------------------------------------------------------------------------

matrix_size = 64;  % Reconstruction grid size for quick inspection

[xrms, x0] = check_orient_xrms(baseFolder, matrix_size);

bmImage(xrms);
bmImage(x0);


%% === 2. Coil Sensitivity Estimation =====================================
% This function will guide you through selecting:
%
%   - Prescan sequence file
%   - Body coil raw data
%   - Head coil raw data
%   - Main measurement sequence (for orientation check)
%   - Main measurement raw data
%
% It also suggests a reconstruction directory (reconDir).
% -------------------------------------------------------------------------

[C1, reconDir, ~, ~, ~] = coilsense_script(baseFolder);

% If preferred, you may override reconDir manually:
% reconDir = 'path/to/your/result/reconDir/';


%% === 3. Save Coil Sensitivity Maps ======================================
% Save the estimated coil sensitivity maps (C1) to disk.
% They will be required for final image reconstruction.
% -------------------------------------------------------------------------

saveCDir = fullfile(reconDir, 'C');

if ~exist(saveCDir, 'dir')
    mkdir(saveCDir);
end

CfileName = 'C.mat';
CfilePath = fullfile(saveCDir, CfileName);

% Save as variable C1 (recommended to keep consistent naming)
save(CfilePath, 'C1');

disp(' ');
disp('====================================================');
disp([' Coil sensitivity saved at: ', CfilePath]);
disp('====================================================');
disp(' ');


%% === 4. Prepare Data (Mitosius Wrapper) =================================
% Mitosius wraps:
%   - Gridded k-space data
%   - Trajectory
%   - Volume elements
%
% into cell format for reconstruction.
%
% This is a preparation step for both:
%   - Standard gridded reconstruction
%   - Compressed sensing reconstruction
% -------------------------------------------------------------------------

mitosius_nobin_ui(baseFolder);


%% === 5. Final Image Reconstruction ======================================
% Select the mitosius file generated in Step 4.
%
% Choose reconstruction mode:
%
%   flagCS = 0 → Standard gridded reconstruction
%   flagCS = 1 → Compressed sensing reconstruction
%
% Outputs:
%   x   : Final reconstructed image
%   x0  : Initial reconstruction (before regularization)
% -------------------------------------------------------------------------

flagCS = 0;   % Set to 1 to enable compressed sensing

[x, x0] = recon_nobin_ui([], [], flagCS);

disp(' ');
disp('====================================================');
disp(' Reconstruction completed successfully.');
disp('====================================================');
disp(' ');