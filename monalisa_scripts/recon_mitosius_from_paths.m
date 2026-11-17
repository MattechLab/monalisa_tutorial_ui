function [mDir,reconstructionCFile] = recon_mitosius_from_paths(opts)
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
[~,base] = fileparts(opts.measureFile);
token = regexp(base, 'MID\d+', 'match');
assert(~isempty(token), 'Could not extract MID\d+ from measurement file path/name.');
measMID = token{1};
disp(['measMID = ', measMID]);

% ---- Output recon directory ----
reconDir = fullfile(opts.reconRoot, [measMID '_recon']);
disp(['reconDir = ', reconDir]);
mDir = fullfile(reconDir,'T1_LIBRE_woBinning','mitosius');
assert(~isfolder(mDir), ...
    'Mitosius output exists (possibly partial); choose a new task/output root.');
if ~isfield(opts,'rovirOptions'), opts.rovirOptions=struct(); end
if ~isfield(opts,'requireProvenance'), opts.requireProvenance=false; end
preparationTimer=tic;
reconstructionCFile = opts.CfilePath;

if ~isfolder(reconDir)
    mkdir(reconDir);
    disp(['Directory created: ', reconDir]);
else
    disp(['Directory already exists: ', reconDir]);
end
previousDir = pwd;
directoryCleanup = onCleanup(@() cd(previousDir));
cd(reconDir);

otherDir = fullfile(reconDir, 'T1_LIBRE_woBinning', 'other');
if ~isfolder(otherDir)
    mkdir(otherDir);
end

% ---- Resolve physical grouping from the sequence before loading samples ----
seqParams = [];
if isfield(opts,'sequenceParams') && ...
        strcmp(opts.sequenceParams.headerprovenance.sourceFile,opts.seqFile)
    seqParams = opts.sequenceParams;
end
if isfield(opts,'nShot_off') && ~isempty(opts.nShot_off)
    nShotOff = opts.nShot_off;
else
    nShotOff = 50;
end
[y_tot,acq,rawInfo,geometry] = load_siemens_main( ...
    opts.measureFile,opts.seqFile,nShotOff,seqParams);
seqParams = rawInfo.sequenceParams;
t_tot = bmTraj(acq);
assert(numel(t_tot)/3==acq.N*size(y_tot,3),'Trajectory/readout layout mismatch.');

sourceC = load(opts.CfilePath,'C','mainIdentity');
assert(isfield(sourceC,'C'),'C file does not contain C.');
C = sourceC.C;
assert(isequal([size(C,1),size(C,2),size(C,3)],[48 48 48]) && ...
    size(C,4)==acq.nCh && all(isfinite(C(:))) && any(abs(C(:))>0), ...
    'C must be finite physical-coil 48^3 maps matching the raw data.');
currentIdentity = main_scan_identity(rawInfo);
if isfield(sourceC,'mainIdentity')
    assert(isequal(sourceC.mainIdentity,currentIdentity),'C belongs to another main acquisition or inputs changed.');
else
    assert(~opts.requireProvenance,'C lacks task provenance; save it through the task calibration/reuse cells.');
    warning('Tutorial:MissingCProvenance','C provenance unavailable; only dimensions and values checked.');
end
disp(['C loaded from: ', opts.CfilePath]);

% ---- Reconstruction parameters ----
if isfield(opts, 'matrix_size') && ~isempty(opts.matrix_size)
    matrix_size = opts.matrix_size;
else
    matrix_size = 120;
end
validateattributes(matrix_size,{'numeric'},{'scalar','integer','positive','even','>=',4});

N_u = [matrix_size, matrix_size, matrix_size];
dK_u  = [1 1 1] ./ geometry.FoV;

Cnormalization = bmImResize(C, [48, 48, 48], N_u);
x_tot = bmMathilda(y_tot, t_tot, bmVolumeElement(t_tot,'voronoi_full_radial3'), Cnormalization, N_u, N_u, dK_u);
clear Cnormalization;

temp_im = x_tot(round(matrix_size/4):round(matrix_size/4*3), ...
                round(matrix_size/4):round(matrix_size/4*3), ...
                round(matrix_size/2));
normalize_val = mean(temp_im(:));
disp(['normalize_val = ', num2str(normalize_val)]);

% NOTE: your original "if real(y_tot)<1" is ambiguous (y_tot array).
% I keep the intent but make it deterministic:
testval = real(y_tot(1));  % choose one sample
normalization = struct('matrixSize',matrix_size,'value',normalize_val, ...
    'firstSampleReal',testval,'applied',testval<1,'method','legacy_complex_central_slice_mean');
if normalization.applied
    assert(isfinite(normalize_val) && abs(normalize_val)>0,'Normalization divisor is zero or nonfinite.');
    y_tot = y_tot / normalize_val;
end
assert(all(isfinite(y_tot(:))),'Normalization produced nonfinite data.');
clear x_tot temp_im;

% ---- ROVir: same QR/conjugation as HPC method 3, after normalization ----
if isfield(opts.rovirOptions,'weightsFile') && ~isempty(opts.rovirOptions.weightsFile)
    weights = load(opts.rovirOptions.weightsFile,'W','rovirProvenance');
    assert(isfield(weights,'W'),'ROVir weights file must contain W.');
    if isfield(weights,'rovirProvenance')
        rp=weights.rovirProvenance;
        assert(isequal(rp.mainIdentity,currentIdentity) && rp.nShotOff==nShotOff, ...
            'ROVir weights do not match this acquisition/shot exclusion.');
        assert(strcmp(rp.sourceCMD5,file_md5(opts.CfilePath)), ...
            'ROVir source C differs; review weights for this calibration.');
    else
        assert(~opts.requireProvenance,'ROVir weights lack task provenance; rerun the ROVir cell.');
        warning('Tutorial:MissingROVirProvenance','ROVir provenance unavailable; only dimensions checked.');
    end
    [y_tot,C,ccInfo] = apply_rovir_compression(y_tot,C,weights.W,opts.rovirOptions.nVirtualCoils);
    ccInfo.weightsFile = opts.rovirOptions.weightsFile;
    ccInfo.sourceCFile = opts.CfilePath;
    if ~isfolder(mDir), mkdir(mDir); end
    reconstructionCFile = fullfile(mDir,sprintf('C_compressed_%d.mat',size(C,4)));
    save(reconstructionCFile,'C','-v7.3');
    save(fullfile(mDir,'rovir_compression.mat'),'ccInfo','-v7.3');
    fprintf('ROVir: %d physical coils -> %d virtual coils. C saved: %s\n', ...
        ccInfo.nPhysicalCoils,ccInfo.nVirtualCoils,reconstructionCFile);
end

% ---- eMask ----
eMask = true(1, seqParams.nshot*seqParams.nseg);
eMask(1:acq.nShot_off*seqParams.nseg) = false;

eMaskFilePath = fullfile(otherDir, 'eMask_woBin.mat');
save(eMaskFilePath, 'eMask');
disp(['eMask saved: ', eMaskFilePath]);

% Select exactly the retained imaging readouts, in original acquisition order.
eyeMask = eMask;
eyeMask = reshape(eyeMask, [numel(eyeMask)/ (seqParams.nseg*seqParams.nshot), seqParams.nseg, seqParams.nshot]); 
if geometry.excludeSI, eyeMask(:,1,:) = []; end
eyeMask(:, :, 1:acq.nShot_off) = [];
eyeMask = bmPointReshape(eyeMask);

% ---- mitosius saving ----
mDir = fullfile(reconDir, 'T1_LIBRE_woBinning', 'mitosius');
if ~isfolder(mDir); mkdir(mDir); end
disp(['mDir = ', mDir]);

[y, t] = bmMitosis(y_tot, t_tot, eyeMask);
y = bmPermuteToCol(y);
ve  = bmVolumeElement(t, 'voronoi_full_radial3');

bmMitosius_create(mDir, y, t, ve);
header = seqParams.header; provenance = seqParams.headerprovenance;
save(fullfile(mDir,'sequence_header.mat'),'header','provenance');
save(fullfile(mDir,'recon_geometry.mat'),'geometry');
save(fullfile(mDir,'raw_info.mat'),'rawInfo');
preparation = struct('options',opts,'normalization',normalization, ...
    'sourceCMD5',file_md5(opts.CfilePath),'reconstructionCFile',reconstructionCFile, ...
    'reconstructionCMD5',file_md5(reconstructionCFile),'mainIdentity',currentIdentity, ...
    'elapsedSeconds',toc(preparationTimer),'matlabVersion',version);
save(fullfile(mDir,'preparation.mat'),'preparation');
complete = true;
save(fullfile(mDir,'complete.mat'),'complete');
disp(['Mitosius files saved: ', mDir]);
fprintf('Matching sensitivity file for reconstruction: %s\n',reconstructionCFile);

end
