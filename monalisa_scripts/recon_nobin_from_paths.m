function [x, x0, outDir] = recon_nobin_from_paths(opts)
%RECON_NOBIN_FROM_PATHS Reconstruct prepared Mitosius with its matching C.
% Required opts: reconDir, CfilePath, Matrix_size, nIter, delta, rho, flagCS.
% Optional mitosiusDir selects the exact export returned by preparation.

arguments
    opts struct
end
validateattributes(opts.Matrix_size,{'numeric'},{'scalar','integer','positive','even'});
validateattributes(opts.nIter,{'numeric'},{'scalar','integer','positive'});
validateattributes(opts.delta,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(opts.rho,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(opts.flagCS,{'numeric','logical'},{'scalar','binary'});
runTimer=tic;

reconDir = opts.reconDir;
disp(['reconDir: ', reconDir]);
previousDir = pwd;
directoryCleanup = onCleanup(@() cd(previousDir));
cd(reconDir);

mDir = fullfile(reconDir, 'T1_LIBRE_woBinning', 'mitosius');
if isfield(opts,'mitosiusDir') && ~isempty(opts.mitosiusDir), mDir = opts.mitosiusDir; end
disp(['mDir: ', mDir]);

assert(isfolder(mDir), 'Mitosius folder not found: %s', mDir);

CfilePath = opts.CfilePath;
disp(['CfilePath: ', CfilePath]);
assert(isfile(CfilePath), 'C.mat not found: %s', CfilePath);
preparationFile = fullfile(mDir,'preparation.mat');
if isfile(preparationFile)
    assert(isfile(fullfile(mDir,'complete.mat')),'Preparation is incomplete.');
    P = load(preparationFile,'preparation');
    assert(strcmp(file_md5(CfilePath),P.preparation.reconstructionCMD5), ...
        'C differs from the exact C paired with this prepared export.');
else
    warning('Tutorial:LegacyPreparation','Legacy export: exact C identity cannot be verified.');
end

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
geometryFile = fullfile(mDir,'recon_geometry.mat');
if isfile(geometryFile)
    savedGeometry = load(geometryFile,'geometry');
    ReconFov = savedGeometry.geometry.FoV;
else
    answer = inputdlg({'Reconstruction FOV [mm] (including readout oversampling):'}, ...
        'Missing Mitosius geometry',[1 80],{''});
    if isempty(answer), error('Tutorial:SequenceInputCancelled','Reconstruction FOV entry cancelled.'); end
    ReconFov = str2double(answer{1});
end
validateattributes(ReconFov,{'numeric'},{'scalar','real','finite','positive'});
N_u   = [Matrix_size, Matrix_size, Matrix_size];
n_u   = [Matrix_size, Matrix_size, Matrix_size];
dK_u  = [1, 1, 1] ./ ReconFov;
fprintf('Reconstruction FOV: %g mm; matrix: %d; grid spacing: %.6g mm\n', ...
    ReconFov,Matrix_size,ReconFov/Matrix_size);

nFr = size(y, 1);

% ---- Load / resize C ----
S = load(CfilePath);
assert(isfield(S, 'C'), 'C.mat does not contain variable C.');
C = S.C;
assert(isequal([size(C,1),size(C,2),size(C,3)],[48 48 48]) && ...
    all(isfinite(C(:))) && any(abs(C(:))>0),'Invalid C values or spatial grid.');
assert(all(cellfun(@(frame) size(frame,2)==size(C,4),y)), ...
    'Mitosius and C have different coil counts. Select the matching C from preparation.');

C = bmImResize(C, [48, 48, 48], N_u);

% ---- Compute x0 ----
x0 = cell(nFr, 1);
for i = 1:nFr
    x0{i} = bmMathilda(y{i}, t{i}, ve{i}, C, N_u, n_u, dK_u, [], [], [], []);
end

% ---- Save x0 ----
outputRoot = fullfile(reconDir,'T1_LIBRE_woBinning','output');
if isfield(opts,'outputRoot') && ~isempty(opts.outputRoot), outputRoot=opts.outputRoot; end
outDir = fullfile(outputRoot, ...
    ['run_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
assert(~isfolder(outDir),'Reconstruction run already exists.');
mkdir(outDir);
runInfo = struct('options',opts,'CFileMD5',file_md5(CfilePath), ...
    'mitosiusDir',mDir,'FoV',ReconFov,'gridSpacingMm',ReconFov/Matrix_size, ...
    'matlabVersion',version,'complete',false);
save(fullfile(outDir,'run_info.mat'),'runInfo');

x0Path = fullfile(outDir, ['x0_' num2str(Matrix_size) '.mat']);
save(x0Path, 'x0', '-v7.3');
disp(['x0 saved: ', x0Path]);

% ---- Sparse matrices ----
    if opts.flagCS
        [~, Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);
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
runInfo.elapsedSeconds=toc(runTimer); runInfo.complete=true;
save(fullfile(outDir,'run_info.mat'),'runInfo');
fprintf('Reconstruction outputs saved in: %s\n',outDir);
end
