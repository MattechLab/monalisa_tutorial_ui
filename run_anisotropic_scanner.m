function plan = run_anisotropic_scanner(datasetFolder,kernelSigmas,nShotOff,outDir,showImages,checkOnly)
%RUN_ANISOTROPIC_SCANNER Independent seq712-714 scanner xRMS pipeline.
% plan=run_anisotropic_scanner(datasetFolder,[.61 .32 .20 .10],20)
% Optional outDir chooses a NEW batch folder; showImages defaults false.
% checkOnly=true validates inputs/output collisions without creating files,
% initializing MONALISA, reading raw signal, or reconstructing.
if nargin<2 || isempty(kernelSigmas), kernelSigmas=[.61 .32 .20 .10]; end
if nargin<3 || isempty(nShotOff), nShotOff=20; end
if nargin<4 || isempty(outDir)
    outDir=fullfile(datasetFolder,'recon_results','anisotropic_scanner', ...
        ['run_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
end
if nargin<5, showImages=false; end
if nargin<6, checkOnly=false; end
validateattributes(kernelSigmas,{'numeric'},{'vector','real','finite','positive'});
validateattributes(nShotOff,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(showImages,{'logical'},{'scalar'});
validateattributes(checkOnly,{'logical'},{'scalar'});
assert(isfolder(datasetFolder),'Dataset folder is unavailable: %s',datasetFolder);
assert(~isfolder(outDir) && ~isfile(outDir),'Choose a new batch output folder: %s',outDir);
kernelSigmas=kernelSigmas(:).';
assert(numel(unique(kernelSigmas))==numel(kernelSigmas),'Kernel widths must be unique.');
tokens=arrayfun(@(s)strrep(sprintf('%.15g',s),'.','p'),kernelSigmas,'UniformOutput',false);
assert(numel(unique(tokens))==numel(tokens),'Kernel widths create colliding output labels.');
plan=struct('pipeline','anisotropic_scanner_xrms','datasetFolder',datasetFolder, ...
    'outDir',outDir,'sequenceIds',712:714,'kernelSigmas',kernelSigmas, ...
    'nShotOff',nShotOff,'kernelWindow',3,'jobs',[]);
% Resolve EVERY input before writing outputs or launching reconstruction.
for id=712:714
    seq=dir(fullfile(datasetFolder,sprintf('yj_seq%d_*.seq',id)));
    raw=dir(fullfile(datasetFolder,sprintf('meas_*seq%d.dat',id)));
    assert(isscalar(seq) && isscalar(raw), ...
        'Expected exactly one matching .seq and .dat for seq%d.',id);
    for k=1:numel(kernelSigmas)
        job=struct('sequenceId',id,'rawFile',fullfile(raw.folder,raw.name), ...
            'seqFile',fullfile(seq.folder,seq.name),'sigma',kernelSigmas(k), ...
            'outDir',fullfile(outDir,sprintf('seq%d_sigma%s',id,tokens{k})));
        plan.jobs=[plan.jobs job];
    end
end
fprintf('Scanner pipeline: %d jobs; sigma=%s; nShotOff=%d\n', ...
    numel(plan.jobs),mat2str(kernelSigmas),nShotOff);
fprintf('New batch folder: %s\n',outDir);
if checkOnly, return; end
repoRoot=fileparts(mfilename('fullpath')); addpath(repoRoot);
config=tutorial_config_7TTerraVE12U();
config.baseFolder=datasetFolder;
config.outputRoot=fullfile(datasetFolder,'recon_results');
setup_tutorial(config);
mkdir(outDir); save(fullfile(outDir,'batch_plan.mat'),'plan');
for j=1:numel(plan.jobs)
    job=plan.jobs(j);
    [xrms,~,~]=recon_anisotropic_xrms(job.rawFile,job.seqFile,job.outDir,nShotOff,job.sigma);
    encodingFile=fullfile(outDir,'encoding',sprintf('seq%d.mat',job.sequenceId));
    export_anisotropic_nifti(fullfile(job.outDir,'reconstruction.mat'),job.outDir,encodingFile);
    % The exporter verifies the NIfTI before this explicit space-saving step.
    delete(fullfile(job.outDir,'reconstruction.mat'));
    if showImages, bmImage(xrms); end
end
fprintf('Anisotropic scanner batch complete: %s\n',outDir);
end
