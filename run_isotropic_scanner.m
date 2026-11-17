function plan = run_isotropic_scanner(datasetFolder,outDir,nShotOff,scanSet)
%RUN_ISOTROPIC_SCANNER Physical-coil scanner xRMS on a 240-cubed grid.
% scanSet: 'seq701_705' (default), 'seq702_705to709' (six selected scans),
% 'seq700_norm_shimming', or 'seq700_linear_shimming' (GRE/LIBRE pairs).
% Uses the 7T main-scan load_siemens_lowres/bmMathilda path, including its
% Pulseq trajectory normalization and Cartesian k-space selection.
% Saves verified max-normalized float32 NIfTI, JSON, logs and previews.
if nargin<4 || isempty(scanSet), scanSet='seq701_705'; end
scanSet=validatestring(scanSet,{'seq701_705','seq702_705to709', ...
    'seq700_norm_shimming','seq700_linear_shimming'});
if nargin<2 || isempty(outDir)
    outDir=fullfile(datasetFolder,'recon_results','isotropic_scanner', ...
        ['run_' scanSet '_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
end
if nargin<3 || isempty(nShotOff), nShotOff=50; end
validateattributes(nShotOff,{'numeric'},{'scalar','integer','nonnegative','<',500});
assert(isfolder(datasetFolder),'Dataset is unavailable: %s',datasetFolder);
assert(~isfolder(outDir) && ~isfile(outDir),'Choose a new output folder: %s',outDir);
repoRoot=fileparts(mfilename('fullpath')); addpath(repoRoot);
config=tutorial_config_7TTerraVE12U();
config.baseFolder=datasetFolder; config.outputRoot=fileparts(outDir);
setup_tutorial(config);
plan=struct('pipeline','isotropic_scanner_xrms','datasetFolder',datasetFolder, ...
    'outDir',outDir,'scanSet',scanSet,'matrixSize',[240 240 240],'nShotOff',nShotOff, ...
    'kernel','gauss','kernelWindow',3,'kernelParameters',[.61 10], ...
    'status','planned','jobs',[]);
if any(strcmp(scanSet,{'seq701_705','seq702_705to709'}))
    if strcmp(scanSet,'seq701_705')
        ids=701:705;
    else
        ids=[702 706 707 705 708 709]; % Acquisition order in the user's selection.
    end
    expectedShots=500;
    labels=arrayfun(@(id)sprintf('seq%d',id),ids,'UniformOutput',false);
    rawPatterns=arrayfun(@(id)sprintf('meas_*seq%d.dat',id),ids,'UniformOutput',false);
    seqPatterns=arrayfun(@(id)sprintf('yj_seq%d_*.seq',id),ids,'UniformOutput',false);
else
    ids=[700 700]; expectedShots=199;
    if strcmp(scanSet,'seq700_norm_shimming'), shimLabel='norm'; else, shimLabel='linear'; end
    labels={['seq700_gre_' shimLabel '_shimming'],['seq700_libre_' shimLabel '_shimming']};
    rawPatterns={['meas_*' labels{1} '.dat'],['meas_*' labels{2} '.dat']};
    seqPatterns={'yj_seq700_t1w_gre_*.seq','yj_seq700_t1w_libre_*.seq'};
end
% Verify all source pairings and dimensions before starting any reconstruction.
for k=1:numel(ids)
    id=ids(k);
    raw=dir(fullfile(datasetFolder,rawPatterns{k}));
    seq=dir(fullfile(datasetFolder,seqPatterns{k}));
    assert(isscalar(raw) && isscalar(seq),'Expected one .dat/.seq pair for %s.',labels{k});
    rawFile=fullfile(raw.folder,raw.name); seqFile=fullfile(seq.folder,seq.name);
    [~,acq,info,geometry]=siemens_sequence_metadata(rawFile,seqFile,nShotOff);
    assert(isfield(info,'B0_T') && info.B0_T>6 && info.B0_T<8,'Expected 7T data.');
    assert(acq.FoV==240 && acq.N==480 && acq.nShot==expectedShots && acq.nSeg==44 && ...
        geometry.excludeSI,'Unexpected acquisition geometry for %s.',labels{k});
    job=struct('sequenceId',id,'label',labels{k},'rawFile',rawFile,'seqFile',seqFile, ...
        'outDir',fullfile(outDir,labels{k}),'status','planned');
    plan.jobs=[plan.jobs job]; %#ok<AGROW>
end
mkdir(outDir); plan.status='running'; save(fullfile(outDir,'batch_plan.mat'),'plan');
for j=1:numel(plan.jobs)
    reconstruct_one(plan.jobs(j),plan);
    plan.jobs(j).status='complete';
    save(fullfile(outDir,'batch_plan.mat'),'plan');
end
plan.status='complete'; save(fullfile(outDir,'batch_plan.mat'),'plan');
fprintf('ISOTROPIC SCANNER BATCH COMPLETE: %s\n',outDir);
end

function reconstruct_one(job,plan)
mkdir(job.outDir); diary(fullfile(job.outDir,'reconstruction.log'));
logCleanup=onCleanup(@()diary('off')); %#ok<NASGU>
started=tic; seqHash=file_md5(job.seqFile);
grid=plan.matrixSize;
[y,t,ve,acq,metadata]=load_siemens_lowres(job.rawFile,job.seqFile,grid,plan.nShotOff);
assert(all(isfinite(y),'all') && all(ve>0,'all'));
sumSquares=zeros(grid,'single');
% Combine incrementally so full complex coil volumes are never stored.
for channel=1:acq.nCh
    coilImage=bmMathilda(y(:,channel),t,ve,[],grid,grid,[1 1 1]/acq.FoV, ...
        plan.kernel,plan.kernelWindow,plan.kernelParameters,'MATLAB',true);
    assert(isequal(size(coilImage),grid) && all(isfinite(coilImage),'all'));
    sumSquares=sumSquares+real(coilImage.*conj(coilImage));
    fprintf('%s coils reconstructed %d/%d\n',job.label,channel,acq.nCh);
end
clear coilImage;
xrms=sqrt(sumSquares/acq.nCh);
assert(all(isfinite(xrms),'all') && max(xrms(:))>0);
metadata.sequenceId=job.sequenceId; metadata.sequenceMD5=seqHash;
metadata.sequenceLabel=job.label; metadata.scanSet=plan.scanSet;
% Preserve the source definition and separately record the user's LIBRE convention.
header=metadata.sequenceParams.header;
metadata.sequenceFileFlipAngleDeg=header.rf_flip_deg;
if header.rfmode==2
    metadata.reportedFlipAngleDeg=2*header.rf_flip_deg;
    metadata.reportedFlipAngleBasis='User convention: twice sequence-file rf_flip_deg for LIBRE; not a measured local flip angle';
else
    metadata.reportedFlipAngleDeg=header.rf_flip_deg;
    metadata.reportedFlipAngleBasis='Sequence-file rf_flip_deg';
end
metadata.matrixSize=grid; metadata.fovMm=repmat(acq.FoV,1,3);
metadata.voxelSizeMm=metadata.fovMm./grid;
metadata.acquisitionNominalVoxelSizeMm=metadata.sequenceParams.header.res_mm;
metadata.dK_cycles_per_mm=[1 1 1]/acq.FoV;
metadata.trajectoryConvention='7T main-scan bmTraj/mlTrajFromPulseq: normalize max radius to N/(2*FoV); then bmLowRes Cartesian grid selection';
metadata.method='MONALISA gridding; sqrt(mean(abs(physicalCoilImage).^2))';
metadata.densityCompensation='voronoi_full_radial3 on retained imaging spokes, before bmLowRes selection';
metadata.retainedShots=acq.nShot-plan.nShotOff;
metadata.imagingReadouts=(acq.nSeg-acq.selfNav_flag)*metadata.retainedShots;
metadata.samplesUsed=size(y,1); metadata.physicalCoils=acq.nCh;
metadata.kernel=plan.kernel; metadata.kernelWidth=plan.kernelWindow;
metadata.kernelParameters=plan.kernelParameters;
metadata.signalConjugated=false; metadata.additionalADCPhaseCorrection=false;
metadata.gradientDelayCorrection=false; metadata.coilCompression=false;
metadata.matlabVersion=version; metadata.runnerMD5=file_md5([mfilename('fullpath') '.m']);
metadata.loaderMD5=file_md5(which('load_siemens_lowres'));
metadata.trajectoryFunctionMD5=file_md5(which('mlTrajFromPulseq'));
metadata.griddingMEX=which('bmGridder_n2u_leight3_mex');
metadata.griddingMEXMD5=file_md5(metadata.griddingMEX);
metadata.elapsedSeconds=toc(started);
after=dir(job.rawFile);
assert(after.bytes==metadata.rawStat.bytes && after.datenum==metadata.rawStat.datenum, ...
    'Raw source changed during reconstruction.');
assert(strcmp(seqHash,file_md5(job.seqFile)),'Sequence changed during reconstruction.');
clear y t ve sumSquares;
x_mm=(-grid(1)/2:grid(1)/2-1)*metadata.voxelSizeMm(1);
y_mm=(-grid(2)/2:grid(2)/2-1)*metadata.voxelSizeMm(2);
z_mm=(-grid(3)/2:grid(3)/2-1)*metadata.voxelSizeMm(3);
% Only the magnitude is staged for the existing verified NIfTI exporter.
staged=fullfile(job.outDir,'reconstruction.mat');
save(staged,'xrms','metadata','x_mm','y_mm','z_mm','-v7.3');
report=export_anisotropic_nifti(staged,job.outDir);
assert(report.normalizedArrayExactReadback && report.relativeRawIntensityRestorationError<1e-7);
delete(staged);
f=figure('Visible','off','Color','w','Position',[100 100 1150 420]);
reference=max(xrms(:));
subplot(1,3,1); imagesc(x_mm,y_mm,max(xrms,[],3)'/reference); axis image xy;
xlabel('Sequence x (mm)'); ylabel('Sequence y (mm)'); title('XY maximum projection');
subplot(1,3,2); imagesc(x_mm,z_mm,squeeze(max(xrms,[],2))'/reference); axis image xy;
xlabel('Sequence x (mm)'); ylabel('Sequence z (mm)'); title('XZ maximum projection');
subplot(1,3,3); imagesc(y_mm,z_mm,squeeze(max(xrms,[],1))'/reference); axis image xy;
xlabel('Sequence y (mm)'); ylabel('Sequence z (mm)'); title('YZ maximum projection');
colormap(gray); drawnow;
exportgraphics(f,fullfile(job.outDir,'xrms_preview.png'),'Resolution',150); close(f);
fprintf('%s VERIFIED 240-CUBED XRMS COMPLETE: %s\n',job.label,report.nifti);
end
