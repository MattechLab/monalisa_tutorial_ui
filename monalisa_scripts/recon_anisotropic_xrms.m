function [xrms,x0,metadata] = recon_anisotropic_xrms(rawFile,seqFile,outDir,nShotOff,kernelSigma)
%RECON_ANISOTROPIC_XRMS Reconstruct physical scanner coils on the declared FOV.
% Run setup_tutorial first. Pulseq coordinates are converted from cycles/m to
% cycles/mm once. The FOV is in mm and voxel size comes from sequence res_mm.
if nargin<4, nShotOff=20; end
if nargin<5, kernelSigma=0.61; end
validateattributes(kernelSigma,{'numeric'},{'scalar','real','finite','positive'});
kernelParameters=[kernelSigma 10];
assert(~isfile(fullfile(outDir,'reconstruction.mat')), ...
    'Reconstruction exists; choose a new output directory.');
if ~isfolder(outDir), mkdir(outDir); end
diary(fullfile(outDir,'reconstruction.log'));
logCleanup=onCleanup(@()diary('off'));
started=tic;
[image,acq,info,geometry,rawLines] = siemens_sequence_metadata( ...
    rawFile,seqFile,nShotOff,[],false,true);
assert(all(abs(geometry.voxelSizeMm-geometry.voxelSizeMm(1))<1e-9), ...
    'This spherical radial reconstruction expects isotropic nominal voxels.');
assert(isfield(info.sequenceParams,'requiredextensions') && ...
    contains(info.sequenceParams.requiredextensions,'ROTATIONS'), ...
    'Expected the anisotropic sequence with physical ROTATIONS.');
grid=geometry.matrixSize;
fovMm=geometry.FoV;
dK=1./fovMm;
seqHash=file_md5(seqFile);
seq=mr.Sequence(); seq.read(seqFile);
[kAll,tAll]=seq.calculateKspacePP();
assert(size(kAll,1)==3 && size(kAll,2)==acq.N*geometry.nSeg*geometry.nShot, ...
    'Pulseq ADC count differs from raw-data count.');
% Select exactly the same readouts for trajectory and complex scanner data.
kLines=reshape(kAll,3,acq.N,[]);
t=double(kLines(:,:,rawLines))/1000;
adcTimes=reshape(tAll,acq.N,[]); adcTimes=adcTimes(:,rawLines);
assert(all(isfinite(t),'all') && all(diff(adcTimes(:))>0));
directions=squeeze(t(:,end,:)-t(:,1,:));
directions=directions./vecnorm(directions);
radius=squeeze(sum(t.*reshape(directions,3,1,[]),1));
orthogonal=t-reshape(directions,3,1,[]).*reshape(radius,1,acq.N,[]);
assert(max(vecnorm(orthogonal),[],'all')<1e-5,'ADC lines are not radial about k=0.');
assert(max(vecnorm(t(:,acq.N/2+1,:)),[],'all')<1e-5, ...
    'Radial DCF requires a centre sample at N/2+1.');
ve=bmVolumeElement(t,'voronoi_full_radial3');
assert(numel(ve)==acq.N*numel(rawLines) && all(isfinite(ve)) && all(ve>0));
raw=read_siemens_blocks(image,acq,rawLines);
assert(isequal(size(raw),[acq.nCh acq.N numel(rawLines)]));
x0=complex(zeros([grid acq.nCh],'single'));
sumSquares=zeros(grid,'single');
for channel=1:acq.nCh
    % Keep scanner complex phase and readout order as returned by TWIX.
    signal=reshape(raw(channel,:,:),1,[]);
    coilImage=bmMathilda(signal,t,ve,[],grid,grid,dK,'gauss',3,kernelParameters,'MATLAB',true);
    assert(isequal(size(coilImage),grid) && all(isfinite(coilImage),'all'));
    x0(:,:,:,channel)=coilImage;
    sumSquares=sumSquares+abs(coilImage).^2;
    fprintf('seq coils reconstructed %d/%d\n',channel,acq.nCh);
end
clear raw signal;
xrms=sqrt(sumSquares/acq.nCh);
assert(all(isfinite(xrms),'all') && max(xrms(:))>0);
metadata=info;
metadata.matrixSize=grid;
metadata.voxelSizeMm=geometry.voxelSizeMm;
metadata.fovMm=fovMm;
metadata.dK_cycles_per_mm=dK;
metadata.trajectoryUnits='cycles/mm; physical calculateKspacePP coordinates / 1000';
metadata.sequenceMD5=seqHash;
metadata.keptReadoutIndices=rawLines;
metadata.retainedShots=acq.nShot-nShotOff;
metadata.imagingReadouts=numel(rawLines);
metadata.samplesUsed=acq.N*numel(rawLines);
metadata.method='MONALISA gridding; sqrt(mean(abs(physicalCoilImage).^2))';
metadata.densityCompensation='voronoi_full_radial3, recalculated after SI/shot removal';
metadata.kernel='gauss'; metadata.kernelWidth=3; metadata.kernelParameters=kernelParameters;
metadata.signalConjugated=false; metadata.additionalADCPhaseCorrection=false;
metadata.gradientDelayCorrection=false;
metadata.matlabVersion=version;
metadata.runnerMD5=file_md5([mfilename('fullpath') '.m']);
metadata.griddingMEX=which('bmGridder_n2u_leight3_mex');
metadata.griddingMEXMD5=file_md5(metadata.griddingMEX);
metadata.elapsedSeconds=toc(started);
x_mm=(-grid(1)/2:grid(1)/2-1)*geometry.voxelSizeMm(1);
y_mm=(-grid(2)/2:grid(2)/2-1)*geometry.voxelSizeMm(2);
z_mm=(-grid(3)/2:grid(3)/2-1)*geometry.voxelSizeMm(3);
% Write a temporary MAT file and check before publishing the final result.
partial=fullfile(outDir,'reconstruction.partial.mat');
save(partial,'xrms','x0','metadata','t','ve','adcTimes','x_mm','y_mm','z_mm','-v7.3');
saved=load(partial,'xrms','metadata');
assert(isequal(saved.xrms,xrms) && isequal(saved.metadata.matrixSize,grid));
rmsCheck=sqrt(mean(abs(x0).^2,4));
relativeRmsError=norm(double(rmsCheck(:))-double(xrms(:)))/norm(double(xrms(:)));
assert(relativeRmsError<1e-6,'Incremental xRMS differs from direct physical-coil RMS.');
assert(strcmp(seqHash,file_md5(seqFile)),'Sequence changed during reconstruction.');
after=dir(rawFile);
assert(after.bytes==info.rawStat.bytes && after.datenum==info.rawStat.datenum, ...
    'Raw source changed during reconstruction.');
movefile(partial,fullfile(outDir,'reconstruction.mat'));
f=figure('Visible','off','Color','w','Position',[100 100 1150 420]);
reference=max(xrms(:));
subplot(1,3,1); imagesc(x_mm,y_mm,max(xrms,[],3)'/reference); axis image xy;
xlabel('Sequence x (mm)'); ylabel('Sequence y (mm)'); title('XY maximum projection');
subplot(1,3,2); imagesc(x_mm,z_mm,squeeze(max(xrms,[],2))'/reference); axis image xy;
xlabel('Sequence x (mm)'); ylabel('Sequence z (mm)'); title('XZ maximum projection');
subplot(1,3,3); imagesc(y_mm,z_mm,squeeze(max(xrms,[],1))'/reference); axis image xy;
xlabel('Sequence y (mm)'); ylabel('Sequence z (mm)'); title('YZ maximum projection');
colormap(gray); exportgraphics(f,fullfile(outDir,'xrms_preview.png'),'Resolution',150); close(f);
validation=struct('status','complete','all_checks_pass',true, ...
    'matrix',grid,'fov_mm',fovMm,'voxel_mm',geometry.voxelSizeMm, ...
    'samples_per_readout',acq.N,'physical_coils',acq.nCh,'initial_shots_excluded',nShotOff, ...
    'retained_shots',metadata.retainedShots,'imaging_readouts',numel(rawLines), ...
    'samples_per_coil',metadata.samplesUsed,'relative_rms_combination_error',relativeRmsError, ...
    'raw_bytes',after.bytes,'sequence_md5',seqHash,'source_files_unchanged',true, ...
    'elapsed_seconds',metadata.elapsedSeconds,'gaussian_sigma',kernelSigma, ...
    'kernel_window_width',3,'kernel_parameters',kernelParameters);
if isfield(metadata,'B0_T'), validation.B0_T=metadata.B0_T; end
fid=fopen(fullfile(outDir,'validation.json'),'w'); assert(fid>=0);
jsonCleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(validation,PrettyPrint=true));
fprintf('ANISOTROPIC SCANNER XRMS COMPLETE: %s\n',outDir);
end
