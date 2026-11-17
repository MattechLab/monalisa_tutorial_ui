function report = export_anisotropic_nifti(inputMat,outDir,encodingFile)
%EXPORT_ANISOTROPIC_NIFTI Verified compact magnitude export; source stays intact.
% max normalization is reversible using the JSON divisor (float32 precision).
% Coordinates are the nominal sequence grid, not verified patient anatomy.
repo=fileparts(which('tutorial_config_7TTerraVE12U'));
assert(~isempty(repo),'Add the MONALISA tutorial repository to the MATLAB path.'); addpath(repo);
addpath(fullfile(repo,'pipeline_utils'));
assert(isfile(inputMat),'Reconstruction MAT is missing.');
if nargin<2 || isempty(outDir), outDir=fileparts(inputMat); end
if nargin<3, encodingFile=''; end
if ~isfolder(outDir), mkdir(outDir); end
niiPath=fullfile(outDir,'xrms_normalized.nii.gz');
jsonPath=fullfile(outDir,'reconstruction.json');
assert(~isfile(niiPath) && ~isfile(jsonPath),'Compact export already exists.');
before=dir(inputMat);
r=load(inputMat,'xrms','metadata','x_mm','y_mm','z_mm');
assert(isa(r.xrms,'single') && all(isfinite(r.xrms),'all') && all(r.xrms>=0,'all'));
normalizationDivisor=max(double(r.xrms(:)));assert(normalizationDivisor>0);
normalized=single(double(r.xrms)/normalizationDivisor);
xyz={double(r.x_mm(:)),double(r.y_mm(:)),double(r.z_mm(:))};
spacing=cellfun(@(q)q(2)-q(1),xyz);
for k=1:3
    assert(all(abs(diff(xyz{k})-spacing(k))<1e-10) && spacing(k)>0);
end
assert(isequal(size(normalized),r.metadata.matrixSize));
assert(all(abs(spacing-r.metadata.voxelSizeMm)<1e-10));
% Native MATLAB writer supplies NIfTI layout; only its transform code is fixed.
stem=fullfile(outDir,'xrms_normalized'); uncompressed=[stem '.nii'];
assert(~isfile(uncompressed),'Uncompressed temporary output already exists.');
niftiwrite(normalized,stem);
info=niftiinfo(uncompressed);info.PixelDimensions=spacing;info.SpaceUnits='Millimeter';
info.Description='Max-normalized xRMS; nominal sequence grid; anatomy unvalidated';
origin=cellfun(@(q)q(1),xyz);
info.Transform=affine3d([diag(spacing) zeros(3,1);origin 1]);
info.TransformName='Sform';niftiwrite(normalized,stem,info);
% MATLAB writes sform_code=1. Code 2 explicitly denotes an aligned nominal grid.
% NIfTI-1 little-endian sform_code is the int16 at zero-based byte offset 254.
fid=fopen(uncompressed,'r+','ieee-le');assert(fid>=0);
assert(fseek(fid,254,'bof')==0 && fwrite(fid,int16(2),'int16')==1);fclose(fid);
gzip(uncompressed);
savedInfo=niftiinfo(niiPath);saved=niftiread(savedInfo);
assert(isequal(saved,normalized),'NIfTI saved values differ from normalized xRMS.');
assert(isequal(savedInfo.Transform.T,info.Transform.T));
assert(strcmp(savedInfo.SpaceUnits,'Millimeter') && savedInfo.raw.sform_code==2 && savedInfo.raw.qform_code==0);
assert(all(abs(savedInfo.PixelDimensions-spacing)<1e-10));
relativeError=norm(double(saved(:))*normalizationDivisor-double(r.xrms(:)))/norm(double(r.xrms(:)));
assert(relativeError<1e-7,'Normalization restoration exceeds float32 tolerance.');
delete(uncompressed);
if ~isempty(encodingFile)
    e=load(inputMat,'t','ve','adcTimes','metadata','x_mm','y_mm','z_mm');
    if isfile(encodingFile)
        old=load(encodingFile,'t','ve','adcTimes','x_mm','y_mm','z_mm');
        for key={'t','ve','adcTimes','x_mm','y_mm','z_mm'}
            assert(isequal(old.(key{1}),e.(key{1})),['Encoding differs: ' key{1}]);
        end
    else
        parent=fileparts(encodingFile);if ~isfolder(parent),mkdir(parent);end
        save(encodingFile,'-struct','e','-v7.3');
        old=load(encodingFile,'t','ve','adcTimes');
        assert(isequal(old.t,e.t) && isequal(old.ve,e.ve) && isequal(old.adcTimes,e.adcTimes));
    end
end
report=struct('format','compressed NIfTI-1 float32', ...
    'inputMat',inputMat,'inputMatMD5',file_md5(inputMat),'nifti',niiPath, ...
    'niftiMD5',file_md5(niiPath),'normalization','divide each magnitude volume by its own maximum; no clipping', ...
    'normalizationDivisor',normalizationDivisor,'restoreRawIntensity','double(niftiread(nifti))*normalizationDivisor', ...
    'normalizedRange',[double(min(saved(:))) double(max(saved(:)))], ...
    'matrixSize',size(saved),'spacingMm',spacing,'originMm',origin, ...
    'coordinateFrame','Nominal centred sequence x/y/z in mm; not validated patient anatomical orientation', ...
    'sformCode',2,'qformCode',0,'noResamplingOrAxisPermutation',true, ...
    'normalizedArrayExactReadback',true,'relativeRawIntensityRestorationError',relativeError, ...
    'sharedEncodingFile',encodingFile,'sourceMATRetainedByExporter',true,'metadata',r.metadata);
fid=fopen(jsonPath,'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));fclose(fid);
after=dir(inputMat);assert(after.bytes==before.bytes && after.datenum==before.datenum,'Input MAT changed.');
fprintf('VERIFIED NIFTI: %s; source MAT still intact\n',niiPath);
end
