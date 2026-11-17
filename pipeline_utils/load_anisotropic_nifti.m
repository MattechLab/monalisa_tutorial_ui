function r = load_anisotropic_nifti(folder,restoreRawIntensity)
%LOAD_ANISOTROPIC_NIFTI Read compact xRMS and optionally restore its raw scale.
if nargin<2, restoreRawIntensity=true; end
report=jsondecode(fileread(fullfile(folder,'reconstruction.json')));
info=niftiinfo(fullfile(folder,'xrms_normalized.nii.gz'));
r.xrms=niftiread(info);
assert(isequal(size(r.xrms),double(report.matrixSize(:).')));
assert(info.raw.sform_code==2 && info.raw.qform_code==0);
if restoreRawIntensity,r.xrms=double(r.xrms)*report.normalizationDivisor;end
r.metadata=report.metadata;r.export=report;
coords=cell(1,3);
for k=1:3
 coords{k}=report.originMm(k)+(0:report.matrixSize(k)-1)*report.spacingMm(k);
end
r.x_mm=coords{1};r.y_mm=coords{2};r.z_mm=coords{3};
end
