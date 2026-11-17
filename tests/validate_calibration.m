function result = validate_calibration(config,files,maskFile,referenceFile)
%VALIDATE_CALIBRATION Run the production path on real files and a saved mask.
% referenceFile is a previous C.mat from the same acquisition and mask.
config=setup_tutorial(config);
options=struct('files',files,'outputRoot',fullfile(config.outputRoot,'validation'), ...
    'maskFile',maskFile,'showFigures',false);
[C,reconDir]=coilsense_script(config.baseFolder,options);
reference=load(referenceFile,'C');
assert(isequal(size(C),size(reference.C)) && all(isfinite(C(:))));
relativeError=norm(double(C(:))-double(reference.C(:)))/norm(double(reference.C(:)));
assert(relativeError<1e-5,'Calibration differs from the approved reference.');
save(fullfile(reconDir,'C.mat'),'C','-v7.3');
result=struct('status','passed','relativeError',relativeError,'size',size(C), ...
    'referenceFile',referenceFile,'output',reconDir,'matlab',version);
fid=fopen(fullfile(reconDir,'validation.json'),'w');
fwrite(fid,jsonencode(result,'PrettyPrint',true)); fclose(fid);
disp(result);
end
