function result = validate_main_preview_cell()
%VALIDATE_MAIN_PREVIEW_CELL Exercise the actual inspection cell without MRI data.
% Mock only reconstruction/display/identity, retaining the script's control flow.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
source = fileread(fullfile(repoRoot,'main_script_3TPrismaXA60.m'));
cells = regexp(source,'(?m)^%%','split');
index = find(cellfun(@(c) startsWith(c,' Optional: Main-data inspection'),cells));
assert(isscalar(index),'Could not locate the task preview cell.');
cellSource = ['%' cells{index}];
task = struct('previewMatrix',60,'mainShotOff',50,'files',struct('mainRaw','fixture'));
baseFolder = 'fixture';
bmImage = @(x) assert(~isempty(x));
main_scan_identity = @(info) info.identity;

% Reproduce the reported failure: setup exists, calibration variables do not.
check_orient_xrms = @(folder,matrix,shots,files,params) mock_preview( ...
    folder,matrix,shots,files,params,[]);
eval(cellSource);
assert(isequal(mainRms,1) && isequal(coilImages,{1}));
assert(exist('mainInfo','var')==0 && exist('mainIdentity','var')==0, ...
    'Preview must not create calibration provenance.');

% Calibration metadata must still be reused and its identity protected.
mainInfo = struct('sequenceParams',struct('fixture',true));
mainIdentity = 'same-scan';
check_orient_xrms = @(folder,matrix,shots,files,params) mock_preview( ...
    folder,matrix,shots,files,params,mainInfo.sequenceParams);
eval(cellSource);
mainIdentity = 'different-scan';
rejected = false;
try
    eval(cellSource);
catch exception
    rejected = strcmp(exception.message,'Main scan changed.');
end
assert(rejected,'Inspection accepted a changed acquisition after calibration.');
result = struct('status','passed','previewBeforeCalibration',true, ...
    'calibrationMetadataReused',true,'changedAcquisitionRejected',true);
disp(result);
end

function [rms,coils,acq,info] = mock_preview(folder,matrix,shots,files,params,expected)
assert(strcmp(folder,'fixture') && matrix==60 && shots==50);
assert(strcmp(files.mainRaw,'fixture'));
assert(isequal(params,expected),'Incorrect sequence metadata passed to preview.');
rms = 1; coils = {1}; acq = []; info = struct('identity','same-scan');
end
