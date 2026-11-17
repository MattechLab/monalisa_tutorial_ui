function [C,reconDir,xrms,x0,rawInfo] = coilsense_main_script(baseFolder,outputRoot,nShotOff,options)
%COILSENSE_MAIN_SCRIPT Direct main-scan relative C, matching the 7T method.
% No BC/HC reference; common object phase remains. This is not ESPIRiT.
% Returns physical-coil images for the existing ROVir cell.
gridSize = [48 48 48];
if nargin<4, options=struct(); end
defaults=struct('files',struct(),'maskFile','','supportThreshold',0.1, ...
    'smoothSigma',1,'laplaceIterations',1000,'laplaceTolerance',1e-4);
keys=fieldnames(options);
assert(all(ismember(keys,fieldnames(defaults))),'Unknown main-calibration option.');
for k=1:numel(keys), defaults.(keys{k})=options.(keys{k}); end
options=defaults;
timer=tic;
[xrms,coilCells,acq,rawInfo] = check_orient_xrms(baseFolder,gridSize(1),nShotOff,options.files);
x0 = cat(4,coilCells{:});
xrss = sqrt(sum(abs(x0).^2,4));
supportThreshold = options.supportThreshold;
if ~isempty(options.maskFile)
    saved=load(options.maskFile,'mask'); mask=saved.mask;
    assert(islogical(mask) && isequal(size(mask),gridSize) && any(mask(:)) && ~all(mask(:)), ...
        'Saved mask must be nonempty, non-full logical 48^3.');
    bmImage({xrms,single(mask)});
    choice=questdlg('Use this saved mask for the selected main scan?', ...
        'Review saved mask','Use mask','Cancel','Cancel');
    assert(strcmp(choice,'Use mask'),'Mask not approved.');
else
while true
    answer = inputdlg({'Signal-support threshold (fraction of maximum RSS, 0 < value < 1):'}, ...
        'Main-scan calibration mask',[1 80],{num2str(supportThreshold)});
    assert(~isempty(answer),'Mask selection cancelled.');
    supportThreshold = str2double(answer{1});
    validateattributes(supportThreshold,{'numeric'},{'scalar','real','finite','>',0,'<',1});
    mask = xrss>supportThreshold*max(xrss(:));
    bmImage({xrms,single(mask)});
    choice = questdlg('Use this signal-support mask for main-scan C estimation?', ...
        'Review calibration mask','Use mask','Adjust','Cancel','Cancel');
    if strcmp(choice,'Use mask'), break; end
    assert(strcmp(choice,'Adjust'),'Mask not approved.');
end
end
smoothSigma = options.smoothSigma; laplaceIterations = options.laplaceIterations;
laplaceTolerance = options.laplaceTolerance;
[C,calibrationDetails] = estimate_relative_sensitivities( ...
    x0,mask,smoothSigma,laplaceIterations,laplaceTolerance);
rawFile = rawInfo.rawFile; seqFile = rawInfo.seqFile;
systemInfo = read_twix_system_info(rawFile);
[~,stem] = fileparts(rawFile);
token = regexp(stem,'MID\d+','match','once');
if ~isempty(token), stem = token; end
runName = [stem '_recon_C_main_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))];
reconDir = fullfile(resolve_output_root(baseFolder,outputRoot),runName);
assert(~isfolder(reconDir),'Output folder already exists.');
mkdir(reconDir);
elapsedSeconds=toc(timer);
save(fullfile(reconDir,'mask.mat'),'mask','supportThreshold');
save(fullfile(reconDir,'calibration.mat'),'systemInfo','rawInfo','rawFile', ...
    'seqFile','gridSize','nShotOff','calibrationDetails','acq','options','elapsedSeconds');
header = rawInfo.sequenceParams.header;
provenance = rawInfo.sequenceParams.headerprovenance;
save(fullfile(reconDir,'main_sequence_header.mat'),'header','provenance');
bmImage(C);
fprintf('Main-scan relative C estimated. Run the save cell next.\n');
end
