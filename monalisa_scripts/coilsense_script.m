function [C1,reconDir,prescan_seqParams,xsos_all,ncc_all,infoMain] = coilsense_script(baseFolder,options)
%COILSENSE_SCRIPT Select prescan/main files, review a mask, and estimate C.
% setup_tutorial must be run first (main_script_3TPrismaXA60 cell 1).
% Optional second argument is a struct with:
%   files: prescanSeq, bodyRaw, headRaw, mainSeq, mainRaw (full paths)
%   outputRoot: result root; maskFile: saved MAT containing logical mask
%   showFigures: true by default; false requires maskFile
%   prescanShotOff/mainShotOff: 14/50
% Saved masks allow an approved calibration to be repeated without new choices.
if nargin<2, options=struct(); end
defaults=struct('files',struct(),'outputRoot','recon_results', ...
    'maskFile','','showFigures',true,'prescanShotOff',14,'mainShotOff',50, ...
    'nIter',5,'orientation',struct());
fields=fieldnames(options);
assert(all(ismember(fields,fieldnames(defaults))),'Unknown coil-sensitivity option.');
for k=1:numel(fields), defaults.(fields{k})=options.(fields{k}); end
options=defaults;
validateattributes(options.nIter,{'numeric'},{'scalar','integer','positive'});
stageTimer=tic; timings=struct();
assert(options.showFigures || ~isempty(options.maskFile), ...
    'A saved maskFile is required when showFigures=false.');
if nargin<1 || isempty(baseFolder), baseFolder=pwd; end
options.outputRoot=resolve_output_root(baseFolder,options.outputRoot);
files=options.files;
keys={'prescanSeq','bodyRaw','headRaw','mainSeq','mainRaw'};
titles={'Select prescan sequence','Select body-coil prescan raw data', ...
    'Select head-coil prescan raw data','Select main sequence','Select main raw data'};
instructions = { ...
    'Select the Pulseq sequence (.seq) used for the BC and HC PRESCANS.', ...
    'Select the BODY-COIL (BC) prescan raw-data file (.dat), often named bc.dat.', ...
    'Select the HEAD-COIL (HC) prescan raw-data file (.dat), often named hc.dat.', ...
    'Select the MAIN scan Pulseq sequence (.seq) for the acquisition to reconstruct.', ...
    'Select the MAIN raw-data file (.dat) matching the main sequence just selected.'};
for k=1:numel(keys)
    if ~isfield(files,keys{k}) || isempty(files.(keys{k}))
        if endsWith(keys{k},'Seq'), filter='*.seq'; else, filter='*.dat'; end
        selection_reminder(sprintf('Coil calibration %d/5: %s',k,titles{k}),instructions{k});
        [name,folder]=uigetfile(filter,titles{k},baseFolder);
        assert(~isequal(name,0),'File selection cancelled.');
        files.(keys{k})=fullfile(folder,name);
    end
    assert(isfile(files.(keys{k})),'File not found (check symbolic links): %s',files.(keys{k}));
end
[~,stem]=fileparts(files.headRaw);
token=regexp(stem,'MID\d+','match','once');
if ~isempty(token), stem=token; end
runName=[stem '_recon_C_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))];
reconDir=fullfile(options.outputRoot,runName);
assert(~isfolder(reconDir),'Output folder already exists.');
mkdir(reconDir);
diary(fullfile(reconDir,'calibration.log'));
logCleanup=onCleanup(@() diary('off'));
fprintf('Coil sensitivity output: %s\n',reconDir);
N_u=[48 48 48]; previewGrid=[64 64 64]; xsos_all=cell(3,1);

fprintf('Preparing body-coil prescan...\n');
[y,t,ve,acqBC,infoBC]=load_siemens_lowres(files.bodyRaw,files.prescanSeq,previewGrid,options.prescanShotOff);
prescan_seqParams=infoBC.sequenceParams;
xsos_all{1}=prescan_sos(y,t,ve,acqBC,previewGrid);
dK_u=[1 1 1]./acqBC.FoV;
[y_body,t_body,ve_body]=bmLowRes(y.',t,ve,N_u,dK_u); y_body=y_body.';
clear y t ve
timings.bodyPreviewSeconds=toc(stageTimer); stageTimer=tic;

fprintf('Preparing head-coil prescan...\n');
[y,t,ve,acqHC,infoHC]=load_siemens_lowres(files.headRaw,files.prescanSeq,previewGrid,options.prescanShotOff,prescan_seqParams);
assert(acqBC.FoV==acqHC.FoV,'BC/HC fields of view differ.');
xsos_all{2}=prescan_sos(y,t,ve,acqHC,previewGrid);
[y_surface,t_surface]=bmLowRes(y.',t,ve,N_u,dK_u); y_surface=y_surface.';
assert(isequal(t_body,t_surface),'BC/HC trajectories differ.');
clear y t ve t_surface
timings.headPreviewSeconds=toc(stageTimer); stageTimer=tic;

fprintf('Preparing main-scan orientation preview...\n');
[y,t,ve,acqMain,infoMain]=load_siemens_lowres(files.mainRaw,files.mainSeq,previewGrid,options.mainShotOff);
assert(acqHC.nCh==acqMain.nCh && acqHC.FoV==acqMain.FoV,'HC/main channel count or FoV differs.');
xsos_all{3}=prescan_sos(y,t,ve,acqMain,previewGrid);
clear y t ve
timings.mainPreviewSeconds=toc(stageTimer); stageTimer=tic;
[orientation_ok,ncc_all,orientationDetails]=checkPrescanOrientation(xsos_all,options.showFigures,options.orientation);
metadataComparison=struct('headCoilSelectionEqual',[], 'headMainSliceGeometryEqual',[]);
if ~isempty(infoHC.coilSelection) && ~isempty(infoMain.coilSelection)
    metadataComparison.headCoilSelectionEqual=isequal(infoHC.coilSelection,infoMain.coilSelection);
    fprintf('HC/main coil-selection metadata equal: %d (review any differences).\n',metadataComparison.headCoilSelectionEqual);
else
    fprintf('Coil-selection metadata unavailable: channel ordering is not verified.\n');
end
if ~isempty(infoHC.sliceGeometry) && ~isempty(infoMain.sliceGeometry)
    metadataComparison.headMainSliceGeometryEqual=isequal(infoHC.sliceGeometry,infoMain.sliceGeometry);
    fprintf('HC/main slice-geometry metadata equal: %d (review previews).\n',metadataComparison.headMainSliceGeometryEqual);
end
save(fullfile(reconDir,'orientation.mat'),'xsos_all','ncc_all','orientation_ok','orientationDetails','metadataComparison');
assert(orientation_ok,'Orientation check failed. No sensitivity maps were calculated.');
if options.showFigures
    bmImage(xsos_all);
    answer=questdlg('Accept BC/HC/main alignment after reviewing the matched previews? This is not a head-motion registration.', ...
        'Review calibration alignment','Accept','Cancel','Cancel');
    assert(strcmp(answer,'Accept'),'Calibration alignment not approved.');
end

fprintf('Constructing sensitivity gridding operators...\n');
[Gn,Gu,Gut]=bmTraj2SparseMat(t_body,ve_body,N_u,dK_u);
if isempty(options.maskFile)
    mask=bmCoilSense_nonCart_mask_automatic(y_body,Gn,false);
else
    saved=load(options.maskFile,'mask');
    assert(isfield(saved,'mask'),'maskFile must contain variable mask.');
    mask=saved.mask;
end
assert(islogical(mask) && isequal(size(mask),N_u) && any(mask(:)) && ~all(mask(:)), ...
    'The mask must be a nonempty, non-full logical 48^3 array.');
if options.showFigures
    bmImage(single(mask));
    answer=questdlg('Use this final calibration mask?', ...
        'Review calibration mask','Use mask','Cancel','Cancel');
    assert(strcmp(answer,'Use mask'),'Calibration mask not approved.');
end
save(fullfile(reconDir,'mask.mat'),'mask');
timings.reviewAndGriddingSeconds=toc(stageTimer); stageTimer=tic;

fprintf('Estimating body-coil reference...\n');
[y_ref,C_ref]=bmCoilSense_nonCart_ref(y_body,Gn,mask,[]);
timings.referenceSeconds=toc(stageTimer); stageTimer=tic;
fprintf('Estimating %d head-coil sensitivity maps...\n',acqHC.nCh);
C_array_prime=bmCoilSense_nonCart_primary(y_surface,y_ref,C_ref,Gn,ve_body,mask);
timings.primarySeconds=toc(stageTimer); stageTimer=tic;
fprintf('Refining sensitivity maps (%d iterations)...\n',options.nIter);
[C1,~]=bmCoilSense_nonCart_secondary(y_surface,C_array_prime,y_ref,C_ref,Gn,Gu,Gut,ve_body,options.nIter,false);
timings.refinementSeconds=toc(stageTimer);
disp(timings);
assert(isequal(size(C1),[N_u acqHC.nCh]) && all(isfinite(C1(:))) && any(abs(C1(:))>0), ...
    'Invalid sensitivity output.');
settings=struct('gridSize',N_u,'nIter',options.nIter,'prescanShotOff',options.prescanShotOff, ...
    'mainShotOff',options.mainShotOff,'maskFile',options.maskFile);
save(fullfile(reconDir,'calibration.mat'),'files','settings','infoBC','infoHC','infoMain','prescan_seqParams','timings');
header=prescan_seqParams.header; provenance=prescan_seqParams.headerprovenance;
save(fullfile(reconDir,'prescan_sequence_header.mat'),'header','provenance');
header=infoMain.sequenceParams.header; provenance=infoMain.sequenceParams.headerprovenance;
save(fullfile(reconDir,'main_sequence_header.mat'),'header','provenance');
if options.showFigures, bmImage(C1); end
fprintf('C estimation complete: %s. Run the save cell next.\n',mat2str(size(C1)));
end
