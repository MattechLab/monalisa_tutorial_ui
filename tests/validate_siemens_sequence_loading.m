function result = validate_siemens_sequence_loading(config)
%VALIDATE_SIEMENS_SEQUENCE_LOADING Verify grouping before any raw sample read.
% Synthetic TWIX has NSeg=1 but physical grouping nSeg=3, nShot=4.
% Mitosius confirms with sequence defaults; missing fields use their input UI.
if nargin<1, config = tutorial_config(); end
repoRoot = fileparts(fileparts(mfilename('fullpath')));
config.baseFolder = repoRoot;
setup_tutorial(config);
workDir = tempname; mkdir(workDir);
fileCleanup = onCleanup(@() rmdir(workDir,'s'));
fixtureDir = fullfile(repoRoot,'tests','fixtures','siemens_sequence_loading');
key = 'monalisa_sequence_loading_test';
assert(~isappdata(0,key),'Another sequence-loading validation is active.');
addpath(fixtureDir,'-begin');
clear mapVBVD_JH_for_monalisa inputdlg createRawDataReader;
fixtureCleanup = onCleanup(@() restore_fixtures(fixtureDir,key));
nSample = 4; nCh = 2; nSeg = 3; nShot = 4;
values = single(reshape(1:nSample*nCh*nSeg*nShot,nSample,nCh,nSeg*nShot));
source = complex(values,-values/3);
sequence = mr.Sequence();
sequence.setDefinition('FOV',[0.1 0.1 0.1]);
sequence.setDefinition('ReadoutOversamplingFactor',2);
sequence.setDefinition('nSam',nSample);
sequence.setDefinition('nSeg',nSeg); sequence.setDefinition('nShot',nShot);
sequence.setDefinition('flagSelfNav',1);
for line = 1:nSeg*nShot
    sequence.addBlock(mr.makeAdc(nSample,'Dwell',1e-5));
end
seqFile = fullfile(workDir,'main.seq'); sequence.write(seqFile,true);
rawFile = fullfile(workDir,'meas_MID00001.dat');
fid = fopen(rawFile,'wb'); fclose(fid);
image = struct('NCol',nSample,'NCha',nCh,'NAcq',nSeg*nShot, ...
    'NSeg',1,'NLin',nSeg*nShot,'NEco',1,'IsReflected',false, ...
    'flagRemoveOS',false,'flagRampSampRegrid',false,'flagDoRawDataCorrect',false);
twix = struct('image',image,'hdr',struct('Dicom',struct('tSequenceVariant',sequence.signatureValue)));
context = struct('twix',twix,'source',source,'sampleReads',0,'dialogs',0, ...
    'confirmations',0,'allowMissing',false,'allowConfirmation',false,'nSeg',nSeg,'nShot',nShot);
setappdata(0,key,context);
[loadedImage,acq,info,~,lines] = siemens_sequence_metadata(rawFile,seqFile,1);
context = getappdata(0,key);
assert(context.dialogs==0 && context.sampleReads==0,'Complete sequence grouping must be resolved before reading samples.');
assert(acq.nSeg==3 && acq.nShot==4 && isequal(lines(:)',[5 6 8 9 11 12]));
context.allowConfirmation = true; setappdata(0,key,context);
[actual,~,~,~] = load_siemens_main(rawFile,seqFile,1,info.sequenceParams);
expected = permute(source(:,:,lines),[2 1 3]);
assert(isequal(actual,expected),'Blocked data differ from original reshape/SI/shot filtering.');
context = getappdata(0,key);
assert(context.dialogs==1 && context.confirmations==1,'Mitosius should show one prefilled confirmation.');
context.sampleReads = 0; setappdata(0,key,context);
blocked = read_siemens_blocks(loadedImage,acq,lines,2);
context = getappdata(0,key);
assert(isequal(blocked,expected) && context.sampleReads==3,'Multi-block reading changed raw order or complex values.');
bad = info.sequenceParams; bad.nshot = nShot+1;
context.sampleReads = 0; setappdata(0,key,context); rejected = false;
try
    siemens_sequence_metadata(rawFile,seqFile,1,bad);
catch failure
    rejected = contains(failure.message,'does not match acquired readouts');
end
context = getappdata(0,key);
assert(rejected && context.sampleReads==0,'Count mismatch must stop before reading samples.');
retryCalls = 0;
confirmed = confirm_sequence_grouping(info.sequenceParams,nSeg*nShot,@retry_counts);
assert(retryCalls==2 && confirmed.nseg==nSeg && confirmed.nshot==nShot);
cancelled = false;
try
    confirm_sequence_grouping(info.sequenceParams,nSeg*nShot,@(varargin) []);
catch failure
    cancelled = strcmp(failure.identifier,'Tutorial:SequenceInputCancelled');
end
assert(cancelled);
% The prompt is used only after counts are actually removed from the header.
remove(sequence.definitions,{'nSeg','nShot'});
missingFile = fullfile(workDir,'missing_counts.seq'); sequence.write(missingFile,true);
twix.hdr.Dicom.tSequenceVariant = sequence.signatureValue;
context.twix = twix; context.allowMissing = true; context.allowConfirmation = false;
setappdata(0,key,context);
[~,resolved,~,~,~] = siemens_sequence_metadata(rawFile,missingFile,0);
context = getappdata(0,key);
assert(context.dialogs==2 && resolved.nSeg==3 && resolved.nShot==4);
result = struct('status','passed','sequenceCountsPrefillConfirmation',true, ...
    'rawNSegCounterIgnored',true,'missingCountsPrompted',true, ...
    'countMismatchStopsBeforeRead',true,'blockedReadMatchesOriginal',true);
disp(result);
clear fixtureCleanup;

    function answer = retry_counts(~,~,~,defaults)
        retryCalls = retryCalls+1;
        if retryCalls==1
            assert(isequal(defaults,{'3','4'}));
            answer = {'0','4'};
        else
            answer = {'3','4'};
        end
    end
end

function restore_fixtures(folder,key)
rmpath(folder);
clear mapVBVD_JH_for_monalisa inputdlg createRawDataReader;
if isvalid(groot) && isappdata(0,key), rmappdata(0,key); end
end
