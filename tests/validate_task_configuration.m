function result = validate_task_configuration()
%VALIDATE_TASK_CONFIGURATION No raw samples, figures, or MEX calls needed.
repoRoot=fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot,fullfile(repoRoot,'pipeline_utils'));
work=tempname; mkdir(work); cleanup=onCleanup(@() rmdir(work,'s'));
task=tutorial_task_config_3T();
keys=fieldnames(task.files);
for k=1:numel(keys)
    name=[keys{k} '.fixture']; fid=fopen(fullfile(work,name),'w');
    fprintf(fid,'%s',keys{k}); fclose(fid); task.files.(keys{k})=name;
end
resolved=resolve_task_inputs(task,work);
assert(strcmp(resolved.files.mainRaw,fullfile(work,task.files.mainRaw)));
assert(resolved.mainShotOff==50 && resolved.normalizationMatrix==120);
bad=task; bad.name='../other'; rejected=false;
try, resolve_task_inputs(bad,work); catch, rejected=true; end
assert(rejected,'Task label allowed an output path escape.');
bad=task; bad.preparedCFile=task.files.mainRaw; rejected=false;
try, resolve_task_inputs(bad,work); catch, rejected=true; end
assert(rejected,'Unpaired prepared paths accepted.');
main=task; main.calibrationMethod='main';
main.files.bodyRaw='does_not_exist'; main.files.headRaw='does_not_exist';
resolve_task_inputs(main,work); % No prescan requirement or picker.
calDir=fullfile(work,'MID00001_recon_C_20261101_120000');
mkdir(fullfile(calDir,'C')); C=single(1); save(fullfile(calDir,'C','C.mat'),'C');
assert(strcmp(choose_calibration_file(work),fullfile(calDir,'C','C.mat')));
% Uneven anatomy gives an unambiguous identity, and a flip must fail.
rng(7); reference=rand(8,8,8,'single');
[ok,~,details]=checkPrescanOrientation({reference,reference,reference},false);
assert(ok && details(1).identityScore>0.999 && details(1).margin>0);
[ok,~,~]=checkPrescanOrientation({reference,flip(reference,1),reference},false);
assert(~ok,'A flipped HC preview was accepted.');
stat=dir(resolved.files.mainRaw);
rawInfo=struct('rawFile',resolved.files.mainRaw,'rawStat',stat, ...
    'seqFile',resolved.files.mainSeq,'nCh',2,'N',8,'nSeg',3,'nShot',100, ...
    'FoV',200,'headerSignature','fixture','sequenceParams',struct('flagselfnav',1));
identity=main_scan_identity(rawInfo);
rawInfo.sequenceParams.flagselfnav=0;
assert(~isequal(identity,main_scan_identity(rawInfo)),'SI exclusion missing from identity.');
result=struct('status','passed','explicitInputs',true,'invalidTasksRejected',true, ...
    'timestampedCalibrationFound',true,'orientationFlipRejected',true,'siIdentityChecked',true);
disp(result);
end
