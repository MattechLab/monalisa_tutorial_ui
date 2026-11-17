function result = validate_sequence_headers()
%VALIDATE_SEQUENCE_HEADERS Normalize both styles; test required-field dialogs.
% All dialogs use scripted answers on temporary fixtures, not real choices.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'seq_utils'));
workDir = tempname; mkdir(workDir);
fileCleanup = onCleanup(@() rmdir(workDir,'s'));
required = {'FOV','ReadoutOversamplingFactor','nSam','nSeg','nShot','flagSelfNav'};
file = fullfile(workDir,'colleague.seq');
write_fixture(file,sprintf(['FOV 0.1 0.1 0.1\nReadoutOversamplingFactor 2\n', ...
    'TE_s 0.003755\nTR_s 0.00848\nnSeg 44\nnShot 1055\nmodality t1\n']));
[header,~] = read_sequence_header(file);
assert(isfield(header,'flagSelfNav') && isempty(header.flagSelfNav));
assert(isempty(header.rfmode) && strcmp(header.modality,'t1'));
calls = 0;
[params,completed,provenance] = extract_seq_params(file,required,struct('nSam',500),@reply);
assert(calls==2,'Invalid input should prompt again.');
assert(completed.nSam==500 && completed.flagSelfNav==1);
assert(isequal(sort(provenance.suppliedFields),sort({'nSam','flagSelfNav'})));
geometry = sequence_recon_geometry(params,500,44*1055);
assert(geometry.FoV==200 && geometry.excludeSI);
assert_rejected(@() sequence_recon_geometry(params,499,44*1055),'');
assert_rejected(@() sequence_recon_geometry(params,500,44*1055-1),'');
assert_rejected(@() extract_seq_params(file,required,struct(),@(varargin) []), ...
    'Tutorial:SequenceInputCancelled');

reference = fullfile(workDir,'reference.seq');
write_fixture(reference,sprintf(['FOV 0.12 0.12 0.12\nReadoutOversamplingFactor 1\n', ...
    'nSam 240\nnSeg 44\nnShot 267\nflagSelfNav 1\nRequiredExtensions ROTATIONS\n']));
[complete,styled,source] = extract_seq_params(reference,required,struct(),@unexpected_dialog);
assert(isempty(source.suppliedFields) && strcmp(styled.RequiredExtensions,'ROTATIONS'));
referenceGeometry = sequence_recon_geometry(complete,240,44*267);
assert(referenceGeometry.FoV==120);
% No-SI acquisitions preserve the first readout of each shot.
complete.flagselfnav = 0;
noSI = sequence_recon_geometry(complete,240,44*267);
assert(~noSI.excludeSI);
% Explicit physical-coordinate geometry supports anisotropic image grids.
aniso = complete; aniso.fov = [.24 .24 .18]; aniso.res_mm = [2 2 2];
aniso.readoutoversamplingfactor = 2; % Physical FOV stays as declared.
aniso.flagselfnav = 1;
g = sequence_recon_geometry(aniso,240,44*267,true);
assert(isequal(g.FoV,[240 240 180]) && isequal(g.matrixSize,[120 120 90]));
assert(isequal(g.voxelSizeMm,[2 2 2]) && g.excludeSI);
assert_rejected(@() sequence_recon_geometry(aniso,240,44*267),'');
badResolution = aniso; badResolution.res_mm=[2 2 2.1];
assert_rejected(@() sequence_recon_geometry(badResolution,240,44*267,true),'');
% Reading/completing headers must not alter the source definitions.
[again,~] = read_sequence_header(file);
assert(isequaln(header,again));
result = struct('status','passed','missingFieldsCompleted',true, ...
    'invalidInputRetried',true,'cancellationChecked',true, ...
    'oversampling2FovMm',geometry.FoV,'oversampling1FovMm',referenceGeometry.FoV);
disp(result);

    function answers = reply(prompts,~,~,defaults)
        calls = calls+1;
        answers = cell(size(prompts));
        for k = 1:numel(prompts)
            if startsWith(prompts{k},'nSam:')
                assert(strcmp(defaults{k},'500'),'Raw sample count should prefill the dialog.');
                answers{k} = '500';
            elseif calls==1
                answers{k} = 'not-a-number';
            else
                answers{k} = '1';
            end
        end
    end
end

function write_fixture(file,definitions)
fid = fopen(file,'wt');
fprintf(fid,'[VERSION]\nmajor 1\nminor 5\nrevision 1\n\n[DEFINITIONS]\n%s\n[BLOCKS]\n',definitions);
fclose(fid);
end

function answer = unexpected_dialog(varargin) %#ok<STOUT>
error('Test:UnexpectedDialog','A complete header must not prompt.');
end

function assert_rejected(action,id)
rejected = false;
try
    action();
catch failure
    rejected = isempty(id) || strcmp(failure.identifier,id);
end
assert(rejected,'Expected invalid input or cancellation to be rejected.');
end
