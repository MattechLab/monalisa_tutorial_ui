function isMatch = check_hash(rawDataFile,seqFile,headerSignature)
%CHECK_HASH Verify Pulseq file content against the Siemens header signature.
% The optional signature avoids rereading an already parsed raw-data header.
% Also recognizes nSeg/nShot lines added to [DEFINITIONS] after signing, but
% only if removing those lines restores the exact signed payload MD5.
if nargin < 3
    twix = mapVBVD_JH_for_monalisa(rawDataFile);
    if ~iscell(twix), twix = {twix}; end
    headerSignature = '';
    for k = numel(twix):-1:1
        if isfield(twix{k}.hdr,'Dicom') && isfield(twix{k}.hdr.Dicom,'tSequenceVariant')
            headerSignature = strtrim(char(twix{k}.hdr.Dicom.tSequenceVariant));
            if ~isempty(headerSignature), break; end
        end
    end
end
assert(~isempty(headerSignature),'Raw file has no Pulseq signature: %s',rawDataFile);
assert(isfile(seqFile),'Sequence file not found (check symbolic links): %s',seqFile);
fid = fopen(seqFile,'rb'); assert(fid>=0,'Cannot read %s',seqFile);
bytes = fread(fid,Inf,'*uint8').'; fclose(fid);
text = char(bytes);
pos = regexp(text,'(?m)^\[SIGNATURE\]','once');
if ~isempty(pos)
    stored = regexp(text(pos:end),'(?m)^Hash ([0-9a-fA-F]{32})','tokens','once');
    last = pos-1;
    if last>0 && bytes(last)==10, last=last-1; end
    if last>0 && bytes(last)==13, last=last-1; end
    bytes = bytes(1:last);
else
    stored = {};
end
temporary = tempname;
cleanup = onCleanup(@() delete(temporary));
actual = payload_md5(bytes,temporary);
verified = actual;
fprintf('Pulseq signature in TWIX:       %s\n',strtrim(char(headerSignature)));
if ~isempty(stored)
    fprintf('Pulseq [SIGNATURE] Hash:        %s (matches TWIX: %d)\n', ...
        stored{1},strcmpi(strtrim(char(headerSignature)),stored{1}));
end
fprintf('Recomputed sequence-content MD5: %s\n',actual);
if ~isempty(stored)
    if ~strcmpi(actual,stored{1})
        [candidate,removedNames] = without_grouping_definitions(bytes);
        if ~isempty(removedNames)
            candidateHash = payload_md5(candidate,temporary);
            fprintf('MD5 excluding added %s definitions: %s\n',strjoin(removedNames,', '),candidateHash);
            if strcmpi(candidateHash,stored{1})
                verified = candidateHash;
                fprintf('Verified exact signed payload; only added grouping definitions differ.\n');
            end
        end
        if ~strcmpi(verified,stored{1})
            error('Tutorial:SequenceContentMismatch', ...
                ['Sequence content differs from its stored signature: %s\n', ...
                 'TWIX: %s\nStored [SIGNATURE] Hash: %s\nCurrent content MD5: %s\n', ...
                 'The stored Hash was read successfully; this is a file-content verification failure.'], ...
                seqFile,strtrim(char(headerSignature)),stored{1},actual);
        end
    end
end
isMatch = strcmpi(strtrim(char(headerSignature)),verified);
if ~isMatch, warning('Tutorial:SequenceMismatch','Sequence and raw data do not match.'); end
end

function digest = payload_md5(bytes,path)
fid = fopen(path,'wb'); assert(fid>=0,'Cannot create a temporary hash input.');
fwrite(fid,bytes,'uint8'); fclose(fid);
digest = file_md5(path);
end

function [candidate,names] = without_grouping_definitions(bytes)
% Restrict the reconstruction to integer-valued grouping metadata definitions.
text = char(bytes); candidate = bytes; names = {};
first = regexp(text,'(?m)^\[DEFINITIONS\][ \t]*\r?\n','end','once');
if isempty(first), return; end
next = regexp(text(first+1:end),'(?m)^\[','once');
last = numel(text);
if ~isempty(next), last = first+next-1; end
definitions = text(first+1:last);
pattern = '(?m)^(nSeg|nShot)[ \t]+[1-9][0-9]*[ \t]*\r?\n';
found = regexp(definitions,pattern,'tokens');
if isempty(found), return; end
names = cellfun(@(token) token{1},found,'UniformOutput',false);
% Duplicated keys are ambiguous and must not be normalized away.
if numel(unique(names))~=numel(names), names = {}; return; end
candidate = uint8([text(1:first),regexprep(definitions,pattern,''),text(last+1:end)]);
end
