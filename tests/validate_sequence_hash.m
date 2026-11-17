function result = validate_sequence_hash(config)
%VALIDATE_SEQUENCE_HASH Check interoperability with the actual Pulseq writer.
% No scanner data, MEX routines, or reconstruction are needed.
if nargin<1, config = tutorial_config(); end
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(config.pulseqRoot,'matlab'), ...
    fullfile(repoRoot,'pipeline_utils'),fullfile(repoRoot,'seq_utils'));
workDir = tempname; mkdir(workDir);
fileCleanup = onCleanup(@() rmdir(workDir,'s'));
seqFile = fullfile(workDir,'signed.seq');
sequence = mr.Sequence();
sequence.addBlock(mr.makeDelay(1e-3));
sequence.write(seqFile,true);
storedHash = sequence.signatureValue;
assert(check_hash('',seqFile,storedHash),'Pulseq-generated signature was rejected.');
readBack = mr.Sequence(); readBack.read(seqFile);
assert(strcmpi(readBack.signatureValue,storedHash),'Pulseq reader disagrees with the stored hash.');
fid = fopen(seqFile,'rb'); bytes = fread(fid,Inf,'*uint8'); fclose(fid);
withGrouping = strrep(char(bytes.'),sprintf('[DEFINITIONS]\n'), ...
    sprintf('[DEFINITIONS]\nnSeg 44\nnShot 1055\n'));
assert(~isequal(uint8(withGrouping(:)),bytes),'Fixture must contain a definitions section.');
fid = fopen(seqFile,'wb'); fwrite(fid,uint8(withGrouping),'uint8'); fclose(fid);
assert(check_hash('',seqFile,storedHash),'Exact signed payload with added grouping metadata was rejected.');
bytes = uint8(withGrouping(:));
bytes(1) = uint8('!'); % Change content without changing its stored signature.
fid = fopen(seqFile,'wb'); fwrite(fid,bytes,'uint8'); fclose(fid);
rejected = false;
try
    check_hash('',seqFile,storedHash);
catch failure
    rejected = strcmp(failure.identifier,'Tutorial:SequenceContentMismatch');
end
assert(rejected,'Changed content with the old stored hash was accepted.');
result = struct('status','passed','pulseqWriterSignature',storedHash, ...
    'pulseqReaderAgrees',true,'addedGroupingVerified',true,'changedContentRejected',true);
disp(result);
end
