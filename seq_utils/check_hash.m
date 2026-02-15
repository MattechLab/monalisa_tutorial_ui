function isMatch = check_hash(rawDataFile, seqFile, mapVBVD_path, pulseq_path)
% Yiwei Jia, ChatGPT
%CHECK_HASH  Compare the MD5 “signature” embedded in a Siemens TWIX file
%            with the signature re-computed from its corresponding Pulseq
%            *.seq file.
%
%   isMatch = CHECK_HASH(rawDataFile, seqFile)
%   isMatch = CHECK_HASH(rawDataFile, seqFile, mapVBVD_path, pulseq_path)
%
%   INPUTS
%     rawDataFile   – full path to the Siemens raw data file (*.dat / *.twix)
%     seqFile       – full path to the Pulseq sequence file (*.seq)
%     mapVBVD_path  – (optional) folder containing the mapVBVD toolbox
%                     default = '/Users/cag/Documents/forclone/mapVBVD_Jaime'
%     pulseq_path   – (optional) folder containing the Pulseq repository
%                     default = '/Users/cag/Documents/forclone/pulseq_v15/pulseq'
%
%   OUTPUT
%     isMatch       – logical scalar, true if the 32-byte hash read from the
%                     TWIX header equals the hash re-calculated from seqFile
%
%   The function issues a warning when the two hashes differ.
%
%   Example
%     ok = check_hash('meas_MID00042_FID12345.dat','mySeq.seq');
%
%  
% -------------------------------------------------------------------------

% --------- handle optional paths ----------------------------------------
if nargin < 3 || isempty(mapVBVD_path)
    mapVBVD_path = '/Users/cag/Documents/forclone/mapVBVD_Jaime';
end

if nargin < 4 || isempty(pulseq_path)
    pulseq_path  = '/Users/cag/Documents/forclone/pulseq_v15/pulseq';
end

addpath(genpath(mapVBVD_path));   % mapVBVD + helpers
addpath(genpath(pulseq_path));    % for completeness (not used here)

% --------- read TWIX header & extract 32-char MD5 ------------------------
twixObj        = mapVBVD_JB(rawDataFile);


twix_idx = 1;
seqHash_twix   = char(twixObj{twix_idx}.hdr.Dicom.tSequenceVariant);  % already 32 chars
while isempty(seqHash_twix)
    warning('Sequence hash in TWIX %d header is empty – aborting.', ...
        twix_idx);
    twix_idx = twix_idx+1;
    seqHash_twix   = char(twixObj{twix_idx}.hdr.Dicom.tSequenceVariant);  % already 32 chars
end
fprintf('TWIX file contains Pulseq signature: %s\n', seqHash_twix);

% --------- load .seq file (strip [SIGNATURE] block) ----------------------
rawSeq = fileread(seqFile);
sigPos = strfind(rawSeq, '[SIGNATURE]');
if ~isempty(sigPos)
    rawSeq = rawSeq(1, 1:sigPos(1)-2);     % drop CR/LF before the tag
end

% --------- re-calculate MD5 ---------------------------------------------
md      = java.security.MessageDigest.getInstance('MD5');
md.update(uint8(rawSeq));
hashRaw = md.digest();                             % int8 Java array
hashStr = lower( reshape( dec2hex(typecast(hashRaw,'uint8'))', 1, [] ) );

fprintf('Re-calculated MD5 hash           : %s\n', hashStr);

% --------- compare -------------------------------------------------------
isMatch = isequal(seqHash_twix, hashStr);

if ~isMatch
    warning(['Hash mismatch – the Pulseq file does NOT correspond to ' ...
             'the supplied raw data.']);
else
    disp('Hash match!')
end

end