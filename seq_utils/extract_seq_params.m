% =====================================================
% Helper function to extract sequence definitions
% =====================================================
function params = extract_seq_params(seqFile)
    params = struct();
    if ~isfile(seqFile), return; end
    try
        seq = mr.Sequence();
        seq.read(seqFile);
        defs = seq.definitions;
        keysList = keys(defs);
        for i = 1:numel(keysList)
            key = keysList{i};
            val = defs(key);
            cleanKey = regexprep(lower(key), '[^a-z0-9_]', '');
            params.(cleanKey) = val;
        end
    catch ME
        warning('Failed to read seq params from %s: %s', seqFile, ME.message);
    end
end
% Sequence parameters
% seqFile = seqFolder + "/" + seqName_list{1};  % main sequence
% seqParams = extract_seq_params(seqFile);