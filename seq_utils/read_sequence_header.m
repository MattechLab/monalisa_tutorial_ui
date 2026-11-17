function [header,provenance] = read_sequence_header(seqFile)
%READ_SEQUENCE_HEADER Normalize Pulseq definitions to Yiwei's header style.
% Reads definitions only. Missing optional values remain []; no file is edited.
% FOV is in metres, TR_s/TE_s in seconds, res_mm in millimetres.
order = {'AdcRasterTime','BlockDurationRasterTime','FOV','FOVCenter', ...
    'GradientRasterTime','Implementation','RadiofrequencyRasterTime', ...
    'ReadoutOversamplingFactor','ReceiverGainHigh','RequiredExtensions', ...
    'TE_s','TR_s','TotalDuration','bw_px','fidNav','flagSelfNav','grad_mode', ...
    'modality','nSam','nSeg','nShot','overlapping_k_center','phase_sign', ...
    'res_mm','rf_flip_deg','rf_gap_us','rf_norm_block','rf_offset_hz', ...
    'rf_spoiling_increment_deg','rf_subpulse_us','rfdelay_ratio','rfmode', ...
    'seqSuffix','seq_mode','swapGrad','traj_design','useGradSpoil', ...
    'useRewinder','use_rfspoiler'};
header = cell2struct(repmat({[]},size(order)),order,2);
fid = fopen(seqFile,'rt'); assert(fid>=0,'Cannot open sequence: %s',seqFile);
fileCleanup = onCleanup(@() fclose(fid));
inDefinitions = false; original = struct(); seen = {};
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);
    if isempty(line) || startsWith(line,'#'), continue; end
    if startsWith(line,'[')
        if inDefinitions, break; end
        inDefinitions = strcmp(line,'[DEFINITIONS]');
        continue;
    end
    if ~inDefinitions, continue; end
    entry = regexp(line,'^(\S+)\s+(.+)$','tokens','once');
    assert(~isempty(entry),'Malformed sequence definition: %s',line);
    key = entry{1}; valueText = strtrim(entry{2});
    assert(isvarname(key),'Unsupported sequence definition name: %s',key);
    match = find(strcmpi(order,key),1);
    if isempty(match), name = key; else, name = order{match}; end
    assert(~any(strcmpi(seen,name)),'Duplicate sequence definition: %s',name);
    seen{end+1} = name; %#ok<AGROW>
    values = str2double(strsplit(valueText));
    if all(isfinite(values)), value = values; else, value = valueText; end
    original.(key) = value;
    header.(name) = value;
end
assert(inDefinitions,'Sequence has no [DEFINITIONS] section: %s',seqFile);
provenance = struct('sourceFile',seqFile,'originalDefinitions',original, ...
    'suppliedFields',{{}},'missingOptionalFields',{{}});
end
