function [params,header,provenance] = extract_seq_params(seqFile,requiredFields,suggestions,inputFunction)
%EXTRACT_SEQ_PARAMS Normalize definitions and request missing reconstruction data.
% header preserves Yiwei-style names; params retains lowercase field names for
% existing callers. Only requiredFields are prompted. Suggestions prefill the
% input dialog but are never accepted without user confirmation.
% inputFunction defaults to inputdlg; injection permits unattended UI tests.
params = struct();
assert(isfile(seqFile),'Sequence file not found (check symbolic links): %s',seqFile);
if nargin<2, requiredFields = {}; end
if nargin<3, suggestions = struct(); end
if nargin<4, inputFunction = @inputdlg; end
[header,provenance] = read_sequence_header(seqFile);
requiredFields = cellstr(requiredFields);
missing = {};
for k = 1:numel(requiredFields)
    field = requiredFields{k};
    assert(isfield(header,field),'Unknown required header field: %s',field);
    if ~valid_value(field,header.(field)), missing{end+1} = field; end %#ok<AGROW>
end
while ~isempty(missing)
    fprintf('Sequence %s needs user input for: %s\n',seqFile,strjoin(missing,', '));
    prompts = cell(size(missing)); defaults = cell(size(missing));
    for k = 1:numel(missing)
        field = missing{k}; prompts{k} = field_prompt(field);
        value = header.(field);
        if isempty(value) && isfield(suggestions,field), value = suggestions.(field); end
        if isempty(value), defaults{k} = ''; elseif isnumeric(value), defaults{k} = num2str(value); else, defaults{k} = value; end
    end
    answer = inputFunction(prompts,'Missing/invalid sequence reconstruction fields',[1 85],defaults);
    if isempty(answer), error('Tutorial:SequenceInputCancelled','Sequence-parameter entry cancelled.'); end
    invalid = {};
    for k = 1:numel(missing)
        field = missing{k}; value = str2double(strsplit(strtrim(answer{k})));
        header.(field) = value;
        if valid_value(field,value)
            provenance.suppliedFields{end+1} = field;
        else
            invalid{end+1} = field; %#ok<AGROW>
        end
    end
    missing = invalid;
end
if ~isempty(header.FOV) && isnumeric(header.FOV)
    if isscalar(header.FOV), header.FOV = repmat(header.FOV,1,3); end
end
keysList = fieldnames(header);
for i = 1:numel(keysList)
    key = keysList{i};
    if isempty(header.(key)), provenance.missingOptionalFields{end+1} = key; continue; end
    cleanKey = regexprep(lower(key),'[^a-z0-9_]','');
    params.(cleanKey) = header.(key);
end
params.header = header;
params.headerprovenance = provenance;
end

function valid = valid_value(field,value)
valid = isnumeric(value) && isreal(value) && ~isempty(value) && all(isfinite(value(:)));
if ~valid, return; end
switch field
    case 'FOV'
        valid = (isscalar(value) || numel(value)==3) && all(value(:)>0);
    case {'nSeg','nShot','nSam'}
        valid = isscalar(value) && value>=1 && value==fix(value);
    case 'flagSelfNav'
        valid = isscalar(value) && any(value==[0 1]);
    case 'ReadoutOversamplingFactor'
        valid = isscalar(value) && value>=1;
    otherwise
        error('Tutorial:UnknownRequiredField','No reconstruction validator for %s.',field);
end
end

function prompt = field_prompt(field)
switch field
    case 'FOV', prompt = 'FOV [m]: one isotropic value or three values, e.g. 0.12 0.12 0.12';
    case 'nSeg', prompt = 'nSeg: physical readouts per shot, INCLUDING SI if present (positive integer)';
    case 'nShot', prompt = 'nShot: physical number of shots, BEFORE discarding initial shots (positive integer)';
    case 'nSam', prompt = 'nSam: acquired imaging samples per readout, INCLUDING readout oversampling';
    case 'flagSelfNav', prompt = 'flagSelfNav: 1 if the FIRST readout of every shot is SI; 0 if there is no SI readout';
    case 'ReadoutOversamplingFactor', prompt = 'ReadoutOversamplingFactor: acquired readout oversampling, e.g. 1 or 2';
    otherwise, error('Tutorial:UnknownRequiredField','No prompt for %s.',field);
end
end
