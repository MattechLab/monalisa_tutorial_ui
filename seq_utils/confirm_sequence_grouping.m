function params = confirm_sequence_grouping(params,nReadouts,inputFunction)
%CONFIRM_SEQUENCE_GROUPING Prefill physical counts from resolved sequence data.
% The raw count includes all shots/readouts before SI or startup exclusions.
if nargin<3, inputFunction = @inputdlg; end
original = [params.nseg params.nshot];
defaults = {num2str(params.nseg),num2str(params.nshot)};
titleText = 'Confirm sequence shot/segment counts';
while true
    prompts = {sprintf(['Raw acquired readouts: %d. Counts below come from the sequence.\n', ...
        'nSeg: readouts per shot, INCLUDING SI if present:'],nReadouts), ...
        'nShot: total shots, BEFORE discarding initial shots:'};
    answer = inputFunction(prompts,titleText,[1 95],defaults);
    if isempty(answer)
        error('Tutorial:SequenceInputCancelled','Shot/segment confirmation cancelled; no samples were read.');
    end
    counts = [str2double(answer{1}),str2double(answer{2})];
    if all(isfinite(counts) & counts>=1 & counts==fix(counts)) && prod(counts)==nReadouts
        break;
    end
    defaults = answer;
    titleText = sprintf('Invalid counts: nSeg times nShot must equal %d',nReadouts);
end
params.nseg = counts(1); params.nshot = counts(2);
params.header.nSeg = counts(1); params.header.nShot = counts(2);
params.headerprovenance.groupingConfirmation = struct( ...
    'sequenceDefaults',original,'confirmedCounts',counts,'rawReadouts',nReadouts);
names = {'nSeg','nShot'};
for k = 1:2
    if counts(k)~=original(k)
        params.headerprovenance.suppliedFields = unique( ...
            [params.headerprovenance.suppliedFields,names(k)],'stable');
    end
end
end
