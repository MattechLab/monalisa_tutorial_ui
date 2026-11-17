function answer = inputdlg(varargin)
% Test fixture: scripted responses; unexpected prompts fail the test.
key = 'monalisa_sequence_loading_test';
context = getappdata(0,key);
prompts = varargin{1}; titleText = varargin{2}; defaults = varargin{4};
context.dialogs = context.dialogs+1;
if strcmp(titleText,'Confirm sequence shot/segment counts')
    assert(context.allowConfirmation,'Unexpected confirmation dialog.');
    assert(isequal(defaults,{num2str(context.nSeg),num2str(context.nShot)}), ...
        'Confirmation counts were not prefilled from the sequence.');
    context.confirmations = context.confirmations+1;
    answer = defaults;
else
    assert(context.allowMissing,'Unexpected missing-field dialog.');
    answer = cell(size(prompts));
    for k = 1:numel(prompts)
        if startsWith(prompts{k},'nSeg:')
            answer{k} = num2str(context.nSeg);
        elseif startsWith(prompts{k},'nShot:')
            answer{k} = num2str(context.nShot);
        else
            error('Test:UnexpectedPrompt','Unexpected field prompt: %s',prompts{k});
        end
    end
end
setappdata(0,key,context);
end
