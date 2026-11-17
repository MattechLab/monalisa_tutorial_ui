function task = resolve_task_inputs(task,baseFolder)
%RESOLVE_TASK_INPUTS Resolve explicit paths; select missing inputs once.
assert(~isempty(regexp(task.name,'^[A-Za-z0-9][A-Za-z0-9_-]*$','once')), ...
    'task.name must contain only letters, digits, underscores, or hyphens.');
assert(ismember(task.calibrationMethod,{'prescan','main','reuse'}),'Unknown calibration method.');
validateattributes(task.mainShotOff,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(task.prescanShotOff,{'numeric'},{'scalar','integer','nonnegative'});
assert(ismember(task.normalizationMatrix,[60 120]),'Choose a normalization matrix of 60 or 120.');
validateattributes(task.previewMatrix,{'numeric'},{'scalar','integer','positive','even'});
validateattributes(task.nVirtualCoils,{'numeric'},{'scalar','integer','positive'});
validateattributes(task.prescanIterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(task.reconstruction.Matrix_size,{'numeric'},{'scalar','integer','positive','even'});
validateattributes(task.reconstruction.nIter,{'numeric'},{'scalar','integer','positive'});
validateattributes(task.reconstruction.delta,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(task.reconstruction.rho,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(task.reconstruction.flagCS,{'numeric','logical'},{'scalar','binary'});
validateattributes(task.useROVir,{'numeric','logical'},{'scalar','binary'});
keys = {'mainRaw','mainSeq'};
if strcmp(task.calibrationMethod,'prescan')
    keys = [keys {'prescanSeq','bodyRaw','headRaw'}];
end
for k=1:numel(keys)
    key = keys{k};
    if isempty(task.files.(key))
        if endsWith(key,'Seq'), filter='*.seq'; else, filter='*.dat'; end
        selection_reminder(sprintf('Task %s: input %d/%d',task.name,k,numel(keys)), ...
            sprintf('Select %s (%s) for this task.\nDataset: %s',key,filter,baseFolder));
        [name,folder] = uigetfile(filter,['Select ' key],baseFolder);
        assert(~isequal(name,0),'Task input selection cancelled.');
        task.files.(key) = fullfile(folder,name);
    else
        task.files.(key) = resolve_output_root(baseFolder,task.files.(key));
    end
    assert(isfile(task.files.(key)),'Missing %s: %s',key,task.files.(key));
end
pathFields = {'maskFile','existingCFile','preparedMitosiusDir','preparedCFile'};
for k=1:numel(pathFields)
    key = pathFields{k};
    if ~isempty(task.(key))
        task.(key) = resolve_output_root(baseFolder,task.(key));
        assert(isfile(task.(key)) || isfolder(task.(key)),'Missing task path: %s',task.(key));
    end
end
if strcmp(task.calibrationMethod,'reuse')
    if isempty(task.existingCFile)
        selection_reminder('Task: existing coil sensitivities', ...
            'Select the physical-coil C.mat matching this main acquisition.');
        [name,folder] = uigetfile('*.mat','Select existing physical-coil C',baseFolder);
        assert(~isequal(name,0),'Calibration file selection cancelled.');
        task.existingCFile = fullfile(folder,name);
    end
    assert(isfile(task.existingCFile),'Missing existing physical-coil C file.');
end
assert(isempty(task.preparedMitosiusDir)==isempty(task.preparedCFile), ...
    'Supply both preparedMitosiusDir and its matching preparedCFile.');
end
