function runtime = setup_monalisa_runtime(monalisaRoot,cacheRoot)
%SETUP_MONALISA_RUNTIME Use locally built MEX files; unify OpenMP on macOS.
% No MONALISA source, installed binary, or global MATLAB setting is changed.
names = {'bmSparseMat_cC_oBlock_omp_mex','bmImLaplaceEquationSolver3_omp_mex'};
relative = {fullfile('sparseMat','mex','bmSparseMat_vec','bmSparseMat_cC_oBlock_omp'), ...
    fullfile('image3','mex','bmImLaplaceEquationSolver3_omp')};
sources = cell(size(names));
for k = 1:numel(names)
    sources{k} = fullfile(monalisaRoot,'src',relative{k},[names{k} '.' mexext]);
    assert(isfile(sources{k}), ...
        'Build MONALISA for this MATLAB/platform first (compile_mex_for_monalisa), then restart MATLAB. Missing: %s',sources{k});
end
if ~ismac
    runtime = 'MONALISA native platform binaries';
    return;
end
library = fullfile(matlabroot,'bin',computer('arch'),'libomp.dylib');
assert(isfile(library),'This macOS setup requires MATLAB''s bundled libomp: %s',library);
% A foreign runtime can stay loaded even after changing the MATLAB path.
[~,loaded] = inmem('-completenames');
for k = 1:numel(loaded)
    if contains(loaded{k},'bm') && endsWith(loaded{k},['.' mexext])
        assert(isfile(loaded{k}),'A loaded MEX file was moved or removed. Restart MATLAB before setup.');
        deps = dependencies(loaded{k});
        foreign = regexp(deps,'(?m)^\s*([^\r\n]*?lib(?:i)?omp[^\r\n]*?\.dylib)[ \t]+\(compatibility version','tokens');
        for j = 1:numel(foreign)
            assert(strcmp(foreign{j}{1},library), ...
                'A MONALISA MEX already loaded a different OpenMP runtime. Restart MATLAB and run setup_tutorial first.');
        end
    end
end
cache = fullfile(cacheRoot,computer('arch'),version('-release'));
if ~isfolder(cache), mkdir(cache); end
libraryHash = file_md5(library);
for k = 1:numel(names)
    target = fullfile(cache,[names{k} '.' mexext]);
    stampFile = fullfile(cache,[names{k} '_build.mat']);
    stamp = struct('source',sources{k},'sourceHash',file_md5(sources{k}), ...
        'matlabRoot',matlabroot,'library',library,'libraryHash',libraryHash,'version',version);
    rebuild = true;
    if isfile(target) && isfile(stampFile)
        previous = load(stampFile,'stamp');
        rebuild = ~isequal(previous.stamp,stamp) || ~contains(dependencies(target),[library ' (']);
    end
    if rebuild
        fprintf('Creating/refreshing MEX cache copy:\n  source: %s\n  target: %s\n  OpenMP: %s\n', ...
            sources{k},target,library);
        assert(~any(strcmp(loaded,target)), ...
            'A cached MEX changed since it was loaded. Restart MATLAB before refreshing it.');
        copyfile(sources{k},target,'f');
        fileattrib(target,'+w','u');
        deps = dependencies(target);
        old = regexp(deps,'(?m)^\s*([^\r\n]*?lib(?:i)?omp[^\r\n]*?\.dylib)[ \t]+\(compatibility version','tokens');
        assert(isscalar(old),'Expected one OpenMP dependency in %s',sources{k});
        if ~strcmp(old{1}{1},library)
            run_command(['install_name_tool -change ' shell_quote(old{1}{1}) ' ' ...
                shell_quote(library) ' ' shell_quote(target)]);
            run_command(['codesign --force --sign - ' shell_quote(target)]);
        end
    end
    deps = dependencies(target);
    assert(contains(deps,[library ' (']),'MEX does not link to MATLAB libomp: %s',target);
    fprintf('MEX source: %s\nMEX cache:  %s\n',sources{k},target);
    if rebuild, save(stampFile,'stamp'); end
    % Clear a previously loaded, correctly linked diagnostic copy if needed.
    if ~strcmp(which(names{k}),target), clear(names{k}); end
end
addpath(cache,'-begin');
for k = 1:numel(names)
    assert(strcmp(which(names{k}),fullfile(cache,[names{k} '.' mexext])), ...
        'Another MEX shadows the generated local copy: %s. Change out of its folder.',names{k});
end
runtime = library;
end

function result = dependencies(path)
[status,result] = system(['otool -L ' shell_quote(path)]);
assert(status==0,'Cannot inspect MEX dependencies: %s',result);
end

function run_command(command)
[status,result] = system(command);
assert(status==0,'MEX runtime setup failed: %s',result);
end
