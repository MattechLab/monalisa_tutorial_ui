function outputRoot = resolve_output_root(baseFolder,outputRoot)
%RESOLVE_OUTPUT_ROOT Put relative result folders under the raw dataset folder.
% An explicit absolute output path is retained. No directories are created.
if nargin<2 || isempty(outputRoot), outputRoot = 'recon_results'; end
outputRoot = char(outputRoot);
isAbsolute = startsWith(outputRoot,'/') || startsWith(outputRoot,'\') || ...
    ~isempty(regexp(outputRoot,'^[A-Za-z]:[\\/]','once'));
if ~isAbsolute, outputRoot = fullfile(baseFolder,outputRoot); end
end
