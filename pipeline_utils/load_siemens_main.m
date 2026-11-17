function [y,acq,info,geometry] = load_siemens_main(rawFile,seqFile,nShotOff,params)
%LOAD_SIEMENS_MAIN Read main data with sequence grouping and no legacy UI.
% Missing parameters use the shared input dialogs. Mitosius confirms nSeg/nShot
% with prefilled sequence values before raw dimensions are validated/read.
if nargin<4, params = []; end
[image,acq,info,geometry,rawLines] = siemens_sequence_metadata(rawFile,seqFile,nShotOff,params,true);
fprintf('Reading %d retained lines, %d samples/line, %d physical coils...\n', ...
    numel(rawLines),acq.N,acq.nCh);
y = read_siemens_blocks(image,acq,rawLines);
fprintf('Main raw data loaded: %s\n',mat2str(size(y)));
end
