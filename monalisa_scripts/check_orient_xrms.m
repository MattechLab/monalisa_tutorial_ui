function [xrms,x0,acq,rawInfo] = check_orient_xrms(baseFolder,matrix_size,nShotOff,files,params)
%CHECK_ORIENT_XRMS Optional main-data preview using bounded raw-data reads.
% files.mainRaw/mainSeq bypass pickers; params reuses resolved sequence fields.
if nargin<2 || isempty(matrix_size), matrix_size=48; end
if nargin<3 || isempty(nShotOff), nShotOff=50; end
if nargin<4, files=struct(); end
if nargin<5, params=[]; end
if isfield(files,'mainSeq') && ~isempty(files.mainSeq)
    seqFile=files.mainSeq;
else
selection_reminder('Preview 1/2: main sequence', ...
    'Select the MAIN scan Pulseq sequence (.seq) for the preview.');
[name,folder]=uigetfile('*.seq','Select main measurement sequence',baseFolder);
assert(~isequal(name,0),'Sequence selection cancelled.');
seqFile=fullfile(folder,name);
end
if isfield(files,'mainRaw') && ~isempty(files.mainRaw)
    rawFile=files.mainRaw;
else
selection_reminder('Preview 2/2: main raw data', ...
    sprintf('Select the MAIN raw-data file (.dat) matching this sequence:\n%s',seqFile));
[name,folder]=uigetfile('*.dat','Select main raw data',baseFolder);
assert(~isequal(name,0),'Raw-data selection cancelled.');
rawFile=fullfile(folder,name);
end
grid=[matrix_size matrix_size matrix_size];
[y,t,ve,acq,rawInfo]=load_siemens_lowres(rawFile,seqFile,grid,nShotOff,params);

x0=cell(acq.nCh,1); sumSquares=zeros(grid,'single');
for channel=1:acq.nCh
    x0{channel}=bmMathilda(y(:,channel),t,ve,[],grid,grid,[1 1 1]./(acq.FoV));
    % warning('hiiiii change it back later!');
    sumSquares=sumSquares+real(x0{channel}.*conj(x0{channel}));
end
xrms=sqrt(sumSquares/acq.nCh);
end
