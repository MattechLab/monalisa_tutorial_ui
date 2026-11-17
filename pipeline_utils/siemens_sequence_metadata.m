function [image,acq,info,geometry,rawLines] = siemens_sequence_metadata(rawFile,seqFile,nShotOff,params,confirmGrouping,allowAnisotropic)
%SIEMENS_SEQUENCE_METADATA Resolve sequence parameters before reading raw data.
% Uses physical nSeg/nShot from the sequence, not the Siemens segment counter.
% A supplied resolved header is reused only for the matching sequence path.
if nargin<4, params = []; end
if nargin<5, confirmGrouping = false; end
if nargin<6, allowAnisotropic = false; end
twix = mapVBVD_JH_for_monalisa(rawFile);
if iscell(twix), twix = twix{end}; end
image = twix.image;
headerSignature = char(twix.hdr.Dicom.tSequenceVariant);
assert(check_hash(rawFile,seqFile,headerSignature),'Sequence/header mismatch.');
if isempty(params)
    required = {'FOV','ReadoutOversamplingFactor','nSam','nSeg','nShot','flagSelfNav'};
    params = extract_seq_params(seqFile,required,struct('nSam',double(image.NCol)));
else
    assert(strcmp(params.headerprovenance.sourceFile,seqFile),'Resolved parameters belong to another sequence.');
end
if confirmGrouping
    params = confirm_sequence_grouping(params,double(image.NAcq));
end
geometry = sequence_recon_geometry(params,double(image.NCol),double(image.NAcq),allowAnisotropic);
nSeg = geometry.nSeg; nShot = geometry.nShot;
excludeSI = geometry.excludeSI;
validateattributes(nShotOff,{'numeric'},{'scalar','integer','nonnegative','<',nShot});
assert(double(image.NEco)==1,'Only single-echo Siemens data are supported here.');
assert(~any(image.IsReflected),'Blocked reading requires non-reflected readouts.');
assert(~image.flagRemoveOS && ~image.flagRampSampRegrid && ...
    ~image.flagDoRawDataCorrect,'Unexpected raw-data correction flags.');
acq = bmMriAcquisitionParam([]);
acq.N = double(image.NCol); acq.nCh = double(image.NCha);
acq.nSeg = nSeg; acq.nShot = nShot; acq.nEcho = 1;
acq.nShot_off = nShotOff; acq.selfNav_flag = excludeSI;
% Use the declared readout factor rather than assuming a fixed factor of 2.
acq.FoV = geometry.FoV;
acq.traj_type = 'pulseq'; acq.pulseqTrajFile_name = seqFile;
info = struct('rawFile',rawFile,'seqFile',seqFile,'rawStat',dir(rawFile), ...
    'headerSignature',headerSignature, ...
    'N',acq.N,'nCh',acq.nCh,'nSeg',nSeg,'nShot',nShot, ...
    'nShotOff',nShotOff,'FoV',acq.FoV,'sequenceParams',params);
geometry.nShotOff = nShotOff;
geometry.nPhysicalCoils = acq.nCh;
info.coilSelection = []; info.sliceGeometry = [];
if isfield(twix.hdr.Dicom,'flMagneticFieldStrength')
    info.B0_T = double(twix.hdr.Dicom.flMagneticFieldStrength);
end
if isfield(twix.hdr,'MeasYaps')
    if isfield(twix.hdr.MeasYaps,'sCoilSelectMeas')
        info.coilSelection = twix.hdr.MeasYaps.sCoilSelectMeas;
    end
    if isfield(twix.hdr.MeasYaps,'sSliceArray')
        info.sliceGeometry = twix.hdr.MeasYaps.sSliceArray;
    end
end
rawLines = reshape(1:geometry.nSeg*geometry.nShot,geometry.nSeg,geometry.nShot);
if geometry.excludeSI, rawLines(1,:) = []; end
rawLines(:,1:nShotOff) = []; rawLines = rawLines(:);
fprintf('Sequence grouping: nSeg=%d, nShot=%d; verified against %d acquired readouts.\n', ...
    acq.nSeg,acq.nShot,double(image.NAcq));
end
