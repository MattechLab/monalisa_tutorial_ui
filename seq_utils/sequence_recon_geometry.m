function geometry = sequence_recon_geometry(params,nSample,nReadouts,allowAnisotropic)
%SEQUENCE_RECON_GEOMETRY Validate required header values against raw dimensions.
% FOV definition is in metres before readout oversampling, as in both supplied
% examples: 100 mm * 2 -> 200 mm; 120 mm * 1 -> 120 mm.
% The explicit anisotropic path uses physical FOV and res_mm, with ADC
% oversampling represented by the actual coordinates rather than FOV scaling.
if nargin<4, allowAnisotropic=false; end
required = {'fov','nseg','nshot','nsam','flagselfnav','readoutoversamplingfactor'};
assert(all(isfield(params,required)),'Required reconstruction header fields are unresolved.');
validateattributes(params.nseg,{'numeric'},{'scalar','integer','positive'});
validateattributes(params.nshot,{'numeric'},{'scalar','integer','positive'});
validateattributes(params.nsam,{'numeric'},{'scalar','integer','positive'});
validateattributes(params.flagselfnav,{'numeric'},{'scalar','integer','>=',0,'<=',1});
validateattributes(params.readoutoversamplingfactor,{'numeric'},{'scalar','real','finite','>=',1});
validateattributes(params.fov,{'numeric'},{'vector','real','finite','positive'});
assert(params.nsam==nSample,'nSam=%g does not match raw samples per readout=%g.',params.nsam,nSample);
assert(params.nseg*params.nshot==nReadouts, ...
    'nSeg*nShot=%g does not match acquired readouts=%g.',params.nseg*params.nshot,nReadouts);
assert(params.nseg>params.flagselfnav,'There must be imaging readouts after excluding SI.');
fov = double(params.fov(:)')*1000;
if isscalar(fov), fov = repmat(fov,1,3); end
assert(numel(fov)==3,'FOV must have one or three dimensions.');
if allowAnisotropic
    assert(isfield(params,'res_mm'),'Anisotropic reconstruction requires res_mm.');
    resolution = double(params.res_mm(:)');
    validateattributes(resolution,{'numeric'},{'numel',3,'positive','finite','real'});
    matrixSize = fov./resolution;
    assert(all(abs(matrixSize-round(matrixSize))<1e-6) && all(mod(round(matrixSize),2)==0), ...
        'FOV/res_mm must give three even integer matrix dimensions.');
    geometry = struct('FoV',fov,'matrixSize',round(matrixSize),'voxelSizeMm',resolution);
else
    % Preserve the established isotropic reader's oversampling convention.
    fov = fov*params.readoutoversamplingfactor;
    assert(all(abs(fov-fov(1))<1e-6), ...
        'Use the anisotropic scanner xRMS cell for anisotropic FOV data.');
    geometry = struct('FoV',fov(1));
end
geometry.nSeg = params.nseg;
geometry.nShot = params.nshot;
geometry.nSample = params.nsam;
geometry.excludeSI = logical(params.flagselfnav);
end
