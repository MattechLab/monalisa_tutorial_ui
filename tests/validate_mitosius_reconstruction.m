function results = validate_mitosius_reconstruction(config)
%VALIDATE_MITOSIUS_RECONSTRUCTION Test both coil modes and both recon modes.
% Synthetic data only; use actual Mitosius I/O, gridding, and bmSteva.
if nargin<1, config = tutorial_config(); end
repoRoot = fileparts(fileparts(mfilename('fullpath')));
config.baseFolder = repoRoot; % This validation needs no external dataset.
setup_tutorial(config);
workDir = tempname; mkdir(workDir);
fileCleanup = onCleanup(@() rmdir(workDir,'s'));
N = 12; fovMm = 200; nCh = 4; nKeep = 2;
[a,b,c] = ndgrid(single(-N/2:N/2-1));
object = exp(-(a.^2+b.^2+c.^2)/8);
signal = fftshift(fftn(ifftshift(object)));
coefficients = single([1,0.8*exp(0.3i),0.6*exp(-0.2i),0.4*exp(0.5i)]);
coefficients = coefficients/norm(coefficients);
physicalC = repmat(reshape(coefficients,1,1,1,nCh),[48 48 48 1]);
D = signal(:).*coefficients;
yRaw = permute(reshape(D,N,N*N,nCh),[3 1 2]);
tRaw = reshape([a(:),b(:),c(:)].'/fovMm,3,N,N*N);
veRaw = ones(1,N,N*N,'single')/fovMm^3;
W = [1 1i 1 -1; 1i 1 -1 1; -1 1 1i 1; 1 -1 1 1i];
results = struct('compressionEnabled',{},'channels',{}, ...
    'griddedPassed',{},'compressedSensingPassed',{});
for compressed = 0:1
    if compressed
        [data,C] = apply_rovir_compression(yRaw,physicalC,W,nKeep);
        expectedChannels = nKeep;
    else
        data = yRaw; C = physicalC; expectedChannels = nCh;
    end
    reconDir = fullfile(workDir,sprintf('MID%05d_recon',compressed+1));
    mDir = fullfile(reconDir,'T1_LIBRE_woBinning','mitosius');
    [yc,tc,vc] = bmMitosis(data,tRaw,veRaw,true(1,N*N));
    yc = bmPermuteToCol(yc);
    bmMitosius_create(mDir,yc,tc,vc);
    geometry = struct('FoV',fovMm);
    save(fullfile(mDir,'recon_geometry.mat'),'geometry');
    cFile = fullfile(mDir,'C.mat'); save(cFile,'C');
    preparation=struct('reconstructionCMD5',file_md5(cFile));
    save(fullfile(mDir,'preparation.mat'),'preparation');
    complete=true; save(fullfile(mDir,'complete.mat'),'complete');
    assert(size(yc{1},2)==expectedChannels && size(C,4)==expectedChannels);
    opts = struct('reconDir',reconDir,'mitosiusDir',mDir,'CfilePath',cFile, ...
        'Matrix_size',N,'nIter',2,'delta',0.001,'rho',0.01,'flagCS',0);
    [x,initial,griddedDir] = recon_nobin_from_paths(opts);
    assert(isempty(x) && isscalar(initial) && numel(initial{1})==N^3);
    assert(all(isfinite(initial{1}(:))) && any(abs(initial{1}(:))>0));
    assert(isfile(fullfile(griddedDir,'x0_12.mat')));
    opts.flagCS = 1;
    [x,~,csDir] = recon_nobin_from_paths(opts);
    assert(~strcmp(griddedDir,csDir) && isfile(fullfile(griddedDir,'x0_12.mat')));
    assert(numel(x)==N^3 && all(isfinite(x(:))) && any(abs(x(:))>0));
    csFile = fullfile(csDir, ...
        'x_nIter2_delta_0.001_rho_0.010_12.mat');
    assert(isfile(csFile),'Compressed-sensing reconstruction was not saved.');
    run=load(fullfile(csDir,'run_info.mat'),'runInfo');
    assert(run.runInfo.complete && run.runInfo.options.delta==opts.delta);
    % Same dimensions are insufficient: reject a C with a different basis.
    wrongFile=fullfile(mDir,'wrong_C.mat'); C=C*exp(0.3i); save(wrongFile,'C');
    badOpts=opts; badOpts.CfilePath=wrongFile;
    rejected=false;
    try
        recon_nobin_from_paths(badOpts);
    catch failure
        rejected=contains(failure.message,'exact C paired');
    end
    assert(rejected,'Same-size C from another basis was accepted.');
    results(end+1) = struct('compressionEnabled',logical(compressed), ...
        'channels',expectedChannels,'griddedPassed',true,'compressedSensingPassed',true); %#ok<AGROW>
end
disp(results);
end
