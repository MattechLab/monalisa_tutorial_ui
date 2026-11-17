function [yCompressed,CCompressed,info] = apply_rovir_compression(y,C,W,nVirtualCoils)
%APPLY_ROVIR_COMPRESSION Apply the existing HPC method-3 coil mixing.
% y: [physical coils, samples, readouts]; C: [Nx,Ny,Nz,physical coils].
% W columns are ordered ROVir vectors. QR preserves their selected subspace.
% Matches Recon_scripts/hpc/scripts/pipeline_fixed/helper/
% mitosius_write_from_mask.m: compress_rovir, including complex conjugation.
nCh = size(y,1);
assert(isfloat(y) && ~isempty(y) && ndims(y)<=3,'y must be a nonempty coil-first numeric array.');
assert(isfloat(C) && size(C,4)==nCh,'C and raw data must have the same physical coils.');
assert(isfloat(W) && ismatrix(W) && size(W,1)==nCh && all(isfinite(W(:))), ...
    'W must be a finite matrix with one row per physical coil.');
assert(isnumeric(nVirtualCoils) && isscalar(nVirtualCoils) && ...
    isfinite(nVirtualCoils) && nVirtualCoils==fix(nVirtualCoils) && ...
    nVirtualCoils>=1 && nVirtualCoils<=min(size(W)), ...
    'nVirtualCoils must be an integer between 1 and min(size(W)).');
Wk = W(:,1:nVirtualCoils);
assert(rank(Wk)==nVirtualCoils,'Selected ROVir vectors are linearly dependent.');
[Q,~] = qr(Wk,0);
% Q' on coil-first data equals Draw*conj(Q) in the original HPC code.
yCompressed = reshape(cast(Q','like',y)*reshape(y,nCh,[]), ...
    nVirtualCoils,size(y,2),size(y,3));
CCompressed = reshape(reshape(C,[],nCh)*cast(conj(Q),'like',C), ...
    size(C,1),size(C,2),size(C,3),nVirtualCoils);
info = struct('method','ROVir','nPhysicalCoils',nCh, ...
    'nVirtualCoils',nVirtualCoils,'orthonormalWeights',Q);
end
