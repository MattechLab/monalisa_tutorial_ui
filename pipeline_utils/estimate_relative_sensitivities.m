function [C,details] = estimate_relative_sensitivities(coilImages,mask,smoothSigma,laplaceIterations,laplaceTolerance)
%ESTIMATE_RELATIVE_SENSITIVITIES Direct main-scan relative sensitivity maps.
% Smooth complex low-resolution coil images, divide by their RSS within the
% reviewed signal mask, and harmonically extend outside using MONALISA.
% A common object phase remains in C: this is a relative receive basis,
% without a body-coil intensity reference. This is not ESPIRiT.
assert(isfloat(coilImages) && ndims(coilImages)==4 && all(isfinite(coilImages(:))), ...
    'Expected finite complex coil images [Nx,Ny,Nz,nCoils].');
grid = size(coilImages); grid = grid(1:3);
assert(islogical(mask) && isequal(size(mask),grid) && any(mask(:)) && ~all(mask(:)), ...
    'Provide a nonempty, non-full logical signal mask matching the images.');
validateattributes(smoothSigma,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(laplaceIterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(laplaceTolerance,{'numeric'},{'scalar','real','finite','positive'});
C = single(coilImages);
for channel = 1:size(C,4)
    if smoothSigma>0
        image = C(:,:,:,channel);
        C(:,:,:,channel) = complex(imgaussfilt3(real(image),smoothSigma), ...
            imgaussfilt3(imag(image),smoothSigma));
    end
end
rss = sqrt(sum(abs(C).^2,4));
assert(all(rss(mask)>0),'The signal mask contains zero-RSS voxels.');
denominator = rss; denominator(~mask) = 1;
for channel = 1:size(C,4)
    image = C(:,:,:,channel)./denominator;
    image(~mask) = 0;
    C(:,:,:,channel) = bmImLaplaceEquationSolver(image,mask, ...
        laplaceIterations,laplaceTolerance,'omp');
end
power = sqrt(sum(abs(C).^2,4));
assert(all(isfinite(C(:))) && all(power(:)>0), ...
    'Sensitivity extension produced invalid values or a zero coil basis.');
C = C./power;
details = struct('method','direct_smoothed_RSS_relative', ...
    'smoothSigmaVoxels',smoothSigma,'laplaceIterations',laplaceIterations, ...
    'laplaceTolerance',laplaceTolerance,'maskVoxels',nnz(mask), ...
    'normalization','unit root-sum-of-squares across coils', ...
    'bodyCoilReference',false);
end
