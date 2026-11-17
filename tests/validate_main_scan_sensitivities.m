function result = validate_main_scan_sensitivities(config)
%VALIDATE_MAIN_SCAN_SENSITIVITIES Check relative maps and MONALISA combination.
% Synthetic signal only; does not establish performance on a real Terra scan.
if nargin<1, setup_tutorial(); else, setup_tutorial(config); end
[x,y,z] = ndgrid(single(linspace(-1,1,12)));
coilImages = complex(zeros(12,12,12,3,'single'));
object = exp(-2*(x.^2+y.^2+z.^2));
coilImages(:,:,:,1) = object.*(1+0.1*x).*exp(1i*0.2*y);
coilImages(:,:,:,2) = object.*(0.7+0.1*y).*exp(1i*(0.3+0.2*z));
coilImages(:,:,:,3) = object.*(0.8+0.1*z).*exp(1i*(-0.4+0.2*x));
mask = object>0.15;
[C,details] = estimate_relative_sensitivities(coilImages,mask,0,1000,1e-4);
rss = sqrt(sum(abs(coilImages).^2,4));
expected = coilImages./rss;
inside = repmat(mask,[1 1 1 3]);
mapError = norm(double(C(inside)-expected(inside)))/norm(double(expected(inside)));
assert(mapError<1e-6 && all(isfinite(C(:))));
unitError = double(max(abs(sum(abs(C).^2,4)-1),[],'all'));
assert(unitError<1e-6 && isa(C,'single'));
combined = bmCoilSense_pinv(C,coilImages,[12 12 12]);
combinationError = norm(double(combined(mask)-rss(mask)))/norm(double(rss(mask)));
assert(combinationError<1e-6,'MONALISA combination disagrees with the relative-map model.');
[smoothC,~] = estimate_relative_sensitivities(coilImages,mask,1,1000,1e-4);
assert(all(isfinite(smoothC(:))) && max(abs(sum(abs(smoothC).^2,4)-1),[],'all')<1e-6);
assert(~details.bodyCoilReference);
result = struct('status','passed','insideMapRelativeError',mapError, ...
    'unitPowerMaximumError',unitError,'combinationRelativeError',combinationError);
disp(result);
end
