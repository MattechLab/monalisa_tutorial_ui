function result = validate_rovir_compression(config)
%VALIDATE_ROVIR_COMPRESSION Check complex mixing and actual Mitosius round-trip.
% Synthetic inputs only; requires the tutorial's installed MONALISA dependency.
if nargin<1, setup_tutorial(); else, setup_tutorial(config); end
previousRng = rng; rngCleanup = onCleanup(@() rng(previousRng));
rng(31);
nCh = 4; nKeep = 2; nSample = 8; nLine = 8;
C = complex(single(randn(4,4,4,nCh)),single(randn(4,4,4,nCh)));
object = complex(single(randn(64,1)),single(randn(64,1)));
D = reshape(C,[],nCh).*object;
y = permute(reshape(D,nSample,nLine,nCh),[3 1 2]);
W = complex(randn(nCh),randn(nCh));
[actualY,actualC,info] = apply_rovir_compression(y,C,W,nKeep);

% Independent original HPC formula: sample-first reshaping and conj(Q).
[Q,~] = qr(W(:,1:nKeep),0);
Draw = reshape(permute(y,[2 3 1]),[],nCh);
expectedD = double(Draw)*conj(Q);
actualD = reshape(permute(actualY,[2 3 1]),[],nKeep);
dataError = norm(double(actualD)-expectedD,'fro')/norm(expectedD,'fro');
forwardD = reshape(actualC,[],nKeep).*object;
forwardError = norm(double(actualD-forwardD),'fro')/norm(double(actualD),'fro');
assert(dataError<1e-6 && forwardError<1e-6,'Complex mixing differs between raw data and C.');
assert(isa(actualY,'single') && isa(actualC,'single'),'Single precision was not preserved.');
assert(norm(info.orthonormalWeights'*info.orthonormalWeights-eye(nKeep),'fro')<1e-12);
[fullY,fullC] = apply_rovir_compression(y,C,W,nCh);
fullD = reshape(permute(fullY,[2 3 1]),[],nCh);
fullCmat = reshape(fullC,[],nCh);
recovered = sum(conj(fullCmat).*fullD,2)./sum(abs(fullCmat).^2,2);
recoveryError = norm(double(recovered-object))/norm(double(object));
assert(recoveryError<1e-6,'Full-rank coil mixing changed the combined image.');
assert_rejected(@() apply_rovir_compression(y,C,W,nCh+1));
assert_rejected(@() apply_rovir_compression(y,C,W(1:3,:),nKeep));
assert_rejected(@() apply_rovir_compression(y,C,ones(nCh),nKeep));

workDir = tempname; mkdir(workDir);
fileCleanup = onCleanup(@() rmdir(workDir,'s'));
t = single(reshape(1:3*nSample*nLine,3,nSample,nLine));
ve = ones(1,nSample,nLine,'single');
mask = [mod(1:nLine,2)==1; mod(1:nLine,2)==0];
[yc,tc,vc] = bmMitosis(actualY,t,ve,mask);
yc = bmPermuteToCol(yc);
bmMitosius_create(workDir,yc,tc,vc);
loadedY = bmMitosius_load(workDir,'y');
loadedT = bmMitosius_load(workDir,'t');
loadedV = bmMitosius_load(workDir,'ve');
assert(isequal(loadedY,yc) && isequal(loadedT,tc) && isequal(loadedV,vc), ...
    'Mitosius serialization changed the compressed data, trajectory, or weights.');
assert(all(cellfun(@(frame) size(frame,2)==nKeep,loadedY)));
result = struct('status','passed','dataRelativeError',dataError, ...
    'forwardModelRelativeError',forwardError,'fullRankRecoveryError',recoveryError);
disp(result);
end

function assert_rejected(action)
rejected = false;
try
    action();
catch
    rejected = true;
end
assert(rejected,'Invalid compression inputs were accepted.');
end
