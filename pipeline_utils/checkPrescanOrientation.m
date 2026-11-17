function [orientation_ok,ncc_all,details] = checkPrescanOrientation(xsos_all,showFigures,options)
%CHECKPRESCANORIENTATION Compare HC/main with BC over 48 permutations/flips.
% xsos_all contains the normalized 64^3 SoS previews {BC,HC,main}.
% This checks array axes and flips, not rigid-body registration or head motion.
if nargin<2, showFigures=true; end
if nargin<3, options=struct(); end
if ~isfield(options,'minimumCorrelation'), options.minimumCorrelation=[]; end
if ~isfield(options,'minimumMargin'), options.minimumMargin=[]; end
if ~isempty(options.minimumCorrelation)
    validateattributes(options.minimumCorrelation,{'numeric'},{'scalar','real','finite','>=',-1,'<=',1});
end
if ~isempty(options.minimumMargin)
    validateattributes(options.minimumMargin,{'numeric'},{'scalar','real','finite','nonnegative'});
end
details=struct('identityScore',{},'bestScore',{},'margin',{},'bestPermutation',{},'bestFlip',{},'passed',{});
ref = xsos_all{1}; pp = perms(1:3); ff = dec2bin(0:7)-'0';
orientation_ok = true; ncc_all = cell(2,1);
for scan = 2:3
    assert(isequal(size(ref),size(xsos_all{scan})),'Preview grids differ.');
    bestScore=-Inf; index=0; values=zeros(48,1);
    for p=1:6
        for f=1:8
            index=index+1; tmp=permute(xsos_all{scan},pp(p,:));
            for axisIndex=1:3
                if ff(f,axisIndex), tmp=flip(tmp,axisIndex); end
            end
            values(index)=corr(ref(:),tmp(:));
            if isequal(pp(p,:),[1 2 3]) && isequal(ff(f,:),[0 0 0])
                identityIndex=index;
            end
            if values(index)>bestScore
                bestScore=values(index); bestPerm=pp(p,:); bestFlip=ff(f,:);
            end
        end
    end
    assert(isfinite(bestScore),'Invalid orientation correlation.');
    ncc_all{scan-1}=values;
    aligned=isequal(bestPerm,[1 2 3]) && isequal(bestFlip,[0 0 0]);
    ordered=sort(values,'descend'); margin=ordered(1)-ordered(2);
    if ~isempty(options.minimumCorrelation), aligned=aligned && bestScore>=options.minimumCorrelation; end
    if ~isempty(options.minimumMargin), aligned=aligned && margin>=options.minimumMargin; end
    details(scan-1)=struct('identityScore',values(identityIndex),'bestScore',bestScore, ...
        'margin',margin,'bestPermutation',bestPerm,'bestFlip',bestFlip,'passed',aligned);
    orientation_ok=orientation_ok && aligned;
    fprintf('Scan %d: axes=%s flips=%s correlation=%.6f\n', ...
        scan,mat2str(bestPerm),mat2str(bestFlip),bestScore);
    fprintf('Identity correlation=%.6f; best-to-runner-up margin=%.6f\n',values(identityIndex),margin);
    if showFigures
        figure; bar(values); xlabel('Permutation/flip combination'); ylabel('Correlation');
        title(sprintf('Scan %d relative to BC',scan));
    end
end
end
