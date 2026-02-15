function [orientation_ok, xsos_all, ncc_all] = checkPrescanOrientation(readerBC, readerHC, readerMain, frameSize)
% checkPrescanOrientation - Verify orientation consistency of prescans
%
% This function reconstructs sum-of-squares (SoS) images from prescan raw data
% using provided reader objects, and checks whether the coil prescans
% (HC prescan and main scan) are aligned with the body coil prescan (BC).
%
% It computes full-volume pixel-wise correlation for all 48 possible 
% permutations and flips, prints the best matching orientation, and
% optionally plots bar charts and central slice comparisons.
%
% Usage:
%   orientation_ok = checkPrescanOrientation(readerBC, readerHC, readerMain)
%   [orientation_ok, xsos_all, ncc_all] = checkPrescanOrientation(readerBC, readerHC, readerMain, 64)
%
% Inputs:
%   readerBC    - bmRawDataReader for body coil prescan (reference)
%   readerHC    - bmRawDataReader for high-res prescan
%   readerMain  - bmRawDataReader for main scan
%   frameSize   - optional, integer specifying reconstruction size (default 64)
%
% Outputs:
%   orientation_ok - boolean flag: true if both HC and main scan are aligned to BC
%   xsos_all       - cell array of SoS images {BC, HC, Main}
%   ncc_all        - cell array of correlation values for scan2 & scan3

if nargin < 4
    frameSize = 64;
end

FrameSize = [frameSize, frameSize, frameSize];
readers = {readerBC, readerHC, readerMain};
numScans = 3;

xsos_all = cell(numScans,1);
ncc_all = cell(numScans-1,1); % only for scan2 and scan3
orientation_ok = true;         % assume good orientation initially

%% Step 1: Reconstruct SoS images
for i = 1:numScans
    reader = readers{i};
    
    % Compute coil images
    [y, t, v] = bmCoilSense_nonCart_data(reader, FrameSize);
    x = bmMathilda(y, t, v, [], FrameSize, FrameSize, [1 1 1]./reader.acquisitionParams.FoV);
    
    % Sum-of-squares
    [nx, ny, nz, numCoils] = size(x);
    xsos = zeros(nx, ny, nz, 'single');
    for c = 1:numCoils
        xsos = xsos + abs(x(:,:,:,c)).^2;
    end
    xsos = sqrt(xsos);
    
    % Log + percentile normalization
    xsos = log1p(xsos);
    lo = prctile(xsos(:),1);
    hi = prctile(xsos(:),99);
    xsos = (xsos - lo)/(hi-lo);
    xsos = max(min(xsos,1),0);
    
    xsos_all{i} = xsos;
end

%% Step 2: Orientation check (HC & Main relative to BC)
ref = xsos_all{1};
permsXYZ = perms(1:3);
flips = dec2bin(0:7) - '0';

for i = 2:numScans
    mov = xsos_all{i};
    bestScore = -Inf;
    bestPerm  = [];
    bestFlip  = [];
    bestTmp   = [];
    
    allScores = zeros(size(permsXYZ,1)*size(flips,1),1);
    labels = strings(length(allScores),1);
    idx = 0;
    
    for p = 1:size(permsXYZ,1)
        for f = 1:size(flips,1)
            idx = idx + 1;
            tmp = permute(mov, permsXYZ(p,:));
            if flips(f,1), tmp = flip(tmp,1); end
            if flips(f,2), tmp = flip(tmp,2); end
            if flips(f,3), tmp = flip(tmp,3); end
            
            % Pixel-wise correlation
            score = corr(ref(:), tmp(:));
            allScores(idx) = score;
            labels(idx) = sprintf('perm[%d%d%d]-flip[%d%d%d]', permsXYZ(p,:), flips(f,:));
            
            if score > bestScore
                bestScore = score;
                bestPerm  = permsXYZ(p,:);
                bestFlip  = flips(f,:);
                bestTmp   = tmp;
            end
        end
    end
    
    ncc_all{i-1} = allScores; % store all correlations
    
    % Plot bar chart of correlations
    figure;
    bar(allScores);
    xticks(1:length(allScores));
    xticklabels(labels);
    xtickangle(90);
    ylabel('Pixel-wise correlation');
    title(sprintf('Scan %d: correlations for all 48 orientations', i));
    grid on;
    
    % Warn & set orientation flag only if axes differ
    if ~isequal(bestPerm, [1 2 3]) || ~isequal(bestFlip, [0 0 0])
        warning('Scan %d orientation mismatch! Best perm [%d %d %d], flip [%d %d %d], corr=%.3f', ...
            i, bestPerm, bestFlip, bestScore);
        orientation_ok = false;
        
        % Show central slice comparison
        k = round(size(ref,3)/2);
        figure;
        bmImage(cat(2, ref(:,:,k), bestTmp(:,:,k)));
        title(sprintf('Reference vs Scan %d (best orientation)', i));
    else
        fprintf('Scan %d aligned (perm [1 2 3], flip [0 0 0], corr=%.3f)\n', i, bestScore);
    end
end

end
