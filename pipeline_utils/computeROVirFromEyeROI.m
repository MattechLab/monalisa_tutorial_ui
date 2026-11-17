function [W, lambda, masksOut] = computeROVirFromEyeROI(x0, mask, irMask2D, xrms, varargin)
debugFlag = 0;
% computeROVirFromEyeROI
%   Compute ROVir (GEVD) virtual coil weights from an eye ROI and an interference region.
%
% Options (name-value):
%   'IRMode'          : 'manual' (default) | 'complement' | 'complement_object'
%   'DilateROIPixels' : nonnegative integer (default 1)
%   'ObjectThresh'    : in [0,1], relative threshold on abs(xrms) for object mask (default 0.01)
%   'ZRange'          : [zStart zEnd] slice bounds for 3D mask extrusion (default full volume)
%   'SelectZRange'    : true to interactively choose ZRange from a sagittal view of xrms
%   'IRGapPixels'     : 3D safety margin used when IR is auto-derived from ROI complement (default 1)
%   'CheckMasks'      : true to open an overlay viewer for ROI/IR masks on xrms (default false)
%   'EpsFactor'       : regularization factor multiplier (default 1e-6)
% Usage:
% [W, lambda] = computeROVirFromEyeROI(x0, mask, irMask2D);

% [W, lambda] = computeROVirFromEyeROI(x0, mask, [], xrms, ...
    % 'IRMode','complement', 'DilateROIPixels', 6);

% [W, lambda] = computeROVirFromEyeROI(x0, mask, [], xrms, ...
    % 'IRMode','complement_object', 'DilateROIPixels', 6, 'ObjectThresh', 0.05);
% ---------------- Options ----------------
p = inputParser;
p.addParameter('IRMode', 'manual', @(s) ischar(s) || isstring(s));
p.addParameter('DilateROIPixels', 1, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('ObjectThresh', 0.01, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('ZRange', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.addParameter('SelectZRange', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('IRGapPixels', 1, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('CheckMasks', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('EpsFactor', 1e-6, @(x) isnumeric(x) && isscalar(x) && x>0);
p.parse(varargin{:});
opt = p.Results;

opt.IRMode = lower(string(opt.IRMode));

% ---------------- Load x0 if it is a file path ----------------
if ischar(x0) || isstring(x0)
    x0path = char(x0);
    if ~exist(x0path, 'file')
        error('computeROVirFromEyeROI:FileNotFound', 'x0 path not found: %s', x0path);
    end
    s = load(x0path);
    if ~isfield(s, 'x0')
        error('computeROVirFromEyeROI:MissingVar', 'MAT file does not contain variable ''x0'': %s', x0path);
    end
    x0 = s.x0;
end

% ---------------- Convert x0 -> g_cal 4D ----------------
g_cal = x0;
if iscell(g_cal)
    try
        g_cal = cat(4, g_cal{:});   % [Nx,Ny,Nz,nCh]
    catch ME
        error('computeROVirFromEyeROI:BadCell', ...
            'Failed to cat(4, x0{:}). Ensure each cell contains same-size [Nx,Ny,Nz].\n%s', ME.message);
    end
end
if ndims(g_cal) ~= 4
    error('computeROVirFromEyeROI:BadX0', 'x0 must be cell of 3D or 4D [Nx,Ny,Nz,nCh].');
end

[Nx, Ny, Nz, nCh] = size(g_cal);

% ---------------- Get a 2D view of xrms for interactive drawing / object mask ----------------
xrms2D = [];
xrms3D = [];
sliceIndex = 1;
if ~isempty(xrms)
    if ndims(xrms) == 3
        xrms3D = xrms;
        sliceIndex = round(size(xrms,3)/2);
        xrms2D = xrms(:,:,sliceIndex);
    else
        xrms2D = xrms;
    end
end

% ---------------- If ROI missing: interactive draw (needs xrms) ----------------
if isempty(mask)
    if isempty(xrms2D)
        error(['mask is empty but xrms is not provided. ', ...
               'Tip: run check_orient_xrms to generate xrms, then call this function again.']);
    end
    figure; imagesc(abs(xrms2D)); axis image off; colormap gray;
    title('Draw ROI around the eyes, double-click to finish');
    mask = logical(roipoly);
    close(gcf);
    warning('computeROVirFromEyeROI:InteractiveROI', ...
        'ROI mask was empty. Used interactive roipoly (sliceIndex=%d).', sliceIndex);
end

% ---------------- Validate ROI mask size ----------------
if ~isequal(size(mask), [Nx, Ny])
    error('computeROVirFromEyeROI:MaskSizeMismatch', ...
        'ROI mask size must match [Nx,Ny]=[%d,%d], but got [%d,%d].', ...
        Nx, Ny, size(mask,1), size(mask,2));
end

% ---------------- Optionally dilate ROI (2D) ----------------
roi2D = logical(mask);
if opt.DilateROIPixels > 0
    se = strel('disk', opt.DilateROIPixels, 0);
    roi2D = imdilate(roi2D, se);
end

% ---------------- Select z-range for 3D extrusion ----------------
if isempty(opt.ZRange)
    zRange = [1, Nz];
    if opt.SelectZRange
        if isempty(xrms3D)
            error(['SelectZRange=true requires xrms as a 3D volume. ', ...
                   'Tip: run check_orient_xrms to generate xrms, then call this function again.']);
        end
        zRange = selectZRangeFromSagittal(xrms3D, roi2D);
    end
else
    zRange = round(double(opt.ZRange(:)'));
    zRange = sort(zRange);
    if numel(zRange) ~= 2 || any(~isfinite(zRange))
        error('computeROVirFromEyeROI:BadZRange', ...
            'ZRange must be a finite 2-element numeric vector [zStart zEnd].');
    end
    zRange(1) = max(1, zRange(1));
    zRange(2) = min(Nz, zRange(2));
    if zRange(1) > zRange(2)
        error('computeROVirFromEyeROI:BadZRange', ...
            'ZRange must satisfy zStart <= zEnd after clamping to [1, %d].', Nz);
    end
end

% ---------------- Build IR mask according to IRMode ----------------
irUsesFullComplement3D = false;
switch opt.IRMode
    case "manual"
        if isempty(irMask2D)
            ir2D = ~roi2D;
            irUsesFullComplement3D = true;
            warning('computeROVirFromEyeROI:AutoIRComplement', ...
                ['irMask2D is empty. Using IR as the full 3D complement of the ROI ', ...
                 'instead of drawing a separate IR mask.']);
        else
            if ~isequal(size(irMask2D), [Nx, Ny])
                error('computeROVirFromEyeROI:IRMaskSizeMismatch', ...
                    'IR mask size must match [Nx,Ny]=[%d,%d], but got [%d,%d].', ...
                    Nx, Ny, size(irMask2D,1), size(irMask2D,2));
            end
            ir2D = logical(irMask2D);
        end

    case "complement"
        % IR = everything except (dilated) ROI
        ir2D = ~roi2D;
        irUsesFullComplement3D = true;

    case "complement_object"
        % IR = object mask (from xrms) minus (dilated) ROI
        if isempty(xrms2D)
            error(['IRMode=complement_object requires xrms. ', ...
                   'Tip: run check_orient_xrms to generate xrms, then call this function again.']);
        end
        mag = abs(xrms2D);
        thr = opt.ObjectThresh * max(mag(:));
        obj2D = mag > thr;
        % optional cleanup
        obj2D = imfill(obj2D, 'holes');
        obj2D = bwareaopen(obj2D, 50);
        ir2D = obj2D & ~roi2D;

    otherwise
        error('computeROVirFromEyeROI:BadIRMode', ...
            'Unknown IRMode: %s. Use ''manual'', ''complement'', or ''complement_object''.', opt.IRMode);
end

% ---------------- Expand 2D masks to 3D ----------------
zMask = false(1, 1, Nz);
zMask(1, 1, zRange(1):zRange(2)) = true;
maskROI = repmat(roi2D, [1 1 Nz]) & zMask;
if irUsesFullComplement3D
    if opt.IRGapPixels > 0
        seIRGap = strel('sphere', opt.IRGapPixels);
        roiExclusion3D = imdilate(maskROI, seIRGap);
    else
        roiExclusion3D = maskROI;
    end
    maskIR = ~roiExclusion3D;
else
    maskIR = repmat(ir2D, [1 1 Nz]) & zMask;
end

% Remove overlap if any
if any(maskROI(:) & maskIR(:))
    warning('computeROVirFromEyeROI:MaskOverlap', 'ROI and IR overlap. Removing overlap from IR.');
    maskIR(maskROI) = false;
end

idxO = maskROI(:);
idxG = maskIR(:);

if nnz(idxO) == 0, error('computeROVirFromEyeROI:EmptyROI', 'ROI mask has 0 voxels.'); end
if nnz(idxG) == 0, error('computeROVirFromEyeROI:EmptyIR',  'IR mask has 0 voxels.'); end

% ---------------- Form A and B ----------------
G = reshape(g_cal, [], nCh);   % [Nvox, nCh]
G_O = G(idxO, :);
G_G = G(idxG, :);

A = (G_O' * G_O) / nnz(idxO);
B = (G_G' * G_G) / nnz(idxG);

% enforce Hermitian
A = (A + A') / 2;
B = (B + B') / 2;
if debugFlag
    % Diagnostics
    fprintf('Nroi=%d, Nir=%d\n', nnz(idxO), nnz(idxG));
    fprintf('trace(A)=%.3e, trace(B)=%.3e\n', real(trace(A)), real(trace(B)));
    fprintf('rcond(B)=%.3e\n', rcond(B));
    
    eA = eig(A); eB = eig(B);
    fprintf('eig(A) min/median/max = %.3e / %.3e / %.3e\n', min(real(eA)), median(real(eA)), max(real(eA)));
    fprintf('eig(B) min/median/max = %.3e / %.3e / %.3e\n', min(real(eB)), median(real(eB)), max(real(eB)));
end
% ---------------- Regularize B ----------------
tB = real(trace(B));
if ~isfinite(tB) || tB <= 0
    tB = norm(B, 'fro');
end
epsReg = opt.EpsFactor * (tB / nCh);
if ~isfinite(epsReg) || epsReg == 0
    epsReg = opt.EpsFactor;
end
Breg = B + epsReg * eye(nCh);

% ---------------- GEVD ----------------
[W, d] = eig(A, Breg, 'vector');
lambda = real(d);
[lambda, order] = sort(lambda, 'descend');
W = W(:, order);
W = W ./ vecnorm(W, 2, 1);

% ---------------- Return masks (optional) ----------------
masksOut = struct();
masksOut.roi2D = roi2D;
masksOut.ir2D  = ir2D;
masksOut.zRange = zRange;
masksOut.zMask = zMask;
masksOut.roi3D = maskROI;
masksOut.ir3D  = maskIR;

if opt.CheckMasks
    if isempty(xrms3D)
        warning('computeROVirFromEyeROI:CheckMasksNeeds3D', ...
            'CheckMasks=true requires xrms as a 3D volume. Skipping overlay viewer.');
    else
        visualizeMaskOverlay3D(xrms3D, maskROI, maskIR);
    end
end

end

function zRange = selectZRangeFromSagittal(xrms3D, roi2D)
magVol = abs(xrms3D);
[Nx, Ny, Nz] = size(magVol);

ySupport = find(any(roi2D, 1));
if isempty(ySupport)
    ySupport = 1:Ny;
end

xSupport = find(any(roi2D, 2));
if isempty(xSupport)
    xSupport = 1:Nx;
end

% Collapse over ROI-bearing y-columns so the eyeball remains visible in the sagittal view.
sagittalRef = squeeze(max(magVol(:, ySupport, :), [], 2));

f = figure('Name', 'Select Z range from sagittal view', 'NumberTitle', 'off');
imagesc(1:Nz, 1:Nx, sagittalRef);
axis image;
axis xy;
colormap gray;
xlabel('z slice');
ylabel('x');
title({
    'Draw a box around the eyeball in sagittal view'
    'Only the horizontal extent is used as ZRange; double-click inside the box when done'
});

defaultPos = [1, min(xSupport), Nz, max(1, numel(xSupport))];
h = drawrectangle(gca, 'Position', defaultPos, 'Color', 'g');
wait(h);

pos = round(h.Position);
zStart = max(1, min(Nz, pos(1)));
zEnd = max(1, min(Nz, pos(1) + pos(3) - 1));
zRange = sort([zStart, zEnd]);

close(f);
end

function visualizeMaskOverlay3D(xrms3D, maskROI, maskIR)
magVol = abs(xrms3D);
[Nx, Ny, Nz] = size(magVol);

ix = max(1, min(Nx, round(Nx / 2)));
iy = max(1, min(Ny, round(Ny / 2)));
iz = max(1, min(Nz, round(Nz / 2)));

f = figure('Name', 'ROI/IR mask overlay', 'NumberTitle', 'off', ...
    'Position', [100 100 1400 700]);

ax1 = subplot(1,3,1, 'Parent', f);
ax2 = subplot(1,3,2, 'Parent', f);
ax3 = subplot(1,3,3, 'Parent', f);

sldX = uicontrol(f, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.10 0.06 0.22 0.035], 'Min', 1, 'Max', Nx, 'Value', ix, ...
    'SliderStep', [1/max(Nx-1,1), min(10/max(Nx-1,1),1)]);
sldY = uicontrol(f, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.39 0.06 0.22 0.035], 'Min', 1, 'Max', Ny, 'Value', iy, ...
    'SliderStep', [1/max(Ny-1,1), min(10/max(Ny-1,1),1)]);
sldZ = uicontrol(f, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.68 0.06 0.22 0.035], 'Min', 1, 'Max', Nz, 'Value', iz, ...
    'SliderStep', [1/max(Nz-1,1), min(10/max(Nz-1,1),1)]);

    function updateViews(~, ~)
        ix = round(get(sldX, 'Value'));
        iy = round(get(sldY, 'Value'));
        iz = round(get(sldZ, 'Value'));

        drawOverlay(ax1, magVol(:, :, iz), maskROI(:, :, iz), maskIR(:, :, iz), ...
            sprintf('Axial z=%d', iz));
        drawOverlay(ax2, squeeze(magVol(:, iy, :)), squeeze(maskROI(:, iy, :)), squeeze(maskIR(:, iy, :)), ...
            sprintf('Sagittal y=%d', iy));
        drawOverlay(ax3, squeeze(magVol(ix, :, :))', squeeze(maskROI(ix, :, :))', squeeze(maskIR(ix, :, :))', ...
            sprintf('Coronal x=%d', ix));
    end

set(sldX, 'Callback', @updateViews);
set(sldY, 'Callback', @updateViews);
set(sldZ, 'Callback', @updateViews);
updateViews();
end

function drawOverlay(ax, baseSlice, roiSlice, irSlice, plotTitle)
imagesc(ax, baseSlice);
axis(ax, 'image');
axis(ax, 'off');
colormap(ax, gray);
hold(ax, 'on');

roiRgb = cat(3, zeros(size(roiSlice)), ones(size(roiSlice)), zeros(size(roiSlice)));
irRgb  = cat(3, ones(size(irSlice)), zeros(size(irSlice)), zeros(size(irSlice)));

hROI = imagesc(ax, roiRgb);
set(hROI, 'AlphaData', 0.5 * double(roiSlice));
hIR = imagesc(ax, irRgb);
set(hIR, 'AlphaData', 0.5 * double(irSlice));

title(ax, sprintf('%s | ROI=green IR=red', plotTitle));
hold(ax, 'off');
end
