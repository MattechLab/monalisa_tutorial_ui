function [y,t,ve,acq,info] = load_siemens_lowres(rawFile,seqFile,gridSize,nShotOff,params)
%LOAD_SIEMENS_LOWRES Read Pulseq radial data without loading the full raw file.
% Retains the original bmTraj, Voronoi weights and bmLowRes selection. Physical
% shot/segment counts come from the sequence, not Siemens counter meanings.
validateattributes(gridSize,{'numeric'},{'vector','numel',3,'integer','positive'});
assert(all(mod(gridSize,2)==0),'Grid dimensions must be even.');
if nargin<5, params = []; end
[image,acq,info,geometry,rawLines] = siemens_sequence_metadata(rawFile,seqFile,nShotOff,params);
nSeg = geometry.nSeg; nShot = geometry.nShot; excludeSI = geometry.excludeSI;
fprintf('Loading %d channels: nSeg=%d nShot=%d, discard=%d shots\n',acq.nCh,nSeg,nShot,nShotOff);
t = bmPointReshape(bmTraj(acq));
assert(size(t,2)==acq.N*(nSeg-excludeSI)*(nShot-nShotOff),'Trajectory/readout layout mismatch.');
ve = bmVolumeElement(t,'voronoi_full_radial3');
[selected,t,ve] = bmLowRes(1:size(t,2),t,ve,gridSize,[1 1 1]./acq.FoV);
assert(all(isfinite(t(:))) && all(isfinite(ve(:))),'Nonfinite trajectory/weights.');
y = complex(zeros(numel(selected),acq.nCh,'single'));
nextOut = 1; blockSize = 128;
for first = 1:blockSize:numel(rawLines)
    last = min(first+blockSize-1,numel(rawLines));
    lastOut = find(selected<=last*acq.N,1,'last');
    if isempty(lastOut) || lastOut<nextOut, continue; end
    block = image.unsorted(rawLines(first:last));
    block = reshape(permute(block,[2 1 3]),acq.nCh,[]);
    y(nextOut:lastOut,:) = block(:,selected(nextOut:lastOut)-(first-1)*acq.N).';
    nextOut = lastOut+1;
end
assert(nextOut==numel(selected)+1 && all(isfinite(y(:))),'Invalid raw-data output.');
fprintf('Low-resolution data: %s, grid %s\n',mat2str(size(y)),mat2str(gridSize));
end
