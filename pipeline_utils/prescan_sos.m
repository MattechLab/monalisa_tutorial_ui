function image = prescan_sos(y,t,ve,acq,gridSize)
%PRESCAN_SOS Original orientation-preview normalization, in coil order.
image = zeros(gridSize,'single');
for first = 1:4:acq.nCh
    last = min(first+3,acq.nCh);
    x = bmMathilda(y(:,first:last),t,ve,[],gridSize,gridSize,[1 1 1]./acq.FoV);
    for channel = 1:size(x,4), image = image+abs(x(:,:,:,channel)).^2; end
end
image = log1p(sqrt(image));
lo = prctile(image(:),1); hi = prctile(image(:),99);
assert(hi>lo,'The orientation preview has no intensity range.');
image = max(min((image-lo)/(hi-lo),1),0);
end
