function y = read_siemens_blocks(image,acq,rawLines,blockSize)
%READ_SIEMENS_BLOCKS Read selected TWIX lines in original order, coil first.
% Output: [nCh,nSample,nRetainedLines]. Initial shots/SI are selected by the
% sequence-aware metadata loader before reading, without a full raw-data copy.
if nargin<4, blockSize = 128; end
validateattributes(blockSize,{'numeric'},{'scalar','integer','positive'});
y = complex(zeros(acq.nCh,acq.N,numel(rawLines),'single'));
for first = 1:blockSize:numel(rawLines)
    last = min(first+blockSize-1,numel(rawLines));
    block = image.unsorted(rawLines(first:last));
    block = reshape(block,acq.N,acq.nCh,last-first+1);
    assert(all(isfinite(block(:))),'TWIX block contains nonfinite raw samples.');
    y(:,:,first:last) = permute(block,[2 1 3]);
end
end
