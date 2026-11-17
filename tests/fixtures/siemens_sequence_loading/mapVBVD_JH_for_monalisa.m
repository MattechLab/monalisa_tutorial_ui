function twix = mapVBVD_JH_for_monalisa(varargin) %#ok<INUSD>
% Test fixture: return synthetic TWIX metadata/data without a scanner file.
context = getappdata(0,'monalisa_sequence_loading_test');
twix = context.twix;
twix.image.unsorted = @read_test_samples;
end

function block = read_test_samples(indices)
context = getappdata(0,'monalisa_sequence_loading_test');
context.sampleReads = context.sampleReads+1;
setappdata(0,'monalisa_sequence_loading_test',context);
block = context.source(:,:,indices);
end
