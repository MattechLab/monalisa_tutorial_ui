% Single-pulse frequency response of the actual exported seq700 RF events.
% Uses Pulseq's native mr.simRf. No relaxation, RF train, measured B0/B1 map,
% or scanner validation is included. Frequency is a model axis, not a B0 map.
repoRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(fileparts(repoRoot),'pulseq','matlab'));
outDir=fileparts(mfilename('fullpath'));
dataset='/Volumes/Extreme SSD/261110';
for kind={'gre','libre'}
    name=kind{1}; files=dir(fullfile(dataset,['yj_seq700_t1w_' name '_*.seq']));
    assert(isscalar(files));
    seqFile=fullfile(files.folder,files.name);
    seq=mr.Sequence();seq.read(seqFile);
    first=find(cellfun(@(events)events(2)>0,seq.blockEvents),1);
    assert(~isempty(first));block=seq.getBlock(first);rf=block.rf;
    [Mz,Mxy,F]=mr.simRf(rf,0,0);
    assert(all(isfinite(Mz)) && all(isfinite(Mxy)) && all(isfinite(F)));
    assert(max(abs(abs(Mxy).^2+Mz.^2-1))<1e-8);
    response=table(F(:),abs(Mxy(:)),Mz(:),'VariableNames',{'frequency_Hz','Mxy_magnitude','Mz'});
    writetable(response,fullfile(outDir,[name '_single_pulse_response.csv']));
    save(fullfile(outDir,[name '_exported_rf.mat']),'rf','seqFile','first','-v7');
    fprintf('%s: RF duration %.3f ms; peak %.6g Hz; response %d frequency points\n', ...
        name,rf.shape_dur*1000,max(abs(rf.signal)),numel(F));
end
