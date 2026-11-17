repo='/Users/cag/forclone/monalisa_tutorial_ui';addpath(repo);addpath(fullfile(repo,'pipeline_utils'));
sim='/Users/cag/forclone/pulseq4mreye_rectfov/simulation/koma_water';addpath(sim);
plan=jsondecode(fileread('/Users/cag/.codex/.chatgpt-projects/g-p-6a82c561fb88819194c50cc5d389246c/work/pipeline_cleanup_execute_20261111/execution_plan.json'));
newRoot=[tempname '_readonly'];
p=run_anisotropic_scanner('/Volumes/Extreme SSD/261110',[],[],newRoot,false,true);
assert(isequal(p.kernelSigmas,[.61 .32 .20 .10]) && numel(p.jobs)==12 && ~isfolder(newRoot));
for j=1:numel(plan.jobs)
 job=plan.jobs(j);r=load_anisotropic_nifti(job.destination);
 old=load(job.source,'xrms','x_mm','y_mm','z_mm');
 err=norm(r.xrms(:)-double(old.xrms(:)))/norm(double(old.xrms(:)));
 assert(err<1e-7);
 assert(isequal(r.x_mm,old.x_mm(:).') && isequal(r.y_mm,old.y_mm(:).') && isequal(r.z_mm,old.z_mm(:).'));
 n=load_anisotropic_nifti(job.destination,false);assert(isa(n.xrms,'single') && max(n.xrms(:))==1);
end
for id=712:714
 e=load(fullfile(plan.retainedRoot,'encoding',sprintf('seq%d.mat',id)));
 assert(~isfield(e,'xrms') && ~isfield(e,'x0') && isfield(e,'t') && isfield(e,'ve'));
 assert(e.metadata.imagingReadouts==7697);
end
assert(contains(fileread(fullfile(sim,'inspect_water_acquisition.m')),'''output'',''seq712_714'',''seq712_water'''));
assert(contains(fileread(fullfile(sim,'reconstruct_water_xrms.m')),'''output'',''seq712_714'',''seq712_water'''));
for file={fullfile(repo,'run_anisotropic_scanner.m'),fullfile(repo,'pipeline_utils/export_anisotropic_nifti.m'),fullfile(repo,'pipeline_utils/load_anisotropic_nifti.m')}
 findings=checkcode(file{1},'-id');
 for k=1:numel(findings),fprintf('ANALYZER %s:%d %s %s\n',file{1},findings(k).line,findings(k).id,findings(k).message);end
 assert(isempty(findings),'Code analyzer findings need review.');
end
report=struct('default_four_width_plan',true,'twelve_nifti_loader_roundtrips',true, ...
 'normalized_float32_and_reversible_scale',true,'shared_encodings_valid',true, ...
 'koma_default_points_to_retained_batch',true,'matlab_code_analysis',true, ...
 'full_new_simulation_or_reconstruction_run',false);
fid=fopen('/Users/cag/.codex/.chatgpt-projects/g-p-6a82c561fb88819194c50cc5d389246c/work/pipeline_cleanup_execute_20261111/pipeline_validation.json','w');
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));fclose(fid);
disp('UPDATED PIPELINES AND ALL TWELVE NIFTI LOADERS PASSED');
