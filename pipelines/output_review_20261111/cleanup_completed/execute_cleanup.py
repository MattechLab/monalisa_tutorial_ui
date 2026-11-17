from pathlib import Path
import csv,hashlib,json,os,re,shutil
stage=Path(__file__).resolve().parent
plan=json.loads((stage/'execution_plan.json').read_text());preserve=json.loads((stage/'preservation_manifest.json').read_text())
assert json.loads((stage/'pipeline_validation.json').read_text())['twelve_nifti_loader_roundtrips']
oldInventory=Path('/Users/cag/forclone/monalisa_tutorial_ui/pipelines/output_review_20261111/output_files.csv')
with oldInventory.open() as f:snapshot={r['path']:r for r in csv.DictReader(f)}
updates=json.loads((stage/'updates_plan.json').read_text())
changedDest={i['destination']:i['source'] for i in updates}

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
 return h.hexdigest()

def tree_stats(root):
 stats=dict(logical_bytes=0,allocated_bytes=0,files=0)
 for p in Path(root).rglob('*'):
  if p.is_file():
   st=p.stat();stats['logical_bytes']+=st.st_size;stats['allocated_bytes']+=st.st_blocks*512;stats['files']+=1
 return stats

# Recheck all verified replacements and protected input stats before deletion.
for rec in preserve['nifti_validation']:
 assert sha(rec['nifti'])==rec['nifti_sha256']
 assert sha(rec['source'])==rec['source_mat_sha256']
for rec in preserve['copies']:
 expected=sha(changedDest[rec['destination']]) if rec['destination'] in changedDest else rec['sha256']
 assert sha(rec['destination'])==expected,rec['destination']
for path,r in snapshot.items():
 if r['group'].startswith('P0_'):
  p=Path(path);assert p.exists() and p.stat().st_size==int(r['logical_bytes']) and p.stat().st_mtime_ns==int(r['mtime_ns']),path
base,sweep=map(Path,plan['sourceResultRoots'])
assert {p.name for p in base.iterdir() if p.is_dir()}=={f'seq{i}' for i in (712,713,714)}
assert {p.name for p in sweep.iterdir() if p.is_dir()}=={f'seq{i}_sigma{s}' for i in (712,713,714) for s in ('050','040','032','027','020','015','010')}
roots=[base,sweep,sweep.parent/'box_phantom_kernel_curve_20261111_105211_216']+list(map(Path,plan['discardOtherRoots']))
expectedRoots={
 '/Volumes/Extreme SSD/261110/recon_results/anisotropic_xrms_20261111_102743_031',
 '/Volumes/Extreme SSD/261110/recon_results/anisotropic_kernel_sweep_20261111_103727_973',
 '/Volumes/Extreme SSD/261110/recon_results/box_phantom_kernel_curve_20261111_105211_216',
 '/Users/cag/forclone/pulseq4mreye_rectfov/simulation/koma_water/output/seq712_water',
 '/Users/cag/forclone/pulseq4mreye_rectfov/simulation/koma_water/output/seq712_714/shape_phantom',
 '/Users/cag/forclone/pulseq4mreye_rectfov/libre_3d_radial/output/rectfov'}
assert {str(p) for p in roots}==expectedRoots
selected=Path('/Users/cag/forclone/pulseq4mreye_rectfov/simulation/koma_water/output/seq712_714')
for id in (712,713,714):
 assert (selected/f'seq{id}_water/acquisition.mat').is_file()
 assert (selected/f'seq{id}_water/recon_xrms/reconstruction.mat').is_file()
assert (selected/'fixed_shape_phantom/results.json').is_file()
records=[]
for root in roots:
 assert root.is_dir() and not root.is_symlink()
 for p in root.rglob('*'):
  assert not p.is_symlink(),p
  if not p.is_file():continue
  st=p.stat();old=snapshot.get(str(p))
  if not p.name.startswith('._') and p.name!='.DS_Store':
   assert old is not None,f'Unreviewed file: {p}'
   assert st.st_size==int(old['logical_bytes']) and st.st_mtime_ns==int(old['mtime_ns']),f'Changed file: {p}'
  records.append(dict(path=str(p),bytes=st.st_size,allocated_bytes=st.st_blocks*512,mtime_ns=st.st_mtime_ns))
# All planned retained source MATs are inside the reviewed roots, replaced by NIfTI.
assert sum(Path(r['path']).name=='reconstruction.mat' and any(Path(r['path']).is_relative_to(q) for q in (base,sweep)) for r in records)==24
beforeExternal={k:sum(int(r[k]) for r in snapshot.values() if r['root']=='scanner_results') for k in ('logical_bytes','allocated_bytes')}
(stage/'authorized_deletion_manifest.json').write_text(json.dumps(dict(authorization=plan['authorized'],roots=[str(p) for p in roots],files=records,retained_nifti_verified=12,protected_inputs_verified=True),indent=2)+'\n')
for root in roots:
 shutil.rmtree(root)
 sibling=root.parent/('._'+root.name)
 if sibling.is_file():sibling.unlink()
 assert not root.exists();print('CLEANED:',root,flush=True)
# Preserve current diagnostics, sources and all compact outputs.
for rec in preserve['nifti_validation']:assert Path(rec['nifti']).is_file()
for rec in preserve['copies']:assert Path(rec['destination']).is_file()
for id in (712,713,714):assert (selected/f'seq{id}_water/acquisition.mat').is_file()
report=dict(status='cleanup_complete',removed_roots=[str(p) for p in roots],
 full_scanner_reconstruction_MATs_removed=24,retained_nifti_volumes=12,
 retained_sigmas=[.61,.32,.20,.10],removed_sigmas=[.50,.40,.27,.15],
 removed_file_count=len(records),removed_content_bytes=sum(r['bytes'] for r in records),
 original_external_results=beforeExternal,current_external_results=tree_stats(sweep.parent),
 compact_scanner_batch=tree_stats(plan['retainedRoot']),
 normalized_nifti_bytes=sum(Path(r['nifti']).stat().st_size for r in preserve['nifti_validation']),
 protected_scanner_and_Koma_inputs_unchanged=True,selected_Koma_batch_and_fixed_phantom_preserved=True)
(stage/'cleanup_result.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
