from pathlib import Path
import csv,gzip,hashlib,json,struct
import h5py,numpy as np
stage=Path(__file__).resolve().parent
plan=json.loads((stage/'execution_plan.json').read_text());ret=Path(plan['retainedRoot'])
proof=json.loads((stage/'preservation_manifest.json').read_text())
updates=json.loads((stage/'updates_plan.json').read_text());changed={i['destination']:i['source'] for i in updates}
mono=Path('/Users/cag/forclone/monalisa_tutorial_ui');sim=Path('/Users/cag/forclone/pulseq4mreye_rectfov/simulation/koma_water')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for p in plan['sourceResultRoots']+plan['discardOtherRoots']:assert not Path(p).exists(),p
assert len(list(ret.glob('seq*/xrms_normalized.nii.gz')))==12
assert not list(ret.glob('seq*/reconstruction.mat'))
seen=set()
for rec in proof['nifti_validation']:
 p=Path(rec['nifti']);assert sha(p)==rec['nifti_sha256']
 m=json.loads((p.parent/'reconstruction.json').read_text())
 sigma=float(m['metadata']['kernelParameters'][0]);assert sigma in (.61,.32,.20,.10)
 seen.add((rec['sequence'],sigma))
 data=gzip.open(p,'rb').read();assert struct.unpack_from('<hh',data,252)==(0,2)
 assert np.isfinite(m['normalizationDivisor']) and m['normalizationDivisor']>0
assert len(seen)==12
for sid in (712,713,714):
 with h5py.File(ret/'encoding'/f'seq{sid}.mat') as f:
  assert {'t','ve','adcTimes','metadata','x_mm','y_mm','z_mm'}<=set(f)
  assert 'x0' not in f and 'xrms' not in f
for rec in proof['copies']:
 expected=sha(changed[rec['destination']]) if rec['destination'] in changed else rec['sha256']
 assert sha(rec['destination'])==expected,rec['destination']
with (mono/'pipelines/output_review_20261111/output_files.csv').open() as f:
 for r in csv.DictReader(f):
  if r['group'].startswith('P0_') or r['group'] in ('S1_selected_Koma_batch','S4_current_fixed_phantom_and_shots','S8_Koma_validation_checks'):
   p=Path(r['path']);assert p.exists(),p
   if p.name=='review_batch.m':continue # intentional source-path repair, documented
   st=p.stat();assert st.st_size==int(r['logical_bytes']) and st.st_mtime_ns==int(r['mtime_ns']),p
for name in ('run_directory.txt','kernel_sweep_directory.txt','extended_kernel_sweep_directory.txt','lower_kernel_sweep_directory.txt'):
 assert Path((mono/'results/scanner_anisotropic_261110'/name).read_text().strip())==ret
for p in (sim/'HANDOVER_KOMA_WATER_20261111.md',mono/'pipelines/anisotropic_scanner/HANDOVER_ANISOTROPIC_RECON_20261111.md'):
 assert p.is_file() and p.stat().st_size>3000
for n in ('test_phantom_sigma.m','test_phantom_extended.m'):
 s=(ret.parents[1]/'water_phantom_diagnostics/kernel_sweep_20261111'/n).read_text()
 assert str(ret/'encoding') in s and "sprintf('seq%d.mat',id)" in s
checkdir=Path('/private/tmp/koma_cleanup_final_check');assert not checkdir.exists()
result=dict(status='final_checks_passed',nifti_count=12,retained_widths=[.61,.32,.20,.10],
 raw_inputs_and_kept_simulation_arrays_unchanged=True,preserved_reports_and_diagnostics_sha256_valid=True,
 retired_paths_absent=True,full_scanner_MATs_absent=True,result_pointers_valid=True,
 both_handover_files_present=True,nominal_geometry_and_divisors_valid=True)
(stage/'final_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
