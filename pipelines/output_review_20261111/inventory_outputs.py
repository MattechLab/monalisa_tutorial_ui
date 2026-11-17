from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
import csv,hashlib,json,os

work=Path(__file__).resolve().parent
repo=Path('/Users/cag/forclone/pulseq4mreye_rectfov')
mono=Path('/Users/cag/forclone/monalisa_tutorial_ui')
workspace=Path('/Users/cag/.codex/.chatgpt-projects/g-p-6a82c561fb88819194c50cc5d389246c')
ssd=Path('/Volumes/Extreme SSD/261110')
roots={
 'koma_output':repo/'simulation/koma_water/output',
 'koma_checks':repo/'simulation/koma_water/checks',
 'koma_inputs':repo/'simulation/koma_water/input',
 'sequence_outputs':repo/'libre_3d_radial/output',
 'scanner_results':ssd/'recon_results',
 'monalisa_outputs':mono/'outputs',
 'monalisa_pointers':mono/'results',
 'workspace_outputs':workspace/'outputs',
 'visualizations':Path('/Users/cag/.codex/visualizations/2026/11/06/01a1b12a-d3db-7af1-8466-2a072fdb7fbe')}
for p in sorted(Path('/private/tmp').glob('kernel*'))+sorted(Path('/private/tmp').glob('orbital*'))+sorted(Path('/private/tmp').glob('waterfig*')):
 if p.is_dir():roots['temporary_'+p.name]=p

def group_for(name,rel):
 parts=rel.parts
 if name=='koma_output':
  if parts[0]=='seq712_water':return 'S2_earlier_seq712_run'
  if len(parts)>1 and parts[1]=='shape_phantom':return 'S3_superseded_separate_spheres'
  if len(parts)>1 and parts[1]=='fixed_shape_phantom':return 'S4_current_fixed_phantom_and_shots'
  return 'S1_selected_Koma_batch'
 if name=='scanner_results':
  if parts[0].startswith('box_phantom_kernel_curve'):return 'S5_box_phantom_kernel_curve'
  if parts[0].startswith('anisotropic_xrms_'):return 'R1_scanner_baseline_061'
  if len(parts)>2 and parts[1].startswith('seq') and '_sigma' in parts[1]:return 'R2_scanner_sigma'+parts[1].split('_sigma')[-1]
  if len(parts)>1 and parts[1].startswith('phantom_'):return 'S6_three_shape_phantom_kernel_sweep'
  if len(parts)>1 and (parts[1].startswith('test_phantom') or 'phantom' in parts[1]):return 'S6_three_shape_phantom_kernel_sweep'
  return 'R3_scanner_comparisons_and_scripts'
 if name=='workspace_outputs':
  top=parts[0]
  if top.startswith(('radial_investigation','trajectory_coverage','gaussian_comparison','window_comparison')):return 'S7_static_investigation_'+top
  if top.startswith('orbital_kernel_comparison'):return 'R4_local_orbital_figures'
  return 'X1_unrelated_workspace_outputs'
 if name=='sequence_outputs':
  if parts[0]=='anisotropic_ptp':return 'P1_current_sequence_exports'
  if parts[0]=='anisotropic_auto_sampling_validation':return 'P2_sampling_validation'
  if parts[0]=='rectfov':return 'P3_earlier_cone_ring_prototype'
  return 'X2_other_sequence_outputs'
 if name=='koma_checks':return 'S8_Koma_validation_checks'
 if name=='koma_inputs':return 'P0_frozen_Koma_sequence_inputs'
 if name=='monalisa_outputs':return 'X3_other_MONALISA_outputs'
 if name=='monalisa_pointers':return 'R5_result_path_pointers'
 if name=='visualizations':return 'S9_earlier_visualization'
 return 'T1_temporary_working_files'

files=[];errors=[]
for name,root in roots.items():
 if not root.exists():errors.append(f'Missing root: {root}');continue
 def walkerror(e):errors.append(str(e))
 for folder,dirs,names in os.walk(root,followlinks=False,onerror=walkerror):
  for filename in names:
   p=Path(folder)/filename
   try:
    st=p.lstat();rel=p.relative_to(root)
    files.append(dict(root=name,group=group_for(name,rel),path=str(p),relative_path=str(rel),
      logical_bytes=st.st_size,allocated_bytes=st.st_blocks*512,mtime_ns=st.st_mtime_ns,
      is_symlink=p.is_symlink(),extension=p.suffix.lower()))
   except OSError as e:errors.append(str(e))
# Source acquisitions are inventoried as protected inputs, without reading their contents.
for p in ssd.iterdir():
 if p.is_file() and p.suffix in ('.dat','.seq') and not p.name.startswith('._'):
  st=p.stat();files.append(dict(root='scanner_inputs',group='P0_protected_scanner_inputs',path=str(p),relative_path=p.name,
   logical_bytes=st.st_size,allocated_bytes=st.st_blocks*512,mtime_ns=st.st_mtime_ns,is_symlink=False,extension=p.suffix))
files.sort(key=lambda r:(r['group'],r['path']))
groups=defaultdict(lambda:dict(files=0,logical_bytes=0,allocated_bytes=0,paths=[]))
for r in files:
 g=groups[r['group']];g['files']+=1;g['logical_bytes']+=r['logical_bytes'];g['allocated_bytes']+=r['allocated_bytes'];g['paths'].append(r['path'])
# Exact duplicate search is restricted to small reports/scripts/sequence exports,
# not large MAT/NPZ arrays. Similar runs are not declared identical.
by_size=defaultdict(list)
for r in files:
 p=Path(r['path'])
 if not p.name.startswith('._') and 0<r['logical_bytes']<12*1024**2 and p.suffix.lower() in ('.png','.pdf','.json','.m','.jl','.py','.seq','.txt','.md'):
  by_size[r['logical_bytes']].append(r)
duplicates=[]
for size,rows in by_size.items():
 if len(rows)<2:continue
 matches=defaultdict(list)
 for r in rows:matches[hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()].append(r['path'])
 for sha,paths in matches.items():
  if len(paths)>1:duplicates.append(dict(sha256=sha,bytes_each=size,paths=paths))
with (work/'output_files.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(files[0]));w.writeheader();w.writerows(files)
manifest=dict(scanned_utc=datetime.now(timezone.utc).isoformat(),roots={k:str(v) for k,v in roots.items()},
 scope='Both experiment repositories, scanner result root, local output roots and task temporary folders. Other projects flagged X; scanner .dat/.seq protected.',
 deletion_performed=False,errors=errors,groups=dict(groups),exact_small_file_duplicates=duplicates,
 duplicate_check_limit='Only selected report/script/sequence file types below 12 MiB. No duplicate claim for acquisitions or large MAT/NPZ files.')
(work/'output_inventory.json').write_text(json.dumps(manifest,indent=2)+'\n')
lines=['# Output inventory — 11 November 2026','',
 'No existing output was moved, rewritten or deleted. Sizes below are measured, not estimates. Logical size is file content; allocated size is filesystem space.',
 '', '| Group | Files | Logical MiB | Allocated MiB |','|---|---:|---:|---:|']
for name,g in sorted(groups.items()):lines.append(f"| {name} | {g['files']} | {g['logical_bytes']/1024**2:.2f} | {g['allocated_bytes']/1024**2:.2f} |")
lines+=['','Complete absolute paths: [output_files.csv](output_files.csv). Group paths and exact small-file duplicate hashes: [output_inventory.json](output_inventory.json).',
 '', 'X groups belong to other experiments and are excluded from this cleanup. P0 groups are protected source inputs. No acquired scanner raw file is a discard candidate.',
 '', 'Duplicate checking covers only selected small reports/scripts/sequence exports, not full acquisitions or large MAT/NPZ arrays. The earlier seq712 run is superseded but is not claimed byte-identical to the selected batch.']
(work/'OUTPUT_INVENTORY.md').write_text('\n'.join(lines)+'\n')
for k,g in sorted(groups.items()):print(f"{k}: {g['files']} files, {g['logical_bytes']/1024**2:.2f} logical MiB, {g['allocated_bytes']/1024**2:.2f} allocated MiB")
print('Exact small-file duplicate sets:',len(duplicates),'Errors:',errors)
