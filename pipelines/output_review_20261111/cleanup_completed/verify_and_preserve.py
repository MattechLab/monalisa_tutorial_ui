from pathlib import Path
import csv,gzip,hashlib,json,shutil,struct
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
stage=Path(__file__).resolve().parent
plan=json.loads((stage/'execution_plan.json').read_text())
retained=Path(plan['retainedRoot']);base,sweep=map(Path,plan['sourceResultRoots'])
oldCurve=sweep.parent/'box_phantom_kernel_curve_20261111_105211_216'
staticRoot=sweep.parent/'water_phantom_diagnostics'
comparisons=retained/'comparisons';history=comparisons/'historical_scripts'
static=staticRoot/'kernel_sweep_20261111';curve=staticRoot/'box_kernel_curve_20261111'
for p in (comparisons,history,static,curve):p.mkdir(parents=True,exist_ok=True)

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()

def read_known_nifti(p):
 # Independent verification of the native MATLAB writer's NIfTI-1 export.
 data=gzip.open(p,'rb').read();assert struct.unpack_from('<i',data,0)[0]==348
 dim=struct.unpack_from('<8h',data,40);assert dim[0]==3 and tuple(dim[1:4])==(120,120,90)
 assert struct.unpack_from('<h',data,70)[0]==16 and struct.unpack_from('<h',data,72)[0]==32
 assert struct.unpack_from('<hh',data,252)==(0,2)
 assert data[123]&7==2
 offset=round(struct.unpack_from('<f',data,108)[0])
 a=np.frombuffer(data,dtype='<f4',offset=offset,count=120*120*90).reshape((120,120,90),order='F').copy()
 return a

validation=[];volumes={}
for job in plan['jobs']:
 source=Path(job['source']);dest=Path(job['destination'])
 meta=json.loads((dest/'reconstruction.json').read_text())
 a=read_known_nifti(dest/'xrms_normalized.nii.gz')
 with h5py.File(source,'r') as f:
  original=f['xrms'][:].transpose()
  expected=(original.astype(float)/meta['normalizationDivisor']).astype(np.float32)
  assert np.array_equal(a,expected)
  coords=[f[k][:].ravel() for k in ('x_mm','y_mm','z_mm')]
  for k,c in enumerate(coords):
   assert np.array_equal(c,meta['originMm'][k]+np.arange(len(c))*meta['spacingMm'][k])
  np.testing.assert_allclose(f['metadata/kernelParameters'][:].ravel(),[job['sigma'],10],rtol=0,atol=1e-14)
 assert np.isfinite(a).all() and a.min()>=0 and a.max()==1
 assert meta['relativeRawIntensityRestorationError']<1e-7
 assert hashlib.md5((dest/'xrms_normalized.nii.gz').read_bytes()).hexdigest()==meta['niftiMD5']
 volumes[job['sequence'],job['sigma']]=a.astype(float)*meta['normalizationDivisor']
 validation.append(dict(sequence=job['sequence'],sigma=job['sigma'],source=str(source),
  nifti=str(dest/'xrms_normalized.nii.gz'),exact_normalized_values=True,
  coordinate_vectors_exact=True,nifti_sha256=sha(dest/'xrms_normalized.nii.gz'),source_mat_sha256=sha(source)))
 print('Independent NIfTI check:',job['sequence'],job['sigma'],flush=True)

copies=[]
def copy_verified(src,dst):
 assert not dst.exists(),f'Destination exists: {dst}'
 dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
 checksum=sha(src);assert checksum==sha(dst)
 copies.append(dict(source=str(src),destination=str(dst),sha256=checksum,bytes=src.stat().st_size))
# Preserve all small sweep-level reports/settings/scripts, while separating analytic phantom work.
for src in sweep.iterdir():
 if src.is_dir() or src.name.startswith('._') or src.name=='.DS_Store':continue
 n=src.name
 isstatic=n.startswith(('phantom','plot_phantom','test_phantom'))
 dst=(static if isstatic else history if src.suffix in ('.py','.m','.mat') else comparisons)/n
 copy_verified(src,dst)
for src in oldCurve.iterdir():
 if src.is_file() and not src.name.startswith('._') and src.name!='.DS_Store':copy_verified(src,curve/src.name)
for src in base.iterdir():
 if src.is_file() and not src.name.startswith('._') and src.name!='.DS_Store':copy_verified(src,history/('baseline_'+src.name))
# Current four-width figures restore the original scale before applying one shared display range.
x,y,z=coords
xi=abs(x)<=55;yi=(y>=-86)&(y<=-46);planes=(-24,-20);sigmas=(.61,.32,.20,.10)
reference=np.concatenate([volumes[sid,.20][:,:,np.argmin(abs(z-zz))][np.ix_(xi,yi)].ravel() for sid in (712,713,714) for zz in planes])
vmax=float(np.percentile(reference,99.5))
for zz in planes:
 fig,axs=plt.subplots(3,4,figsize=(15,7));fig.subplots_adjust(left=.12,right=.91,bottom=.10,top=.80,hspace=.35,wspace=.12)
 for row,sid in enumerate((712,713,714)):
  for col,sigma in enumerate(sigmas):
   ax=axs[row,col];a=volumes[sid,sigma][:,:,np.argmin(abs(z-zz))][np.ix_(xi,yi)].T
   im=ax.imshow(a,origin='lower',extent=(x[xi][0]-1,x[xi][-1]+1,y[yi][0]-1,y[yi][-1]+1),cmap='gray',norm=Normalize(0,vmax),interpolation='nearest')
   if row==0:ax.set_title(f'sigma = {sigma:.2f}',fontsize=14)
   if col==0:ax.set_ylabel(f'seq{sid}\ny (mm)',fontsize=12)
   else:ax.set_yticklabels([])
   if row==2:ax.set_xlabel('Sequence x (mm)')
   else:ax.set_xticklabels([])
 cax=fig.add_axes([.94,.21,.012,.48]);fig.colorbar(im,cax=cax,label='xRMS on restored original intensity scale')
 fig.suptitle(f'Retained orbital reconstructions: seq712–714 | z = {zz} mm',fontsize=18,y=.97)
 fig.text(.5,.89,'0.61 baseline · 0.32 · 0.20 · 0.10 | same crop and intensity scale across all panels and both planes',ha='center')
 for suffix in ('png','pdf'):fig.savefig(comparisons/f'orbital_retained_zm{abs(zz):02d}.{suffix}',dpi=210)
 plt.close(fig)
(stage/'preservation_manifest.json').write_text(json.dumps(dict(nifti_validation=validation,copies=copies,
 staticRoot=str(staticRoot),staticSweep=str(static),boxCurve=str(curve),comparisons=str(comparisons),
 copied_files_all_sha256_verified=True,all_12_nifti_independently_verified=True),indent=2)+'\n')
print('All NIfTI arrays verified; reports and static diagnostics preserved with SHA-256 checks.',flush=True)
