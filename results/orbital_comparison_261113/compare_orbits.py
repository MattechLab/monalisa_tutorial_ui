"""Descriptive six-scan orbital comparison using existing verified exports.

Project standard: prior compare_orbits.py/review_lens.py readers and figures.
Raw and reconstruction inputs stay intact. Native coordinates, no registration.
"""
from pathlib import Path
import gzip, hashlib, json, re, struct
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path('/Volumes/Extreme SSD/261113/recon_results/isotropic_scanner/run_seq702_705to709_20261113_155007_367')
OUT=Path(__file__).resolve().parent
IDS=[702,706,707,705,708,709]
coord=np.arange(-120,120)
idx=[np.flatnonzero((coord>=lo)&(coord<=hi)) for lo,hi in [(-65,65),(-110,-35),(-50,15)]]
xyz=[coord[i] for i in idx]
reports={sid:json.loads((ROOT/f'seq{sid}'/'reconstruction.json').read_text()) for sid in IDS}
cache=OUT/'orbital_crops.npz'
if not cache.exists():
    crops={}
    for sid,r in reports.items():
        path=ROOT/f'seq{sid}'/'xrms_normalized.nii.gz'
        assert hashlib.md5(path.read_bytes()).hexdigest()==r['niftiMD5']
        with gzip.open(path,'rb') as f:
            h=f.read(352)
            assert struct.unpack_from('<4h',h,40)==(3,240,240,240)
            assert struct.unpack_from('<h',h,70)[0]==16
            for offset,expected in [(280,(1.,0.,0.,-120.)),(296,(0.,1.,0.,-120.)),(312,(0.,0.,1.,-120.))]:
                assert struct.unpack_from('<4f',h,offset)==expected
            v=np.frombuffer(f.read(),dtype='<f4').reshape((240,240,240),order='F')
        crops[str(sid)]=v[np.ix_(*idx)].copy()*r['normalizationDivisor']
    np.savez_compressed(cache,**crops)
with np.load(cache) as f:data={sid:f[str(sid)] for sid in IDS}

header_pattern=re.compile(r'^(sGRADSPEC\.asGPAData\[0\]\.lOffset[XYZ]|sGRADSPEC\.alShimCurrent\[\d+\]|sTXSPEC\.asNucleusInfo\[0\]\.(?:lFrequency|flReferenceAmplitude|flAmplitudeCorrection)|sTXSPEC\.aRFPULSE\[0\]\.flAmplitude|sAdjData\.uiAdjShimMode)\s*=\s*([^\r\n]+)',re.M)
headers={}
for sid,r in reports.items():
    with Path(r['metadata']['rawFile']).open('rb') as f:text=f.read(1500000).decode('latin1')
    headers[sid]={k:float(v.strip()) for k,v in header_pattern.findall(text)}
common_keys=['matrixSize','voxelSizeMm','fovMm','kernel','kernelWidth','kernelParameters','nShotOff','retainedShots','imagingReadouts','physicalCoils','loaderMD5','trajectoryFunctionMD5','griddingMEXMD5']
for sid in IDS:
    for k in common_keys:assert reports[sid]['metadata'][k]==reports[702]['metadata'][k],(sid,k)
(OUT/'acquisition_audit.json').write_text(json.dumps({'processing_keys_matched':common_keys,'raw_headers':headers,
    'scans':[{'sequence':sid,'TR_ms':r['metadata']['sequenceParams']['header']['TR_s']*1000,
              'TE_ms':r['metadata']['sequenceParams']['header']['TE_s']*1000,
              'FA_deg':r['metadata']['reportedFlipAngleDeg'],
              'duration_s':r['metadata']['sequenceParams']['header']['TotalDuration'],
              'source_nifti_MD5':r['niftiMD5']} for sid,r in reports.items()]},indent=2)+'\n')

def label(sid):
    h=reports[sid]['metadata']['sequenceParams']['header']
    return f"seq{sid} | FA {reports[sid]['metadata']['reportedFlipAngleDeg']}°\nTR/TE {h['TR_s']*1000:g}/{h['TE_s']*1000:g} ms"

def panel(sid,axis,pos,limits=None,scale=1e-8):
    i=np.flatnonzero(xyz[axis]==pos);assert len(i)==1
    a=np.take(data[sid],i[0],axis=axis)
    remaining=[k for k in range(3) if k!=axis]
    u,v=[xyz[k] for k in remaining]
    if limits:
        keep=[(q>=limits[k][0])&(q<=limits[k][1]) for q,k in zip([u,v],remaining)]
        a=a[np.ix_(*keep)];u=u[keep[0]];v=v[keep[1]]
    return a.T/scale,[u[0]-.5,u[-1]+.5,v[0]-.5,v[-1]+.5]

def atlas(axis,positions,filename,limits=None,scales=None,vmax=5.):
    fig,axes=plt.subplots(len(positions),6,figsize=(18,2.3*len(positions)),layout='constrained')
    for row,p in enumerate(positions):
        for col,sid in enumerate(IDS):
            a,extent=panel(sid,axis,p,limits,1e-8 if scales is None else scales[sid])
            ax=axes[row,col]
            ax.imshow(a,origin='lower',extent=extent,cmap='gray',vmin=0,vmax=vmax,interpolation='nearest')
            if axis==2:ax.invert_yaxis()
            if row==0:ax.set_title(label(sid),fontsize=10)
            if col==0:ax.set_ylabel(f'Native {"xyz"[axis]}={p} mm',fontsize=9)
            ax.tick_params(labelsize=7)
    fig.suptitle(('Common restored magnitude window 0–5e-8' if scales is None else 'Vitreous-matched display; identical window across scans')+
                 '\nNative orbital slices; no interpolation or registration',fontsize=13)
    fig.savefig(OUT/filename,dpi=145);plt.close(fig)

atlas(2,[-36,-32,-28,-24,-20,-16,-12,-8],'axial_raw_atlas.png',limits={0:(-60,60),1:(-105,-40),2:(-50,15)})
atlas(0,[28,32,36,40,44],'positive_sagittal_raw_atlas.png',limits={1:(-100,-45),2:(-43,5)})
atlas(0,[-44,-40,-36,-32,-28],'negative_sagittal_raw_atlas.png',limits={1:(-105,-45),2:(-43,5)})
atlas(1,[-94,-90,-86,-82,-78],'anterior_coronal_raw_atlas.png',limits={0:(-60,60),2:(-43,5)})
region_file=OUT/'region_comparison.json'
if region_file.exists():
    metrics=json.loads(region_file.read_text())['metrics']
    for eye,positions in [('positive_x',[28,32,36,40,44]),('negative_x',[-44,-40,-36,-32,-28])]:
        scales={r['sequence']:r['vitreous_median'] for r in metrics if r['eye']==eye}
        atlas(0,positions,eye+'_sagittal_vitreous_matched_atlas.png',
              limits={1:(-100,-45),2:(-43,5)},scales=scales,vmax=1.8)
print('Saved source-verified orbital crop cache, header audit and native raw-scale atlases.')
