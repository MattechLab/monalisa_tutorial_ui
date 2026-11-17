"""Read-only native-volume orbital review; exploratory anatomical measurements.

Follows the existing compare_orbits.py NIfTI reader and plotting convention.
Sources are never rewritten. No registration or tissue-model fit is applied.
"""
from pathlib import Path
import gzip
import json
import struct
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('/Volumes/Extreme SSD/261110/recon_results/isotropic_scanner/run_20261111_174401_195')
OUT = Path(__file__).resolve().parent
IDS = [701, 702, 704, 703, 705]  # Increasing FA within the common-TR series.
FA_FACTOR = 2  # User correction: report 12/14/16/20 deg for sequence 6/7/8/10.
coord = np.arange(-120, 120)
indices = [np.flatnonzero((coord >= low) & (coord <= high))
           for low, high in [(-60, 60), (-98, -44), (-43, 3)]]
xyz = [coord[i] for i in indices]
cache = OUT / 'orbital_crops.npz'
reports = {sid: json.loads((ROOT/f'seq{sid}'/'reconstruction.json').read_text()) for sid in IDS}
if not cache.exists():
    crops = {}
    for sid in IDS:
        with gzip.open(ROOT/f'seq{sid}'/'xrms_normalized.nii.gz', 'rb') as stream:
            header = stream.read(352)
            assert struct.unpack_from('<4h', header, 40) == (3,240,240,240)
            assert struct.unpack_from('<h', header, 70)[0] == 16
            assert struct.unpack_from('<f', header, 108)[0] == 352
            for offset, expected in [(280,(1,0,0,-120)),(296,(0,1,0,-120)),(312,(0,0,1,-120))]:
                assert struct.unpack_from('<4f', header, offset) == expected
            volume = np.frombuffer(stream.read(), dtype='<f4').reshape((240,240,240), order='F')
        crops[str(sid)] = volume[np.ix_(*indices)].copy() * reports[sid]['normalizationDivisor']
    np.savez_compressed(cache, **crops)
with np.load(cache) as loaded:
    data = {sid: loaded[str(sid)] for sid in IDS}

X,Y,Z = np.meshgrid(*xyz, indexing='ij')
eye_centres = {'positive_x': (36,-70,-21), 'negative_x': (-32,-73,-21)}
vitreous_masks = {name: ((X-c[0])**2+(Y-c[1])**2+(Z-c[2])**2 <= 7**2)
                 for name,c in eye_centres.items()}
scales = {(sid,name): float(np.median(data[sid][mask]))
          for name,mask in vitreous_masks.items() for sid in IDS}

def title(sid):
    h = reports[sid]['metadata']['sequenceParams']['header']
    return f"seq{sid} | {FA_FACTOR*h['rf_flip_deg']}°\nTR {h['TR_s']*1000:g} ms"

def slice_data(sid, axis, pos, eye=None):
    index = np.flatnonzero(xyz[axis] == pos)
    assert len(index) == 1
    a = np.take(data[sid], index[0], axis=axis)
    remaining = [k for k in range(3) if k != axis]
    u,v = [xyz[k] for k in remaining]
    if eye:
        cx = eye_centres[eye][0]
        limits = {0:(cx-17,cx+17),1:(-94,-51),2:(-40,-5)}
        masks = [(q>=limits[k][0]) & (q<=limits[k][1]) for q,k in zip([u,v],remaining)]
        a = a[np.ix_(*masks)]; u=u[masks[0]]; v=v[masks[1]]
        a = a/scales[sid,eye]
    else:
        a = a/1e-8
    return a.T, [u[0]-.5,u[-1]+.5,v[0]-.5,v[-1]+.5], remaining

def atlas(axis, positions, eye, filename):
    fig,axes=plt.subplots(len(positions),5,figsize=(15,2.3*len(positions)),layout='constrained')
    for row,pos in enumerate(positions):
        for col,sid in enumerate(IDS):
            a,extent,remaining = slice_data(sid,axis,pos,eye)
            ax=axes[row,col]
            ax.imshow(a,origin='lower',extent=extent,cmap='gray',vmin=0,vmax=1.7 if eye else 4,interpolation='nearest')
            if axis==2: ax.invert_yaxis()
            if row==0: ax.set_title(title(sid),fontsize=10)
            if col==0: ax.set_ylabel(f'Native {"xyz"[axis]}={pos} mm',fontsize=9)
            ax.tick_params(labelsize=7)
    fig.suptitle(f'{eye or "Both eyes"}: native orbital slices; '+
                 ('each eye scaled by its own posterior-vitreous median' if eye else 'shared xRMS scale 0–4e-8')+
                 '\nNo resampling or registration; columns reordered 701, 702, 704, 703, 705',fontsize=12)
    fig.savefig(OUT/filename,dpi=145); plt.close(fig)

atlas(2,[-30,-27,-24,-21,-18,-15,-12],None,'axial_orbit_atlas.png')
atlas(0,[28,31,34,37,40,43],'positive_x','positive_eye_sagittal_atlas.png')
atlas(0,[-40,-37,-34,-31,-28,-25],'negative_x','negative_eye_sagittal_atlas.png')
atlas(1,[-89,-86,-83,-80,-77],None,'anterior_coronal_atlas.png')

metrics=[]
for name,mask in vitreous_masks.items():
    for sid in IDS:
        a=data[sid][mask].astype(float); median=np.median(a)
        smooth=gaussian_filter(data[sid].astype(float),sigma=2,mode='nearest')[mask]
        residual=a-smooth
        metrics.append({'eye':name,'sequence':sid,'mask_centre_mm':eye_centres[name],
                        'mask_radius_mm':7,'n_voxels':int(mask.sum()),
                        'median_xrms':float(median),'IQR_over_median':float(np.subtract(*np.percentile(a,[75,25]))/median),
                        'P90_minus_P10_over_median':float(np.subtract(*np.percentile(a,[90,10]))/median),
                        'smooth_P90_minus_P10_over_median':float(np.subtract(*np.percentile(smooth,[90,10]))/median),
                        'residual_MAD_over_median':float(np.median(abs(residual-np.median(residual)))/median)})
(OUT/'posterior_vitreous_variation.json').write_text(json.dumps({'interpretation':
    'Exploratory spatial variation; not calibrated noise, SNR, or intrinsic tissue heterogeneity. Fixed 7-mm-radius interior spheres; no motion registration.',
    'measurements':metrics},indent=2)+'\n')
print('Saved full orbital atlases and exploratory posterior-vitreous variation metrics.')
