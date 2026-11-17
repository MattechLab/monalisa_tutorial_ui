"""Descriptive phantom comparison, using the established NIfTI/JSON contract.

Restore original xRMS scale, use identical native planes, and compare each
sequence with its own normal-shimming acquisition. No voxelwise inference,
registration, B0-map estimation, or source modification is performed.
"""
from pathlib import Path
import gzip, hashlib, json, re, struct
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
DATASET=Path('/Volumes/Extreme SSD/261110')
ROOT=DATASET/'recon_results'/'isotropic_scanner'
BATCHES={'normal':ROOT/'run_seq700_norm_shimming_20261112_094850_901',
         'linear':ROOT/'run_seq700_linear_shimming_20261112_100256_121'}
IDS=[('gre','normal'),('gre','linear'),('libre','normal'),('libre','linear')]
coord=np.arange(-120,120)
volumes={};reports={};headers={}
header_pattern=re.compile(r'^(sGRADSPEC\.asGPAData\[0\]\.lOffset[XYZ]|sGRADSPEC\.alShimCurrent\[\d+\]|sTXSPEC\.asNucleusInfo\[0\]\.(?:lFrequency|flReferenceAmplitude|flAmplitudeCorrection)|sTXSPEC\.aRFPULSE\[0\]\.flAmplitude|sAdjData\.uiAdjShimMode)\s*=\s*([^\r\n]+)',re.M)
for kind,shim in IDS:
    suffix='norm' if shim=='normal' else 'linear'
    folder=BATCHES[shim]/f'seq700_{kind}_{suffix}_shimming'
    report=json.loads((folder/'reconstruction.json').read_text());reports[kind,shim]=report
    path=folder/'xrms_normalized.nii.gz'
    assert hashlib.md5(path.read_bytes()).hexdigest()==report['niftiMD5']
    with gzip.open(path,'rb') as stream:
        h=stream.read(352)
        assert struct.unpack_from('<4h',h,40)==(3,240,240,240)
        assert struct.unpack_from('<h',h,70)[0]==16
        for offset,expected in [(280,(1.,0.,0.,-120.)),(296,(0.,1.,0.,-120.)),(312,(0.,0.,1.,-120.))]:
            assert struct.unpack_from('<4f',h,offset)==expected
        v=np.frombuffer(stream.read(),dtype='<f4').reshape((240,240,240),order='F').copy()
    v*=report['normalizationDivisor'];volumes[kind,shim]=v
    raw=Path(report['metadata']['rawFile'])
    with raw.open('rb') as stream:text=stream.read(1500000).decode('latin1')
    headers[f'{kind}_{shim}']={k:float(value.strip()) for k,value in header_pattern.findall(text)}
    assert report['metadata']['rawStat']['bytes']==raw.stat().st_size

keys=['matrixSize','voxelSizeMm','fovMm','kernel','kernelWidth','kernelParameters','nShotOff','retainedShots','imagingReadouts','physicalCoils','samplesUsed','sequenceMD5','loaderMD5','trajectoryFunctionMD5','griddingMEXMD5']
for kind in ['gre','libre']:
    for key in keys:assert reports[kind,'normal']['metadata'][key]==reports[kind,'linear']['metadata'][key],(kind,key)
    for shim in ['normal','linear']:
        assert reports[kind,shim]['metadata']['nShot']==199

inside=np.flatnonzero((coord>=-90)&(coord<=90))
extent=[-90.5,90.5,-90.5,90.5]
fig,axes=plt.subplots(2,4,figsize=(15.5,8),layout='constrained')
slice_vmax=np.percentile(np.concatenate([v[30:211,30:211,120].ravel() for v in volumes.values()]),99.8)
for row,(axis,label) in enumerate([(2,'Axial z=0 mm'),(1,'Coronal y=0 mm')]):
    for col,key in enumerate(IDS):
        a=np.take(volumes[key],120,axis=axis)[np.ix_(inside,inside)].T
        ax=axes[row,col];im=ax.imshow(a,origin='lower',extent=extent,cmap='gray',vmin=0,vmax=slice_vmax,interpolation='nearest')
        if row==0:ax.set_title(f'{key[0].upper()} | {key[1]} shimming',fontsize=11)
        if col==0:ax.set_ylabel(label+'\nSequence '+('y' if axis==2 else 'z')+' (mm)')
        ax.set_xlabel('Sequence x (mm)');ax.tick_params(labelsize=8)
fig.colorbar(im,ax=axes,shrink=.7,pad=.015,label='Restored xRMS; same window in all panels')
fig.suptitle('Seq700 phantom: matched native slices on a shared reconstructed-magnitude scale\n240³, 1 mm; same exclusions, DCF and gridding settings; no registration',fontsize=13)
fig.savefig(OUT/'matched_slices.png',dpi=160);plt.close(fig)

# Ratios within each sequence reduce the common transmit/receive shading.
# They are not B0, flip-angle, or tissue-property maps.
xx,yy=np.meshgrid(coord,coord,indexing='ij')
display_mask=(xx**2+yy**2)<=70**2
profile_x=np.arange(-60,61)
disk=(yy**2+xx**2)<=40**2
ratios={};profiles={}
fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained')
for row,kind in enumerate(['gre','libre']):
    denominator=volumes[kind,'normal']
    ratio=np.divide(volumes[kind,'linear'],denominator,out=np.full_like(denominator,np.nan),where=denominator>0)
    ratios[kind]=ratio
    profiles[kind]=np.array([np.median(ratio[x+120,:,:][disk]) for x in profile_x])
    for col,axis in enumerate([2,1]):
        a=np.take(ratio,120,axis=axis).T
        a=np.ma.masked_where(~display_mask.T,a)
        ax=axes[row,col]
        im=ax.imshow(a,origin='lower',extent=[-120.5,119.5,-120.5,119.5],vmin=0,vmax=1.5,cmap='viridis',interpolation='nearest')
        ax.set_xlim(-85,85);ax.set_ylim(-85,85)
        ax.set_title(f'{kind.upper()}: linear / normal\n'+('z=0 mm' if axis==2 else 'y=0 mm'),fontsize=10)
        ax.set_xlabel('Sequence x (mm)')
    axes[row,2].plot(profile_x,profiles[kind],color='#1662a1' if kind=='gre' else '#c34c22',lw=2)
    axes[row,2].axhline(1,color='gray',ls='--');axes[row,2].set_ylim(0,1.6)
    axes[row,2].set_title(f'{kind.upper()}: median linear / normal\nacross central yz disk, radius 40 mm',fontsize=10)
    axes[row,2].set_xlabel('Sequence x (mm)');axes[row,2].grid(alpha=.2)
fig.colorbar(im,ax=axes[:,:2],shrink=.65,pad=.01,label='Magnitude ratio; 1 = unchanged')
fig.suptitle('Shim effect within GRE and LIBRE: native-position ratios, without gain fitting',fontsize=13)
fig.savefig(OUT/'shim_effect_ratios.png',dpi=160);plt.close(fig)

band_x=int(profile_x[np.argmin(profiles['libre'])])
band_indices=np.flatnonzero(abs(coord-band_x)<=2)
summary={}
for kind in ['gre','libre']:
    band=ratios[kind][band_indices,:,:][:,disk]
    summary[kind]={'band_slab_median_linear_over_normal':float(np.median(band)),
                   'band_slab_P10_P90':np.percentile(band,[10,90]).tolist(),
                   'central_cylinder_median_linear_over_normal':float(np.median(ratios[kind][profile_x+120,:,:][:,disk]))}
profile_table=np.c_[profile_x,profiles['gre'],profiles['libre'],profiles['libre']/profiles['gre']]
np.savetxt(OUT/'ratio_profiles.csv',profile_table,delimiter=',',header='x_mm,GRE_linear_over_normal,LIBRE_linear_over_normal,double_ratio',comments='')
result={'status':'descriptive comparison; no independent repeat acquisitions or voxelwise inference',
        'band_slab_definition':'Data-selected minimum LIBRE median-ratio x, then x +/-2 mm and central yz disk radius 40 mm',
        'band_centre_x_mm':band_x,'processing_match_within_each_sequence':True,
        'raw_scale_restored':True,'registration':False,'headers':headers,'summary':summary,
        'no_B0_map':'Recorded shim offsets are scanner units, not a measured spatial frequency map.'}
(OUT/'comparison_summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'band_centre_x_mm':band_x,'summary':summary},indent=2))
