"""Exploratory native-grid region comparison; anatomical labels remain provisional.

No segmentation, registration, inferred SNR, T1/T2* fitting or voxelwise tests.
"""
from pathlib import Path
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d, gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
IDS=[702,706,707,705,708,709]
xyz=[np.arange(-65,66),np.arange(-110,-34),np.arange(-50,16)]
X,Y,Z=np.meshgrid(*xyz,indexing='ij')
with np.load(OUT/'orbital_crops.npz') as f:data={sid:f[str(sid)].astype(float) for sid in IDS}
audit=json.loads((OUT/'acquisition_audit.json').read_text())
settings={r['sequence']:r for r in audit['scans']}
regions={
 'positive_x':{'vitreous':(37,-66,-18),'dark_crescent':(39,-73,-26),'bright_cap':(37,-77,-13)},
 'negative_x':{'vitreous':(-32,-68,-16),'dark_crescent':(-32,-74,-26),'bright_cap':(-32,-78,-12)}}

def sphere(c,r):return (X-c[0])**2+(Y-c[1])**2+(Z-c[2])**2<=r*r
def title(sid):
 s=settings[sid]
 return f"seq{sid} | FA {s['FA_deg']}°\nTR/TE {s['TR_ms']:g}/{s['TE_ms']:g} ms"

metrics=[];scales={};smooth={sid:gaussian_filter(v,2,mode='nearest') for sid,v in data.items()}
for eye,points in regions.items():
 for sid in IDS:
    vitmask=sphere(points['vitreous'],4); v=float(np.median(data[sid][vitmask]));scales[sid,eye]=v
    d=float(np.median(data[sid][sphere(points['dark_crescent'],2)]))
    a=float(np.median(data[sid][sphere(points['bright_cap'],2)]))
    residual=(data[sid]-smooth[sid])[vitmask]
    metrics.append({'sequence':sid,'eye':eye,'vitreous_median':v,'dark_region_median':d,'bright_cap_median':a,
                    'dark_over_vitreous':d/v,'cap_over_vitreous':a/v,
                    'signed_dark_contrast':(v-d)/(v+d),'signed_cap_contrast':(a-v)/(a+v),
                    'relative_fine_variation_MAD':float(np.median(abs(residual-np.median(residual)))/v)})

# Matched raw and vitreous-relative sagittal views; points make sampling reviewable.
for mode in ['raw','vitreous_matched']:
 fig,axes=plt.subplots(2,6,figsize=(18,6.5),layout='constrained')
 for row,(eye,points) in enumerate(regions.items()):
  pos=39 if eye=='positive_x' else -32
  ix=np.flatnonzero(xyz[0]==pos)[0]
  iy=(xyz[1]>=-99)&(xyz[1]<=-50);iz=(xyz[2]>=-40)&(xyz[2]<=3)
  yy,zz=xyz[1][iy],xyz[2][iz]
  for col,sid in enumerate(IDS):
   ax=axes[row,col];scale=1e-8 if mode=='raw' else scales[sid,eye]
   ax.imshow(data[sid][ix,:,:][np.ix_(iy,iz)].T/scale,origin='lower',
       extent=[yy[0]-.5,yy[-1]+.5,zz[0]-.5,zz[-1]+.5],cmap='gray',vmin=0,
       vmax=5 if mode=='raw' else 1.8,interpolation='nearest')
   if row==0:ax.set_title(title(sid),fontsize=10)
   if col==0:ax.set_ylabel(eye.replace('_',' ')+f'\nx={pos} mm; z (mm)',fontsize=9)
   ax.set_xlabel('Sequence y (mm)',fontsize=8);ax.tick_params(labelsize=8)
 fig.suptitle(('Shared restored magnitude window 0–5e-8' if mode=='raw' else 'Each eye matched to its posterior-vitreous median; shared window 0–1.8')+
              '\nNative representative sagittal planes; no registration; see atlases for adjacent planes',fontsize=13)
 fig.savefig(OUT/f'sagittal_{mode}.png',dpi=160);plt.close(fig)

# ROI verification in the sagittal plane of each sphere's own x coordinate.
for eye,points in regions.items():
 fig,axes=plt.subplots(3,6,figsize=(17,9),layout='constrained')
 for row,(region,c) in enumerate(points.items()):
  ix=np.flatnonzero(xyz[0]==c[0])[0];iy=(xyz[1]>=-96)&(xyz[1]<=-52);iz=(xyz[2]>=-39)&(xyz[2]<=0)
  yy,zz=xyz[1][iy],xyz[2][iz]
  for col,sid in enumerate(IDS):
   ax=axes[row,col]
   ax.imshow(data[sid][ix,:,:][np.ix_(iy,iz)].T/scales[sid,eye],origin='lower',
      extent=[yy[0]-.5,yy[-1]+.5,zz[0]-.5,zz[-1]+.5],cmap='gray',vmin=0,vmax=1.8,interpolation='nearest')
   ax.add_patch(plt.Circle((c[1],c[2]),4 if region=='vitreous' else 2,fill=False,color='#35d4d0',linewidth=1.2))
   if row==0:ax.set_title(title(sid),fontsize=9)
   if col==0:ax.set_ylabel(region.replace('_',' ')+f'\nx={c[0]} mm',fontsize=9)
   ax.tick_params(labelsize=7)
 fig.suptitle(eye+': fixed exploratory sampling regions; no validated lens segmentation',fontsize=12)
 fig.savefig(OUT/f'{eye}_roi_check.png',dpi=145);plt.close(fig)

pairs=[('TR at 12°',702,706),('TR at 20°',705,708),('TE at 12°',706,707),('TE at 20°',708,709)]
pair_results=[]
for eye in regions:
 lookup={m['sequence']:m for m in metrics if m['eye']==eye}
 for name,base,test in pairs:
  b,a=lookup[base],lookup[test]
  pair_results.append({'comparison':name,'eye':eye,'baseline':base,'changed':test,
   'vitreous_signal_ratio':a['vitreous_median']/b['vitreous_median'],
   'dark_region_signal_ratio':a['dark_region_median']/b['dark_region_median'],
   'cap_signal_ratio':a['bright_cap_median']/b['bright_cap_median'],
   'dark_over_vitreous_before':b['dark_over_vitreous'],'dark_over_vitreous_after':a['dark_over_vitreous'],
   'cap_over_vitreous_before':b['cap_over_vitreous'],'cap_over_vitreous_after':a['cap_over_vitreous']})

# Line profiles through the anterior crescent and posterior globe boundary.
fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
colors=['#1f77b4','#2ca02c','#17becf','#d62728','#9467bd','#e377c2']
profiles={}
for row,(eye,points) in enumerate(regions.items()):
 for col,region in enumerate(['dark_crescent','vitreous']):
  c=points[region];xi=np.flatnonzero(abs(xyz[0]-c[0])<=1);zi=np.flatnonzero(abs(xyz[2]-c[2])<=1)
  yi=np.flatnonzero((xyz[1]>=-94)&(xyz[1]<=-50))
  for sid,color in zip(IDS,colors):
   p=np.median(data[sid][np.ix_(xi,yi,zi)],axis=(0,2))/scales[sid,eye]
   profiles[eye,region,sid]=p
   axes[row,col].plot(xyz[1][yi],p,label=f'seq{sid}',color=color,lw=1.7)
  axes[row,col].set_title(eye.replace('_',' ')+': '+('anterior-region profile' if col==0 else 'posterior-boundary profile'),fontsize=11)
  axes[row,col].set_xlabel('Sequence y (mm)');axes[row,col].set_ylabel('Magnitude / posterior-vitreous median')
  axes[row,col].grid(alpha=.2);axes[row,col].legend(fontsize=8,ncol=2)
fig.suptitle('Matched local magnitude profiles; no physical-resolution or layer-separation claim',fontsize=13)
fig.savefig(OUT/'relative_profiles.png',dpi=155);plt.close(fig)

result={'status':'ANALYZED; exploratory descriptive assessment','region_centres_mm':regions,
        'radii_mm':{'vitreous':4,'dark_crescent':2,'bright_cap':2},
        'anatomical_caveat':'Anterior cap/dark crescent are image-region labels; not validated lens/aqueous masks.',
        'no_registration':True,'metrics':metrics,'parameter_pairs':pair_results}
(OUT/'region_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
for row in pair_results:print(row)
