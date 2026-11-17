"""Exploratory localization review, not a validated lens segmentation.

Locate a smoothed low-signal anterior focus within manually bounded inner
globe regions. Orthogonal overlays must be inspected before interpreting it.
"""
from pathlib import Path
import json
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
IDS=[701,702,704,703,705]
FA={701:12,702:12,704:14,703:16,705:20}
xyz=[np.arange(-60,61),np.arange(-98,-43),np.arange(-43,4)]
X,Y,Z=np.meshgrid(*xyz,indexing='ij')
specs={'positive_x':((36,-72,-21),((30,44),(-85,-77),(-28,-17))),
       'negative_x':((-32,-76,-22),((-39,-25),(-90,-80),(-28,-17)))}
variation=json.loads((OUT/'posterior_vitreous_variation.json').read_text())['measurements']
scales={(r['sequence'],r['eye']):r['median_xrms'] for r in variation}
results=[]
with np.load(OUT/'orbital_crops.npz') as source:
 for name,(c,bounds) in specs.items():
    mask=((X-c[0])/12)**2+((Y-c[1])/13)**2+((Z-c[2])/13)**2<.84**2
    for q,(lo,hi) in zip([X,Y,Z],bounds):mask &= (q>=lo)&(q<=hi)
    fig,axes=plt.subplots(3,5,figsize=(15,9),layout='constrained')
    for col,sid in enumerate(IDS):
        data=source[str(sid)]
        smoothed=gaussian_filter(data,.7)
        point=np.unravel_index(np.argmin(np.where(mask,smoothed,np.inf)),data.shape)
        centre=[int(q[k]) for q,k in zip(xyz,point)]
        # Diagnostic only: these samples are selected by a local minimum.
        values={}
        for radius in [1,1.5,2]:
            core=(X-centre[0])**2+(Y-centre[1])**2+(Z-centre[2])**2<=radius**2
            values[str(radius)]=float(np.median(data[core])/scales[sid,name])
        results.append({'eye':name,'sequence':sid,'candidate_mm':centre,
                        'candidate_selected_by':'minimum after 0.7-voxel Gaussian smoothing within constrained search region',
                        'native_core_to_vitreous_median_by_radius_mm':values})
        for row,axis in enumerate([2,0,1]):
            remaining=[k for k in range(3) if k!=axis]
            a=np.take(data,point[axis],axis=axis)
            u,v=[xyz[k] for k in remaining]
            keep=[abs(q-centre[k])<=13 for q,k in zip([u,v],remaining)]
            a=a[np.ix_(*keep)].T/scales[sid,name]
            uu,vv=u[keep[0]],v[keep[1]]
            ax=axes[row,col]
            ax.imshow(a,origin='lower',extent=[uu[0]-.5,uu[-1]+.5,vv[0]-.5,vv[-1]+.5],cmap='gray',vmin=0,vmax=1.7,interpolation='nearest')
            ax.add_patch(plt.Circle([centre[k] for k in remaining],1.5,fill=False,color='#35d4d0',linewidth=1.3))
            if axis==2:ax.invert_yaxis()
            if row==0:ax.set_title(f'seq{sid} | FA {FA[sid]}°\ncandidate {tuple(centre)} mm',fontsize=9)
            if col==0:ax.set_ylabel(f'{"xyz"[axis]}={centre[axis]} mm\n{"xyz"[remaining[1]]} (mm)',fontsize=9)
            ax.set_xlabel(f'{"xyz"[remaining[0]]} (mm)',fontsize=8);ax.tick_params(labelsize=7)
    fig.suptitle(f'{name}: candidate anterior low-signal focus, orthogonal native planes\nExploratory localization; circles are not validated lens masks; vitreous-matched display',fontsize=12)
    fig.savefig(OUT/f'{name}_candidate_focus.png',dpi=155);plt.close(fig)
(OUT/'candidate_focus_measurements.json').write_text(json.dumps({'status':'exploratory localization only; no motion estimate or validated lens contrast',
    'measurements':results},indent=2)+'\n')
for r in results:print(r)
