"""Native planes through manually selected representative lens regions.

Points are approximate visual localizations, not validated lens centres or
segmentations; source atlases are saved alongside this script for inspection.
No quantitative lens contrast or displacement is inferred from these points.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
IDS=[701,702,704,703,705]
FA={701:12,702:12,704:14,703:16,705:20}
xyz=[np.arange(-60,61),np.arange(-98,-43),np.arange(-43,4)]
points={
 'positive_x':{701:(40,-79,-23),702:(39,-79,-22),704:(36,-82,-23),703:(36,-82,-22),705:(36,-82,-23)},
 'negative_x':{701:(-32,-85,-24),702:(-32,-85,-24),704:(-33,-85,-24),703:(-32,-85,-23),705:(-33,-85,-24)}}
metrics=json.loads((OUT/'posterior_vitreous_variation.json').read_text())['measurements']
scales={(r['sequence'],r['eye']):r['median_xrms'] for r in metrics}
with np.load(OUT/'orbital_crops.npz') as source:
 fig,axes=plt.subplots(2,5,figsize=(15,6),layout='constrained')
 for row,(eye,locations) in enumerate(points.items()):
  for col,sid in enumerate(IDS):
   centre=locations[sid]; ix=np.flatnonzero(xyz[0]==centre[0])[0]
   iy=(xyz[1]>=-94)&(xyz[1]<=-51); iz=(xyz[2]>=-40)&(xyz[2]<=-5)
   a=source[str(sid)][ix,:,:][np.ix_(iy,iz)].T/scales[sid,eye]
   yy,zz=xyz[1][iy],xyz[2][iz];ax=axes[row,col]
   ax.imshow(a,origin='lower',extent=[yy[0]-.5,yy[-1]+.5,zz[0]-.5,zz[-1]+.5],cmap='gray',vmin=0,vmax=1.7,interpolation='nearest')
   ax.add_patch(plt.Circle((centre[1],centre[2]),1.5,fill=False,color='#35d4d0',linewidth=1.2))
   if row==0:ax.set_title(f'seq{sid} | FA {FA[sid]}°',fontsize=11)
   ax.set_xlabel(f'Native x={centre[0]} mm; sequence y (mm)',fontsize=8)
   if col==0:ax.set_ylabel(eye.replace('_',' ')+'\nSequence z (mm)',fontsize=9)
   ax.tick_params(labelsize=8)
 fig.suptitle('Representative lens-region native sagittal planes, vitreous-matched display\nCyan: approximate manually selected lens-region location; no registration or validated segmentation',fontsize=12)
 fig.savefig(OUT/'lens_region_sagittal.png',dpi=160);plt.close(fig)
(OUT/'lens_region_localizations.json').write_text(json.dumps({
 'method':'approximate manual positions selected by visual review of native orbital atlases',
 'positions_mm':points,'reported_FA_deg':FA,'sequence_file_FA_factor':2,
 'limitations':['These are representative regions, not validated lens centres or masks.',
                'No lens displacement or gaze angle is inferred from these points.',
                'Automatic darkest-focus searches reached the anterior/inferior boundary in some scans and were rejected for lens-contrast quantification.']},indent=2)+'\n')
