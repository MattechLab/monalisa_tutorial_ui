"""Interpret Pulseq's single-pulse response without fitting an unmeasured B0 map."""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
curves={};nulls=[]
for kind in ['gre','libre']:
 a=np.genfromtxt(OUT/f'{kind}_single_pulse_response.csv',delimiter=',',names=True)
 f=a['frequency_Hz'];s=a['Mxy_magnitude'];s0=float(np.interp(0,f,s))
 curves[kind]=(f,s/s0)
 if kind=='libre':
  indices,_=find_peaks(-s,prominence=.001)
  nulls=[{'frequency_Hz':float(f[i]),'Mxy_relative_to_on_resonance':float(s[i]/s0)}
         for i in indices if 500<abs(f[i])<1400]
rf=loadmat(OUT/'libre_exported_rf.mat',simplify_cells=True)['rf']
signal=rf['signal'];t=rf['t'];dt=float(np.median(np.diff(t)))
active=np.flatnonzero(abs(signal)>.5*max(abs(signal)))
segments=np.split(active,np.flatnonzero(np.diff(active)>1)+1)
assert len(segments)==2 and all(len(q)==500 for q in segments)
separation=float(t[segments[1]].mean()-t[segments[0]].mean())
simple_null=1/(2*separation)
assert abs(simple_null-943.3962264150944)<1e-6
assert len(nulls)==2 and all(abs(abs(q['frequency_Hz'])-simple_null)<10 for q in nulls)

# Independent small-tip Fourier check of the same exported complex waveform.
# Only zero locations are checked; this is not a steady-state image model.
sta_checks=[]
for null in nulls:
 fgrid=np.arange(null['frequency_Hz']-20,null['frequency_Hz']+20.01,.25)
 response=abs(np.exp(-2j*np.pi*fgrid[:,None]*t[None,:]) @ signal)*2*np.pi*dt
 fsta=float(fgrid[np.argmin(response)])
 sta_checks.append({'native_rotation_simulation_null_Hz':null['frequency_Hz'],'small_tip_Fourier_null_Hz':fsta})
 assert abs(fsta-null['frequency_Hz'])<10

fig,ax=plt.subplots(figsize=(10.5,4.7),layout='constrained')
for kind,color in [('gre','#1662a1'),('libre','#c34c22')]:
 f,s=curves[kind];ax.plot(f,s,label=kind.upper(),color=color,lw=2)
for q in nulls:
 ax.axvline(q['frequency_Hz'],color='#c34c22',ls=':',alpha=.8)
 ax.text(q['frequency_Hz'],.12,f"{q['frequency_Hz']:.0f} Hz",ha='center',fontsize=9)
ax.axvline(0,color='gray',ls='--',alpha=.5)
ax.set_xlim(-1600,1600);ax.set_ylim(0,2)
ax.set_xlabel('Local water frequency offset from nominal water resonance (Hz)')
ax.set_ylabel('Single-pulse transverse response / own response at 0 Hz')
ax.set_title('Actual exported seq700 RF pulses: frequency-selective LIBRE cancellation\nPulseq rotation simulation; equilibrium initial magnetization; no relaxation or measured field map',fontsize=11)
ax.grid(alpha=.2);ax.legend()
fig.savefig(OUT/'rf_frequency_response.png',dpi=170);plt.close(fig)

profiles=np.genfromtxt(OUT/'ratio_profiles.csv',delimiter=',',names=True)
indices,_=find_peaks(-profiles['LIBRE_linear_over_normal'],prominence=.15,distance=15)
bands=[{'x_mm':float(profiles['x_mm'][i]),'GRE_linear_over_normal':float(profiles['GRE_linear_over_normal'][i]),
        'LIBRE_linear_over_normal':float(profiles['LIBRE_linear_over_normal'][i])} for i in indices]
estimate=(nulls[1]['frequency_Hz']-nulls[0]['frequency_Hz'])/(bands[1]['x_mm']-bands[0]['x_mm']) if len(bands)==2 else None
report={'simulation':'Native Pulseq mr.simRf on the first actual exported RF pulse',
        'subpulse_centre_separation_us':separation*1e6,'simple_interference_null_magnitude_Hz':simple_null,
        'native_simulation_nulls':nulls,'independent_small_tip_null_checks':sta_checks,
        'observed_spatial_ratio_minima':bands,
        'conditional_frequency_gradient_Hz_per_mm':estimate,
        'conditional_gradient_caveat':'Inferred by assigning the two observed spatial bands to the two simulated RF nulls; not an independent B0 measurement.',
        'limitations':['No measured spatial B0 or B1 map was available.',
                       'Single-pulse simulation does not predict the full steady-state intensity pattern or readout dephasing.',
                       'RF amplitudes are used exactly as exported; sequence/user flip-angle labels do not substitute for a measured local net rotation.',
                       'Header gradient-offset units were not converted to a physical gradient without scanner calibration.']}
(OUT/'rf_response_interpretation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
