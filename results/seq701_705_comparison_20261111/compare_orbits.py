"""Descriptive comparison of existing reconstructions; source files stay intact.

Coordinate convention: default orthogonal Slicer offset S=z and L=-x.
Half-voxel planes are linearly interpolated only for these display figures.
No registration, denoising, anatomical reorientation or inferential tests.
"""
from pathlib import Path
import gzip
import json
import struct
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('/Volumes/Extreme SSD/261110/recon_results/isotropic_scanner/run_20261111_174401_195')
OUT = Path(__file__).resolve().parent
IDS = list(range(701, 706))
FA_FACTOR = 2  # User correction, 12 Nov 2026: sequence 6/7/8/10 -> 12/14/16/20 deg.
VOLUMES, REPORTS = {}, {}
for sid in IDS:
    folder = ROOT / f'seq{sid}'
    report = json.loads((folder / 'reconstruction.json').read_text())
    with gzip.open(folder / 'xrms_normalized.nii.gz', 'rb') as stream:
        header = stream.read(352)
        assert struct.unpack_from('<4h', header, 40) == (3, 240, 240, 240)
        assert struct.unpack_from('<h', header, 70)[0] == 16
        assert struct.unpack_from('<f', header, 108)[0] == 352
        affine = np.array([struct.unpack_from('<4f', header, p) for p in [280, 296, 312]])
        np.testing.assert_array_equal(affine, np.c_[np.eye(3), [-120, -120, -120]])
        data = np.frombuffer(stream.read(), dtype='<f4').reshape((240, 240, 240), order='F')
    VOLUMES[sid] = data
    REPORTS[sid] = report

COORDS = np.arange(240) - 120
SPECS = [(2, -19.5, 'Axial S = -19.5 mm'), (2, -21.5, 'Axial S = -21.5 mm'),
         (0, 36.5, 'Sagittal L = -36.5 mm (x = +36.5)'),
         (0, -28.5, 'Sagittal L = +28.5 mm (x = -28.5)')]

def plane(data, axis, position):
    index = position + 120
    lo = int(np.floor(index)); weight = index - lo
    return ((1 - weight) * np.take(data, lo, axis=axis) +
            weight * np.take(data, lo + 1, axis=axis))

def cropped(data, axis, position):
    a = plane(data, axis, position)
    if axis == 2:
        i = (COORDS >= -65) & (COORDS <= 65)
        j = (COORDS >= -98) & (COORDS <= -38)
    else:
        i = (COORDS >= -98) & (COORDS <= -38)
        j = (COORDS >= -46) & (COORDS <= 5)
    return a[np.ix_(i, j)].T, (COORDS[i][0]-.5, COORDS[i][-1]+.5,
                                     COORDS[j][0]-.5, COORDS[j][-1]+.5)

def title(sid):
    h = REPORTS[sid]['metadata']['sequenceParams']['header']
    return f"seq{sid}\nTR {1000*h['TR_s']:g} ms | FA {FA_FACTOR*h['rf_flip_deg']:g}°"

for mode in ['common_raw', 'max_normalized']:
    fig, axes = plt.subplots(4, 5, figsize=(17, 10.5), layout='constrained')
    for row, (axis, pos, label) in enumerate(SPECS):
        panels = []
        for sid in IDS:
            data, extent = cropped(VOLUMES[sid], axis, pos)
            if mode == 'common_raw':
                data = data * REPORTS[sid]['normalizationDivisor']
            panels.append(data)
        vmax = np.percentile(np.concatenate([a.ravel() for a in panels]), 99.5)
        for col, (sid, a) in enumerate(zip(IDS, panels)):
            ax = axes[row, col]
            im = ax.imshow(a, origin='lower', extent=extent, cmap='gray',
                           vmin=0, vmax=vmax, interpolation='nearest')
            if axis == 2: ax.invert_yaxis()
            ax.set_xlabel('Sequence x (mm)' if axis == 2 else 'Sequence y (mm)', fontsize=9)
            if col == 0: ax.set_ylabel(label + '\n' + ('Sequence y (mm)' if axis == 2 else 'Sequence z (mm)'), fontsize=9)
            if row == 0: ax.set_title(title(sid), fontsize=11)
            ax.tick_params(labelsize=8)
        fig.colorbar(im, ax=axes[row, :], shrink=.75, pad=.01,
                     label='Restored xRMS' if mode == 'common_raw' else 'xRMS / own maximum')
    fig.suptitle(('Shared reconstructed magnitude scale' if mode == 'common_raw' else 'Existing max-normalized volumes') +
                 '\nIdentical window across five scans within each row; half-voxel planes linearly interpolated; no registration', fontsize=13)
    fig.savefig(OUT / f'{mode}.png', dpi=155)
    plt.close(fig)

# Native neighbouring planes show whether a half-voxel reslice drives a finding.
for axis, pos, label in SPECS:
    fig, axes = plt.subplots(2, 5, figsize=(16, 5.5), layout='constrained')
    native_positions = [np.floor(pos), np.ceil(pos)]
    panels = {(r, sid): cropped(VOLUMES[sid], axis, p)[0] * REPORTS[sid]['normalizationDivisor']
              for r, p in enumerate(native_positions) for sid in IDS}
    vmax = np.percentile(np.concatenate([a.ravel() for a in panels.values()]), 99.5)
    for r, p in enumerate(native_positions):
        for c, sid in enumerate(IDS):
            _, extent = cropped(VOLUMES[sid], axis, p)
            axes[r, c].imshow(panels[r, sid], origin='lower', extent=extent,
                              cmap='gray', vmin=0, vmax=vmax, interpolation='nearest')
            if axis == 2: axes[r, c].invert_yaxis()
            if r == 0: axes[r, c].set_title(title(sid), fontsize=10)
            if c == 0: axes[r, c].set_ylabel(f'Native {"z" if axis == 2 else "x"}={p:g} mm')
            axes[r, c].tick_params(labelsize=8)
    fig.suptitle(label + ' — adjacent native slices, shared magnitude window', fontsize=13)
    fig.savefig(OUT / f'native_axis{axis}_at{pos:g}.png', dpi=155)
    plt.close(fig)

manifest = {'source': str(ROOT), 'coordinate_mapping': 'default Slicer orthogonal offsets: S=z, L=-x',
            'half_voxel_display': 'linear interpolation of two adjacent native planes',
            'registration': False, 'purpose': 'descriptive image comparison, not validated IQM ranking',
            'flip_angle_factor': FA_FACTOR, 'flip_angle_basis': 'User correction on 12 Nov 2026: double sequence-file rf_flip_deg',
            'scans': []}
for sid in IDS:
    h = REPORTS[sid]['metadata']['sequenceParams']['header']
    manifest['scans'].append({'sequence': sid, 'TR_ms': h['TR_s']*1000, 'TE_ms': h['TE_s']*1000,
                              'FA_deg': FA_FACTOR*h['rf_flip_deg'], 'sequence_file_FA_deg': h['rf_flip_deg'], 'duration_s': h['TotalDuration'],
                              'maximum_xrms': REPORTS[sid]['normalizationDivisor']})
(OUT / 'comparison_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')

# Exploratory interior-vitreous samples, away from the lens and globe wall.
# Same native coordinates in every scan; no segmentation or noise calibration.
roi_centres = {'positive_x_eye': (36, -72, -21), 'negative_x_eye': (-31, -75, -21)}
roi_radius = 3
measurements = []
fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
for ax, (name, centre) in zip(axes, roi_centres.items()):
    ref = None
    for sid in IDS:
        ix = [np.flatnonzero(abs(COORDS - c) <= roi_radius) for c in centre]
        xx, yy, zz = np.meshgrid(*(COORDS[i] for i in ix), indexing='ij')
        mask = ((xx-centre[0])**2 + (yy-centre[1])**2 + (zz-centre[2])**2) <= roi_radius**2
        a = VOLUMES[sid][np.ix_(*ix)][mask].astype(float) * REPORTS[sid]['normalizationDivisor']
        median = np.median(a)
        if ref is None: ref = median
        measurements.append({'sequence': sid, 'region': name, 'centre_sequence_mm': centre,
                             'radius_mm': roi_radius, 'n_voxels': int(mask.sum()),
                             'median_restored_xrms': float(median),
                             'median_relative_to_seq701': float(median/ref),
                             'IQR_over_median': float(np.subtract(*np.percentile(a, [75,25]))/median)})
        # Posterior-boundary profiles averaged over a 5x5-mm native patch.
        xp = np.flatnonzero(abs(COORDS - centre[0]) <= 2)
        zp = np.flatnonzero(abs(COORDS - centre[2]) <= 2)
        yp = np.flatnonzero((COORDS >= -94) & (COORDS <= -50))
        profile = VOLUMES[sid][np.ix_(xp, yp, zp)].mean(axis=(0,2)) * REPORTS[sid]['normalizationDivisor']
        ax.plot(COORDS[yp], profile/median, label=f'seq{sid}', linewidth=1.8)
    ax.set_title(name.replace('_', ' ')); ax.set_xlabel('Sequence y (mm)')
    ax.set_ylabel('Profile / own interior-vitreous median'); ax.grid(alpha=.25); ax.legend()
fig.suptitle('Exploratory anterior–posterior profiles: relative contrast, not a resolution or SNR measurement')
fig.savefig(OUT / 'vitreous_relative_profiles.png', dpi=160); plt.close(fig)
(OUT / 'exploratory_vitreous_samples.json').write_text(json.dumps({
    'interpretation': 'Descriptive unregistered interior samples; IQR includes spatial signal variation/artifacts and is not a noise or SNR estimate.',
    'measurements': measurements}, indent=2)+'\n')

# Visually verify that the exploratory samples fall inside vitreous.
fig, axes = plt.subplots(2, 5, figsize=(15, 6), layout='constrained')
for row, (name, centre) in enumerate(roi_centres.items()):
    panels = [cropped(VOLUMES[sid], 0, centre[0])[0] * REPORTS[sid]['normalizationDivisor'] for sid in IDS]
    vmax = np.percentile(np.concatenate([a.ravel() for a in panels]), 99.5)
    for col, sid in enumerate(IDS):
        _, extent = cropped(VOLUMES[sid], 0, centre[0])
        ax = axes[row, col]
        ax.imshow(panels[col], origin='lower', extent=extent, cmap='gray', vmin=0, vmax=vmax, interpolation='nearest')
        ax.add_patch(plt.Circle((centre[1], centre[2]), roi_radius, fill=False, edgecolor='#35d4d0', linewidth=1.4))
        if row == 0: ax.set_title(title(sid), fontsize=10)
        if col == 0: ax.set_ylabel(name.replace('_', ' ') + f'\nx={centre[0]} mm')
        ax.tick_params(labelsize=8)
fig.suptitle('Exploratory 3-mm-radius interior-vitreous samples; unchanged masks across scans')
fig.savefig(OUT / 'vitreous_roi_check.png', dpi=155); plt.close(fig)

fig, axes = plt.subplots(2, 5, figsize=(15, 6), layout='constrained')
for row, ((name, _), (_, position, label)) in enumerate(zip(roi_centres.items(), SPECS[2:])):
    for col, sid in enumerate(IDS):
        sample = next(m for m in measurements if m['sequence'] == sid and m['region'] == name)
        a, extent = cropped(VOLUMES[sid], 0, position)
        a = a * REPORTS[sid]['normalizationDivisor'] / sample['median_restored_xrms']
        ax = axes[row, col]
        ax.imshow(a, origin='lower', extent=extent, cmap='gray', vmin=0, vmax=1.7, interpolation='nearest')
        if row == 0: ax.set_title(title(sid), fontsize=10)
        if col == 0: ax.set_ylabel(label + '\nSequence z (mm)', fontsize=9)
        ax.set_xlabel('Sequence y (mm)', fontsize=9); ax.tick_params(labelsize=8)
fig.suptitle('Vitreous-matched display: compare relative tissue contrast\nEach eye divided by its own interior-vitreous median; identical window 0–1.7 for every panel', fontsize=13)
fig.savefig(OUT / 'vitreous_matched_sagittal.png', dpi=155); plt.close(fig)

(OUT / 'analysis_validation.json').write_text(json.dumps({
    'verification_status': 'ANALYZED',
    'comparison_type': 'exploratory descriptive assessment of five scans in one session',
    'coordinate_assumption': 'default orthogonal Slicer views, without additional transforms',
    'sources_modified': False, 'inferential_tests': 0,
    'limitations': ['No repeated acquisitions or independent subjects',
                    'No cross-scan rigid or eye-motion registration',
                    'ROI locations are manually selected exploratory samples, not validated segmentations',
                    'Spatial variation inside magnitude images is not a calibrated noise estimate',
                    '1-mm output grid does not establish true resolution or separation of ocular-wall layers'],
    'fallacy_scan_coverage': '11/11 checked',
    'fallacy_scan_notes': {
        'Simpson': 'Both eyes reported separately; no reversal in main signal direction',
        'ecological': 'No population inference; voxels are not independent replicates',
        'Berkson': 'No selected-population association estimated',
        'collider': 'No covariate-conditioned model',
        'base_rate': 'No diagnostic accuracy claim',
        'regression_to_mean': 'No selection of subjects by extreme baseline',
        'survivorship': 'All five completed scans included',
        'look_elsewhere': 'No hypothesis tests or selective significance reporting',
        'forking_paths': 'Exploratory ROI selection disclosed, native neighbouring planes inspected',
        'causality': 'Apparent parameter trends remain confounded by scan-to-scan motion and field/receive effects',
        'reverse_causality': 'No reverse-causal association estimated'}}, indent=2)+'\n')
for row in measurements: print(row)
print('Saved matched requested-plane panels and adjacent native-slice checks:', OUT)
