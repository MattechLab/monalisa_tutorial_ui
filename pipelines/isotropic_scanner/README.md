# Isotropic scanner xRMS

`run_isotropic_scanner` follows the main-scan reconstruction in
`main_script_7TTerraVE12U.m`: `load_siemens_lowres`, then per-physical-coil
`bmMathilda`, and `sqrt(mean(abs(coilImage).^2))`. No calibration or ROVir is
required. Coil images are combined incrementally and are not saved.

```matlab
addpath('/Users/cag/forclone/monalisa_tutorial_ui');
plan=run_isotropic_scanner('/Volumes/Extreme SSD/261110');
```

Arguments are `(datasetFolder,outDir,nShotOff,scanSet)`. Defaults are a new timestamped
folder under `<dataset>/recon_results/isotropic_scanner`, and 50 discarded
initial shots, matching the 7T main-scan script. Existing output folders are
refused. All five raw/sequence pairings are checked before reconstruction.

For the six scans selected from the 13 November dataset:

```matlab
plan=run_isotropic_scanner('/Volumes/Extreme SSD/261113',[],50,'seq702_705to709');
```

This selects only seq702, seq706, seq707, seq705, seq708 and seq709, in that
acquisition order, using each sequence's unique matching raw file. Settings
remain 240³, 1-mm output spacing, 50 excluded shots, physical-coil xRMS and
verified max-normalized float32 NIfTI. Each sidecar retains the sequence-file
flip angle and separately records the user's ×2 LIBRE reporting convention;
the reported angle is not a measurement of local excitation.

For the seq700 normal-shimming GRE/LIBRE pair, use the same reconstruction:

```matlab
plan=run_isotropic_scanner('/Volumes/Extreme SSD/261110',[],50,'seq700_norm_shimming');
% The linear-shimming pair uses the same settings and a separate new batch:
plan=run_isotropic_scanner('/Volumes/Extreme SSD/261110',[],50,'seq700_linear_shimming');
```

This selects only `meas_*seq700_gre_norm_shimming.dat` and
`meas_*seq700_libre_norm_shimming.dat`, verifies each against its matching
GRE/LIBRE sequence, and saves separate labelled folders. Each source has 199
shots; discarding 50 initial shots and each SI leaves 149 × 43 = 6,407 imaging
spokes. The matrix, FOV, kernel, coil combination and verified compact export
are identical to the seq701–705 path. Original sequence flip-angle definitions
are retained in reconstruction metadata.

`seq700_linear_shimming` selects only the two `*_linear_shimming.dat` files,
with the same GRE/LIBRE sequence matching and output checks. Its separate
labels prevent mixing normal- and linear-shimming outputs.

The reconstruction grid is 240 × 240 × 240 across the declared 240-mm cubic
FOV: 1-mm isotropic output voxels. Acquisition `res_mm` remains 0.5 mm in the
source metadata. Each scan has 500 × 44 readouts, including one SI per shot;
after exclusions, 450 × 43 = 19,350 imaging spokes remain.

The existing isotropic loader uses Pulseq ROTATIONS through `bmTraj`, scales
the maximum trajectory radius to `N/(2*FoV)`, computes spherical Voronoi DCF
after shot/SI exclusion, then selects samples inside the 240³ Cartesian
k-space grid using `bmLowRes`. Gridding uses the 7T defaults: Gaussian window
3, parameters `[0.61 10]`, normalized gridding, MATLAB IFFT and deapodization.
No added phase or gradient-delay correction is applied.

Each `seq701`–`seq705` folder saves `xrms_normalized.nii.gz`,
`reconstruction.json`, `reconstruction.log`, and `xrms_preview.png`.
The existing verified exporter is shared with the anisotropic pipeline.
Normalization divides each volume by its own maximum, without clipping.
The JSON divisor restores the original magnitude scale at float32 precision.
NIfTI read-back must match exactly, with relative restoration error <1e-7,
before deleting the temporary magnitude MAT. No full coil-image or encoding
MAT is retained. Original raw/sequence files are preserved.

Coordinates are the nominal centred sequence grid, with 1-mm spacing and
first voxel centre (−120, −120, −120) mm; patient anatomical orientation has
not been established. The export performs no axis permutation or resampling.
`batch_plan.mat` records completion per sequence and for the batch.

## Completed batch — 11 November 2026

`/Volumes/Extreme SSD/261110/recon_results/isotropic_scanner/run_20261111_174401_195`

All five scans completed with the defaults above. The five compressed NIfTIs
total 226.73 MiB (45.29–45.45 MiB each). MATLAB verified exact normalized
float32 read-back for every volume; maximum relative intensity-restoration
error was 2.59e-8. An independent Python standard-library check verified all
five NIfTI headers, 240³ payload sizes, 1-mm spacing, float32 datatype, MD5
checksums and gzip integrity. No temporary magnitude MAT or uncompressed
NIfTI remains in the scan output folders. Raw and sequence files were unchanged
by reconstruction. These checks establish numerical/export integrity, not
image-quality or patient-orientation validation.

## Completed seq700 normal-shimming pair — 12 November 2026

`/Volumes/Extreme SSD/261110/recon_results/isotropic_scanner/run_seq700_norm_shimming_20261112_094850_901`

Both GRE and LIBRE completed at 240³ with 1-mm spacing, 50 discarded initial
shots, 6,407 retained imaging spokes and the same `[0.61 10]` Gaussian kernel.
Separate `seq700_gre_norm_shimming` and `seq700_libre_norm_shimming` folders
contain verified normalized float32 NIfTI, reversible scaling JSON, logs and
previews. NIfTI sizes are 46.53 and 46.61 MiB respectively. Exact MATLAB
read-back passed; maximum relative intensity-restoration error was 2.43e-8.
Independent checks passed for MD5, gzip integrity, header geometry/datatype,
payload size and expected acquisition/method metadata. Temporary MAT and
uncompressed NIfTI files were removed after verification. Local verification:
`results/seq700_norm_shimming_20261112_verification.json`.

Both previews depict a spherical phantom, not the orbital subject. These
scans can assess phantom excitation/shimming behaviour but are not direct
anatomical controls for the lens-region finding.

## Completed seq700 linear-shimming pair — 12 November 2026

`/Volumes/Extreme SSD/261110/recon_results/isotropic_scanner/run_seq700_linear_shimming_20261112_100256_121`

Both GRE and LIBRE completed with the same matrix, spacing, exclusion, kernel,
coil combination and export method as the normal-shimming pair. Verified
float32 normalized NIfTI sizes are 46.69 MiB (GRE) and 46.32 MiB (LIBRE).
Maximum relative intensity-restoration error was 2.37e-8. Exact MATLAB
read-back and independent NIfTI/header/payload/MD5/gzip checks passed. The
verification also matched the normal-shimming grid, kernel, retained spokes,
sample count, sequence hash, loader, trajectory function and gridding MEX.
Temporary magnitude MAT and uncompressed NIfTI files were removed only after
export verification. Local audit:
`results/seq700_linear_shimming_20261112_verification.json`.

Both previews show the spherical phantom. The LIBRE linear-shimming preview
contains a prominent dark signal-loss band in the XY/XZ projections; its
physical cause has not been established by these reconstruction checks.

## Seq700 normal/linear phantom comparison — 12 November 2026

The subsequent comparison is saved under
`results/seq700_shimming_comparison_20261112/`. It restores the JSON intensity
divisors, checks identical processing within each GRE/LIBRE pair, and compares
native z=0/y=0 slices plus central-volume linear/normal ratio profiles.
LIBRE profile minima occur near x=−39 and +42 mm (8.1% and 6.2% of normal
signal); GRE retains 98.6% and 89.8% at those same profile positions.
Raw headers show X-gradient offset 273→3759 in scanner units, unchanged
Y/Z offsets and listed higher-order shim currents, and centre frequency
297191177→297191221 Hz. Per-sequence reference/RF amplitudes are unchanged.

Pulseq `mr.simRf` on the actual exported RF pulses predicts LIBRE excitation
nulls near −940/+945 Hz. A Fourier small-tip calculation independently agrees
within 2 Hz. The 530-µs subpulse-centre spacing gives the simple cancellation
scale 1/(2×530 µs)≈943 Hz. Together these results strongly support off-resonance
LIBRE excitation cancellation under the altered X offset. No measured B0/B1
map was available; exact voxelwise frequency assignment and the full
steady-state intensity pattern remain unverified. This phantom result does
not establish the cause of the orbital lens-region crescent.

## Completed selected orbital batch — 13 November 2026

`/Volumes/Extreme SSD/261113/recon_results/isotropic_scanner/run_seq702_705to709_20261113_155007_367`

All six highlighted scans completed in acquisition order: seq702, seq706,
seq707, seq705, seq708 and seq709. The unselected seq720–723 scans were not
included. Each output is 240³ physical-coil xRMS, with 1-mm spacing, 50
discarded initial shots, 19,350 retained imaging spokes, and Gaussian
window 3 / parameters `[0.61 10]`.

| Scan | TR (ms) | TE (ms) | Reported LIBRE FA | NIfTI size (MiB) |
|---|---:|---:|---:|---:|
| 702 | 6.7 | 2.875 | 12° | 45.26 |
| 706 | 13 | 2.875 | 12° | 45.90 |
| 707 | 13 | 4 | 12° | 45.34 |
| 705 | 6.7 | 2.875 | 20° | 45.41 |
| 708 | 13 | 2.875 | 20° | 45.31 |
| 709 | 13 | 4 | 20° | 45.30 |

The reported flip angle follows the user's ×2 LIBRE convention; the original
sequence-file values remain separately recorded in every JSON. All six
max-normalized float32 compressed NIfTIs passed exact MATLAB read-back and
independent source-pairing, header/affine/datatype/payload, finite normalized
values, metadata, MD5 and gzip-integrity checks. Source sequence hashes and
raw sizes also matched. The six NIfTIs total 272.52 MiB. Temporary magnitude
MAT and uncompressed NIfTI files were removed after verification.

Audit and log: `results/selected_261113/verification.json`,
`results/selected_261113/matlab_batch.log`, and `completed.mat`. The independent
verifier is `results/selected_261113/verify_exports.py`. These checks establish
numerical/export integrity, not image quality or patient orientation.
