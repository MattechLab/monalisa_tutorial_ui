# Anisotropic scanner pipeline

Current handover: [HANDOVER_ANISOTROPIC_RECON_20261111.md](HANDOVER_ANISOTROPIC_RECON_20261111.md).

This pipeline reads scanner Siemens `.dat` and matching Pulseq `.seq`; it does not run Koma.

```matlab
addpath('/Users/cag/forclone/monalisa_tutorial_ui');
dataset='/Volumes/Extreme SSD/261110';
sigmas=[.61 .32 .20 .10];
% Read-only plan (no folders, no raw-signal loading):
plan=run_anisotropic_scanner(dataset,sigmas,20,[],false,true);
% Recompute only when intended, into a new timestamped batch:
% plan=run_anisotropic_scanner(dataset,sigmas,20);
```

The independent anisotropic cell in `main_script_7TTerraVE12U.m` uses the same four widths. New runs save verified max-normalized float32 NIfTI plus metadata/logs/previews and one shared encoding MAT per sequence. The transient full physical-coil MAT is deleted only after NIfTI/encoding verification. The low-level `recon_anisotropic_xrms` still saves its full MAT if invoked directly.

## Current saved batch

`/Volumes/Extreme SSD/261110/recon_results/anisotropic_scanner/retained_20261111`

Each `seq712_sigma0p61`, `seq712_sigma0p32`, `seq712_sigma0p2`, `seq712_sigma0p1` directory contains `xrms_normalized.nii.gz`, `reconstruction.json`, and the scanner validation/log/preview. The same four directories exist for seq713 and seq714. Shared `encoding/seq71*.mat` contains `t`, `ve`, ADC times, metadata and coordinate vectors, not coil images.

Load without losing the original comparison scale:

```matlab
addpath('/Users/cag/forclone/monalisa_tutorial_ui/pipeline_utils');
r=load_anisotropic_nifti('/Volumes/Extreme SSD/261110/recon_results/anisotropic_scanner/retained_20261111/seq712_sigma0p32');
% r.xrms is restored to its original magnitude scale (float32 export precision).
% Pass false as argument 2 to read the normalized values instead.
```

Normalization is `single(double(xrms)/max(double(xrms(:))))`; no clipping, interpolation or axis permutation. The JSON divisor reverses the amplitude scaling. Shared-scale comparisons must restore that divisor. NIfTI sform code 2 denotes the nominal sequence grid, not validated patient anatomy.

The scientific reconstruction remains physical Pulseq ROTATIONS, cycles/mm, 240 × 240 × 180 mm FOV, 2-mm isotropic nominal voxels, 20 discarded initial shots, 7,697 retained imaging spokes, spherical Voronoi DCF, Gaussian window 3/parameters `[sigma 10]`, and physical-coil xRMS. No added phase or gradient-delay correction.

Only 0.61/0.32/0.20/0.10 remain. Widths 0.50/0.40/0.27/0.15 and full coil-image MATs were explicitly removed. Old comparison figures/metrics remain as history in `comparisons`; new `orbital_retained_zm24.*` and `zm20.*` show the four retained widths on one restored raw intensity scale. Historical source scripts are under `comparisons/historical_scripts` and refer to pre-cleanup experiments; use the current entry point for new scanner work.

Static phantom diagnostics are separate at `/Volumes/Extreme SSD/261110/recon_results/water_phantom_diagnostics`. Their active scripts use the retained shared encoding. The cleanup audit is in `../output_review_20261111/cleanup_completed/`.
