# Handover — anisotropic scanner reconstruction

Updated 11 November 2026 after completed output cleanup and verified NIfTI conversion. This is the current scanner reconstruction handover.

## Entry point and separation

Repository: `/Users/cag/forclone/monalisa_tutorial_ui`.
Standalone runner: `run_anisotropic_scanner.m`.
7T entry: `main_script_7TTerraVE12U.m`, independent anisotropic cell; earlier calibration cells are not required.

This pipeline reads scanner Siemens `.dat` and matching Pulseq `.seq`, reconstructs physical receive channels, and saves compact normalized xRMS. It never simulates a phantom. The companion [Koma handover](/Users/cag/forclone/pulseq4mreye_rectfov/simulation/koma_water/HANDOVER_KOMA_WATER_20261111.md) covers simulated data and static Fourier diagnostics.

## Current retained results

```text
/Volumes/Extreme SSD/261110/recon_results/anisotropic_scanner/retained_20261111/
  seq712_sigma0p61/  seq712_sigma0p32/  seq712_sigma0p2/  seq712_sigma0p1/
  seq713_sigma0p61/  seq713_sigma0p32/  seq713_sigma0p2/  seq713_sigma0p1/
  seq714_sigma0p61/  seq714_sigma0p32/  seq714_sigma0p2/  seq714_sigma0p1/
  encoding/seq712.mat, seq713.mat, seq714.mat
  comparisons/
  conversion_manifest.json
```

Exactly 12 scanner magnitude images remain. The 0.61 volumes are the baseline. Each image folder contains:

- `xrms_normalized.nii.gz`: float32, compressed, 120 × 120 × 90, 2-mm isotropic nominal voxels.
- `reconstruction.json`: normalization divisor, metadata, source MAT hash, NIfTI hash and export validation.
- Original scanner `validation.json`, `reconstruction.log`, and `xrms_preview.png`.

The three shared encoding MATs store actual retained `t`, `ve`, ADC times, original metadata and x/y/z coordinate vectors. They contain no image/coil arrays. Full complex coil images and individual full `reconstruction.mat` containers are no longer retained. Regenerate from protected raw data if phase or coil images are needed later.

## Normalization, coordinates and safe loading

Normalization is:

```matlab
normalized = single(double(xrms) / max(double(xrms(:))));
```

No clipping, resampling, registration or axis permutation. The saved maximum is `normalizationDivisor`. Each volume has its own divisor; for common-scale comparisons, restore it:

```matlab
addpath('/Users/cag/forclone/monalisa_tutorial_ui/pipeline_utils');
folder='/Volumes/Extreme SSD/261110/recon_results/anisotropic_scanner/retained_20261111/seq712_sigma0p32';
r=load_anisotropic_nifti(folder);  % default: restores original xRMS scale
% r.xrms, r.x_mm, r.y_mm, r.z_mm, r.metadata, r.export
n=load_anisotropic_nifti(folder,false); % normalized single array instead
```

Direct loading: `double(niftiread(file))*jsondecode(fileread(sidecar)).normalizationDivisor`.

NIfTI uses nominal centred sequence x/y/z coordinates, mm spacing, first voxel centre (−120,−120,−90) mm, sform code 2 and qform code 0. Patient anatomical orientation has NOT been established by this export; do not infer verified RAS/LPS labels or register using an assumed scanner affine. Anterior appears at negative sequence y in the current orbital display.

Every normalized NIfTI matched its expected float32 values exactly on MATLAB read-back and an independent Python check against the original MAT. Maximum relative intensity-restoration L2 error was 2.587e-08, below 1e-7. Both loader modes passed for all 12 images. Data arrays were verified before deletion.

## Scanner inputs and reconstruction contract

Dataset: `/Volumes/Extreme SSD/261110`.

| Sequence | Raw file | Support | Samples/readout | Actual TE (ms) |
|---|---|---|---:|---:|
| seq712 | `meas_MID00193_FID15362_seq712.dat` | Box; conservative enclosing angular-design cylinder | 196 | 2.230 |
| seq713 | `meas_MID00194_FID15363_seq713.dat` | Ellipsoid | 120 | 2.035 |
| seq714 | `meas_MID00195_FID15364_seq714.dat` | Cylinder | 152 | 2.145 |

Use the matching `yj_seq71*_*.seq` in that dataset. Scanner B0 measured in the header: 6.980936 T. All three have 199 shots × 44 readouts (8,756), 32 physical receive channels, TR 8 ms, FOV 240 × 240 × 180 mm, spherical k-space boundary and nominal 2-mm isotropic voxels. Spatial support differs; it is not a different reconstruction box.

1. Verify TWIX/Pulseq pairing and grouping via the existing metadata helper.
2. Physical `calculateKspacePP()` ADC coordinates with ROTATIONS; convert cycles/m to cycles/mm once.
3. Discard 20 initial shots and each SI navigator: 179 retained shots × 43 = 7,697 imaging spokes. Match exact data/coordinate readout indices.
4. Recalculate `bmVolumeElement(t,'voronoi_full_radial3')` for that retained trajectory.
5. MONALISA `bmMathilda`, 120 × 120 × 90 grid, `dK=1./[240 240 180]`, Gaussian window 3/parameters `[sigma 10]`, MATLAB IFFT and matching deapodization, normalized gridding enabled.
6. Preserve raw complex phase; no conjugation, second ADC-phase correction, coil compression/sensitivity combination or additional gradient-delay correction.
7. Physical-coil xRMS is `sqrt(mean(abs(x0).^2,4))`.

The low-level engine remains `monalisa_scripts/recon_anisotropic_xrms.m`. Its scientific computation was not changed by cleanup. Direct invocation still saves the full MAT; the batch runner adds the verified compact export and deletes its transient full MAT only after successful export.

## Run new reconstructions

```matlab
addpath('/Users/cag/forclone/monalisa_tutorial_ui');
dataset='/Volumes/Extreme SSD/261110';
sigmas=[.61 .32 .20 .10];
plan=run_anisotropic_scanner(dataset,sigmas,20,[],false,true); % plan only
% plan=run_anisotropic_scanner(dataset,sigmas,20); % starts a NEW 12-job batch
```

Arguments: `(datasetFolder,kernelSigmas,nShotOff,outDir,showImages,checkOnly)`.
Defaults now use the four retained widths, 20 excluded shots, no image windows, and a new `<dataset>/recon_results/anisotropic_scanner/run_<timestamp>` folder. The 7T cell uses the same widths with image display enabled. Input matches must be unique; existing output roots and duplicate widths are refused. Check mode does not load raw signal or create folders.

A newly completed batch saves compact NIfTI/JSON/log/preview and one shared encoding MAT per sequence, with `batch_plan.mat`. `export_anisotropic_nifti` preserves its input MAT by itself; `run_anisotropic_scanner` performs the explicit post-verification deletion. A failed export leaves the MAT intact.

## Comparisons and interpretation

Current four-width common-scale orbital panels: `comparisons/orbital_retained_zm24.png` and `orbital_retained_zm20.png` (also PDF). Bilateral crop x approximately ±55 mm, y −86 to −46 mm; both planes/rows/columns share one intensity scale after restoring divisors.

Older seven-width/six-width/lower-width PNGs and metrics remain as historical evidence in `comparisons`; their deleted-width full arrays are not retained. Historical scanner analysis scripts are under `comparisons/historical_scripts` and refer to the original experiment layout. They are not the current pipeline; use the runner/loader above. Recreating a deleted-width comparison requires reconstructing that width again from raw data.

0.10 looked smoother than 0.20/0.15, and 0.15 was close to 0.20 in the orbital images. These are visual changes, not proof of accuracy. The 0.32 phantom minimum concerns the box-specific interior unit-water RMS metric; it is not a universal scanner optimum or a default selected for all designs. No new simulation/acquisition/full scanner reconstruction was launched during cleanup.

## Static diagnostics and retired paths

Static phantom kernel sweeps moved to `/Volumes/Extreme SSD/261110/recon_results/water_phantom_diagnostics/kernel_sweep_20261111` and `/Volumes/Extreme SSD/261110/recon_results/water_phantom_diagnostics/box_kernel_curve_20261111`. They load shared encoding rather than deleted scanner MATs. Detailed simulation evidence is in the Koma handover.

Removed by explicit user request:

- Scanner widths 0.50, 0.40, 0.27, 0.15 for all three sequences.
- All 24 original full scanner MATs, with the 12 selected magnitudes replaced by verified compact exports.
- Superseded earlier seq712-only Koma run, separate-sphere phantom and cone/ring exports (B/C/D/H review groups).

The old `anisotropic_xrms_20261111_102743_031`, `anisotropic_kernel_sweep_20261111_103727_973` and `box_phantom_kernel_curve_20261111_105211_216` result folders no longer exist. Their retained images/reports/diagnostics were converted or copied and verified before removal. All `.txt` pointers in `results/scanner_anisotropic_261110` now point to the current compact batch or static curve folder. Old manifests retain original source paths as historical provenance.

All scanner raw `.dat`, matching `.seq`, Koma frozen inputs, current Koma batch/fixed phantom and source generators are preserved. Unrelated Halloween/calibration outputs were excluded.

## Storage and verification audit

SSD reconstruction-result file contents decreased from 7.990 GiB to 270.1 MiB at cleanup completion. This includes compact scanner output and retained phantom diagnostics; filesystem allocated space differs. The 12 NIfTI files together occupy 52.2 MiB; encoding, metadata, reports and figures occupy additional space.

Audit: `/Users/cag/forclone/monalisa_tutorial_ui/pipelines/output_review_20261111/cleanup_completed`. Read `CLEANUP_REPORT.md`, `cleanup_result.json`, `authorized_deletion_manifest.json`, `preservation_manifest.json` and validation records for exact scope/hashes. Cleanup matched the user decision; no raw acquisition was deleted. The original pre-cleanup inventory is retained one directory above the audit.
