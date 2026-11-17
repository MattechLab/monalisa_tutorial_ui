# Completed cleanup — 11 November 2026

User decision: remove review groups B/C/D/H; keep only 0.61 baseline, 0.32, 0.20 and 0.10 for seq712–714; replace retained scanner image MATs with normalized compressed NIfTI; update both pipelines and handovers.

## Completed

- 12 normalized float32 NIfTI scanner images retained; all 12 verified exactly against normalized source arrays before deletion.
- Saved per-volume normalization divisors; maximum relative raw-intensity restoration error 2.587e-08.
- 24 full scanner MATs removed: 12 replaced selected images and 12 unwanted-width results.
- B/C/D/H outputs deleted. Selected Koma batch, fixed phantom, static diagnostics, raw inputs and source generators preserved.
- Shared trajectory/DCF stored once per sequence; no full coil-image arrays remain in the compact batch.
- Comparison reports/figures and static diagnostic moves were SHA-256 verified. Old report source paths remain historical provenance.
- Updated MATLAB batch defaults, NIfTI exporter/loader, 7T entry cell, Koma legacy defaults, diagnostic encoding paths, result pointers and both handovers.

## Current locations

- [Scanner handover](/Users/cag/forclone/monalisa_tutorial_ui/pipelines/anisotropic_scanner/HANDOVER_ANISOTROPIC_RECON_20261111.md)
- [Koma handover](/Users/cag/forclone/pulseq4mreye_rectfov/simulation/koma_water/HANDOVER_KOMA_WATER_20261111.md)
- Scanner compact batch: `/Volumes/Extreme SSD/261110/recon_results/anisotropic_scanner/retained_20261111`.
- Static scanner-encoding phantom studies: `/Volumes/Extreme SSD/261110/recon_results/water_phantom_diagnostics`.
- Selected Koma batch: `/Users/cag/forclone/pulseq4mreye_rectfov/simulation/koma_water/output/seq712_714`.

## Measured storage

| Footprint | Before | After |
|---|---:|---:|
| SSD reconstruction results, file contents | 7.990 GiB | 270.1 MiB |
| SSD reconstruction results, allocated blocks | 8.319 GiB | 549.0 MiB |
| 12 normalized NIfTI files alone | — | 52.2 MiB |

Counts/sizes are the measured cleanup snapshot; later audit/index files add a small overhead. Both normalization and conversion preserve spatial sampling; NIfTI geometry is explicitly nominal sequence coordinates, not validated patient anatomy.

`authorized_deletion_manifest.json` lists removed files. `preservation_manifest.json` records original/copy hashes and independent NIfTI verification. `pipeline_validation.json` records the MATLAB tests. No new full simulation or reconstruction was run. New compact output does not retain coil phase; regenerate from protected raw data if needed.
