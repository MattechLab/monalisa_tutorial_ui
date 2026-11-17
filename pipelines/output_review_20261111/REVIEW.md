> Cleanup authorized and completed in the subsequent pass: B/C/D/H removed; A reduced to 0.61/0.32/0.20/0.10 and compact NIfTI. See [cleanup_completed/CLEANUP_REPORT.md](cleanup_completed/CLEANUP_REPORT.md). The inventory below is the historical pre-cleanup snapshot.

# Pipeline separation and output retention review

Prepared 11 November 2026. No existing scientific output was moved, deleted or rewritten. No new Bloch simulation or scanner reconstruction was run.

## Separate entry points

- [Koma water pipeline](/Users/cag/forclone/pulseq4mreye_rectfov/simulation/koma_water/KOMA_PIPELINE.md): `run_water_batch.jl`, then `reconstruct_water_batch.m`.
- [Anisotropic scanner pipeline](/Users/cag/forclone/monalisa_tutorial_ui/pipelines/anisotropic_scanner/README.md): `run_anisotropic_scanner.m`; the 7T entry cell delegates to it.
- [Static diagnostics branch](/Users/cag/forclone/pulseq4mreye_rectfov/simulation/koma_water/diagnostics/README.md): analytic Fourier phantom, kernel studies and direct-sum/PSF controls, distinct from Bloch and scanner data.

Existing runners and their scientific settings are preserved. The 7T entry-cell sigma stays 0.50; the simulation reconstruction sigma stays 0.61. A selected 0.32 box-phantom metric optimum is not installed as a universal scanner default.

## Retention decisions

| ID | Outputs | Content size | Allocated space | Recommendation |
|---|---|---:|---:|---|
| A | Requested scanner widths: 0.61, 0.50, 0.40, 0.32, 0.20, 0.15, 0.10; seq712-714; 21 full volumes | 7067.7 MiB | 7227.0 MiB | KEEP |
| B | Additional scanner width 0.27; seq712-714; three full volumes | 1009.8 MiB | 1032.0 MiB | USER CHOICE; outside latest seven-column comparison |
| C | Earlier seq712-only Koma run; superseded by selected three-sequence batch; not byte-identical | 156.0 MiB | 156.0 MiB | DISCARD CANDIDATE; update legacy default paths first |
| D | Earlier separate-sphere analytic phantom, superseded by fixed body+A+B+cavity phantom | 41.6 MiB | 41.6 MiB | DISCARD CANDIDATE; retain canonical fixed phantom |
| E | Selected three-sequence Koma batch | 404.9 MiB | 405.0 MiB | KEEP |
| F | Current fixed phantom, 199-vs-987 shots and Cartesian reference | 78.6 MiB | 78.7 MiB | KEEP |
| G | Other static diagnostic studies, kernel curves, trajectory previews, numerical checks | 602.5 MiB | 667.1 MiB | KEEP; optional checkpoint slimming only on a separate request |
| H | Earlier rectangular/cone-ring sequence prototype exports | 57.6 MiB | 57.6 MiB | ARCHIVE or KEEP; different earlier design |

## Proposed organization after retention choices

The two repositories remain separate. Source sequence exports and input snapshots stay at their current paths. The original `pulseq4mreye` repository and synced `sources/` are excluded.

New scanner batches already target `<dataset>/recon_results/anisotropic_scanner/run_<timestamp>`. Koma batches target `simulation/koma_water/output/<new_batch>`. Static diagnostics are indexed under `simulation/koma_water/diagnostics`.

Historical scanner baseline/sweep folders remain in place pending decisions. Their mixed `phantom_*.mat` and phantom scripts/metrics belong to the simulation diagnostic branch. A physical move would require repairing historical script paths and result pointers, then verifying copies before removing old locations.

Deleting group C requires changing the legacy default input folders in `inspect_water_acquisition` and `reconstruct_water_xrms`, plus their old README examples, to the retained selected batch. Deleting B removes saved 0.27 volumes used by older plots; those plots/metrics can remain, but recreating them would require reconstructing 0.27 again from protected raw data.

Raw scanner `.dat`, matching `.seq`, frozen Koma input snapshots, receipts, dependency manifests, source code, validation records, and current final comparison figures are recommended KEEP. No exact-duplicate claim is made for the two seq712 Bloch runs.

## Inventory scope and checks

- Checked 1103 files across both repositories' output roots, scanner results, local workspace outputs, the task visualization root and matching temporary folders.
- Other experiment outputs are listed as X groups, including Halloween data and older coil-calibration outputs. They are outside this cleanup.
- Exact small-file duplicate checks include selected reports/scripts/sequence exports below 12 MiB, excluding Apple resource-fork files. Large acquisitions and MAT/NPZ arrays were not exhaustively hashed.
- Both new entry points passed plan checks without creating output folders; duplicate scanner widths, existing output folders and missing simulation inputs were refused.
- Every inventoried existing file retained its recorded size and modification time. Reconstruction engines, phantom definition, the 3T entry script and earlier 7T cells were verified unchanged.

[All file paths and sizes](output_files.csv) · [Grouped inventory and duplicate hashes](output_inventory.json) · [Exact paths for choices A–H](retention_choices.json) · [Validation record](separation_validation.json)

Reply with KEEP / DISCARD / ARCHIVE for B, C, D and H, plus any changes to the recommended KEEP groups. No deletion will occur before those decisions.
