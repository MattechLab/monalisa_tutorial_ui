# MONALISA tutorial UI

**Isotropic seq701–705 scanner xRMS:** see [pipelines/isotropic_scanner/README.md](pipelines/isotropic_scanner/README.md). The independent `run_isotropic_scanner.m` entry point saves verified normalized 240³ NIfTI at 1-mm isotropic spacing.

**Anisotropic scanner pipeline:** see [pipelines/anisotropic_scanner/README.md](pipelines/anisotropic_scanner/README.md) and the independent `run_anisotropic_scanner.m` entry point. Koma water simulation is maintained separately in `/Users/cag/forclone/pulseq4mreye_rectfov/simulation/koma_water/KOMA_PIPELINE.md`.

A MATLAB tutorial for estimating coil sensitivities from body-coil (BC) and
head-coil (HC) prescans, with optional data preparation and reconstruction.
Run **main_script_3TPrismaXA60.m cell by cell**. Running the entire script stops after saving C.
For a **quick raw-data preview**, run only its last cell, **Quick raw-data
preview ONLY**, from a fresh workspace. It sets up dependencies and asks for
the dataset folder and matching main `.seq`/`.dat` pair. It needs no task
setup, prescans, or C estimation and creates no task output folder. Edit
`quickPreviewMatrix` and `quickPreviewShotOff` in that cell as needed;
missing sequence metadata can still prompt for input. Results use
`quickPreview*` variables so an existing task's images and provenance remain intact.

For **7T Terra VE12U**, use the independent `main_script_7TTerraVE12U.m`
and edit `tutorial_config_7TTerraVE12U.m`. This entry estimates relative C
from the main scan, without BC/HC prescans. It reconstructs complex 48³ coil
images, asks you to review a signal-support mask, smooths the coil images,
normalizes by their root-sum-of-squares, and uses MONALISA's Laplace solver
to extend the maps outside the mask. The saved complex `C` has unit RSS
across coils and can be consumed by MONALISA. It retains a common object
phase and has no body-coil intensity reference. Real-scan quality still
requires validation; numerical format compatibility alone does not prove it.

The supplied `meas_MID00014_FID03945_cfa1.dat` header reports Terra,
**syngo MR E12**, baseline **N4_VE12U_LATEST_20181126**, and **6.980936 T**.
Its accompanying sequence has `nSeg 44` and `nShot 1055` definitions added
after signing. Removing only those two lines in memory reproduces the exact
stored/TWIX signature. `check_hash` recognizes and reports this case; other
unverified content changes are rejected. The data loader checks grouping
against the acquired readout count. The script offers
main-scan C estimation/save, optional ROVir weight estimation, Mitosius export,
and gridded or compressed-sensing reconstruction. The FOV comes from the
resolved sequence metadata, not an assumed 240 mm.

In the 7T **Prepare main data for reconstruction (Mitosius)** cell, edit
`useROVir` to enable/disable compression and **`nVirtualCoils`** to set the
retained count. `computeROVirFromEyeROI` still saves the complete basis;
the count is applied to raw data and C during Mitosius preparation. The
normalization-grid default is explicitly 120, and the initial-shot exclusion
is reused from calibration. Both appear in the preparation settings dialog.
Select the same raw scan/sequence and its physical-coil C file. Reconstruction
automatically uses the returned Mitosius/C pair and its FOV metadata. Its
dialog sets the final reconstruction grid, iteration count, delta, and rho.
The reconstruction cell uses separate output variables so that the original
coil images remain available for ROVir review. Full real-scan 7T preparation
and reconstruction still require validation on the intended inputs.

For context, the [BART calibration example](https://mrirecon.github.io/bart/examples.html)
compares direct calibration with ESPIRiT and demonstrates a projection-error
check. The 7T implementation here uses direct relative calibration, not ESPIRiT,
and introduces no BART dependency.

## Install dependencies once

- MATLAB with Image Processing Toolbox and Statistics and Machine Learning Toolbox.
- [MONALISA](https://mattechlab.github.io/monalisa/), including its Siemens reader
  in `third_part/twix_for_monalisa` and Pulseq trajectory support.
- [Pulseq](https://github.com/pulseq/pulseq).
- Build MONALISA's MEX files for your MATLAB release and platform, following its
  [installation guide](https://mattechlab.github.io/monalisa/1-2_installation.html).
  If needed, select a compiler with `mex -setup C++`, then run MONALISA's
  `compile_mex_for_monalisa`. Restart MATLAB after building.

Separate mapVBVD installations and personal paths are not required by this
tutorial: it uses MONALISA's bundled `mapVBVD_JH_for_monalisa` consistently.

## Configure a different 3T task in the same dataset

Edit **`tutorial_task_config_3T.m`** before running setup cell 1. This is the
single place for 3T input files and processing settings. Keep installation
paths and the dataset root in `tutorial_config.m`.
Both 3T and 7T configurations leave `baseFolder` empty by default. Setup asks
you to choose the dataset folder on each new run, unless you explicitly set
its path in the configuration. No previous dataset choice is restored.

| What to change | Configuration field |
| --- | --- |
| Dataset folder | `config.baseFolder` in `tutorial_config.m` |
| Output parent folder | `config.outputRoot` in `tutorial_config.m` |
| Task label | `task.name` (letters, digits, underscore, hyphen) |
| Main data and matching sequence | `task.files.mainRaw`, `task.files.mainSeq` |
| BC/HC inputs | `task.files.prescanSeq`, `bodyRaw`, `headRaw` |
| Default calibration dialog choice | `task.calibrationMethod`: `'prescan'`, `'main'`, `'reuse'` |
| Existing physical-coil C | `task.existingCFile` for `'reuse'`; empty opens a picker |
| Previously approved calibration mask | `task.maskFile` |
| Initial shot exclusions | `task.mainShotOff = 50`, `task.prescanShotOff = 14` |
| Main coil preview grid | `task.previewMatrix = 60` |
| Main-scan C estimator parameters | `task.mainCalibration` |
| BC/HC refinement iterations | `task.prescanIterations = 5` |
| Optional orientation score gates | `task.orientation` (empty disables each gate) |
| ROVir compression and retained channels | `task.useROVir`, `task.nVirtualCoils` |
| Normalization grid | `task.normalizationMatrix`: 60 or 120 (default 120) |
| Final reconstruction grid and CS settings | `task.reconstruction` |
| Reconstruct an existing export | Both `task.preparedMitosiusDir` and `task.preparedCFile` |

Input paths may be absolute or relative to the dataset root. Empty required
scan paths open file pickers once during setup; the resolved files are then
reused by calibration, inspection, ROVir, and preparation. For BC/HC mode the
prescan inputs are also required; main/reuse modes do not request them.
Every native file/folder picker has a separate instruction popup on macOS.

For example, use `task.name = 'static_fullcoil'` with `task.useROVir = false`
for a full-coil task. For another task on the same dataset, change the name
to `'motion_rovir25'`, set its main inputs, and enable 25 virtual coils.
Rerun setup cell 1 to clear previous variables and create a new task folder.
Do not change task inputs midway through a run: start again from cell 1.

Each setup creates:

```text
<dataset>/recon_results/<task.name>_<timestamp>/
    task_config.mat                 resolved inputs and settings
    <MID>_recon_C_<timestamp>/       BC/HC calibration (or main/reused_C)
        C/C.mat                     physical-coil C and task/main provenance
        coilCom/                    reviewed ROVir weights and coil images
    <MID>_recon/
        T1_LIBRE_woBinning/mitosius/ prepared data, matching C, preparation record
    reconstructions/run_<timestamp>/
        run_info.mat                exact settings, C hash, timing, completion
        x0_<matrix>.mat
        x_nIter...mat                when CS is enabled
```

Absolute output overrides remain supported. Completed or partial Mitosius
folders are never overwritten. Each reconstruction gets a new run folder,
so rounding of parameters in a filename cannot overwrite another run.
Preparation writes `complete.mat` only after all data and provenance have
been saved. Reconstruction's `run_info.mat` is marked complete only after
outputs have been saved.

## Run the 3T cells

1. **Setup a new task.** Edit the task configuration and run cell 1. First
   choose **BC/HC prescans**, **Main scan only**, or **Reuse saved C**. Only
   BC/HC mode requests prescan inputs. Inspect the printed settings and
   output path; `task_config.mat` saves the chosen method and resolved inputs.
2. **Estimate/load C.** `'prescan'` retains the BC reference, HC primary maps,
   and refinement estimator. Review BC/HC/main previews and accept the final
   mask. `'main'` uses the existing 7T direct smoothed-RSS/Laplace estimator
   with a reviewed 48³ support mask; it is not ESPIRiT and needs no additional
   dependency. `'reuse'` loads physical-coil C; a legacy C without provenance
   requires explicit review of coil order and spatial correspondence.
3. **Save C.** Writes `C/C.mat` with task and acquisition provenance. Whole
   script execution stops here; run optional cells individually.
4. **Inspect main data.** Required before ROVir in prescan mode; main/reuse
   calibration already returns physical-coil images. This cell also runs
   immediately after setup, before C estimation. It uses the task's main
   inputs, resolving sequence metadata if calibration has not run. After
   calibration, it reuses that metadata and checks the acquisition identity.
5. **ROVir.** Draw eye ROI, choose its z range, and review ROI/IR overlays.
   Rerunning archives the earlier coilCom folder only after the new estimate
   succeeds. Weights save the scan identity, shot exclusion, and source C hash.
6. **Prepare.** Uses the selected main scan and saved physical-coil C directly.
   The existing single-sample decision and complex central-slice normalization
   are retained. The grid is small and independent of the final reconstruction
   grid; the value/decision are saved and a zero/nonfinite divisor is rejected.
   ROVir applies the same QR/conjugated mixing to raw data and C. Preparation
   verifies task provenance and returns its exact Mitosius/C pair.
7. **Reconstruct.** Uses `task.reconstruction` directly. Optionally set both
   prepared paths to reconstruct an existing export. New outputs still go
   under the current task. Reconstruction checks the exact paired C hash
   when the export has a preparation record. Coil images retain separate
   variable names and remain available for ROVir review.

The axes/flip check saves identity/best correlations and the best-to-runner-up
margin. Optional score gates default to disabled because useful thresholds
need acquisition-specific validation. It is not rigid head-motion registration.
HC/main coil-selection metadata are compared when available; missing metadata
is reported. Header equality does not establish channel ordering at every
raw sample. Raw identity uses path, size and modification time, sequence hash,
and acquisition dimensions; it does not hash the entire raw file.

Calibration logs/timings distinguish loading, gridding/review, reference,
primary, and refinement stages. Underlying MONALISA Laplace solvers still use
convergence loops: 1000 is an iteration batch, not a hard total limit. They
have not been changed here.

For standalone legacy 7T/UI calls, pickers and parameter dialogs remain
available. Legacy C/ROVir/export files without provenance are reported as
unverified rather than silently treated as task-verified. The 3T task workflow
requires provenance at preparation.

## Common sequence-header format and missing values

`read_sequence_header` reads `[DEFINITIONS]` into a common structure using
Yiwei's field names and order, including `FOV`, `TR_s`, `TE_s`, `nSam`,
`nSeg`, `nShot`, `flagSelfNav`, and `ReadoutOversamplingFactor`. Optional
fields such as RF-generation parameters remain `[]` when unavailable.
Unknown source definitions are retained. Original signed `.seq` files are
never rewritten or given invented RF/timing values.

The active raw-data loader requires `FOV`, `ReadoutOversamplingFactor`,
`nSam`, `nSeg`, `nShot`, and `flagSelfNav`. `extract_seq_params` opens a text
input dialog for missing/invalid required fields, retries invalid entries,
and stops cleanly on cancellation. The raw sample count prefills missing
`nSam` for confirmation; SI presence and physical grouping are not guessed
from the filename or Siemens `NSeg`. Resolved values are checked against the
actual raw sample/readout counts. Complete headers need no missing-field
dialog; Mitosius still shows the prefilled shot/segment confirmation above.

| Quantity | Common units / meaning |
| --- | --- |
| `FOV` | Metres, before the declared readout oversampling factor |
| `TR_s`, `TE_s` | Seconds; preserved as supplied |
| `nSam` | Acquired imaging samples/readout, including oversampling |
| `nSeg`, `nShot` | Physical readouts/shot and total shots, before exclusions |
| `flagSelfNav` | 1: first readout/shot is SI; 0: no SI readout |
| `ReadoutOversamplingFactor` | Explicit acquisition factor, e.g. 1 or 2 |

Reconstruction FOV in mm is `1000 * FOV * ReadoutOversamplingFactor`.
The two supplied examples therefore give 200 mm (100 mm × 2) and 120 mm
(120 mm × 1). The current radial trajectory path supports isotropic FOV.
Mitosius stores this resolved geometry; final reconstruction reads it back.
An older export without geometry prompts for FOV instead of assuming 240 mm.

Normalized headers and their source/user-entry provenance are saved beside
C as `main_sequence_header.mat` (plus `prescan_sequence_header.mat` for 3T).
Mitosius saves `sequence_header.mat` and `recon_geometry.mat` alongside its data.

## macOS OpenMP setup

Do **not** distribute the diagnostic or generated MEX copies. MEX binaries are
specific to the platform and MATLAB installation.

On macOS, setup copies the two OpenMP MEX files used by sensitivity estimation
from your own MONALISA build into `runtime_cache/mex/<architecture>/<release>`. It links
these local copies to the `libomp.dylib` bundled with the running MATLAB, signs
the changed copies, verifies their dependencies with `otool -L`, and gives them
path priority. MONALISA's installed files are not modified.

See [runtime_cache/README.md](runtime_cache/README.md) for the exact original
MEX locations, generated-file provenance, and rebuild instructions. Setup
prints the full source and cache paths. If the cache is absent, run cell 1
or `config = setup_tutorial();` to recreate it automatically from the compiled
MONALISA originals. If those originals are missing, build MONALISA first.

The cache is refreshed if the source MEX, MATLAB installation/version, or runtime
library changes. Setup refuses to continue if a MONALISA MEX has already loaded
a different OpenMP runtime: restart MATLAB and run the setup cell first.
Rebuilding MONALISA later is fine; restart MATLAB and rerun setup afterwards.

This route requires Apple's `otool`, `install_name_tool`, and `codesign` tools
and a MATLAB installation that bundles `libomp.dylib`. It was validated on
Apple silicon with R2025a Update 1. Windows and Linux use their native MONALISA
MEX builds without this macOS-specific relinking step; those platforms were
not executed during this cleanup.

## Visible settings and repeating a mask

Edit `tutorial_config.m` and rerun cell 1 to change the default paths.
The 3T processing options and input paths are in tutorial_task_config_3T.m.
The task entry prints and saves these settings, then uses them directly.
Standalone UI calls retain their parameter dialogs. Setup never writes settings back or restores choices from a MAT file.

Both configurations default to `config.outputRoot = 'recon_results'`. This
relative path is resolved under the selected raw-data folder, so code and
scan results are separate. Selecting another dataset during setup updates
the destination automatically. An explicit absolute `outputRoot` is retained.
Mitosius also follows the dataset selected in its dialog when using the default.
Its final reconstruction outputs stay under that prepared dataset's directory.
The generated runtime cache remains with the code. Existing results from
earlier runs have not been moved.

For a one-run path override, edit `config` between the two setup calls in cell 1,
or run the following directly (rerunning cell 1 reloads the file's defaults):

```matlab
config = tutorial_config();
config.baseFolder = '/path/to/dataset';
config = setup_tutorial(config);
```

To reuse a reviewed mask, set `coilOptions.maskFile` to a MAT file containing
the logical 48³ variable `mask`. This avoids making new mask choices.
For noninteractive use, `coilsense_script` also accepts `coilOptions.files`
with `prescanSeq`, `bodyRaw`, `headRaw`, `mainSeq`, and `mainRaw` full paths,
and `coilOptions.showFigures=false` when a saved mask is supplied.

## Repository contents and sharing

```text
main_script_3TPrismaXA60.m  cell-by-cell entry point
tutorial_task_config_3T.m  explicit per-task inputs and processing settings
main_script_7TTerraVE12U.m  independent main-scan calibration entry for 7T
tutorial_config.m    visible, editable default paths (no saved overrides)
tutorial_config_7TTerraVE12U.m  independent default paths for 7T
setup_tutorial.m     dependency and runtime setup; prints active settings
monalisa_scripts/    user-facing calibration/reconstruction functions
pipeline_utils/     runtime, bounded data loading, and orientation helpers
seq_utils/          sequence parsing and signature checking
tests/              runtime and real-data validation functions
runtime_cache/      MEX provenance/rebuild README; generated mex/ is ignored
<baseFolder>/recon_results/  scan results, saved beside the raw data
```

Push the source, visible default configuration, documentation, and tests.
Use portable paths in the shared configuration. `.gitignore` excludes `runtime_cache/mex`,
legacy `outputs` folders, raw data, MAT results, and compiled MEX binaries. Do not add the
ignored files with `git add -f`. No installation-specific binary copies are
needed in the shared repository.

## Validation

`validate_task_configuration()` checks explicit/relative paths, invalid task
labels and unpaired exports, main-only input selection, timestamped calibration
discovery, identity-versus-flip orientation, and SI-aware acquisition identity.
`validate_mitosius_reconstruction` also checks separate run folders, saved
completion/settings, and rejection of a different C with identical dimensions.

`validate_siemens_sequence_loading(config)` uses synthetic TWIX metadata with
`NSeg=1` to check automatic sequence grouping and prefilled confirmation,
rejection before sample reads on mismatch, and exact blocked-read equivalence
to the original filtering/reshape operation. It does not open real dialogs.

`validate_sequence_headers()` uses temporary fixtures to test canonical names,
required-field input/retry/cancellation, raw-count checks, SI selection, and
FOV scaling with oversampling factors 1 and 2. It does not open real dialogs
or modify source sequences.

`validate_sequence_hash(config)` checks the hash boundary against a sequence
signed by the installed Pulseq writer, verifies the native reader sees the
same stored signature, and checks rejection after changing the file content.
Hash diagnostics print TWIX, stored `[SIGNATURE]` Hash, and recomputed MD5
separately; matching stored hashes alone do not establish current file integrity.
The test also checks that added integer `nSeg`/`nShot` definitions are accepted
only when removing them restores the original signed payload exactly.

After setup, `addpath('tests'); validate_runtime_setup(config)` checks actual
parallel sparse multiplication against MATLAB, tests a constant-boundary
Laplace solution, and checks that setup reuses an unchanged cache.

`validate_rovir_compression(config)` checks the complex mixing convention
against the original HPC formula, consistency of compressed data and C,
full-rank image recovery, invalid inputs, and a two-bin Mitosius save/load
round-trip using synthetic data.

`validate_mitosius_reconstruction(config)` exports synthetic data with physical
coils and with ROVir compression, then runs the actual gridded and bmSteva
compressed-sensing paths on a small grid and verifies their saved outputs.
`recon_nobin_from_paths` provides the shared reconstruction core without UI;
`recon_nobin_ui` collects settings and calls that same core.

`validate_main_scan_sensitivities(config)` checks direct relative calibration,
Laplace extension, unit RSS, Gaussian smoothing, and combination through the
actual MONALISA `bmCoilSense_pinv` function with synthetic coil images.

`validate_calibration(config,files,maskFile,referenceFile)` runs the production
calibration on real data with an approved mask and compares C with a previous
matching `C.mat`. It writes a separate result directory and does not overwrite
the reference. Tests require your installed dependencies; acquisition files
are not distributed with this repository.

See [MONALISA documentation](https://mattechlab.github.io/monalisa/) for the
reconstruction methods and dependency build instructions.
