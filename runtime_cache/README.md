# Generated MONALISA MEX cache

This visible directory contains generated binaries, not user settings.
Edit **tutorial_config.m** in the repository root to change paths, including
`config.mexCacheRoot`. Edit processing options in each main_script_3TPrismaXA60 cell.

## Where the cache comes from

On macOS, `setup_tutorial` calls `pipeline_utils/setup_monalisa_runtime.m`.
It copies these files from your compiled MONALISA checkout (`config.monalisaRoot`):

| Original file under MONALISA | Generated copy under the default cache |
| --- | --- |
| `src/sparseMat/mex/bmSparseMat_vec/bmSparseMat_cC_oBlock_omp/bmSparseMat_cC_oBlock_omp_mex.<mexext>` | `mex/<architecture>/<release>/bmSparseMat_cC_oBlock_omp_mex.<mexext>` |
| `src/image3/mex/bmImLaplaceEquationSolver3_omp/bmImLaplaceEquationSolver3_omp_mex.<mexext>` | `mex/<architecture>/<release>/bmImLaplaceEquationSolver3_omp_mex.<mexext>` |

For example, Apple silicon R2025a uses `mex/maca64/2025a/` and `.mexmaca64`.
The original binaries remain unchanged. Each cache copy is linked to
`<matlabroot>/bin/<architecture>/libomp.dylib`, signed with `codesign`, and
checked with `otool -L`. Setup prints the exact source, cache, and runtime paths.
Files ending in `_build.mat` record the source/runtime hashes and MATLAB version;
they are generated provenance, not editable configuration.

Windows and Linux use MONALISA's original native binaries without this cache.

## Recreate a missing cache

1. Restart MATLAB if MEX files from this cache have already been loaded.
2. Open the tutorial repository as MATLAB's Current Folder.
3. Run **main_script_3TPrismaXA60 cell 1**, or run:

   ```matlab
   config = setup_tutorial();
   ```

Missing folders, binary copies, and build records are recreated automatically.
There is no separate cache compilation step. Changes to the original binaries,
MATLAB installation/version, or OpenMP runtime also trigger a cache refresh.

If the **original MONALISA MEX files are missing**, first build them in your
MONALISA checkout with `compile_mex_for_monalisa` (configure the C++ compiler
with `mex -setup C++` if needed). Restart MATLAB, then rerun tutorial setup.
The cache cannot be recreated without these compiled originals.

To force a fresh cache, close MATLAB, remove only the generated `mex/` folder,
and rerun setup. Keep this README. Do not commit or distribute `mex/` files;
they depend on the local MATLAB installation and platform.
