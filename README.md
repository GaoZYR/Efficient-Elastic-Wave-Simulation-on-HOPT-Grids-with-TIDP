# Efficient Elastic Wave Simulation on HOPT Grids with TIDP

A MATLAB program for two-dimensional elastic-wave forward modeling on adaptive grids.

HOPT solves the two-dimensional isotropic elastic velocity–stress equations on hierarchical patch-based adaptive grids. The included example uses a basin model, a point source, and a surface receiver array. The program supports a free surface, absorbing boundaries, wavefront-based patch activation, and wavefield snapshots.

This README describes the current code and bundled data. Length and time are expressed in meters and seconds; velocity is in m/s, and density is in kg/m³. Source amplitude follows the program's definition; its relationship to a physical seismic moment is not established here.

## Author

**Zhangyongrui Gao**  
Email: [gaozhangyongrui24@mails.ucas.ac.cn](mailto:gaozhangyongrui24@mails.ucas.ac.cn)

## Directory Structure

```text
project_root/
├─ README.md
├─ src/                         MATLAB scripts and functions
│  ├─ main_basin.m              Basin forward-modeling entry point
│  ├─ plot_HOPT.m               Model and AMR grid visualization
│  ├─ Initialize_wavefield_and_plan.m
│  ├─ read_model_parameters.m
│  ├─ cpml_coeff.m
│  ├─ add_mpml.m
│  ├─ add_cpml.m
│  └─ ...
├─ data/
│  ├─ model/
│  │  ├─ vs.bin                 Original Vs model for visualization
│  │  ├─ level_vp/              Binary Vp data for each level
│  │  ├─ level_vs/              Binary Vs data for each level
│  │  ├─ level_rho/             Binary density data for each level
│  │  └─ patch_vs/              Patch metadata and velocity files
│  ├─ patch/                   Per-level patch geometry for plotting
│  └─ T_first_arrive/           First-arrival time files: T_level*.mat
└─ snapshot_basin/              Wavefield snapshots, possibly from earlier runs
```

## Quick Start

1. Keep the relative locations of `src/` and `data/`, including all model and first-arrival files.
2. In MATLAB, set the current folder to the project's `src` directory.
3. Run:

```matlab
main_basin
```

After the simulation, the Command Window displays the total runtime and timings for individual computational stages. Receiver records and the final patch wavefields remain in the MATLAB workspace.

`main_basin.m` starts with `clear` and `close all`. Save any workspace results you want to retain before running it. Runtime depends on the model, number of time steps, and hardware; no fixed runtime is specified here.

To visualize the original model and AMR grids, run separately:

```matlab
plot_HOPT
```

`plot_HOPT.m` also clears the workspace and closes figures. Save forward-modeling results before running it.

## Current Example Configuration

The main parameters are defined in `src/main_basin.m`.

| Parameter | Current Value | Description |
|---|---:|---|
| `dt` | 0.0005 s | Time step shared by all levels |
| `Nt` | 40001 | Number of main-loop iterations; receiver records cover approximately 20 s |
| `N` | 4 | Finite-difference half-width and ghost-layer count, giving eighth-order staggered spatial differences |
| `fm` | 3 Hz | Dominant frequency of the Ricker wavelet |
| `src_amp` | 0.5 | Source amplitude coefficient used by the program |
| `nx00, nz00` | 272, 96 | Base-level full-domain array dimensions for the bundled model |
| `free_surface` | 1 | Enable the top free-surface treatment |
| `wavefront_flag` | 1 | Enable patch activation based on first-arrival times |
| `t_safe_offset` | 0 s | Amount by which patch activation is advanced |
| `force_retained_pml_tags` | 0 | Identify boundary patches across levels rather than only on the finest level |
| `Lx, Lz` | 1000 m, 500 m | Parameters used to calculate absorbing-layer thicknesses, accounting for model offsets |
| `mpml_eta` | 0.05 | Cross-direction damping ratio |
| `nrec` | 280 | Number of surface receivers, with 50 m horizontal spacing |

The bundled data contain four levels, numbered 0–3. The original model used for visualization has 161×111 samples, with horizontal and vertical spacings of 100 m and 50 m. Exact source and receiver coordinates, including staggered-grid offsets, are specified by `src_x/src_z` and `rec_x/rec_z` in the main script.

## Outputs and Plotting

### Receiver Records

The main script leaves the following results in the workspace. **Receiver records are not automatically saved to a MAT file.**

- `record_Vx`, `record_Vz`: horizontal and vertical velocity receiver records. Both are currently 40000×280 arrays, with rows representing time and columns representing receivers.
- `Vx_patches`, `Vz_patches`: final velocity fields stored by level and patch.
- `Txx_patches`, `Tzz_patches`, `Txz_patches`: final stress fields.
- `time_stat`, `total_time`: computational-stage and total runtimes.

After the simulation, and before clearing the workspace, plot the vertical-component shot record with:

```matlab
t = (1:size(record_Vz, 1)) * dt;
figure;
imagesc(rec_x, t, record_Vz);
xlabel('Distance (m)');
ylabel('Time (s)');
colorbar;
```

To save the receiver records, execute the following in the workspace, for example:

```matlab
out_file = fullfile(base_dir, '..', 'record_basin.mat');
save(out_file, 'record_Vx', 'record_Vz', 'rec_x', 'rec_z', ...
    'src_x', 'src_z', 'dt', 'fm', 'mpml_eta', 'wavefront_flag', '-v7.3');
```

### Wavefield Snapshots

The entire `Save wavefield snapshots` block is currently commented out, so **new snapshots are not written by default**. Existing files in `snapshot_basin/` may belong to earlier runs.

To enable snapshots, uncomment that block. The existing settings cover 0–10 s, with a snapshot every 0.02 s. Filenames follow the pattern `snapshot_000040.mat`. Each file contains `Vx_gather`, `Vz_gather`, the snapshot time, the iteration number, and grid-positioning parameters.

Saving a snapshot overwrites an existing file with the same name. Back up earlier results or select a new snapshot directory if you need to retain them.

## Main Scripts and Functions

| File | Purpose |
|---|---|
| `main_basin.m` | Parameter setup, initialization, time stepping, receiver recording, and output control |
| `levels_patches_info.m` | Read level counts, patch lists, and geometric metadata |
| `Initialize_wavefield_and_plan.m` | Initialize wavefields and build copy, interpolation, and restriction plans |
| `apply_copy_plan_*.m` | Copy data between patches on the same level |
| `apply_interp_plan_*.m` | Supply required data from coarser to finer levels |
| `apply_restriction_plan_*.m` | Transfer data from finer to coarser levels |
| `read_model_parameters.m` | Read Vp, Vs, and density; construct elastic parameters and staggered averages |
| `build_source_receivers.m` | Locate the source and receivers; build the wavelet, source distribution, and receiver interpolation weights |
| `find_wavefront_active_patches.m` | Determine patch activation steps from first-arrival times |
| `cpml_coeff.m` | Generate the original CPML coefficients |
| `add_mpml.m` | Add cross-direction damping and update absorbing-layer memory variables and derivatives |
| `add_cpml.m` | Apply the original CPML update; called when `mpml_eta=0` |
| `make_memory.m` | Allocate memory variables for absorbing-layer patches |
| `record_receivers.m` | Interpolate and record receiver wavefields |
| `gather_wavefield_patch.m` | Gather patch wavefields onto the visualization grid |
| `plot_HOPT.m` | Display the model and hierarchical grid layout |
