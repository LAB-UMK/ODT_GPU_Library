# Usage Guide  
ODT GPU Library — Compact Workflow Overview

This document provides a concise description of the two supported processing workflows. It is intended as a quick, practical reference for users integrating the library into their own software.

---

## 1. Preprocessed-Data Workflow (from sinograms)

The first and most basic mode is used in the post-processing application. It operates on externally prepared sinograms (amplitude and phase) and provides access to intermediate reconstruction data.

### Typical Call Sequence

```c
// Set reconstruction and sampling parameters
HL01_setParams(K_xy, K_z, dx, n_imm, n_proj, Nx, Ny,
               dkP, dkPz, NA, lambda_all, kxp, kyp,
               approxBornNotRytov);

// Upload sinograms and generate K-space
HL02_sendDataAndGenerateKO(sinoAmp, sinoPh, FpmaskLogical);

// Tomographic reconstruction (DI or GP)
HL03_setParamsAndStartDIandGP(nGPi, do_TC, do_NNC, relaxM, betaM, kn_mean,
                             objShift, objSupport);

// Retrieve reconstruction and free GPU memory
HL04_takeReconstructionAndFreeMemory(n_rec);

```
Each exported function is documented directly in the header file (ODT_GPU.h).

### Key Notes
- Input data: sinogram amplitude and phase (`float`).
- This workflow does not support reference holograms and operates on already preprocessed sinograms.
- This workflow allows **retrieval of intermediate data** such as K-space.
- Recommended for **offline processing, debugging, and validation**.
---


## 2. Raw-Data Workflow (from holograms)

This workflow starts directly from raw detector holograms and performs preprocessing, K-space generation, and reconstruction on the GPU.

### Typical Call Sequence

```c
// (Optional) Upload reference hologram
HL_addReference(hologram_ref, X, Y, nproj, NA, lambda, cam_pix, M, n_imm, do_NNC, fftWindowScale);

// Raw holograms -> preprocessing -> K-space
HL00to02_FromPreprocToGenKO(hologram, X, Y, nproj, NA, lambda, cam_pix, M, n_imm,
                           do_NNC, &Kxy_dim, Kspace_oversampling_z,
                           cosFactor, fftWindowScale, approxBornNotRytov);

// Tomographic reconstruction (DI or GP)
HL03_setParamsAndStartDIandGP(nGPi, do_TC, do_NNC, relaxM, betaM, kn_mean,
                             objShift, objSupport);

// Retrieve reconstruction and free GPU memory
HL04_takeReconstructionAndFreeMemory(n_rec);
```

Each exported function is documented directly in the header file (ODT_GPU.h).


### Key Notes
- Input data: raw holograms (`int16` in MATLAB), dimensions `X × Y × nproj`.
- This workflow supports optional reference holograms for background compensation.
- `nGPi = 0` selects **Direct Inverse (DI)** reconstruction.
- `nGPi > 0` selects **Gerchberg–Papoulis (GP)** iterative reconstruction.

---

## 3. Direct Inverse vs. Gerchberg–Papoulis

- **Direct Inverse (DI):** `nGPi = 0`
  - Fastest reconstruction.
  - Suitable for real-time preview.

- **Gerchberg–Papoulis (GP):** `nGPi > 0`
  - Iterative reconstruction.
  - Higher quality at the cost of longer computation time.

---

## 4. MATLAB Integration Notes

- Use `loadlibrary()` and `calllib()` to access the DLL.
- Convert data to C-compatible types (`int16`, `single`).
- Permute arrays if necessary to match the expected C memory layout.
- Example MATLAB script is provided with the library.

---

## 5. API Reference and Additional Functions

Two versions of the library header are provided:

| File | Contents | Intended use |
|------|----------|--------------|
| `include/ODT_GPU.h` | Core workflow functions (HL00–HL04), optional diagnostic functions, and utility functions | Working with the precompiled DLL |
| `src/ODT_GPU/ODT_GPU.h` | The above **plus** the legacy low-level API | Building from source; integrating with older code |

Both headers document every function using Doxygen-style comments, including parameter descriptions and return codes.

The functions are grouped into the following categories:

- **Core workflow** — the main reconstruction pipeline (`HL00to02_*`, `HL01_*`, `HL02_sendDataAndGenerateKO`, `HL03_*`, `HL04_*`), plus reference-hologram management (`HL_addReference`, `HL_removeReference`) and their backward-compatible `HL00_*` aliases.
- **Optional / diagnostic** — retrieval of intermediate results such as the generated K-space (`HL02_B_optionTakeKO`), the Ewald sphere occupancy map (`HL02_B_optionTakeEW`), and reference sinograms (`HL00_B_optionTakeSinoAmpRef`, `HL00_B_optionTakeSinoPhRef`). Useful for validation and comparison against reference implementations.
- **Utility** — GPU device initialization (`cudaInitDev`) and memory reporting (`memoryInfo`).
- **Legacy / low-level** *(source header only)* — an earlier step-by-step interface that predates the HL00–HL04 functions. Retained for backward compatibility with existing integrations; new projects should use the core workflow instead.

---

## 6. Optional Debug and Profiling Output

The source code includes optional instrumentation that writes intermediate arrays to disk and reports execution times for individual processing steps. This is useful for validating a rebuilt library against a reference implementation, or for profiling performance on new hardware.

The instrumentation is **disabled by default** and excluded from the compiled library. To enable it, uncomment the corresponding directives at the top of `src/ODT_GPU/ODT_GPU.h` and rebuild:

```c
//#define Save_data        // dump intermediate arrays
//#define Save_data2       // dump additional intermediate arrays
//#define Save_timings     // report per-step execution times
//#define Save_timings_GP  // report Gerchberg-Papoulis loop timings
```

The output files are written to the current working directory of the host application.

---

© 2025–2026 ODT GPU Library — Licensed under the GNU General Public License v3.0.
