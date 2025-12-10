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

© 2025 ODT GPU Library — For non-commercial research use only.
