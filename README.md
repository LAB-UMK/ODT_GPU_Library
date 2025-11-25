# ODT_GPU_Library
GPU-accelerated library for real-time limited-angle optical diffraction tomography (LaODT). Enables fast 3D image reconstruction using Direct Inverse and Gerchberg–Papoulis algorithms. Handles very large datasets with adjustable iteration count. Provided as a compiled DLL for easy integration.

---

## Table of Contents
1. [Overview](#overview)
2. [Requirements](#-requirements)
3. [Installation](#installation)
4. [Repository Structure](#repository-structure)
5. [Workflows](#workflows)
6. [Usage Example (MATLAB)](#usage-example-matlab)
7. [Relation to the Collaborating Group’s Library](#relation-to-the-collaborating-groups-library)
8. [License](#license)
9. [Citation](#citation)
10. [Contact](#contact)

---

## Overview

The ODT_GPU_Library was developed as a response to the growing demand for high-speed, high-quality data processing in limited-angle optical diffraction tomography (LaODT).  

This library provides a fully GPU-accelerated implementation of key numerical procedures — including preprocessing, K-space generation, and tomographic reconstruction using Direct Inverse and Gerchberg–Papoulis algorithms.

Here, we propose for the first time, to the best of our knowledge,
a high-speed CUDA implementation of the ODT reconstruction algorithm that enables full and accurate reconstruction through an iterative procedure, without compromising image quality, limiting the measurement volume, the number of angular projections, or requiring a real-to-complex Hermitian Fourier transform.

The library is distributed as a precompiled dynamic-link library (DLL) that can be easily integrated with MATLAB, LabVIEW, Python, or custom C/C++ software, allowing researchers to seamlessly incorporate high-performance tomographic reconstruction into their experimental workflows.

---

## 🧠 Requirements
- **Operating system:** Windows 10 or 11 (64-bit)  
- **GPU:** NVIDIA GPU supporting CUDA Technology
- **CUDA Toolkit:** Supported builds  
  - 10.2  
  - 12.1  
  - 12.6  
  - 13.0  
- **Microsoft Visual C++ Redistributable packages for Visual Studio 2019**
- **MATLAB** (optional, for example scripts)  

---


## Installation

💡 For detailed setup instructions and troubleshooting, see the [Installation Guide](docs/INSTALL.md).

To use the ODT GPU Library, make sure the following components are installed:

- **MATLAB** (R2020a or newer) – required for example scripts and integration.
- **NVIDIA CUDA Toolkit** (matching the version of your selected DLL).
- **Microsoft Visual C++ Redistributable Package** – required for loading the compiled DLL.

### Steps:
1. Install MATLAB
2. Install MATLAB Add-On: MinGW-w64 C/C++ Compiler
3. Install the NVIDIA CUDA Toolkit (version 10.2, 12.1, 12.6, or 13.0)
4. Install Microsoft Visual C++ Redistributable (x64)
5. Download ODT_GPU_Library (using Git LFS or manually download each binary file)
6. Run the example script provided in the `examples/` directory to verify the setup (ensure the library version used in the script is compatible with the installed CUDA Toolkit version)

💡 For detailed setup instructions and troubleshooting, see the [Installation Guide](docs/INSTALL.md).

---

## 📁 Repository Structure

```
ODT_GPU_Library/
│
├── include/                     # Header files
│   └── ODT_GPU.h
│
├── bin/                         # Precompiled DLLs
│   └── v-1_0/
│       ├── cuda-10_2/ODT_GPU.dll
│       ├── cuda-12_1/ODT_GPU.dll
│       ├── cuda-12_6/ODT_GPU.dll
│       └── cuda-13_0/ODT_GPU.dll
│
├── examples/
│   └── matlab/
│       ├── example_reconstruction.m
│       ├── checkError.m
│       ├── convertFloatMatrixForC.m
│       ├── convertIntMatrixForC.m
│       ├── convertMatrixFromCToMatlab.m
│       └── dataAndParams.mat
│
├── docs/
│   ├── INSTALL.md
│   └── USAGE.md
├── LICENSE
└── README.md
```

---

## Workflows

The library supports two main workflows:

1. **Raw-data Workflow** — starting from detector holograms, performs full preprocessing and reconstruction directly on the GPU.
2. **Preprocessed-data Workflow** — starts from preprocessed sinograms, offering fine control and access to intermediate data.

### Example sequence:

**Raw-data Workflow**
```c
HL_addReference(...);                           // optional reference
HL00to02_FromPreprocToGenKO(...);               // hologram -> preprocessing -> K-space
HL03_setParamsAndStartDIandGP(...);       		// reconstruction (Direct Inverse or iterative)
HL04_takeReconstructionAndFreeMemory(...);      // get final reconstruction
```

✅ **Key features:**
- Input: raw holograms (`short int`)
- Automatic preprocessing (FFT, filtering, windowing, normalization)
- Optional reference handling (`HL_addReference()`)
- Fastest and simplest path for real-time or automated processing

For detailed step-by-step examples of both workflows, see the [Usage Guide](docs/USAGE.md).

---

## ⚙️ Usage Example (MATLAB)

```matlab
% Load the DLL and header
loadlibrary(dllPath, headerPath);

% Preprocessing: from raw hologram to K-space generation
err = calllib('ODT_GPU','HL00to02_FromPreprocToGenKO', hologramInt16, Nx, Ny, nproj, ...
              NA, lambda, cam_pix, M, n_imm, do_NNC, ...
              Kxy_dim_p, Kspace_oversampling_z, cos_factor, fftWindowScale, ...
              strcmp(Approx,'Born'));

% Run reconstruction
err = calllib('ODT_GPU','HL03_setParamsAndStartDIandGP', nGPi, do_TC, do_NNC, ...
              relaxM_cuda, betaM, -1, objshift, convertFloatMatrixForC(object_support));

% Retrieve reconstruction
rec = libpointer('singlePtr', zeros(Nx, Ny, Nz, 'single'));
err = calllib('ODT_GPU','HL04_takeReconstructionAndFreeMemory', pNrec);

% Unload the library
unloadlibrary('ODT_GPU');
```

For a complete MATLAB example, see  
[`examples/matlab/example_reconstruction.m`](examples/matlab/example_reconstruction.m).

---

## Relation to the Collaborating Group’s Library

A MATLAB-based implementation of the LaODT reconstruction algorithm, developed and published by a collaborating research group, is available at [EWALD](https://github.com/biopto/EWALD).  
The GPU-accelerated library presented here was developed in close collaboration with that group but implemented independently and entirely from scratch, with the goal of achieving real-time performance and enhanced computational efficiency on modern graphics hardware.  

Although the numerical methods and data processing strategies differ fundamentally from the MATLAB implementation, the API structure and naming conventions were intentionally designed to remain compatible with it.  
This approach facilitates direct comparison and validation of reconstruction results between the two environments while maintaining complete algorithmic and implementation independence.

## 🧠 Requirements
- **Operating system:** Windows 10 or 11 (64-bit)  
- **GPU:** NVIDIA GPU with CUDA Compute Capability ≥ 7.5  
- **CUDA Toolkit:** Supported builds  
  - 10.2  
  - 12.1  
  - 12.6  
  - 13.0  
- **MATLAB** (optional, for example scripts)  
- **Microsoft Visual C++ Redistributable packages for Visual Studio 2019** (for C/C++ integration)

---

## 📜 License
This work is distributed under the **CC BY-NC-ND 4.0 License**.  
Non-commercial use only. Modifications are not permitted. 
For full license text, see the [LICENSE](LICENSE) file or visit  
👉 [https://creativecommons.org/licenses/by-nc-nd/4.0/](https://creativecommons.org/licenses/by-nc-nd/4.0/) 
If you use this library in your research, please cite the related paper below.

---

## 📚 Citation
If you use this library in academic work, please cite:

> Marcin Sylwestrzak, Paweł Ossowski, Wojciech Krauze, Maria Baczewska, Szymon Tamborski, Arkadiusz Kuś, Małgorzata Kujawińska, Maciej Szkulmowski, *Wide-field, real-time limited-angle optical diffraction tomography using massively parallel data processing*, 2025.  
> DOI: [link to paper]

---

## 📬 Contact
For questions or collaboration inquiries, please contact:  
📧 [mars@umk.pl]
