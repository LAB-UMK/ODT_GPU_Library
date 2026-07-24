# ODT_GPU_Library
GPU-accelerated library for real-time limited-angle optical diffraction tomography (LaODT). Enables fast 3D image reconstruction using Direct Inverse and Gerchberg–Papoulis algorithms. Handles very large datasets with adjustable iteration count. Provided as a compiled DLL for easy integration, together with the complete CUDA source code.

---

## Table of Contents
1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Repository Structure](#repository-structure)
5. [Building from Source](#building-from-source)
6. [Workflows](#workflows)
7. [Usage Example (MATLAB)](#usage-example-matlab)
8. [Relation to the Collaborating Group’s Library](#relation-to-the-collaborating-groups-library)
9. [License](#license)
10. [Citation](#citation)
11. [Contact](#contact)

---
<a id="overview"></a>
## Overview

The ODT_GPU_Library was developed as a response to the growing demand for high-speed, high-quality data processing in limited-angle optical diffraction tomography (LaODT).  

This library provides a fully GPU-accelerated implementation of key numerical procedures — including preprocessing, K-space generation, and tomographic reconstruction using Direct Inverse and Gerchberg–Papoulis algorithms.

Here, we propose for the first time, to the best of our knowledge,
a high-speed CUDA implementation of the ODT reconstruction algorithm that enables full and accurate reconstruction through an iterative procedure, without compromising image quality, limiting the measurement volume, the number of angular projections, or requiring a real-to-complex Hermitian Fourier transform.

The library is distributed as a precompiled dynamic-link library (DLL) that can be easily integrated with MATLAB, LabVIEW, Python, or custom C/C++ software, allowing researchers to seamlessly incorporate high-performance tomographic reconstruction into their experimental workflows. The complete CUDA source code is provided in the `src/` directory, so the library can also be inspected, modified, and rebuilt.

---
<a id="requirements"></a>
## 🧠 Requirements
- **Operating system:** Windows 10 or 11 (64-bit)  
- **GPU:** NVIDIA GPU supporting CUDA Technology
- **CUDA Toolkit:** Supported builds  
  - 10.2  
  - 12.1  
  - 12.6  
  - 13.0  
- **Microsoft Visual C++ Redistributable (x64)**
- **MATLAB** (optional, for example scripts)  
- **Visual Studio 2022** with the *Desktop development with C++* workload (only if building from source)

---

<a id="installation"></a>
## Installation

💡 For detailed setup instructions see the [Installation Guide](docs/INSTALL.md).

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

💡 For detailed setup instructions see the [Installation Guide](docs/INSTALL.md).

---
<a id="repository-structure"></a>
## 📁 Repository Structure

```
ODT_GPU_Library/
│
├── include/                          # Public headers for the compiled DLL
│   ├── ODT_GPU.h                     # Simplified header (core workflow API)
│   └── helper.h                      # EXPORTED_FUNCTION macro definitions
│
├── src/                              # Complete CUDA source code
│   └── ODT_GPU/                      # Visual Studio 2022 project
│       ├── ODT_GPU.sln               # Solution file
│       ├── ODT_GPU.vcxproj
│       ├── ODT_GPU.vcxproj.filters
│       ├── ODT_GPU.h                 # Full header (core + diagnostic + legacy API)
│       ├── Reconstruction.cu         # DI and Gerchberg-Papoulis reconstruction
│       ├── PreprocessingAndKspace.cu # Preprocessing and K-space generation
│       ├── helper.h
│       ├── pch.h
│       ├── pch.cpp
│       └── dllmain.cpp
│
├── bin/                              # Precompiled DLLs (with import libraries)
│   └── v-1_0/
│       ├── cuda-10_2/                # ODT_GPU.dll + ODT_GPU.lib
│       ├── cuda-12_1/                # ODT_GPU.dll + ODT_GPU.lib
│       ├── cuda-12_6/                # ODT_GPU.dll + ODT_GPU.lib
│       └── cuda-13_0/                # ODT_GPU.dll + ODT_GPU.lib
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
<a id="building-from-source"></a>
## 🔧 Building from Source

The complete CUDA source code is provided in the `src/` directory. Rebuilding the library is optional — precompiled DLLs are available in `bin/` — but may be useful for inspecting the implementation, adapting it to specific needs, or targeting a different CUDA Toolkit version.

### Prerequisites
- **Visual Studio 2022** with the *Desktop development with C++* workload
- **NVIDIA CUDA Toolkit 12.6** (see below for building against a different version)

### Steps
1. Open `src/ODT_GPU/ODT_GPU.sln` in Visual Studio 2022.
2. Select the **Release | x64** configuration.
3. Build the solution (**Build → Build Solution**, or `Ctrl+Shift+B`).

The resulting `ODT_GPU.dll` is placed in the configuration output directory.

### Building against a different CUDA Toolkit version
The project references CUDA 12.6 by default. To build with another installed version, right-click the project in Solution Explorer, choose **Build Dependencies → Build Customizations…**, and select the CUDA version available on your system. Include and library paths are resolved through the `$(CUDA_PATH)` environment variable and do not require manual editing.

### Notes
- The build targets **x64 only**; 32-bit configurations are not supported.
- The header in `src/ODT_GPU/ODT_GPU.h` documents the complete API, including optional diagnostic functions and the legacy low-level interface. The header in `include/ODT_GPU.h` is a simplified version exposing the core workflow functions recommended for typical use.
- Optional debug and profiling output (intermediate data dumps and per-step timings) can be enabled by uncommenting the corresponding `#define` directives at the top of `src/ODT_GPU/ODT_GPU.h`.

---
<a id="workflows"></a>
## Workflows

The library supports two main workflows:

1. **Preprocessed-data Workflow** — starting from preprocessed sinograms, performs reconstruction on the GPU.
2. **Raw-data Workflow** — starting from detector holograms, performs full preprocessing and reconstruction directly on the GPU.

For detailed step-by-step examples of both workflows, see the [Usage Guide](docs/USAGE.md).

---
<a id="usage-example-matlab"></a>
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
err = calllib('ODT_GPU','HL04_takeReconstructionAndFreeMemory', rec);

% Unload the library
unloadlibrary('ODT_GPU');
```

For a complete MATLAB example, see  
[`examples/matlab/example_reconstruction.m`](examples/matlab/example_reconstruction.m).

---
<a id="relation-to-the-collaborating-groups-library"></a>
## Relation to the Collaborating Group’s Library

A MATLAB-based implementation of the LaODT reconstruction algorithm, developed and published by a collaborating research group, is available on GitHub via the [EWALD](https://github.com/biopto/EWALD) repository.  
The GPU-accelerated library presented here was developed in close collaboration with that group but implemented independently and entirely from scratch, with the goal of achieving real-time performance and enhanced computational efficiency on modern graphics hardware.  

Although the numerical methods and data processing strategies differ fundamentally from the MATLAB implementation, the API structure and naming conventions were intentionally designed to remain compatible with it.  
This approach facilitates direct comparison and validation of reconstruction results between the two environments while maintaining complete algorithmic and implementation independence.

---
<a id="license"></a>
## 📜 License
This work is distributed under the **GNU General Public License v3.0 (GPL-3.0)**.  
You are free to use, modify, and redistribute this software, including for commercial purposes, provided that derivative works are distributed under the same license and that the source code remains available.  
For full license text, see the [LICENSE](LICENSE) file or visit  
👉 [https://www.gnu.org/licenses/gpl-3.0.html](https://www.gnu.org/licenses/gpl-3.0.html)  
If you use this library in your research, please cite the related paper below.

---
<a id="citation"></a>
## 📚 Citation
If you use this library in academic work, please cite:

> Marcin Sylwestrzak, Wojciech Krauze, Paweł Ossowski, Maria Baczewska, Szymon Tamborski, Arkadiusz Kuś, Małgorzata Kujawińska, Maciej Szkulmowski, *Wide-field, real-time limited-angle optical diffraction tomography using massively parallel data processing*, Biomedical Optics Express, 2026.  
> DOI: [link to paper]

---
<a id="contact"></a>
## 📬 Contact
For questions or collaboration inquiries, please contact:  
📧 [mars@umk.pl]


---

© 2025–2026 ODT GPU Library — Licensed under the GNU General Public License v3.0.
