# ODT_GPU_Library
GPU-accelerated library for real-time limited-angle optical diffraction tomography (LaODT). Enables fast 3D image reconstruction using Direct Inverse and Gerchberg–Papoulis algorithms. Handles very large datasets with adjustable iteration count. Provided as a compiled DLL for easy integration.

---
## Workflows

The library provides two equivalent reconstruction workflows that differ in how input data are provided and processed.

### **1. Raw-data Workflow**
This workflow starts directly from **raw holograms** (e.g., data from the detector). It performs all preprocessing, including optional reference correction, on the GPU.

#### Example sequence:
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

---

### **2. Preprocessed-data Workflow**
This workflow starts from **preprocessed sinograms** (amplitude and phase). It gives full control over each stage of the reconstruction and allows retrieval of intermediate results.

#### Example sequence:
```c
HL01_setParams(...);
HL02_sendDataAndGenerateKO(...);                // sinograms -> K-space
HL03_setParamsAndStartDIandGP(...);
HL04_takeReconstructionAndFreeMemory(...);
```

✅ **Key features:**
- Input: preprocessed sinograms (`float`)
- User-controlled preprocessing (external or custom)
- No optional reference handling (`HL_addReference()`)

---

### **Backward Compatibility**
The function `HL00_addReference()` remains available as an alias for `HL_addReference()` for backward compatibility.

---


## 📁 Repository Structure

```
ODT_GPU_Library/
│
├── include/                     # Header files
│   └── CUDAprocessing.h
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
│       └── workspace.mat
│
├── LICENSE
└── README.md
```

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
