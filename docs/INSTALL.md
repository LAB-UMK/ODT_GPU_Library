# Installation Guide

This document describes the steps required to correctly install and configure the ODT GPU Library for use with MATLAB and CUDA-based systems.

---

## 1. System Requirements

- **Operating System:** Windows 10 / 11 (64-bit)
- **GPU:** NVIDIA GPU supporting CUDA Technology
- **Software:**
  - **MATLAB R2020a** or newer with MinGW-w64 C/C++ Compiler add-on
  - **NVIDIA CUDA Toolkit** (version matching the chosen DLL build, e.g. 10.2, 12.1, 12.6, 13.0)
  - **Microsoft Visual C++ Redistributable Package (x64)**

---

## 2. Preparing the Environment

### 1. Install MATLAB

### 2. Install MATLAB Add-On: MinGW-w64 C/C++ Compiler
Required to load external DLLs via loadlibrary().

If missing, install it via Add-On Manager:
Home → Add-Ons → Get Add-Ons → “MinGW-w64 C/C++ Compiler”

### 3. Install the NVIDIA CUDA Toolkit
Install the CUDA version matching the DLL folder (10.2, 12.1, 12.6, or 13.0).

### 4. Install Microsoft Visual C++ Redistributable (x64)
Download:
https://aka.ms/vc14/vc_redist.x64.exe

---

## 3. Downloading the Library (with CUDA Compatibility)

⚠️ Important Note on Large Files (Git LFS)
This repository contains large binary files (e.g., .dll, datasets) managed by Git LFS (Large File Storage). Note that using "Download ZIP" will only provide the actual text files and small pointer files for the binaries, which must be replaced. To correctly obtain the full binary files, you must either clone the repository using an installed Git LFS client or manually download each large file from its corresponding page on GitHub.

Navigate to the repository folder:

bin/  
  v-1_0/  
    cuda-10_2/  
    cuda-12_1/  
    cuda-12_6/  
    cuda-13_0/  

Download two files:
- CUDAprocessing.dll  
- CUDAprocessing.h  

---

## 4. MATLAB Integration

Copy the DLL and header file to your MATLAB working directory.

Load the library:

dll = 'path\\CUDAprocessing.dll';  
hdr = 'path\\CUDAprocessing.h';  
loadlibrary(dll, hdr);

Verify available functions:

libfunctions('CUDAprocessing', '-full')

---

## 5. Verifying the Installation

1. Load the library  
2. Run a minimal example workflow  
3. Confirm that the GPU initializes correctly  
4. A reconstructed image slice should appear  

---

© 2025 ODT GPU Library — For non-commercial research use only.