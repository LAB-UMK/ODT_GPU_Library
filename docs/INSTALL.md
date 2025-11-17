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
Ensure MATLAB is installed and functioning correctly.

### 2. Install MATLAB MinGW-w64 C/C++ Compiler
Required to load external DLLs via loadlibrary().

Check installation in MATLAB:
mex -setup cpp

If missing, install it via Add-On Manager:
Home → Add-Ons → Get Add-Ons → “MinGW-w64 C/C++ Compiler”

### 3. Install the NVIDIA CUDA Toolkit
Install the CUDA version matching the DLL folder (10.2, 12.1, 12.6, or 13.0).

Verify the installation:
nvcc --version

### 4. Install Microsoft Visual C++ Redistributable (x64)
Download:
https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist

---

## 3. Downloading the Library (with CUDA Compatibility)

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

### CUDA Compatibility Table

Installed CUDA Toolkit | DLL Folder  
----------------------|----------------  
CUDA 10.2             | cuda-10_2/  
CUDA 12.1             | cuda-12_1/  
CUDA 12.6             | cuda-12_6/  
CUDA 13.0             | cuda-13_0/  

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