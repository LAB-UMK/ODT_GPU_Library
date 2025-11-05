# Installation Guide

This document describes the steps required to correctly install and configure the ODT GPU Library for use with MATLAB and CUDA-based systems.

---

## 1. System Requirements

- **Operating System:** Windows 10 / 11 (64-bit)
- **GPU:** NVIDIA GPU supporting CUDA Technology
- **Software:**
  - **MATLAB R2020a** or newer
  - **NVIDIA CUDA Toolkit** (version matching the chosen DLL build, e.g. 10.2, 12.1, 12.6)
  - **Microsoft Visual C++ Redistributable Package (x64)**

---

## 2. Downloading the Library

1. Navigate to the `bin/` directory in this repository.  
2. Choose the appropriate version:
   ```
   bin/
     v-1_0/
       cuda-10_2/
       cuda-12_1/
       cuda-12_6/
   ```
3. Download:
   - The DLL file (`CUDAprocessing.dll`)
   - The header file (`CUDAprocessing.h`)

---

## 3. MATLAB Integration

1. Copy both files (`CUDAprocessing.dll` and `CUDAprocessing.h`) to your MATLAB working directory.
2. Load the library:
   ```matlab
   fullpathToDll = 'path\to\CUDAprocessing.dll';
   fullpathToHeader = 'path\to\CUDAprocessing.h';
   loadlibrary(fullpathToDll, fullpathToHeader);
   ```
3. Verify the installation by calling:
   ```matlab
   libfunctions('CUDAprocessing', '-full')
   ```
   This should display a list of all exported functions.

---

## 4. CUDA Toolkit Compatibility

Ensure that the version of the CUDA Toolkit installed on your system **matches the DLL version** you downloaded.  
If you are using CUDA 12.1, download the DLL from:
```
bin/v-1_0/cuda-12_1/
```

> ⚠️ If versions do not match, MATLAB may report an error such as:  
> *“The specified module could not be found.”*  
> This typically indicates a mismatch between DLL and runtime toolkit.

---

## 5. Microsoft Visual C++ Redistributable

The library depends on the standard Microsoft Visual C++ runtime.  
Install the latest version from the official Microsoft website:  
👉 https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist

---

## 6. Verifying the Installation

1. Open MATLAB.
2. Load the library as shown above.
3. Run one of the provided example scripts:
   ```matlab
   example_real_time.m
   ```
4. If the GPU is properly initialized, you should see console output similar to:
   ```
   CUDA device initialized successfully.
   ```
5. A test reconstruction will appear in a MATLAB figure window.

---

© 2025 ODT GPU Library — For non-commercial research use only.
