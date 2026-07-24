# Installation Guide  

This document describes the steps required to correctly install and configure the ODT GPU Library for use with MATLAB and CUDA-based systems. Instructions for rebuilding the library from source are provided at the end of this document.

## System Requirements

- **Operating System:** Windows 10 / 11 (64-bit)
- **GPU:** NVIDIA GPU supporting CUDA Technology (see CUDA version compatibility table below)
- **Software:**
  - **MATLAB R2020a** or newer with MinGW-w64 C/C++ Compiler add-on
  - **NVIDIA CUDA Toolkit** (version matching the chosen DLL build, e.g. 10.2, 12.1, 12.6, 13.0)
  - **Microsoft Visual C++ Redistributable Package (x64)**
  - **Visual Studio 2022** with the *Desktop development with C++* workload — only required when building from source

## GPU Compatibility per CUDA Toolkit Version

Each DLL in the `bin/v-1_0/` directory is compiled with a specific CUDA Toolkit.
The available GPU architectures depend on the corresponding CUDA version. Driver version must also be compatible with the selected CUDA Toolkit.

| DLL version folder | CUDA Toolkit | Supported GPU architectures | Example GPUs |
|--------------------|--------------|------------------------------|--------------|
| `cuda-10_2`        | CUDA 10.2    | Kepler (SM 3.x), Maxwell (SM 5.x), Pascal (SM 6.x), Turing (SM 7.5), Ampere (SM 8.6\*) | GTX 750/960/1060, GTX 1650/1660, RTX 2060/2080, RTX 3060/3080 |
| `cuda-12_1`        | CUDA 12.1    | Turing (SM 7.5), Ampere (SM 8.x), Ada Lovelace (SM 8.9) | RTX 2060/2080, RTX 3050–3090, RTX 4060–4090 |
| `cuda-12_6`        | CUDA 12.6    | Same as Toolkit 12.1 (Turing+, Ampere, Ada) | RTX 20xx, 30xx, 40xx |
| `cuda-13_0`        | CUDA 13.0    | Ampere (SM 8.x), Ada (SM 8.9), Blackwell (SM 9.0) | RTX 30xx, 40xx, 50xx |

\* CUDA 10.2 unofficially works on some Ampere GPUs through updated drivers, but not guaranteed.

Older GPUs (GTX 900/1000/1600 series) **are not supported** by CUDA 12.x and 13.x.
If you have one of these cards, you must use the `cuda-10_2` DLL.

Newer GPUs (RTX 20xx/30xx/40xx/50xx) can use any of the 12.x or 13.x DLL versions.

The library has been tested on the following NVIDIA GPUs: RTX 2050, RTX 3070, and RTX Pro 500 (laptops), as well as GTX 1080, RTX 3090, and RTX 4090 (workstations). Other CUDA-compatible GPUs should also work but have not been explicitly tested.

---

## Running the MATLAB Example

Although the example installation below is presented for MATLAB, the ODT GPU Library is a standard Windows DLL and can be used from any environment capable of calling external C functions.
Typical options include:
- C / C++ applications
- Python (via ctypes or cffi)
- LabVIEW
- C# / .NET
- MATLAB (presented here)

The following instructions describe the setup required to run the provided MATLAB examples.

### 1. Install MATLAB

### 2. Install MATLAB Add-On: MinGW-w64 C/C++ Compiler
This compiler is required by MATLAB to work with external DLL libraries via loadlibrary().  

### 3. Install the NVIDIA CUDA Toolkit (version 10.2, 12.1, 12.6, or 13.0 depending on your GPU and the compatibility table above)
Download from: https://developer.nvidia.com/cuda-toolkit-archive

### 4. Install Microsoft Visual C++ Redistributable (x64)  
Download from: https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist

### 5. Download ODT_GPU_Library (using Git LFS or manually download each binary file)
⚠️ Important Note on Large Files (Git LFS)
This repository contains large binary files (e.g., .dll, datasets) managed by Git LFS (Large File Storage). Note that using "Download ZIP" will only provide the actual text files and small pointer files for the binaries, which must be replaced. To correctly obtain the full binary files, you must either clone the repository using an installed Git LFS client or manually download each large file from its corresponding page on GitHub.

Download the `ODT_GPU.dll` file from the appropriate folder and the `ODT_GPU.h` file from `include/` directory.

### 6. Run the example script provided in the `examples/` directory to verify the setup (ensure the library version used in the script is compatible with the installed CUDA Toolkit version)

Before running the example, first add the entire `examples/matlab/` directory (containing the sample script, helper functions, and data) to the MATLAB search path (by right-click on this folder and choosing "Add to Path/Selected Folders"). This ensures that MATLAB can access all required files.

Next, open `matlab/` directory and `example_reconstruction.m` file, place a breakpoint on the `loadlibrary(dllPath, headerPath);` line, and start the script. When execution stops at the breakpoint, verify that the variables
`dllPath` and `headerPath` point to the correct locations of `ODT_GPU.dll` and `ODT_GPU.h`.  
If the files were moved or renamed, adjust these variables accordingly.

Once the paths are confirmed, continue execution of the script.  
If reconstruction appears, installation is correct.

If MATLAB fails to load the DLL with the message “The specified module could not be found”, verify that the CUDA runtime libraries (e.g. `cudart64*.dll`, `cufft64_*.dll`) are visible in the system PATH.

---

## Building the Library from Source

Rebuilding is **optional** — precompiled DLLs are provided in `bin/`. Building from source is useful for inspecting the implementation, modifying the algorithms, or targeting a CUDA Toolkit version other than the ones provided.

### Prerequisites

- **Visual Studio 2022** with the *Desktop development with C++* workload
- **NVIDIA CUDA Toolkit 12.6** (or another version — see below)

Install Visual Studio **before** the CUDA Toolkit. The CUDA installer integrates with Visual Studio and adds the build customizations required to compile `.cu` files; if CUDA is installed first, the integration may be missing and the project will fail to load its CUDA build rules.

### Build Steps

1. Open `src/ODT_GPU/ODT_GPU.sln` in Visual Studio 2022.
2. Select the **Release | x64** configuration from the toolbar.
3. Build the solution (**Build → Build Solution**, or `Ctrl+Shift+B`).

The compiled `ODT_GPU.dll` is written to the configuration output directory (e.g. `src/ODT_GPU/x64/Release/`).

### Building against a different CUDA Toolkit version

The project references **CUDA 12.6** by default. To build against a different installed version:

1. Right-click the **ODT_GPU** project in Solution Explorer.
2. Select **Build Dependencies → Build Customizations…**
3. Uncheck the current CUDA entry, check the version installed on your system, and confirm.

Include and library directories are resolved through the `$(CUDA_PATH)` environment variable set by the CUDA installer, so no manual path editing is required.

### GPU Architectures

The project is configured to generate code for compute capabilities 6.1 through 9.0, plus PTX for forward compatibility with future architectures. To target a narrower set of GPUs and shorten build times, edit **Project Properties → CUDA C/C++ → Device → Code Generation**.

### Optional Debug and Profiling Output

The source includes optional instrumentation for dumping intermediate data to disk and reporting per-step execution times. It is disabled by default. To enable it, uncomment the relevant `#define` directives at the top of `src/ODT_GPU/ODT_GPU.h` and rebuild:

```c
//#define Save_data
//#define Save_data2
//#define Save_timings
//#define Save_timings_GP
```

The output files are written to the current working directory of the host application.

### Troubleshooting

- **“Cannot find CUDA 12.6.props”** — the CUDA Toolkit version referenced by the project is not installed. Either install CUDA 12.6 or switch the build customization as described above.
- **Unresolved external symbols referring to `cufft` or `cudart`** — verify that `CUDA_PATH` is set in the system environment and points to the installed toolkit.
- **Project fails to load `.cu` files** — reinstall the CUDA Toolkit after Visual Studio so that the Visual Studio integration is registered.

---

© 2025–2026 ODT GPU Library — Licensed under the GNU General Public License v3.0.
