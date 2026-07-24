/**
 * @file ODT_GPU.h
 * @brief Header file for the GPU-accelerated library for real-time limited-angle optical diffraction tomography (LaODT).
 *
 * This library provides a complete set of high-level functions (HLxx) and supporting GPU procedures
 * for tomographic reconstruction using Direct Inverse and Gerchberg–Papoulis algorithms.
 * It enables fast 3D imaging with adjustable reconstruction parameters and real-time preview capabilities.
 *
 * This library is distributed as an open-source project under the GPL-3.0 license.
 *
 * @version 1.1
 * @date 2025
 *   - Marcin Sylwestrzak (Nicolaus Copernicus University, Toruń, Poland)
 *   - Wojciech Krauze (Warsaw University of Technology, Warsaw, Poland)
 *   - Paweł Ossowski (Nicolaus Copernicus University, Toruń, Poland)
 *   - Maria Baczewska (Warsaw University of Technology, Warsaw, Poland)
 *   - Szymon Tamborski (Nicolaus Copernicus University, Toruń, Poland)
 *   - Arkadiusz Kuś (Warsaw University of Technology, Warsaw, Poland)
 *   - Małgorzata Kujawińska (Warsaw University of Technology, Warsaw, Poland)
 *   - Maciej Szkulmowski (Nicolaus Copernicus University, Toruń, Poland)
 * @copyright
 * This software is licensed under the GNU General Public License v3.0 (GPL-3.0).
 * See the LICENSE file or visit https://www.gnu.org/licenses/gpl-3.0.html for details.
 * @see The full project repository and documentation are available at:
 *      https://github.com/LAB-UMK/ODT_GPU_Library
 */

#ifndef ODT_GPU_H
#define ODT_GPU_H

#include "helper.h"

 // ================================================================
 //  DEBUG / PROFILING FLAGS
 //  Uncomment to enable saving intermediate data or timing results
 //  to disk. Used for validation and performance analysis.
 //  Corresponding #ifdef blocks are located in the .cu source files.
 // ================================================================
 //#define Save_data
 //#define Save_data2
 //#define Save_timings
 //#define Save_timings_GP

#ifdef  __cplusplus
extern "C" {
#endif

	// =============================================================================================
	//  CORE WORKFLOW FUNCTIONS
	//  Main reconstruction pipeline. Call in the order listed below.
	//  HL00_addReference and HL00_removeReference are backward compatibility aliases.
	// =============================================================================================

	/**
	 * @brief Add a reference hologram to GPU memory for later reuse.
	 *
	 * Transfers a reference hologram (raw data) from host memory to GPU device memory.
	 * The stored reference is reused by subsequent processing steps (e.g., for background/reference
	 * subtraction or normalization) until explicitly removed via HL_removeReference().
	 *
	 * @param hologram             Pointer to the reference hologram stored as int16 (short int).
	 *                             Expected layout: contiguous array with dimensions [_X * _Y * _nproj].
	 *                             Ordering: projection-major or as used by the rest of the pipeline
	 *                             (typically: for each projection, a block of size _X*_Y is stored).
	 * @param _X                   Number of pixels in X (detector width) per projection.
	 * @param _Y                   Number of pixels in Y (detector height) per projection.
	 * @param _nproj               Number of projection angles contained in the hologram dataset.
	 * @param _NA                  Numerical aperture of the imaging objective used during acquisition.
	 * @param lambda               Illumination wavelength (in micrometers).
	 * @param cam_pix              Camera pixel size (micrometers).
	 * @param M                    Transverse magnification of the imaging system.
	 * @param _n_imm               Refractive index of the immersion medium.
	 * @param _do_NNC              Flag controlling the sign:
	 *                             - _do_NNC > 0 : apply **non-negativity** constraint (restrict values >= 0),
	 *                             - _do_NNC < 0 : apply **non-positivity** constraint (restrict values <= 0),
	 *                             - _do_NNC == 0: **no** non-negativity/non-positivity constraint applied.
	 *                             This setting defines how the object refractive index is constrained
	 *                             relative to the immersion medium.
	 * @param _fftWindowsScale     Scaling factor applied to the FFT windowing during preprocessing.
	 *
	 * @return 0 on success; non-zero error code on failure (e.g., invalid parameters, insufficient GPU memory).
	 *
	 * @note The function stores the reference hologram on the GPU until HL_removeReference() is called.
	 *       Ensure that the memory layout of `hologram` matches the expectations of the rest of the pipeline.
	 * @note If calling from MATLAB, convert data to int16 (short) using an appropriate conversion helper
	 *       (e.g., convertIntMatrixForC). Be mindful of MATLAB column-major vs. C row-major ordering if applicable.
	 *
	 * @see HL_removeReference, preproc_sendHologram, HL00to02_FromPreprocToGenKO
	 */
	EXPORTED_FUNCTION int HL_addReference(short int* hologram, int _X, int _Y, int _nproj, float _NA, float lambda,
		float cam_pix, float M, float _n_imm, int _do_NNC, float _fftWindowsScale);

	// Backward compatibility alias:
	// This macro preserves compatibility with older code that used the legacy
	// function name HL00_addReference. The canonical API function name is now
	// HL_addReference.
	EXPORTED_FUNCTION int HL00_addReference(short int* hologram, int _X, int _Y, int _nproj, float _NA, float lambda,
		float cam_pix, float M, float _n_imm, int _do_NNC, float _fftWindowsScale);

	/**
	 * @brief Remove the reference hologram and associated sinogram data from GPU memory.
	 *
	 * This function releases all GPU memory previously allocated for the reference dataset
	 * added via HL_addReference(). It removes both the amplitude and phase components
	 * of the reference sinogram and all intermediate buffers related to reference preprocessing.
	 *
	 * After this call, the reference data is no longer available for subsequent processing steps.
	 * Any attempt to retrieve it (e.g., using HL00_B_optionTakeSinoAmpRef or HL00_B_optionTakeSinoPhRef)
	 * will result in an error.
	 *
	 * @return 0 on success; non-zero error code on failure (e.g., no reference data present in GPU memory).
	 *
	 * @note This function should be called once the reference data is no longer needed,
	 *       for example before loading a new reference hologram or prior to freeing GPU resources.
	 * @note It is safe to call this function even if the reference data has already been removed —
	 *       redundant calls will have no effect beyond returning an appropriate status code.
	 *
	 * @see HL_addReference, HL00_B_optionTakeSinoAmpRef, HL00_B_optionTakeSinoPhRef
	 */
	EXPORTED_FUNCTION int HL_removeReference();

	// Backward compatibility alias for legacy API name:
	EXPORTED_FUNCTION int HL00_removeReference();


	/**
	 * @brief Perform complete preprocessing and K-space generation from raw hologram data
	 * to K-space generation.
	 *
	 * This function executes a combined preprocessing and K-space generation sequence
	 * directly from raw hologram frames acquired by the detector. It consolidates several
	 * processing stages (HL00–HL02) into a single call for maximum efficiency.
	 *
	 * The function performs normalization using the reference hologram previously loaded
	 * with HL_addReference(), applies filtering, performs Fourier transforms, and constructs
	 * the three-dimensional K-space representation of the measured object.
	 *
	 * This combined mode is primarily intended for **real-time imaging** and **fast data pipelines**,
	 * as it minimizes host–GPU data transfers and intermediate allocations.
	 *
	 * @param hologram               Pointer to the raw hologram data (short int*).
	 *                               Expected layout: contiguous array of dimensions `X × Y × _nproj`.
	 *                               Each projection frame must be stored sequentially in memory.
	 * @param X                      Number of pixels in the X dimension (detector width).
	 * @param Y                      Number of pixels in the Y dimension (detector height).
	 * @param _nproj                 Number of projections (illumination angles) in the dataset.
	 * @param NA                     Numerical aperture of the imaging objective used during acquisition.
	 * @param lambda                 Illumination wavelength (in micrometers).
	 * @param cam_pix                Camera pixel size (in micrometers).
	 * @param M                      Optical magnification of the imaging system.
	 * @param _n_imm                 Refractive index of the immersion medium.
	 * @param _do_NNC                Constraint mode flag applied during reconstruction:
	 *                               - `_do_NNC > 0` → apply **non-negativity** constraint (values ≥ 0)
	 *                               - `_do_NNC < 0` → apply **non-positivity** constraint (values ≤ 0)
	 *                               - `_do_NNC == 0` → disable sign constraint
	 *                               This setting defines how the object refractive index is constrained
	 *                               relative to the immersion medium.
	 * @param _K_xy                  Pointer to an integer (int32) that receives the computed lateral size
	 *                               of the generated K-space. The value is updated internally depending
	 *                               on FFT padding and scaling parameters.
	 * @param Kspace_oversampling_z  Oversampling factor applied in the Z dimension during K-space generation.
	 * @param _cosFactor             Parameter defining the **Tukey window** applied to each projection
	 *                               before Fourier transformation.
	 *                               The Tukey window is used to suppress high-frequency noise
	 *                               and edge artifacts originating from sharp intensity discontinuities
	 *                               at the projection boundaries.
	 *                               - `_cosFactor = 0` → rectangular window (no tapering)
	 *                               - `_cosFactor = 1` → full cosine window (strong tapering)
	 *                               Intermediate values (e.g., 0.1–0.3) define the fraction
	 *                               of the window edges that are smoothly tapered.
	 * @param _fftWindowsScale       Scaling factor controlling FFT window amplitude.
	 *                               This parameter allows tuning between noise suppression and resolution.
	 * @param _approxBornNotRytov    Selects the forward scattering model:
	 *                               - `0` → use **Rytov approximation** (phase-based)
	 *                               - non-zero → use **Born approximation** (amplitude-based)
	 *
	 * @return 0 on success; non-zero error code on failure (e.g., invalid parameters, insufficient GPU memory,
	 *         or missing reference data if required).
	 *
	 * @note The function internally performs preprocessing equivalent to HL_addReference(), HL01_setParams
	 *       and HL02_sendDataAndGenerateKO(), and thus does not require separate calls to these functions.
	 * @note The resulting K-space data is stored in GPU memory for subsequent reconstruction steps,
	 *       typically performed using HL03_setParamsAndStartDIandGP() followed by HL04_takeReconstructionAndFreeMemory().
	 * @note The parameter `_K_xy` must point to a valid integer variable; its value is updated on return.
	 *
	 * @see HL_addReference, HL03_setParamsAndStartDIandGP, HL04_takeReconstructionAndFreeMemory
	 */
	EXPORTED_FUNCTION int HL00to02_FromPreprocToGenKO(short int* hologram, int X, int Y, int _nproj, float NA, float lambda,
			float cam_pix, float M, float _n_imm, int _do_NNC, int* _K_xy, int Kspace_oversampling_z, float _cosFactor, float _fftWindowsScale, int _approxBornNotRytov);

	/**
	 * @brief Initialize and allocate GPU resources for tomographic reconstruction.
	 *
	 * This function sets all numerical and physical parameters required for the reconstruction
	 * process and allocates GPU memory buffers accordingly. It defines the geometry of the K-space,
	 * the reconstruction grid, and optical parameters of the imaging system.
	 *
	 * The function must be called before launching the reconstruction algorithm
	 * (HL03_setParamsAndStartDIandGP) when operating in the **manual (multi-step) processing mode**,
	 * where preprocessing and K-space generation are handled separately.
	 *
	 * It initializes all FFT-related structures and prepares memory for in-place transforms
	 * between the object space and the K-space to minimize GPU memory usage.
	 *
	 * @param _K_xy                 Lateral size (in pixels) of the K-space grid (X and Y dimensions).
	 *                              Defines the width and height of the 3D Fourier domain.
	 * @param _K_z                  Axial (Z) dimension of the K-space grid (number of depth samples).
	 * @param _dx                   Lateral pixel size in object space units (µm/pixel).
	 *                              Defines spatial sampling in the reconstruction space.
	 * @param _n_imm                Refractive index of the immersion medium used during imaging.
	 * @param _n_proj               Number of illumination projections (angles) used for reconstruction.
	 * @param _Nx                   Number of pixels in the reconstructed object volume along the X-axis.
	 * @param _Ny                   Number of pixels in the reconstructed object volume along the Y-axis.
	 * @param _dkP                  Step size in K-space along the lateral dimensions (X/Y),
	 *                              corresponding to sampling frequency in Fourier space.
	 * @param _dkPz                 Step size in K-space along the Z dimension.
	 * @param _NA                   Numerical aperture of the imaging objective.
	 * @param _lambda_all           Pointer to an array of illumination wavelengths used in measurement.
	 *                              The array should have `_n_proj` elements if multiple wavelengths
	 *                              are employed, otherwise a single-element array may be provided.
	 * @param _kxp                  Pointer to an array containing the X components of the illumination
	 *                              wave vectors (spatial frequencies) for each projection.
	 *                              Defines the illumination geometry.
	 * @param _kyp                  Pointer to an array containing the Y components of the illumination
	 *                              wave vectors (spatial frequencies) for each projection.
	 * @param _approxBornNotRytov   Selects the forward scattering model used for reconstruction:
	 *                              - `0` → **Rytov approximation** (phase-based)
	 *                              - non-zero → **Born approximation** (amplitude-based)
	 *
	 * @return 0 on success; non-zero error code on failure (e.g., invalid parameters or GPU memory allocation error).
	 *
	 * @note This function initializes both K-space and reconstruction-space buffers.
	 *       The same memory region is reused for in-place FFT and iFFT operations to reduce GPU memory footprint.
	 * @note The function is used only in the full (step-by-step) reconstruction mode.
	 *       In the simplified mode (HL00to02_FromPreprocToGenKO), it is not required.
	 * @note It is good practice to free allocated GPU memory after reconstruction
	 *       using HL04_takeReconstructionAndFreeMemory() or HL02_C_optionFreeGPUmemory().
	 *
	 * @see HL02_sendDataAndGenerateKO, HL03_setParamsAndStartDIandGP, HL04_takeReconstructionAndFreeMemory
	 */
	EXPORTED_FUNCTION int HL01_setParams(int _K_xy, int _K_z, float _dx, float _n_imm, int _n_proj, int _Nx, int _Ny,
		float _dkP, float _dkPz, float _NA, float* _lambda_all, float* _kxp, float* _kyp, int _approxBornNotRytov);


	/**
	 * @brief Send sinogram data to the GPU and generate the K-space representation of the object.
	 *
	 * This function transfers amplitude and phase sinograms from the host to the GPU memory,
	 * applies optional pupil or region-of-interest masks, and performs the generation of the
	 * 3D K-space (Fourier domain) representation of the measured object.
	 *
	 * It should be called after HL01_setParams(), which allocates and initializes all GPU buffers
	 * and defines the reconstruction geometry.
	 *
	 * The function performs all necessary Fourier transforms and normalization steps required
	 * to convert the input sinograms into their spectral (K-space) form, ready for iterative
	 * reconstruction using HL03_setParamsAndStartDIandGP().
	 *
	 * @param _sinoAmp              Pointer to the amplitude sinogram data (float*).
	 *                              The array should contain data of size `_Nx * _Ny * _n_proj`,
	 *                              where each projection is stored sequentially in memory.
	 * @param _sinoPh               Pointer to the phase sinogram data (float*).
	 *                              The layout must match that of `_sinoAmp`.
	 *                              Both amplitude and phase arrays represent the same set
	 *                              of projections acquired at different illumination angles.
	 * @param FpmaskLogical         Pointer to an optional binary (logical) mask array (unsigned char*).
	 *                              This mask defines the effective pupil or region in Fourier space
	 *                              used during K-space generation.
	 *                              - `1` → active region (included in computation)
	 *                              - `0` → masked region (excluded)
	 *                              The mask should have the same lateral size as a single projection (XY).
	 *                              If no mask is needed, pass `NULL`.
	 *
	 * @return 0 on success; non-zero error code on failure (e.g., invalid parameters, missing initialization,
	 *         or insufficient GPU memory).
	 *
	 * @note This function performs the same operation as the internal preprocessing and K-space generation
	 *       performed inside HL00to02_FromPreprocToGenKO(), but allows explicit control over input sinograms.
	 * @note Both amplitude and phase data must be preprocessed prior to calling this function.
	 * @note The generated K-space data remains in GPU memory and is used directly by
	 *       HL03_setParamsAndStartDIandGP() and HL04_takeReconstructionAndFreeMemory().
	 * @note For efficiency, this function reuses buffers allocated in HL01_setParams(),
	 *       performing FFTs in-place.
	 *
	 * @see HL01_setParams, HL03_setParamsAndStartDIandGP, HL04_takeReconstructionAndFreeMemory,
	 *      HL02_B_optionTakeKO, HL02_B_optionTakeEW, HL02_C_optionFreeGPUmemory
	 */
	EXPORTED_FUNCTION int HL02_sendDataAndGenerateKO(float* _sinoAmp, float* _sinoPh, unsigned char* FpmaskLogical);


	/**
	 * @brief Configure reconstruction parameters and start the tomographic reconstruction
	 *        using either the Direct Inverse (DI) or Gerchberg–Papoulis (GP) algorithm.
	 *
	 * This function performs the core tomographic reconstruction from data previously stored
	 * in GPU memory. It uses either the **Direct Inverse** (DI) method (single-step reconstruction)
	 * or the **Gerchberg–Papoulis** (GP) iterative algorithm, depending on the number
	 * of iterations specified in `_nGPi`.
	 *
	 * All data required for reconstruction (K-space, configuration, masks) must already be present
	 * in GPU memory as a result of prior calls to HL00to02_FromPreprocToGenKO() or
	 * HL01_setParams() + HL02_sendDataAndGenerateKO().
	 *
	 * @param _nGPi                 Number of iterations of the Gerchberg–Papoulis algorithm.
	 *                              - `_nGPi = 0` → Direct Inverse (DI) mode (no iteration).
	 *                              - `_nGPi > 0` → perform `_nGPi` iterations of the GP algorithm.
	 *                              Each iteration enforces amplitude and support constraints
	 *                              in both object and K-space domains.
	 * @param _do_TC                Flag enabling or disabling the **transparency constraint**:
	 *                              - `1` → remove the imaginary part of the refractive index
	 *                                (assuming non-absorbing samples).
	 *                              - `0` → skip transparency enforcement.
	 * @param _do_NNC               Flag controlling the **non-negativity / non-positivity constraint**:
	 *                              - `_do_NNC > 0` → enforce **non-negativity** constraint
	 *                                (object refractive index ≥ background).
	 *                              - `_do_NNC < 0` → enforce **non-positivity** constraint
	 *                                (object refractive index ≤ background).
	 *                              - `_do_NNC == 0` → no constraint applied.
	 *                              This constraint is typically used to improve convergence.
	 * @param _relaxM               Relaxation factor (0 < `_relaxM` ≤ 1) controlling convergence speed
	 *                              in the GP algorithm. Lower values yield smoother convergence
	 *                              at the cost of more iterations.
	 * @param _betaM                Parameter weighting the amplitude constraint; defines how strongly
	 *                              the amplitude of the object field is forced to match measured data.
	 * @param _kn_mean              Average wavenumber (`2π * n_imm / λ_mean`) used for scaling and
	 *                              normalization of the refractive index contrast.
	 * @param _objShift             Integer shift value defining the axial position (Z-offset) of the
	 *                              reconstructed object volume relative to the imaging plane.
	 * @param _objSupport           Pointer to a float array defining the **spatial support mask**.
	 *                              Each voxel corresponds to a location in the reconstruction space,
	 *                              where `1` indicates a voxel belonging to the object and `0`
	 *                              indicates background.
	 *                              This mask accelerates convergence by constraining the solution
	 *                              to known physical boundaries.
	 *                              If not used, pass `NULL`.
	 *
	 * @return 0 on success; non-zero error code on failure (e.g., invalid parameters, missing K-space data,
	 *         or insufficient GPU memory).
	 *
	 * @note The function performs in-place FFT and iFFT transformations between object and K-space.
	 *       During iterative reconstruction, multiple constraints are applied in both domains.
	 * @note After the final iteration, the reconstructed object remains in GPU memory and can be
	 *       retrieved using HL04_takeReconstructionAndFreeMemory() or HL04a_takeRecSingleSliceAndFreeMemory().
	 * @note The DI mode (no iteration) is highly optimized for **real-time imaging**, achieving up to
	 *       several frames per second depending on the dataset size and GPU hardware.
	 *
	 * @see HL00to02_FromPreprocToGenKO, HL01_setParams, HL02_sendDataAndGenerateKO,
	 *      HL04_takeReconstructionAndFreeMemory, HL04a_takeRecSingleSliceAndFreeMemory
	 */
	EXPORTED_FUNCTION int HL03_setParamsAndStartDIandGP(int _nGPi, int _do_TC, int _do_NNC, float _relaxM,
		float _betaM, float _kn_mean, int _objShift, float* _objSupport);

	/**
	 * @brief Retrieve the full complex-valued 3D reconstruction from GPU memory and release GPU resources.
	 *
	 * This function transfers the final reconstructed 3D object from GPU device memory to host memory.
	 * The data represent the **complex refractive index contrast** (or scattered field) obtained after
	 * performing the Direct Inverse (DI) or Gerchberg–Papoulis (GP) reconstruction.
	 *
	 * The reconstruction is stored in GPU memory as interleaved complex values, where real and imaginary
	 * components alternate in memory. This function performs a direct copy of that buffer to host memory
	 * and then releases all GPU memory used during reconstruction.
	 *
	 * @param complexData_host      Pointer to a preallocated host array (`float*`) that will receive
	 *                              the reconstructed **complex 3D data**.
	 *                              The array must have a length of `2 * Nx * Ny * Nz` elements, where
	 *                              each voxel is represented by two consecutive floats:
	 *                              `[Re(0), Im(0), Re(1), Im(1), …]`.
	 *
	 * @return 0 on success; non-zero error code on failure (e.g., invalid pointer, missing GPU data,
	 *         or CUDA memory transfer error).
	 *
	 * @note The returned dataset contains both the real and imaginary parts of the reconstruction.
	 *       It does **not** contain the magnitude (absolute value); to retrieve a single magnitude slice,
	 *       use HL04a_takeRecSingleSliceAndFreeMemory().
	 * @note After execution, all GPU buffers related to the reconstruction, K-space, and intermediate
	 *       computations are released.
	 * @note The function must be called only once per reconstruction sequence. Repeated calls without
	 *       reprocessing new data will return an error.
	 *
	 * @par MATLAB usage:
	 * When calling this function from MATLAB, preallocate a `single` array of size `[2 * Nx * Ny * Nz, 1]`
	 * and pass its pointer to receive the interleaved complex data.
	 * Afterward, use the provided helper function:
	 * @code
	 * n_rec = convertMatrixFromCToMatlab(pNrec.Value, Nx, Ny, Nz);
	 * @endcode
	 * This function automatically reorders the data from C-style indexing (z–y–x) to MATLAB-style (y–x–z)
	 * and reconstructs the complex matrix `n_rec = Re + i·Im`.
	 *
	 * @see HL03_setParamsAndStartDIandGP, HL04a_takeRecSingleSliceAndFreeMemory,
	 *      HL02_C_optionFreeGPUmemory
	 */
	EXPORTED_FUNCTION int HL04_takeReconstructionAndFreeMemory(float* complexData_host);

	/**
	 * @brief Retrieve a single axial slice (magnitude) from the reconstructed 3D volume and release GPU memory.
	 *
	 * Transfers one selected Z-slice of the final reconstruction from GPU to host as real-valued magnitude.
	 * After the transfer completes, all GPU buffers allocated during reconstruction are freed.
	 *
	 * This is useful for quick inspection or visualization without downloading the entire 3D dataset.
	 *
	 * @param singleSlice_host      Pointer to a preallocated host array (float*) that will receive the slice
	 *                              **magnitude** values. Expected length: `Nx * Ny` (real-only).
	 * @param noOfSlice             Index (1-based) of the Z-slice to retrieve, where
	 *                              `1` corresponds to the first slice and `Nz` to the last one.
	 *                              Values outside this range will produce an error.
	 *
	 * @return 0 on success; non-zero error code on failure (e.g., invalid pointer, out-of-range slice index,
	 *         or missing GPU data).
	 *
	 * @note The returned data are **real** magnitudes computed on the GPU by a dedicated kernel
	 *       (no real/imag interleaving). If you need the full complex 3D result, use
	 *       HL04_takeReconstructionAndFreeMemory().
	 * @note This call **releases all GPU memory** used by the reconstruction, similar to HL04_* functions.
	 * @note When calling from MATLAB, allocate a `single` array of size `[Ny, Nx]` (or linear length `Nx*Ny`)
	 *       and pass its pointer.
	 *
	 * @see HL03_setParamsAndStartDIandGP, HL04_takeReconstructionAndFreeMemory
	 */
	EXPORTED_FUNCTION int HL04a_takeRecSingleSliceAndFreeMemory(float* singleSlice_host, int noOfSlice);


	// =============================================================================================
	//  OPTIONAL / DIAGNOSTIC FUNCTIONS
	//  Not part of the standard reconstruction workflow.
	//  Primarily used for validation, debugging, and intermediate data inspection.
	// =============================================================================================

	/**
	 * @brief Retrieve the reference sinogram amplitude data from GPU memory.
	 *
	 * This function copies the **amplitude component** of the reference sinogram
	 * (previously generated from the reference hologram via HL_addReference)
	 * from GPU device memory to host memory.
	 *
	 * It is primarily intended for diagnostic or debugging purposes,
	 * or when the user wishes to inspect or reuse the reference sinogram amplitude
	 * in external software (e.g., MATLAB or Python).
	 *
	 * @param sinAmpRef_host       Pointer to a preallocated host array (float*)
	 *                             that will receive the amplitude data.
	 *                             The expected size of the array is determined
	 *                             by the dimensions of the reference dataset:
	 *                             `_X * _Y * _nproj` elements, as specified
	 *                             when calling HL_addReference().
	 *
	 * @return 0 on success; non-zero error code if data transfer fails or
	 *         if the reference sinogram is not available in GPU memory.
	 *
	 * @note The function does not modify or free the reference data stored on the GPU.
	 *       The reference sinogram remains available until HL_removeReference() is called.
	 * @note When calling from MATLAB, provide a single-precision (`single`) preallocated array
	 *       of appropriate size (e.g., using `libpointer('singlePtr', zeros(...,'single'))`).
	 *
	 * @see HL_addReference, HL00_B_optionTakeSinoPhRef, HL_removeReference
	 */
	EXPORTED_FUNCTION int HL00_B_optionTakeSinoAmpRef(float* sinAmpRef_host);


	/**
	 * @brief Retrieve the reference sinogram phase data from GPU memory.
	 *
	 * This function copies the **phase component** of the reference sinogram
	 * (previously generated from the reference hologram via HL_addReference)
	 * from GPU device memory to host memory.
	 *
	 * It is primarily intended for diagnostic or verification purposes,
	 * or when the user wishes to analyze or visualize the reference sinogram phase
	 * in external software (e.g., MATLAB or Python).
	 *
	 * @param sinPhRef_host        Pointer to a preallocated host array (float*)
	 *                             that will receive the phase data.
	 *                             The expected size of the array is determined
	 *                             by the dimensions of the reference dataset:
	 *                             `_X * _Y * _nproj` elements, as specified
	 *                             when calling HL_addReference().
	 *
	 * @return 0 on success; non-zero error code if data transfer fails or
	 *         if the reference sinogram is not available in GPU memory.
	 *
	 * @note The function does not modify or free the reference data stored on the GPU.
	 *       The reference sinogram remains available until HL_removeReference() is called.
	 * @note Ensure that HL_addReference() was successfully executed prior to calling this function.
	 * @note When calling from MATLAB, provide a single-precision (`single`) preallocated array
	 *       of appropriate size (e.g., using `libpointer('singlePtr', zeros(...,'single'))`).
	 *
	 * @see HL_addReference, HL00_B_optionTakeSinoAmpRef, HL_removeReference
	 */
	EXPORTED_FUNCTION int HL00_B_optionTakeSinoPhRef(float* sinPhRef_host);


	/**
	 * @brief Retrieve the generated K-space data from GPU memory to the host.
	 *
	 * This function copies the complex-valued K-space data (Fourier domain representation of the object)
	 * from GPU device memory to the host memory. It allows the user to access the spectral data
	 * after K-space generation (performed by HL02_sendDataAndGenerateKO or HL00to02_FromPreprocToGenKO),
	 * for verification, debugging, or visualization purposes.
	 *
	 * The K-space data represents the 3D Fourier transform of the measured sinograms,
	 * containing both amplitude and phase information stored as interleaved real and imaginary values.
	 *
	 * @param _complexData_host     Pointer to a preallocated host array (float*) that will receive
	 *                              the complex K-space data.
	 *                              The expected array length is `2 * _K_xy * _K_xy * _K_z`,
	 *                              where factor 2 accounts for real and imaginary components
	 *                              stored in alternating order: `[Re0, Im0, Re1, Im1, …]`.
	 *
	 * @return 0 on success; non-zero error code on failure (e.g., missing data in GPU memory
	 *         or invalid pointer).
	 *
	 * @note The K-space data must already exist in GPU memory before calling this function.
	 *       It is generated by either HL02_sendDataAndGenerateKO() (multi-step mode)
	 *       or HL00to02_FromPreprocToGenKO() (simplified mode mode).
	 * @note The function does not modify or free GPU memory — it only performs data transfer.
	 * @note When calling from MATLAB, ensure the host array is of type `single`
	 *       and has sufficient size to hold all real and imaginary components.
	 *
	 * @see HL02_sendDataAndGenerateKO, HL00to02_FromPreprocToGenKO,
	 *      HL02_B_optionTakeEW, HL02_C_optionFreeGPUmemory
	 */
	EXPORTED_FUNCTION int HL02_B_optionTakeKO(float* _complexData_host);


	/**
	 * @brief Retrieve the Ewald-sphere index map or equivalent wave-vector mask from GPU memory.
	 *
	 * This function copies the integer array representing the Ewald-sphere geometry
	 * wave-vector indices (EW array) from GPU device memory to the host. The EW array encodes
	 * the mapping between the K-space coordinates and the illumination geometry used during
	 * optical diffraction tomography (ODT) data acquisition.
	 *
	 * In practical terms, it allows the user to inspect which regions of the 3D K-space are filled
	 * with measured data, and which remain empty due to the limited angular coverage of the system.
	 * Such information is particularly useful for diagnostic visualization or reconstruction validation.
	 *
	 * @param ew_host               Pointer to a preallocated integer array (int*) on the host side,
	 *                              which will receive the EW (Ewald-sphere) index map.
	 *                              The expected array length corresponds to the full size of the
	 *                              generated K-space grid: `_K_xy * _K_xy * _K_z`.
	 *
	 * @return 0 on success; non-zero error code on failure (e.g., EW data not available in GPU memory
	 *         or invalid pointer).
	 *
	 * @note The EW array must already exist in GPU memory prior to calling this function.
	 *       It is generated internally by HL02_sendDataAndGenerateKO() or
	 *       HL00to02_FromPreprocToGenKO() during K-space construction.
	 * @note Each integer value in the EW array typically corresponds to a specific projection index
	 *       (illumination angle) or may encode a binary occupancy flag, depending on implementation.
	 * @note The function only performs a data transfer; it does not modify or free GPU memory.
	 * @note When calling from MATLAB, ensure that the receiving array is of type `int32`.
	 *
	 * @see HL02_sendDataAndGenerateKO, HL02_B_optionTakeKO, HL02_C_optionFreeGPUmemory
	 */
	EXPORTED_FUNCTION int HL02_B_optionTakeEW(int* ew_host);

	/**
	 * @brief Free GPU memory allocated during K-space generation.
	 *
	 * This function releases all GPU memory buffers that were allocated for
	 * K-space generation and related intermediate data structures (e.g., sinograms,
	 * masks, Ewald-sphere maps). It should be called when the user has completed
	 * data processing at the K-space stage and no further reconstruction (DI/GP)
	 * will be performed.
	 *
	 * @return 0 on success; non-zero error code on failure (e.g., invalid device context
	 *         or already released memory).
	 *
	 * @note After calling this function, any attempt to access previously generated
	 *       K-space or Ewald-sphere data (e.g., via HL02_B_optionTakeKO or HL02_B_optionTakeEW)
	 *       will result in an error, as the associated GPU memory has been deallocated.
	 * @note If the reconstruction step (HL03_setParamsAndStartDIandGP) will be executed,
	 *       do not call this function — GPU memory will be automatically reused and released
	 *       by HL04_takeReconstructionAndFreeMemory().
	 * @note It is safe to call this function multiple times; redundant calls will have no effect
	 *       beyond returning an appropriate status code.
	 *
	 * @see HL02_sendDataAndGenerateKO, HL02_B_optionTakeKO, HL02_B_optionTakeEW,
	 *      HL04_takeReconstructionAndFreeMemory
	 */
	EXPORTED_FUNCTION int HL02_C_optionFreeGPUmemory();


	// =============================================================================================
	//  UTILITY FUNCTIONS
	// =============================================================================================

	/**
	 * @brief Query the current GPU memory usage and retrieve the amount of free and total device memory.
	 *
	 * This utility function provides information about the current GPU memory state.
	 * It returns the total and available (free) memory on the active CUDA device in megabytes (MB).
	 * This information is useful for estimating whether the GPU has sufficient memory
	 * to handle upcoming datasets or for diagnostic logging during processing.
	 *
	 * @param freeMB    Pointer to an integer variable that will receive the amount of free GPU memory [MB].
	 * @param totalMB   Pointer to an integer variable that will receive the total GPU memory available [MB].
	 *
	 * @return 0 on success; non-zero error code if the device is not initialized or memory information
	 *         cannot be queried.
	 *
	 * @note The reported memory values are queried directly from the CUDA driver using
	 *       `cudaMemGetInfo()` and converted to megabytes (1 MB = 1024 × 1024 bytes).
	 * @note To ensure valid results, this function should be called only after successful
	 *       CUDA initialization (e.g., via `cudaInitDev()` or any HL00/HL01 function).
	 * @note The values returned may change during execution as GPU memory is allocated or freed
	 *       by other processes or CUDA streams.
	 *
	 * @par MATLAB usage:
	 * @code
	 * freeMB  = libpointer('int32Ptr', 0);
	 * totalMB = libpointer('int32Ptr', 0);
	 * err = calllib('CUDAprocessing', 'memoryInfo', freeMB, totalMB);
	 * fprintf('Free: %d MB / Total: %d MB\n', freeMB.Value, totalMB.Value);
	 * @endcode
	 *
	 * @see cudaInitDev
	 */
	EXPORTED_FUNCTION int memoryInfo(int* freeMB, int* totalMB);


	/**
	 * @brief Initialize the CUDA device and create a GPU computation context.
	 *
	 * Initializes the CUDA runtime environment on the default GPU device (device 0)
	 * and prepares it for further computations. It verifies that a CUDA-capable
	 * device is present and properly configured.
	 *
	 * This function should be called explicitly only when using **low-level or diagnostic functions**
	 * (e.g., `memoryInfo`, `prepareKspaceGeneration`, `sendData_logical`, etc.).
	 * When using any **high-level (HL00–HL04)** reconstruction functions,
	 * CUDA initialization is performed automatically by the library and does not require
	 * manual invocation of this function.
	 *
	 * @return Status code:
	 *         - `0` — success, CUDA successfully initialized.
	 *         - `1` — no compatible CUDA device detected.
	 *         - `2` — initialization failed (e.g., driver or context error).
	 *
	 * @note Safe to call multiple times; repeated calls will return `0` if the device
	 *       has already been initialized.
	 * @note This function only prepares the CUDA execution context — it does not
	 *       allocate or free GPU memory.
	 *
	 * @see memoryInfo
	 */
	EXPORTED_FUNCTION int cudaInitDev();

	// =============================================================================================
	//  LEGACY / LOW-LEVEL API
	//  Lower-level pipeline that predates the HL00-HL04 interface.
	//  Preserved for backward compatibility with existing integrations.
	//  For new projects, use the Core Workflow functions above.
	// =============================================================================================

	// --- From Reconstruction.cu ---

	/**
	 * @brief Set reconstruction volume dimensions and allocate GPU memory.
	 *
	 * Sets the lateral (x, y) and axial (z) dimensions of the reconstruction volume,
	 * allocates GPU memory for the reconstruction buffer and the optional spatial support mask,
	 * and creates a 3D cuFFT plan for in-place complex-to-complex transforms.
	 *
	 * @param x  Lateral size of the reconstruction volume in the X dimension (pixels).
	 * @param y  Lateral size of the reconstruction volume in the Y dimension (pixels).
	 * @param z  Axial depth of the reconstruction volume in the Z dimension (pixels).
	 *
	 * @return 0 on success; non-zero on failure:
	 *         - cuFFT error code if FFT plan creation fails.
	 *         - `3` if GPU memory allocation fails.
	 *
	 * @note Must be called before setParams(), sendData_n_rec(), and cudaProcessing().
	 * @note Uses a shared memory buffer for the reconstruction array and K-space (KOi_dev),
	 *       minimizing GPU memory usage by performing in-place FFT transforms.
	 *
	 * @see setParams, sendData_n_rec, sendData_KO_Re_Im_and_EW8bit, cudaProcessing
	 */
	EXPORTED_FUNCTION int setDims_n_rec(int x, int y, int z);

	/**
	 * @brief Set reconstruction parameters for the legacy GP pipeline.
	 *
	 * Configures all numerical parameters required by the legacy Gerchberg-Papoulis
	 * iterative reconstruction loop executed by cudaProcessing(). Must be called
	 * after setDims_n_rec() and before cudaProcessing().
	 *
	 * @param _do_TC      Transparency constraint flag:
	 *                    - `1` → remove imaginary part of refractive index (non-absorbing sample).
	 *                    - `0` → skip transparency enforcement.
	 * @param _do_NNC     Non-negativity / non-positivity constraint flag:
	 *                    - `> 0` → enforce non-negativity constraint.
	 *                    - `< 0` → enforce non-positivity constraint.
	 *                    - `0`   → no sign constraint applied.
	 * @param _n_imm      Refractive index of the immersion medium.
	 * @param _relaxM     Relaxation factor (0 < _relaxM ≤ 1) controlling convergence speed.
	 * @param _betaM      Amplitude constraint weighting factor.
	 * @param _kn_mean    Mean wavenumber (`n_imm / λ_mean`) used for normalization.
	 * @param _objshift   Integer axial shift of the reconstructed object volume (Z-offset).
	 * @param _dxo        Lateral pixel size in object space (µm/pixel).
	 * @param _nGPi       Number of Gerchberg-Papoulis iterations to perform.
	 *
	 * @return Always returns 0.
	 *
	 * @see setDims_n_rec, cudaProcessing, sendData_KO_Re_Im_and_EW8bit
	 */
	EXPORTED_FUNCTION int setParams(int _do_TC, int _do_NNC, float _n_imm, float _relaxM,
		float _betaM, float _kn_mean, int _objshift, float _dxo, int _nGPi);

	/**
	 * @brief Transfer a complex reconstruction volume from host to GPU memory.
	 *
	 * Copies a complex-valued 3D array representing an initial reconstruction estimate
	 * from host memory to the GPU reconstruction buffer (n_rec_dev).
	 * The array is stored as interleaved real and imaginary floats.
	 *
	 * @param n_rec_host  Pointer to a host array (float*) containing the complex volume.
	 *                    Expected size: `2 * K_xy * K_xy * K_z` floats
	 *                    (`[Re(0), Im(0), Re(1), Im(1), …]`).
	 *
	 * @return 0 on success; non-zero on failure:
	 *         - `4` if the GPU buffer has not been allocated (call setDims_n_rec() first).
	 *         - `5` if the CUDA memory transfer fails.
	 *
	 * @note The GPU buffer must already be allocated via setDims_n_rec().
	 *
	 * @see setDims_n_rec, takeComplexData, cudaProcessing
	 */
	EXPORTED_FUNCTION int sendData_n_rec(float* n_rec_host);

	/**
	 * @brief Transfer a spatial support mask from host to GPU memory.
	 *
	 * Copies a float-valued 3D spatial support mask to GPU memory.
	 * The mask constrains the reconstruction to a known physical region of the object,
	 * which accelerates convergence of the iterative algorithm.
	 * If NULL is passed, the spatial support constraint is disabled.
	 *
	 * @param _objSupport  Pointer to a host float array defining the support mask,
	 *                     where `1.0` marks voxels belonging to the object and
	 *                     `0.0` marks background.
	 *                     Expected size: `K_xy * K_xy * K_z` floats.
	 *                     Pass `NULL` to disable the support constraint.
	 *
	 * @return 0 on success; non-zero on failure:
	 *         - `4` if the GPU buffer has not been allocated.
	 *         - `5` if the CUDA memory transfer fails.
	 *
	 * @see setDims_n_rec, cudaProcessing
	 */
	EXPORTED_FUNCTION int sendData_objectSupport(float* _objSupport);

	/**
	 * @brief Convert a dense K-space and occupancy mask to a sparse representation and send to GPU.
	 *
	 * Takes a dense complex K-space array and an 8-bit binary occupancy mask,
	 * extracts only the non-zero (measured) K-space samples, scales them by `1/dxo³`,
	 * and stores both their values and indices on the GPU as sparse arrays.
	 * This sparse representation is used by the legacy cudaProcessing() loop
	 * to efficiently reinsert measured data at each iteration.
	 *
	 * @param _KO     Pointer to the dense complex K-space array (float*, interleaved Re/Im).
	 *                Expected size: `2 * K_xy * K_xy * K_z` floats.
	 * @param EW8bit  Pointer to an 8-bit binary occupancy mask (unsigned char*).
	 *                Non-zero elements indicate voxels containing measured data.
	 *                Expected size: `K_xy * K_xy * K_z` bytes.
	 * @param n       Number of non-zero elements expected in the sparse representation.
	 *                Must match the actual count of non-zero entries in EW8bit.
	 *
	 * @return 0 on success; non-zero on failure:
	 *         - `4` if host or GPU memory allocation fails.
	 *         - `5` if CUDA memory transfer fails.
	 *         - `16` if the actual number of non-zero elements does not match `n`.
	 *
	 * @note The scaling factor `dxo` must be set beforehand via setParams().
	 * @note The sparse arrays (EKO_dev, EKOind_dev) are used internally by cudaProcessing().
	 *
	 * @see setParams, cudaProcessing, sendData_n_rec
	 */
	EXPORTED_FUNCTION int sendData_KO_Re_Im_and_EW8bit(float* _KO, unsigned char* EW8bit, int n);

	/**
	 * @brief Retrieve the complex reconstruction result from GPU memory and free GPU resources.
	 *
	 * Transfers the final complex-valued 3D reconstruction from GPU memory to host,
	 * then releases all GPU memory allocated during the legacy reconstruction pipeline.
	 *
	 * @param complexData_host  Pointer to a preallocated host array (float*) that will receive
	 *                          the reconstruction. Expected size: `2 * K_xy * K_xy * K_z` floats,
	 *                          stored as interleaved `[Re(0), Im(0), Re(1), Im(1), …]`.
	 *
	 * @return 0 on success; `6` if the CUDA memory transfer fails.
	 *
	 * @note After this call, all GPU buffers are released. A new call to setDims_n_rec()
	 *       is required before processing a new dataset.
	 *
	 * @see setDims_n_rec, cudaProcessing, sendData_n_rec
	 */
	EXPORTED_FUNCTION int takeComplexData(float* complexData_host);

	/**
	 * @brief Execute the legacy Gerchberg-Papoulis iterative reconstruction loop.
	 *
	 * Runs the GP iterative algorithm using the sparse K-space representation
	 * prepared by sendData_KO_Re_Im_and_EW8bit() and the parameters set by setParams().
	 * At each iteration, constraints are applied in the reconstruction space
	 * (transparency, non-negativity/non-positivity, optional spatial support),
	 * followed by a 3D FFT, K-space data reinsertion, and inverse 3D FFT.
	 *
	 * Automatically selects between two execution paths:
	 * - **Fast mode** (fasterMode=1): used when no spatial support is active and both
	 *   transparency and sign constraints are enabled; fuses multiple operations
	 *   into fewer kernel calls for improved performance.
	 * - **Standard mode** (fasterMode=0): applies all constraints as separate steps,
	 *   supporting the full range of parameter combinations.
	 *
	 * @return 0 on success; non-zero CUDA or cuFFT error code on failure.
	 *
	 * @note Requires prior calls to: setDims_n_rec(), setParams(),
	 *       sendData_n_rec(), sendData_KO_Re_Im_and_EW8bit().
	 * @note For new projects, prefer HL03_setParamsAndStartDIandGP(), which provides
	 *       the same functionality through the high-level HL00-HL04 interface.
	 *
	 * @see setParams, sendData_n_rec, sendData_KO_Re_Im_and_EW8bit, takeComplexData,
	 *      HL03_setParamsAndStartDIandGP
	 */
	EXPORTED_FUNCTION int cudaProcessing();

	/**
	 * @brief Execute a single forward 3D FFT on the reconstruction buffer.
	 *
	 * Performs one in-place forward complex-to-complex 3D cuFFT transform,
	 * moving data from the reconstruction space (n_rec_dev) to the K-space
	 * buffer (KOi_dev) using the plan created by setDims_n_rec().
	 *
	 * @return 0 on success; `13` if the cuFFT execution fails.
	 *
	 * @note This function is a low-level building block of the legacy pipeline.
	 *       It is typically called as part of a manual reconstruction sequence,
	 *       not during normal cudaProcessing() execution.
	 *
	 * @see setDims_n_rec, cudaProcessing
	 */
	EXPORTED_FUNCTION int cudaFFT_n_rec_2_KOi();

	// --- From PreprocessingAndKspace.cu ---

	/**
	 * @brief Initialize the K-space generation pipeline and allocate GPU memory.
	 *
	 * Sets all geometric and optical parameters required for K-space generation,
	 * allocates all necessary GPU buffers (K-space, projections, masks, wave vectors),
	 * computes Ewald sphere radii and illumination wave vectors for each projection,
	 * calculates pixel-level lateral shifts (shiftX, shiftY) in K-space,
	 * and creates a batched 2D cuFFT plan for projection Fourier transforms.
	 *
	 * This is the legacy counterpart of HL01_setParams() and must be called before
	 * sendSinograms(), sendFpmask_logical(), and generateKspace().
	 *
	 * @param _K_xy               Lateral size of the 3D K-space grid (X and Y dimensions, pixels).
	 * @param _K_z                Axial size of the 3D K-space grid (Z dimension, pixels).
	 * @param _dx                 Lateral pixel size in object space (µm/pixel).
	 * @param _n_imm              Refractive index of the immersion medium.
	 * @param _n_proj             Number of projection angles.
	 * @param _Nx                 Number of pixels in X per projection (detector width).
	 * @param _Ny                 Number of pixels in Y per projection (detector height).
	 * @param _dkP                K-space sampling step in the lateral (XY) dimensions.
	 * @param _dkPz               K-space sampling step in the axial (Z) dimension.
	 * @param _NA                 Numerical aperture of the imaging objective.
	 * @param _lambda_all         Pointer to an array of illumination wavelengths (µm),
	 *                            one value per projection.
	 * @param _kxp                Pointer to an array of X-components of illumination
	 *                            wave vectors, one value per projection.
	 * @param _kyp                Pointer to an array of Y-components of illumination
	 *                            wave vectors, one value per projection.
	 * @param _approxBornNotRytov Scattering model selector:
	 *                            - `0` → Rytov approximation.
	 *                            - non-zero → Born approximation.
	 *
	 * @return 0 on success; non-zero error code on failure:
	 *         - `51` if memory allocation fails.
	 *         - `52` if CUDA data transfer fails.
	 *         - `55` if cuFFT plan creation fails.
	 *         - `59` if wave vector or shift transfer fails.
	 *
	 * @note The legacy counterpart of HL01_setParams(). For new projects,
	 *       prefer the HL01_setParams() + HL02_sendDataAndGenerateKO() workflow.
	 *
	 * @see sendSinograms, sendFpmask_logical, generateKspace, freeGPUmemory,
	 *      HL01_setParams
	 */
	EXPORTED_FUNCTION int prepareKspaceGeneration(int _K_xy, int _K_z, float _dx, float _n_imm,
		int _n_proj, int _Nx, int _Ny, float _dkP, float _dkPz, float _NA,
		float* _lambda_all, float* _kxp, float* _kyp, int _approxBornNotRytov);

	/**
	 * @brief Transfer amplitude and phase sinograms from host to GPU memory.
	 *
	 * Copies preprocessed sinogram data (amplitude and phase) to the GPU buffers
	 * allocated by prepareKspaceGeneration(). Both arrays must be provided and
	 * have the same layout: `Nx × Ny × n_proj` floats, with projections stored sequentially.
	 *
	 * @param _sinoAmp  Pointer to the amplitude sinogram (float*).
	 *                  Size: `Nx * Ny * n_proj` floats.
	 * @param _sinoPh   Pointer to the phase sinogram (float*).
	 *                  Size: `Nx * Ny * n_proj` floats.
	 *
	 * @return 0 on success; `53` if any pointer is NULL, buffers are not allocated,
	 *         or the CUDA memory transfer fails.
	 *
	 * @note prepareKspaceGeneration() must be called before this function.
	 *
	 * @see prepareKspaceGeneration, sendFpmask_logical, generateKspace
	 */
	EXPORTED_FUNCTION int sendSinograms(float* _sinoAmp, float* _sinoPh);

	/**
	 * @brief Convert a logical Fourier-plane mask to bit-packed format and transfer to GPU.
	 *
	 * Takes a binary (logical) Fourier-plane mask stored as one byte per element
	 * and converts it to a bit-packed representation where each 32-bit unsigned integer
	 * encodes 32 consecutive mask elements. The packed mask is then transferred to GPU memory.
	 * This compact representation reduces memory usage and improves access efficiency
	 * during K-space generation.
	 *
	 * @param Fpmask  Pointer to the logical mask (unsigned char*), with one element per pixel.
	 *                Non-zero values indicate active (included) pixels.
	 *                Size: `Nx * Ny * n_proj` bytes.
	 *
	 * @return 0 on success; `58` if the CUDA memory transfer fails.
	 *
	 * @note prepareKspaceGeneration() must be called before this function
	 *       to allocate the GPU mask buffer.
	 *
	 * @see prepareKspaceGeneration, sendSinograms, generateKspace
	 */
	EXPORTED_FUNCTION int sendFpmask_logical(unsigned char* Fpmask);

	/**
	 * @brief Execute the full K-space generation pipeline on the GPU.
	 *
	 * Processes the sinogram data previously transferred by sendSinograms() and
	 * populates the 3D K-space (KO_dev) by mapping each projection onto its
	 * corresponding Ewald sphere. The pipeline performs the following steps:
	 * applies the Born or Rytov approximation (selected via prepareKspaceGeneration()),
	 * applies 2D FFT shifts, executes batched 2D FFTs on all projections,
	 * scales by `dx²`, performs a parallel scatter operation to map projection
	 * samples onto Ewald spheres using atomic additions, and normalizes the
	 * accumulated K-space values by the per-voxel occupancy count.
	 *
	 * @return 0 on success; non-zero error code on failure.
	 *
	 * @note Requires prior calls to prepareKspaceGeneration(), sendSinograms(),
	 *       and sendFpmask_logical().
	 * @note The generated K-space data can be retrieved via takeDataKspace()
	 *       or used directly for reconstruction via sendData_KO_Re_Im_and_EW8bit().
	 * @note The legacy counterpart of HL02_sendDataAndGenerateKO(). For new projects,
	 *       prefer the HL01–HL02 workflow.
	 *
	 * @see prepareKspaceGeneration, sendSinograms, takeDataKspace, takeDataEW,
	 *      HL02_sendDataAndGenerateKO
	 */
	EXPORTED_FUNCTION int generateKspace();

	/**
	 * @brief Retrieve the complex hologram field from GPU to host.
	 *
	 * Copies the intermediate complex hologram field (complexFieldCSG_device),
	 * computed during GPU preprocessing, from device to host memory.
	 * Useful for diagnostic inspection of the preprocessing intermediate result.
	 *
	 * @param complexData_host  Pointer to a preallocated host array (float*).
	 *                          Size: `2 * D * D * n_proj` floats
	 *                          (interleaved real and imaginary components).
	 *
	 * @return 0 on success; `56` if the CUDA memory transfer fails.
	 *
	 * @note D is the FFT window size computed during preprocessing
	 *       (set internally by preproc_sendHologram()).
	 *
	 * @see preproc_sendHologram, takeProjections
	 */
	EXPORTED_FUNCTION int takeDataHoloCSGdev(float* complexData_host);

	/**
	 * @brief Retrieve the generated 3D K-space data from GPU to host.
	 *
	 * Copies the complex-valued 3D K-space array (KO_dev) from GPU device
	 * memory to host memory after K-space generation has been completed.
	 * Functionally equivalent to HL02_B_optionTakeKO() but part of the legacy API.
	 *
	 * @param complexData_host  Pointer to a preallocated host array (float*).
	 *                          Size: `2 * K_xy * K_xy * K_z` floats
	 *                          (interleaved real and imaginary components).
	 *
	 * @return 0 on success; `56` if the CUDA memory transfer fails.
	 *
	 * @see generateKspace, takeDataEW, HL02_B_optionTakeKO
	 */
	EXPORTED_FUNCTION int takeDataKspace(float* complexData_host);

	/**
	 * @brief Retrieve integer-valued peak position data from GPU to host.
	 *
	 * Copies the integer peak position array (peaksStep2_dev), containing
	 * the pixel-level detected peak coordinates for each projection,
	 * from GPU memory to host. Used for diagnostic inspection of the
	 * peak detection step during preprocessing.
	 *
	 * @param peaks_host  Pointer to a preallocated host array (int*).
	 *                    Size: `2 * n_proj` integers
	 *                    (X and Y peak coordinates for each projection).
	 *
	 * @return 0 on success; `57` if the CUDA memory transfer fails.
	 *
	 * @see takePeaksFloat, preproc_sendHologram
	 */
	EXPORTED_FUNCTION int takePeaksInt(int* peaks_host);

	/**
	 * @brief Retrieve float-valued peak position data from GPU to host.
	 *
	 * Copies the float peak position array (peaksStep1_dev), containing
	 * sub-pixel-resolution detected peak profiles for each projection,
	 * from GPU memory to host. Used for diagnostic inspection of the
	 * sub-pixel peak detection step during preprocessing.
	 *
	 * @param peaks_host  Pointer to a preallocated host array (float*).
	 *                    Size: `2 * D * n_proj` floats.
	 *
	 * @return 0 on success; `57` if the CUDA memory transfer fails.
	 *
	 * @see takePeaksInt, preproc_sendHologram
	 */
	EXPORTED_FUNCTION int takePeaksFloat(float* peaks_host);

	/**
	 * @brief Retrieve the Ewald sphere occupancy map from GPU to host.
	 *
	 * Copies the integer Ewald sphere occupancy/index map (EW_dev) from GPU
	 * memory to host after K-space generation. Each element stores the
	 * accumulated number of projection samples that contributed to the
	 * corresponding K-space voxel. Functionally equivalent to HL02_B_optionTakeEW()
	 * but part of the legacy API.
	 *
	 * @param EW_host  Pointer to a preallocated host integer array.
	 *                 Size: `K_xy * K_xy * K_z` integers.
	 *
	 * @return 0 on success; `57` if the CUDA memory transfer fails.
	 *
	 * @see generateKspace, takeDataKspace, HL02_B_optionTakeEW
	 */
	EXPORTED_FUNCTION int takeDataEW(int* EW_host);

	/**
	 * @brief Retrieve the complex projection data from GPU to host.
	 *
	 * Copies the complex projection array (projections_dev) from GPU memory
	 * to host after preprocessing or K-space generation. Contains the
	 * Fourier-transformed projections before mapping onto Ewald spheres.
	 * Useful for diagnostic inspection at intermediate stages of the pipeline.
	 *
	 * @param complexData_host  Pointer to a preallocated host array (float*).
	 *                          Size: `2 * Nx * Ny * n_proj` floats
	 *                          (interleaved real and imaginary components).
	 *
	 * @return 0 on success; `58` if the CUDA memory transfer fails.
	 *
	 * @see generateKspace, takeDataKspace, preproc_sendHologram
	 */
	EXPORTED_FUNCTION int takeProjections(float* complexData_host);

	/**
	 * @brief Release all GPU memory allocated by the legacy K-space generation pipeline.
	 *
	 * Calls freeDeviceMemory() to release all GPU buffers allocated during
	 * prepareKspaceGeneration(), including K-space, projection, mask,
	 * wave vector, and auxiliary arrays.
	 *
	 * @return Always returns 0.
	 *
	 * @note Should be called after all data has been retrieved from the GPU
	 *       and no further processing will be performed in the current session.
	 * @note The legacy counterpart of HL02_C_optionFreeGPUmemory() and
	 *       HL04_takeReconstructionAndFreeMemory().
	 *
	 * @see prepareKspaceGeneration, HL02_C_optionFreeGPUmemory,
	 *      HL04_takeReconstructionAndFreeMemory
	 */
	EXPORTED_FUNCTION int freeGPUmemory();

	/**
	 * @brief Initiate GPU-based preprocessing from raw hologram data.
	 *
	 * Sets all optical and acquisition parameters and transfers raw detector
	 * hologram frames to the GPU for preprocessing. Enables GPU preprocessing mode,
	 * computes the FFT window size D based on the optical parameters, and prepares
	 * internal buffers for subsequent sinogram extraction and K-space generation.
	 *
	 * This is the legacy counterpart of HL00to02_FromPreprocToGenKO() — it handles
	 * only the parameter setup and data transfer phase, without running
	 * preprocessing or K-space generation. Additional calls to generateKspace()
	 * and related functions are required to complete the pipeline.
	 *
	 * @param hologram          Pointer to raw hologram frames (short int*).
	 *                          Layout: `X × Y × _nproj` int16 samples.
	 * @param _X                Detector width in pixels.
	 * @param _Y                Detector height in pixels.
	 * @param _nproj            Number of projection angles.
	 * @param _NA               Numerical aperture of the imaging objective.
	 * @param lambda            Illumination wavelength (µm).
	 * @param cam_pix           Camera pixel size (µm).
	 * @param M                 Transverse magnification of the imaging system.
	 * @param _n_imm            Refractive index of the immersion medium.
	 * @param _do_NNC           Sign constraint flag (see HL_addReference for details).
	 * @param _cosFactor        Tukey window tapering factor (0 = rectangular, 1 = full cosine).
	 * @param _fftWindowsScale  Scaling factor for the FFT window size.
	 *
	 * @return 0 on success; non-zero error code on failure.
	 *
	 * @note For new projects, prefer HL00to02_FromPreprocToGenKO(), which combines
	 *       preprocessing and K-space generation in a single high-level call.
	 *
	 * @see generateKspace, freeGPUmemory, HL00to02_FromPreprocToGenKO
	 */
	EXPORTED_FUNCTION int preproc_sendHologram(short int* hologram, int _X, int _Y, int _nproj,
		float _NA, float lambda, float cam_pix, float M, float _n_imm,
		int _do_NNC, float _cosFactor, float _fftWindowsScale);

	/**
	 * @brief Retrieve the unwrapped phase sinogram from GPU to host.
	 *
	 * Copies the unwrapped phase sinogram (sinoPh_dev) from GPU memory to host.
	 * This array contains the phase projections after preprocessing (phase unwrapping,
	 * tilt removal, normalization) and before K-space generation.
	 * Useful for diagnostic inspection of the phase preprocessing result.
	 *
	 * @param unwrapped_host  Pointer to a preallocated host array (float*).
	 *                        Size: `D * D * n_proj` floats,
	 *                        where D is the FFT window size.
	 *
	 * @return 0 on success; `56` if the CUDA memory transfer fails.
	 *
	 * @see preproc_sendHologram, generateKspace, takeProjections
	 */
	EXPORTED_FUNCTION int takeUnwrappedFrame(float* unwrapped_host);



#ifdef  __cplusplus
}
#endif

#endif
