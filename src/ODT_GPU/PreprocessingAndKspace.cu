/**
 * @file PreprocessingAndKspace.cu
 * @brief GPU kernels and exported functions for preprocessing and K-space generation.
 *
 * Implements the first two stages of the LaODT processing pipeline: preprocessing
 * of raw off-axis holograms (FFT, extraction of the +1 diffraction order, phase
 * unwrapping, amplitude normalization, phase tilt removal) and generation of the
 * 3D K-space by mapping the processed projections onto Ewald spheres.
 *
 * Part of the ODT_GPU library.
 * Distributed under the GNU General Public License v3.0 (GPL-3.0).
 * See the LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html for details.
 *
 * @see https://github.com/LAB-UMK/ODT_GPU_Library
 */

#define EXPORT_FCNS
#include "helper.h"

#include "pch.h"
#include "ODT_GPU.h"
#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuComplex.h>
#include <stdio.h>
#include <malloc.h>
#include <Windows.h>
#include <time.h>
#include "device_atomic_functions.h"
#include <thrust/complex.h>
#include <thrust/sort.h>

int n_proj; //number of projections to be inserted into the K matrix
float dkP;
float dkPz;
int K_xy;//size of K-space in x and y
int K_z;// size of K-space in z
extern float n_imm;
extern int do_NNC;
extern float kn_mean;

extern int ifftShift_2(int x, int y, int z, cufftComplex* data_complex_dev);
extern int ifftShift_1(int x, int y, int z, cufftComplex* data_complex_dev);

cufftComplex* KO_dev = NULL; //K-space (KO)
cufftComplex* projections_dev = NULL; // projections after the Rytov approximation
unsigned int* fpmask_dev = NULL; // projection mask, bitwise representation
//cufftComplex* projectionsInSizeK_dev = NULL; // projections after the Rytov approximation
int* EW_dev = NULL; //counter of "entries" into K-space
float* rf_dev = NULL; // vector of Ewald sphere radii for the projections
float* kn_dev = NULL; // wave vector values for the projections
float* kz_dev = NULL; // z-shift in K space due to xy shifts, shifts of the center of Ewald spheres for all projections based on kxp and kyp, that are the locations of peaks within aperture
float* kxp_dev = NULL; // x-coordinates of illumination vectors
float* kyp_dev = NULL; // y-coordinates of illumination vectors
float* kxy_dev = NULL; // length of illumination vectors
float NA = 0; //apertura numeryczna obiektywu
//float* lambda_all = NULL; // vector of wavelengths for each projection
float dx = 0; //the x size on the projection

float* sinoAmp_dev = NULL; //singoram amplituda
float* sinoPh_dev = NULL; //singoram faza
cufftComplex* sinoPhComplex_dev = NULL;
float* sinoAmpRef_dev = NULL; //singoram amplituda referencji
float* sinoPhRef_dev = NULL; //singoram faza referencji
cufftComplex* sinoPhRefComplex_dev = NULL; //singoram faza referencji
int Nx = 0; // projection size x
int Ny = 0; // projection size y
int Nx_Fpmask; // x size in the 32-bit matrix for the mask (bitwise values)
float* kxp_host = NULL; // x-coordinates of illumination vectors
float* kyp_host = NULL; // y-coordinates of illumination vectors
int* xShift_dev = NULL; // xshift of the projection accordnig to the x-coordinates of illumination vectors
int* yShift_dev = NULL; // xshift of the projection accordnig to the y-coordinates of illumination vectors

int approxBornNotRytov; // 0 (default) - Rytov, 1 - Born 

cufftHandle planFFTsino = NULL;
cufftHandle planFFTholo = NULL;
cufftHandle planFFTcomplexField = NULL;

// preprocessing
short int* holograms_device = NULL; //hologram data
cufftComplex* hologramsCSG_device = NULL; //hologram data
cufftComplex* croppedPartsHoloCSG_device = NULL; //hologram data
cufftComplex* complexFieldCSG_device = NULL; //hologram data
float* peaksStep1_dev = NULL; // detected peaks in croppedPartsHolo w,x,w,x... for each line (value, x-coordinate...)
int* peaksStep2_dev = NULL; // final coordinates of detected peaks in croppedPartsHolo x,y,x,y,x... for each frame
float* tukeyWin2D_W_dev = NULL;
float* sMask_dev = NULL;
float* mask_dev = NULL;
int X, Y;  //hologram dimensions
int Nx_holo, Ny_holo; //hologram dimmensions after crop

int D = 0; // dimmension of crooped data (square) from hologram
int Dminus = 0; // reduced by 8 pix (4 pix each side) dimmension of crooped data (square) from hologram (for Fast Phase Unwrapping)
bool preProcOnGPU = false; // indicates whether preprocessing is computed on the GPU
bool referenceCalculated = false; //indicates whether the reference has been calculated

#pragma region fast_Phase_Unwrapping

__global__ void kernel_ifftShift_1_float(float* data_complex_dev)
{
	register int index = blockIdx.x * 2 * blockDim.x + threadIdx.x;

	register int offset = blockDim.x;

	float re = data_complex_dev[index];

	data_complex_dev[index] = data_complex_dev[index + offset];

	data_complex_dev[index + offset] = re;
}

int ifftShift_1_float(float* data_float_dev)
{
	int nThreads = D / 2;
	dim3 nBlocks(D, 1, 1);

	kernel_ifftShift_1_float << <nBlocks, nThreads >> > (data_float_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) && (cudaStatus2 != cudaSuccess)) {
		return 11;
	}
	return 0;
}

__global__ void kernel_ifftShift_2_float(float* data_float_dev)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	register int offset = blockDim.x * gridDim.x;

	float re = data_float_dev[index];

	data_float_dev[index] = data_float_dev[index + offset];

	data_float_dev[index + offset] = re;
}

int ifftShift_2_float(float* data_float_dev)
{
	int nThreads = D;
	dim3 nBlocks(D / 2, 1, 1);

	kernel_ifftShift_2_float << <nBlocks, nThreads >> > (data_float_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) && (cudaStatus2 != cudaSuccess)) {
		return 11;
	}
	return 0;
}


__global__ void kernel_calc_matrix_K_2D(float* K_dev)
{
	//k - threadIdx.x
	//j - blockIdx.x
	//i - blockIdx.y

	//pix - blockDim.x
	//ascans - gridDim.x
	//bscans - gridDim.y

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	/// old Eweliny
	/*K_dev[index] = (1.0 / (float(blockDim.x - 1))*(-float(blockDim.x) / 2 + threadIdx.x)) *(1.0 / (float(blockDim.x - 1))*(-float(blockDim.x) / 2 + threadIdx.x))
		+ (1.0 / (float(gridDim.x - 1))*(-float(gridDim.x) / 2 + blockIdx.x))*(1.0 / (float(gridDim.x - 1))*(-float(gridDim.x) / 2 + blockIdx.x)) + 2.0e-16;*/
	int a = threadIdx.x - (blockDim.x >> 1);
	int b = blockIdx.x - (gridDim.x >> 1);
	K_dev[index] = a * a + b * b + 0.01;

}

int calc_matrix_K_2D(float* K_dev)
{
	int nThreads = D;
	dim3 nBlocks(D, 1, 1);

	kernel_calc_matrix_K_2D << <nBlocks, nThreads >> > (K_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) && (cudaStatus2 != cudaSuccess)) {
		return 83;
	}
	return 0;
}

__global__ void P_vec(cufftComplex* P, float* W, cufftComplex* complex_data_in)
{
	register int i = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	register int input = i;


	register float re = complex_data_in[input].x;
	register float im = complex_data_in[input].y;

	//register float ampl = sqrtf(re*re + im * im);
	//amplitude[i] = ampl;

	W[i] = atan2(im, re);

	P[i].x = cos(W[i]);//re / ampl;//
	P[i].y = sin(W[i]);//im / ampl;//
}

__global__ void P_vec_fromPhase(cufftComplex* P, float* W, float* phase_data_in)
{
	register int i = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	register int input = i;

	W[i] = phase_data_in[i];

	P[i].x = cos(W[i]);
	P[i].y = sin(W[i]);
}


__global__ void P_vec_fromPhase2(cufftComplex* P, float* phase_data_in)
{
	register int i = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	register int input = i;

	float W = phase_data_in[i];

	P[i].x = cos(W);
	P[i].y = sin(W);
}

__global__ void product_ComplexFloatD(cufftComplex* C, float* K_dev)
{
	register int i_loc = blockIdx.x * blockDim.x + threadIdx.x;

	register int i_glob = blockIdx.y * blockDim.x * gridDim.x + i_loc;

	C[i_glob].x = C[i_glob].x * K_dev[i_loc];
	C[i_glob].y = C[i_glob].y * K_dev[i_loc];
}

__global__ void signCorrection(cufftComplex* C)
{
	register int i_loc = blockIdx.x * blockDim.x + threadIdx.x;

	register int i_glob = blockIdx.y * blockDim.x * gridDim.x + i_loc;

	register int sign = (threadIdx.x + blockIdx.x) & 0x01;
	sign = 2 * sign - 1;
	sign = -sign;

	C[i_glob].x = C[i_glob].x * sign;
	C[i_glob].y = C[i_glob].y * sign;
}

__global__ void Imag_z1byz2_v2(cufftComplex* A, cufftComplex* B, float allPix_c)
{
	//imag(z1/z2) = im2re1-im1re2; when |z2|=1;
	register int i_loc = blockIdx.x * blockDim.x + threadIdx.x;

	register int i_glob = blockIdx.y * blockDim.x * gridDim.x + i_loc;

	A[i_glob].x = A[i_glob].y * B[i_glob].x / allPix_c - A[i_glob].x * B[i_glob].y / allPix_c;
	A[i_glob].y = 0;
}

__global__ void Imag_z1byz2(cufftComplex* A, float* B, float allPix_c)
{

	//imag(z1/z2) = im2re1-im1re2; when |z2|=1;
	register int i_loc = blockIdx.x * blockDim.x + threadIdx.x;

	register int i_glob = blockIdx.y * blockDim.x * gridDim.x + i_loc;

	A[i_glob].x = A[i_glob].y * cosf(B[i_glob]) / allPix_c - A[i_glob].x * sinf(B[i_glob]) / allPix_c;
	A[i_glob].y = 0;
}

__global__ void div_ComplexFloatD(cufftComplex* C, float* K_dev)
{

	register int i_loc = blockIdx.x * blockDim.x + threadIdx.x;

	register int i_glob = blockIdx.y * blockDim.x * gridDim.x + i_loc;

	C[i_glob].x = C[i_glob].x / K_dev[i_loc];
	C[i_glob].y = C[i_glob].y / K_dev[i_loc];

	//uzasadnione dzieleniem przez ~0?
	if (i_loc == 0) {
		C[i_glob].x = 0;
		C[i_glob].y = 0;
	}
}

__global__ void unwrap(float* W, cufftComplex* PSI, float PI2, float allPix_c)
{
	register int i_loc = blockIdx.x * blockDim.x + threadIdx.x;

	register int i_glob = blockIdx.y * blockDim.x * gridDim.x + i_loc;

	W[i_glob] = W[i_glob] + PI2 * round((PSI[i_glob].x / allPix_c - W[i_glob]) / (PI2));
	//W[i_glob] = PSI[i_glob].x / allPix_c;
}

__global__ void unwrapVolkovCut_v2(float* sinoPhase_out, cufftComplex* P_vec, float PI2, float allPix_c)
{
	register int i_input_loc = blockIdx.x * 2 * blockDim.x + threadIdx.x;
	register int i_input_glob = blockIdx.y * 4 * blockDim.x * gridDim.x + i_input_loc;

	register int i_output_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int i_output_glob = blockIdx.y * blockDim.x * gridDim.x + i_output_loc;

	sinoPhase_out[i_output_glob] = P_vec[i_input_glob].x / allPix_c;
}

__global__ void unwrapVolkovCut(float* sinoPhase_out, float* W, cufftComplex* PSI, float PI2, float allPix_c)
{
	register int i_input_loc = blockIdx.x * 2 * blockDim.x + threadIdx.x;
	register int i_input_glob = blockIdx.y * 4 * blockDim.x * gridDim.x + i_input_loc;

	register int i_output_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int i_output_glob = blockIdx.y * blockDim.x * gridDim.x + i_output_loc;

	sinoPhase_out[i_output_glob] = W[i_input_glob] + PI2 * round((PSI[i_input_glob].x / allPix_c - W[i_input_glob]) / (PI2));
	//W[i_glob] = PSI[i_glob].x / allPix_c;
}

int fastPhaseUnwrapping(cufftComplex* complexFieldCSG_dev, float* sAmp_dev, float* sPhase_dev)
{
	cufftHandle plan1 = NULL, plan2 = NULL;
	int err;
	cudaError_t status;
	cufftResult_t e;

	cufftComplex* P_dev = NULL;

	float* K_dev;
	int e1, e2, e3;
	e1 = cudaMalloc((void**)&K_dev, D * D * sizeof(float));
	e3 = cudaMalloc((void**)&P_dev, D * D * n_proj * sizeof(cufftComplex));

	if ((e1 != 0) || (e3 != 0) || (P_dev == NULL) || (K_dev == NULL))
		return 82;

	/* Create a 3D FFT ComplexToComplex plan. */
	int NN_D[2] = { D,D };

	int idist_D = D * D;
	int odist_D = D * D;
	int inembed_D[] = { D, D };
	int onembed_D[] = { D, D };
	int istride_D = 1;
	int ostride_D = 1;

	e = cufftPlanMany(&plan1, 2, NN_D, inembed_D, istride_D, idist_D, onembed_D, ostride_D, odist_D, CUFFT_C2C, n_proj);
	if (e != cudaSuccess)
		return 81; // FFT plan creation error

	dim3 block(D, 1, 1);
	dim3 grid(D, n_proj, 1);

	calc_matrix_K_2D(K_dev);
	//ifftShift_1_float(K_dev);
	//ifftShift_2_float(K_dev);
												/*float* temp = (float*)malloc(D*D*sizeof(float));
												cudaError_t cudaStatus1 = cudaMemcpy(temp, K_dev, D*D * sizeof(float), cudaMemcpyDeviceToHost);
												FILE* plik8;
												plik8 = fopen("vecK.bin", "wb");
												fwrite(temp, sizeof(float), D*D, plik8);
												fclose(plik8);
												free(temp);*/


	P_vec << < grid, block >> > (P_dev, sPhase_dev, complexFieldCSG_dev);
	cudaDeviceSynchronize();
	err = cudaGetLastError();

	/*temp = (float*)malloc(2*D*D*n_proj * sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, P_dev, 2*D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("phaseBefore1FFT.bin", "wb");
	fwrite(temp, sizeof(float), 2*D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	//1
	e = cufftExecC2C(plan1, P_dev, P_dev, CUFFT_FORWARD);
	cudaDeviceSynchronize();
	ifftShift_1(D, D, n_proj, P_dev);
	ifftShift_2(D, D, n_proj, P_dev);

	/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("phaseAfter1FFT.bin", "wb");
	fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/


	product_ComplexFloatD << < grid, block >> > (P_dev, K_dev);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("afterMult.bin", "wb");
	fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	//2

	ifftShift_1(D, D, n_proj, P_dev);
	ifftShift_2(D, D, n_proj, P_dev);
	e = cufftExecC2C(plan1, P_dev, P_dev, CUFFT_INVERSE);
	cudaDeviceSynchronize();
	//signCorrection << < grid, block >> > (P_dev);


												/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
												cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
												plik8 = fopen("phaseAfter2FFT.bin", "wb");
												fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
												fclose(plik8);
												free(temp);*/

	Imag_z1byz2 << < grid, block >> > (P_dev, sPhase_dev, D * D);
	cudaDeviceSynchronize();
	err = cudaGetLastError();

	/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("phaseBefore3FFT.bin", "wb");
	fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	//3
	e = cufftExecC2C(plan1, P_dev, P_dev, CUFFT_FORWARD);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	ifftShift_1(D, D, n_proj, P_dev);
	ifftShift_2(D, D, n_proj, P_dev);

	/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("phaseAfter3FFT.bin", "wb");
	fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	div_ComplexFloatD << < grid, block >> > (P_dev, K_dev);
	cudaDeviceSynchronize();
	err = cudaGetLastError();

	/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("phaseBefore4FFT.bin", "wb");
	fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	//4
	ifftShift_1(D, D, n_proj, P_dev);
	ifftShift_2(D, D, n_proj, P_dev);
	e = cufftExecC2C(plan1, P_dev, P_dev, CUFFT_INVERSE);
	cudaDeviceSynchronize();

	/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("phaseAfter4FFT.bin", "wb");
	fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	unwrap << < grid, block >> > (sPhase_dev, P_dev, 2 * PI, D * D); //sPhase_dev - out
	cudaDeviceSynchronize();
	err = cudaGetLastError();

	cudaFree(P_dev); P_dev = NULL;
	cudaFree(K_dev); K_dev = NULL;
	cufftDestroy(plan1); plan1 = NULL;
	return 0;
}

__global__ void volkovExpand(float* phaseExpanded_dev, float* phaseIn_dev, int Dplus_D)
{
	register int index_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int index_glob = blockIdx.y * blockDim.x * gridDim.x + index_loc;

	register int in_index_loc = (blockIdx.x + Dplus_D / 2) * (blockDim.x + Dplus_D) + threadIdx.x + Dplus_D / 2;
	register int in_index_glob = blockIdx.y * (blockDim.x + Dplus_D) * (gridDim.x + Dplus_D) + in_index_loc;

	register int i_expanded_lu = blockIdx.x * 2 * blockDim.x + threadIdx.x;                     //left up
	register int i_expanded_ru = blockIdx.x * 2 * blockDim.x + 2 * blockDim.x - threadIdx.x - 1;    //right up
	register int i_expanded_ld = (2 * gridDim.x - 1 - blockIdx.x) * 2 * blockDim.x + threadIdx.x; //left down
	register int i_expanded_rd = (2 * gridDim.x - 1 - blockIdx.x) * 2 * blockDim.x + 2 * blockDim.x - threadIdx.x - 1;//right down
	register int outputFrameOffset = blockIdx.y * 4 * blockDim.x * gridDim.x;

	register float phase = phaseIn_dev[in_index_glob];

	phaseExpanded_dev[i_expanded_lu + outputFrameOffset] = phase;
	phaseExpanded_dev[i_expanded_ru + outputFrameOffset] = phase;
	phaseExpanded_dev[i_expanded_ld + outputFrameOffset] = phase;
	phaseExpanded_dev[i_expanded_rd + outputFrameOffset] = phase;
}


__global__ void volkovExpand2(float* phaseExpanded_dev, float* phaseIn_dev, int Dplus_D)
{
	register int index_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int index_glob = blockIdx.y * blockDim.x * gridDim.x + index_loc;

	register int in_index_loc = (blockIdx.x + Dplus_D / 2) * (blockDim.x + Dplus_D) + threadIdx.x + Dplus_D / 2;
	register int in_index_glob = blockIdx.y * (blockDim.x + Dplus_D) * (gridDim.x + Dplus_D) + in_index_loc;

	register int i_expanded_lu = blockIdx.x * 2 * blockDim.x + threadIdx.x;                     //left up
	register int i_expanded_ru = blockIdx.x * 2 * blockDim.x + 2 * blockDim.x - threadIdx.x - 1;    //right up
	register int i_expanded_ld = (2 * gridDim.x - 1 - blockIdx.x) * 2 * blockDim.x + threadIdx.x; //left down
	register int i_expanded_rd = (2 * gridDim.x - 1 - blockIdx.x) * 2 * blockDim.x + 2 * blockDim.x - threadIdx.x - 1;//right down
	register int outputFrameOffset = blockIdx.y * 4 * blockDim.x * gridDim.x;

	register float phase = phaseIn_dev[in_index_glob];

	phaseExpanded_dev[i_expanded_lu + outputFrameOffset] = phase;
	phaseExpanded_dev[i_expanded_ru + outputFrameOffset] = phase;
	phaseExpanded_dev[i_expanded_ld + outputFrameOffset] = phase;
	phaseExpanded_dev[i_expanded_rd + outputFrameOffset] = phase;
}

__global__ void complexFieldToAmpPhase(cufftComplex* inputData, float* amplitudeOut_dev, float* phaseOut_dev)
{
	register int index_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int index_glob = blockIdx.y * blockDim.x * gridDim.x + index_loc;

	register cufftComplex data = inputData[index_glob];

	register float re = data.x;
	register float im = data.y;

	amplitudeOut_dev[index_glob] = sqrtf(re * re + im * im);

	phaseOut_dev[index_glob] = atan2(im, re);
}

__global__ void kernel_complexPhaseToPhase(cufftComplex* inputData, float* phaseOut_dev)
{
	register int index_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int index_glob = blockIdx.y * blockDim.x * gridDim.x + index_loc;

	register cufftComplex data = inputData[index_glob];

	register float re = data.x;
	register float im = data.y;

	phaseOut_dev[index_glob] = atan2(im, re);
}

int complexPhaseToPhase(cufftComplex* inputData, float* phaseOut_dev)
{
	dim3 block1(D, 1, 1);
	dim3 grid1(D, n_proj, 1);

	kernel_complexPhaseToPhase << < grid1, block1 >> > (inputData, phaseOut_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 116;
	}

	return 0;
}


__global__ void kernel_complexFieldToComplexPhaseAndAmp(cufftComplex* inputData, cufftComplex* complexPhase, float* amplitudeOut_dev)
{
	register int index_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int index_glob = blockIdx.y * blockDim.x * gridDim.x + index_loc;

	register cufftComplex data = inputData[index_glob];

	register float re = data.x;
	register float im = data.y;

	register float modul = sqrtf(re * re + im * im);
	amplitudeOut_dev[index_glob] = modul;

	complexPhase[index_glob].x = re / modul;
	complexPhase[index_glob].y = im / modul;
}


int complexFieldToComplexPhaseAndAmp(cufftComplex* complexFieldCSG_dev, cufftComplex* pahse, float* sinoAmp_dev)
{
	dim3 block1(D, 1, 1);
	dim3 grid1(D, n_proj, 1);

	kernel_complexFieldToComplexPhaseAndAmp << < grid1, block1 >> > (complexFieldCSG_dev, pahse, sinoAmp_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 106;
	}

	return 0;
}

__global__ void kernel_complexFieldReferenceSubtraction(cufftComplex* inOutData, cufftComplex* complexPhaseRef)
{
	register int index_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int index_glob = blockIdx.y * blockDim.x * gridDim.x + index_loc;

	register cufftComplex data = inOutData[index_glob];
	register cufftComplex dataRef = complexPhaseRef[index_glob];

	register float x1 = data.x;
	register float y1 = data.y;

	register float x2Ref = dataRef.x;
	register float y2Ref = -dataRef.y; // complex conjugate

	register cufftComplex result;
	result.x = x1 * x2Ref - y1 * y2Ref;
	result.y = x1 * y2Ref + y1 * x2Ref;

	inOutData[index_glob] = result;
}


int complexFieldReferenceSubtraction(cufftComplex* inOutData, cufftComplex* complexPhaseRef)
{
	dim3 block1(D, 1, 1);
	dim3 grid1(D, n_proj, 1);

	kernel_complexFieldReferenceSubtraction << < grid1, block1 >> > (inOutData, complexPhaseRef);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 114;
	}

	return 0;
}


// assumption: the frame is square DxD; only one line from the top edge and one from the bottom are processed
#define windowSize 9 // in the shape of the letter C, around the point + this point twice
__global__ void medianFiter(float* ioData, int lineFromEdge)
{
	register int frameSize = blockDim.x * blockDim.x;

	float window[windowSize];

	int midWin = int(windowSize / 2);
	int edge = 1;

	if ((blockIdx.x & 0x03) == 0x00) {   //first horizontal line, top of the frame
		register int index_loc = blockDim.x * (lineFromEdge - 1) + threadIdx.x;
		register int index_glob = blockIdx.y * frameSize + index_loc;

		if ((threadIdx.x > (edge - 1)) && (threadIdx.x < (blockDim.x - edge)))
		{

			window[0] = ioData[index_glob];
			window[1] = ioData[index_glob + 1];
			window[2] = ioData[index_glob - 1];
			window[3] = ioData[index_glob + blockDim.x - 1];
			window[4] = ioData[index_glob + blockDim.x];
			window[5] = ioData[index_glob + blockDim.x + 1];
			window[6] = ioData[index_glob + 2 * blockDim.x - 1];
			window[7] = ioData[index_glob + 2 * blockDim.x];
			window[8] = ioData[index_glob + 2 * blockDim.x + 1];
			__syncthreads();


			int n = windowSize;
			int i = 0, j = 0;
			for (i = 0; i < n - 1; i++)
			{
				for (j = 0; j < n - i - 1; j++)
				{
					if (window[j] > window[j + 1])
					{
						float temp = window[j];
						window[j] = window[j + 1];
						window[j + 1] = temp;
					}
				}
			}

			__syncthreads();
			ioData[index_glob] = window[midWin];
		}
	}
	if ((blockIdx.x & 0x03) == 0x01) {   //second horizontal line, bottom of the frame
		register int index_loc = blockDim.x * (blockDim.x - lineFromEdge) + threadIdx.x;
		register int index_glob = blockIdx.y * frameSize + index_loc;

		if ((threadIdx.x > (edge - 1)) && (threadIdx.x < (blockDim.x - edge)))
		{
			window[0] = ioData[index_glob];
			window[1] = ioData[index_glob + 1];
			window[2] = ioData[index_glob - 1];
			window[3] = ioData[index_glob - blockDim.x - 1];
			window[4] = ioData[index_glob - blockDim.x];
			window[5] = ioData[index_glob - blockDim.x + 1];
			window[6] = ioData[index_glob - 2 * blockDim.x - 1];
			window[7] = ioData[index_glob - 2 * blockDim.x];
			window[8] = ioData[index_glob - 2 * blockDim.x + 1];

			int n = windowSize;

			do
			{
				for (int i = 0; i < n - 1; i++)
				{
					if (window[i] > window[i + 1])
					{
						float temp = window[i + 1];
						window[i + 1] = window[i];
						window[i] = temp;
					}
				}
				n = n - 1;
			} while (n > 1);

			__syncthreads();

			ioData[index_glob] = window[midWin];
		}
	}
	if ((blockIdx.x & 0x03) == 0x02) {   //third line, vertical, for smaller indices
		register int index_loc = (lineFromEdge - 1) + threadIdx.x * blockDim.x;
		register int index_glob = blockIdx.y * frameSize + index_loc;

		if ((threadIdx.x > (edge - 1)) && (threadIdx.x < (blockDim.x - edge)))
		{
			window[0] = ioData[index_glob];
			window[1] = ioData[index_glob + blockDim.x];
			window[2] = ioData[index_glob - blockDim.x];
			window[3] = ioData[index_glob + 1 + blockDim.x];
			window[4] = ioData[index_glob + 1];
			window[5] = ioData[index_glob + 1 - blockDim.x];
			window[6] = ioData[index_glob + 2 + blockDim.x];
			window[7] = ioData[index_glob + 2];
			window[8] = ioData[index_glob + 2 - blockDim.x];

			int n = windowSize;

			do
			{
				for (int i = 0; i < n - 1; i++)
				{
					if (window[i] > window[i + 1])
					{
						float temp = window[i + 1];
						window[i + 1] = window[i];
						window[i] = temp;
					}
				}
				n = n - 1;
			} while (n > 1);

			__syncthreads();

			ioData[index_glob] = window[midWin];
		}
	}

	if ((blockIdx.x & 0x03) == 0x03) {   //fourth line, vertical, for larger indices
		register int index_loc = frameSize - lineFromEdge + threadIdx.x * blockDim.x;
		register int index_glob = blockIdx.y * frameSize + index_loc;

		if ((threadIdx.x > (edge - 1)) && (threadIdx.x < (blockDim.x - edge)))
		{
			window[0] = ioData[index_glob];
			window[1] = ioData[index_glob + blockDim.x];
			window[2] = ioData[index_glob - blockDim.x];
			window[3] = ioData[index_glob - 1 + blockDim.x];
			window[4] = ioData[index_glob - 1];
			window[5] = ioData[index_glob - 1 - blockDim.x];
			window[6] = ioData[index_glob - 2 + blockDim.x];
			window[7] = ioData[index_glob - 2];
			window[8] = ioData[index_glob - 2 - blockDim.x];

			int n = windowSize;

			do
			{
				for (int i = 0; i < n - 1; i++)
				{
					if (window[i] > window[i + 1])
					{
						float temp = window[i + 1];
						window[i + 1] = window[i];
						window[i] = temp;
					}
				}
				n = n - 1;
			} while (n > 1);

			__syncthreads();

			ioData[index_glob] = window[midWin];
		}
	}

}

__global__ void coverFiter(float* ioData, int lineFromEdge)
{
	register int frameSize = blockDim.x * blockDim.x;


	int midWin = int(windowSize / 2);
	int edge = 0;

	if ((blockIdx.x & 0x03) == 0x00) {   //first horizontal line, top of the frame
		register int index_loc = blockDim.x * (lineFromEdge - 1) + threadIdx.x;
		register int index_glob = blockIdx.y * frameSize + index_loc;

		if ((threadIdx.x > (edge - 1)) && (threadIdx.x < (blockDim.x - edge)))
		{
			ioData[index_glob] = ioData[index_glob + blockDim.x];
		}
	}
	if ((blockIdx.x & 0x03) == 0x01) {   //second horizontal line, bottom of the frame
		register int index_loc = blockDim.x * (blockDim.x - lineFromEdge) + threadIdx.x;
		register int index_glob = blockIdx.y * frameSize + index_loc;

		if ((threadIdx.x > (edge - 1)) && (threadIdx.x < (blockDim.x - edge)))
		{
			ioData[index_glob] = ioData[index_glob - blockDim.x];
		}
	}
	if ((blockIdx.x & 0x03) == 0x02) {   //third line, vertical, for smaller indices
		register int index_loc = (lineFromEdge - 1) + threadIdx.x * blockDim.x;
		register int index_glob = blockIdx.y * frameSize + index_loc;

		if ((threadIdx.x > (edge - 1)) && (threadIdx.x < (blockDim.x - edge)))
		{

			ioData[index_glob] = ioData[index_glob + 1];
		}
	}

	if ((blockIdx.x & 0x03) == 0x03) {   //fourth line, vertical, for larger indices
		register int index_loc = frameSize - lineFromEdge + threadIdx.x * blockDim.x;
		register int index_glob = blockIdx.y * frameSize + index_loc;

		if ((threadIdx.x > (edge - 1)) && (threadIdx.x < (blockDim.x - edge)))
		{
			ioData[index_glob] = ioData[index_glob - 1];
		}
	}

}

__global__ void addCuttedEdges(float* sinoPhase_out_dev, float* originalPhaseWrapped_dev, float* phaseUnwrapped_dev, int edgeSize)
{
	register int i_input_loc = (blockIdx.x - edgeSize) * (blockDim.x - 2 * edgeSize) + threadIdx.x - edgeSize;
	register int i_input_glob = blockIdx.y * (blockDim.x - 2 * edgeSize) * (gridDim.x - 2 * edgeSize) + i_input_loc;

	register int i_output_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int i_output_glob = blockIdx.y * blockDim.x * gridDim.x + i_output_loc;

	if ((threadIdx.x < edgeSize) || (threadIdx.x > (blockDim.x - edgeSize - 1)) ||
		(blockIdx.x < edgeSize) || (blockIdx.x > (gridDim.x - edgeSize - 1)))
	{
		sinoPhase_out_dev[i_output_glob] = originalPhaseWrapped_dev[i_output_glob]; // edges
	}
	else
	{
		sinoPhase_out_dev[i_output_glob] = phaseUnwrapped_dev[i_input_glob];
	}
}

__device__ void unwrappedPahseValue(float* wrapped, float unwrappedNeighbour)
{
	float diff = *(wrapped)-unwrappedNeighbour;
	int diff_2pi = round(diff / (PIx2));

	float wrapped_middle = *(wrapped)-diff_2pi * PIx2;

	float results[3];
	results[0] = wrapped_middle - PIx2;
	results[1] = wrapped_middle;
	results[2] = wrapped_middle + PIx2;

	float dif[3];
	dif[0] = fabsf(results[0] - unwrappedNeighbour);
	dif[1] = fabsf(results[1] - unwrappedNeighbour);
	dif[2] = fabsf(results[2] - unwrappedNeighbour);

	if (dif[2] < dif[1]) {
		dif[1] = dif[2];
		results[1] = results[2];
	}

	if (dif[1] < dif[0])
		*(wrapped) = results[1];
	else
		*(wrapped) = results[0];

}

__global__ void unwrapEdges(float* ioData, int lineFromEdge)
{
	register int frameSize = blockDim.x * blockDim.x;

	int edge = lineFromEdge;

	if ((blockIdx.x & 0x03) == 0x0) {   //first horizontal line, top of the frame
		register int index_loc = blockDim.x * (lineFromEdge - 1) + threadIdx.x;
		register int index_glob = blockIdx.y * frameSize + index_loc;

		if (threadIdx.x < (edge))
		{
			unwrappedPahseValue(&ioData[index_glob], ioData[index_glob + blockDim.x + edge]);
		}
		else if (threadIdx.x > (blockDim.x - edge - 1))
		{
			unwrappedPahseValue(&ioData[index_glob], ioData[index_glob + blockDim.x - edge - 1]);;
		}
		else
		{
			unwrappedPahseValue(&ioData[index_glob], ioData[index_glob + blockDim.x]);
		}
	}
	if ((blockIdx.x & 0x03) == 0x01) {   //second horizontal line, bottom of the frame
		register int index_loc = blockDim.x * (blockDim.x - lineFromEdge) + threadIdx.x;
		register int index_glob = blockIdx.y * frameSize + index_loc;

		if (threadIdx.x < (edge))
		{
			unwrappedPahseValue(&ioData[index_glob], ioData[index_glob - blockDim.x + edge]);;
		}
		else if (threadIdx.x > (blockDim.x - edge - 1))
		{
			unwrappedPahseValue(&ioData[index_glob], ioData[index_glob - blockDim.x - edge - 1]);;
		}
		else
		{
			unwrappedPahseValue(&ioData[index_glob], ioData[index_glob - blockDim.x]);
		}
	}
	if ((blockIdx.x & 0x03) == 0x02) {   //third line, vertical, for smaller indices
		register int index_loc = (lineFromEdge - 1) + threadIdx.x * blockDim.x;
		register int index_glob = blockIdx.y * frameSize + index_loc;

		if ((threadIdx.x > (edge - 1)) && (threadIdx.x < (blockDim.x - edge)))
		{

			unwrappedPahseValue(&ioData[index_glob], ioData[index_glob + 1]);
		}
	}

	if ((blockIdx.x & 0x03) == 0x03) {   //fourth line, vertical, for larger indices
		register int index_loc = frameSize - lineFromEdge + threadIdx.x * blockDim.x;
		register int index_glob = (blockIdx.y - 1) * frameSize + index_loc;

		if ((threadIdx.x > (edge)) && (threadIdx.x < (blockDim.x - edge + 1)))
		{
			unwrappedPahseValue(&ioData[index_glob], ioData[index_glob - 1]);
		}
	}

}

int fastPhaseUnwrappingWithVolkov(cufftComplex* complexFieldCSG_dev, float* sinoAmp_dev, float* sinoPhase_dev)
{
	cufftHandle plan1 = NULL, plan2 = NULL;
	int err;
	cudaError_t status;
	cufftResult_t e;

	cufftComplex* P_dev = NULL;

	float* K_dev;
	int e1, e2, e3;


	float* tempAmp_dev = NULL, * tempPh_dev = NULL;
	float* phaseExpanded_dev = NULL;

	e1 = cudaMalloc((void**)&tempAmp_dev, D * D * n_proj * sizeof(float));
	e2 = cudaMalloc((void**)&tempPh_dev, D * D * n_proj * sizeof(float));
	e3 = cudaMalloc((void**)&phaseExpanded_dev, 2 * D * 2 * D * n_proj * sizeof(float));
	if ((e1 != 0) || (e2 != 0) || (e3 != 0) || (phaseExpanded_dev == NULL) || (tempAmp_dev == NULL) || (tempPh_dev == NULL))
		return 84;

	dim3 block1(D, 1, 1);
	dim3 grid1(D, n_proj, 1);

	complexFieldToAmpPhase << < grid1, block1 >> > (complexFieldCSG_dev, sinoAmp_dev, tempPh_dev);

	cudaError_t cudaStatus13 = cudaGetLastError();

	dim3 block(Dminus, 1, 1);
	dim3 grid(Dminus, n_proj, 1);
	volkovExpand << < grid, block >> > (phaseExpanded_dev, tempPh_dev, D - Dminus);

	/*float* temp = (float*)malloc(D*D*n_proj * sizeof(float));
	cudaError_t cudaStatus1 = cudaMemcpy(temp, sinoAmp_dev, D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	FILE* plik8;
	plik8 = fopen("amplituda.bin", "wb");
	fwrite(temp, sizeof(float), D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);
	/*float* temp = (float*)malloc(4 * D*D*n_proj *sizeof(float));
	cudaError_t cudaStatus1 = cudaMemcpy(temp, phaseExpanded_dev, 4*D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	FILE* plik8;
	plik8 = fopen("VolkovExpanded_phase.bin", "wb");
	fwrite(temp, sizeof(float), 4*D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	/*float* temp2 = (float*)malloc( D*D*n_proj *sizeof(float));
	cudaError_t cudaStatus101 = cudaMemcpy(temp2, complexFieldCSG_dev, D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	FILE* plik18;
	plik18 = fopen("complexFieldCSG.bin", "wb");
	fwrite(temp2, sizeof(float), D*D*n_proj, plik18);
	fclose(plik18);
	free(temp2);*/
	int Dorig = D;

	D = Dminus * 2;

	float* sPhase_dev = NULL;

	e1 = cudaMalloc((void**)&K_dev, D * D * sizeof(float));
	e3 = cudaMalloc((void**)&P_dev, D * D * n_proj * sizeof(cufftComplex));
	e2 = cudaMalloc((void**)&sPhase_dev, D * D * n_proj * sizeof(float));

	if ((e1 != 0) || (e2 != 0) || (e3 != 0) || (sPhase_dev == NULL) || (P_dev == NULL) || (K_dev == NULL))
		return 82;

	/* Create a 3D FFT ComplexToComplex plan. */
	int NN_D[2] = { D,D };

	int idist_D = D * D;
	int odist_D = D * D;
	int inembed_D[] = { D, D };
	int onembed_D[] = { D, D };
	int istride_D = 1;
	int ostride_D = 1;

	e = cufftPlanMany(&plan1, 2, NN_D, inembed_D, istride_D, idist_D, onembed_D, ostride_D, odist_D, CUFFT_C2C, n_proj);
	if (e != cudaSuccess)
		return 81; // FFT plan creation error

	block = { (unsigned)D, (unsigned)1, (unsigned)1 };
	grid = { (unsigned)D, (unsigned)n_proj, (unsigned)1 };

	calc_matrix_K_2D(K_dev);
	ifftShift_1_float(K_dev);
	ifftShift_2_float(K_dev);
	/*temp = (float*)malloc(D*D*sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, K_dev, D*D * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("vecK.bin", "wb");
	fwrite(temp, sizeof(float), D*D, plik8);
	fclose(plik8);
	free(temp);*/


	P_vec_fromPhase << < grid, block >> > (P_dev, sPhase_dev, phaseExpanded_dev);
	cudaDeviceSynchronize();
	err = cudaGetLastError();

	/*temp = (float*)malloc(2*D*D*n_proj * sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, P_dev, 2*D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("phaseBefore1FFT.bin", "wb");
	fwrite(temp, sizeof(float), 2*D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	//1
	e = cufftExecC2C(plan1, P_dev, P_dev, CUFFT_FORWARD);
	cudaDeviceSynchronize();
	//ifftShift_1(D, D, n_proj, P_dev);
	//ifftShift_2(D, D, n_proj, P_dev);

												/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
												cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
												plik8 = fopen("phaseAfter1FFT.bin", "wb");
												fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
												fclose(plik8);
												free(temp);*/


	product_ComplexFloatD << < grid, block >> > (P_dev, K_dev);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("afterMult.bin", "wb");
	fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	//2

	//ifftShift_1(D, D, n_proj, P_dev);
	//ifftShift_2(D, D, n_proj, P_dev);
	e = cufftExecC2C(plan1, P_dev, P_dev, CUFFT_INVERSE);
	cudaDeviceSynchronize();
	//signCorrection << < grid, block >> > (P_dev);


												/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
												cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
												plik8 = fopen("phaseAfter2FFT.bin", "wb");
												fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
												fclose(plik8);
												free(temp);*/

	Imag_z1byz2 << < grid, block >> > (P_dev, sPhase_dev, D * D);
	cudaDeviceSynchronize();
	err = cudaGetLastError();

	/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("phaseBefore3FFT.bin", "wb");
	fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	//3
	e = cufftExecC2C(plan1, P_dev, P_dev, CUFFT_FORWARD);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	//ifftShift_1(D, D, n_proj, P_dev);
	//ifftShift_2(D, D, n_proj, P_dev);

												/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
												cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
												plik8 = fopen("phaseAfter3FFT.bin", "wb");
												fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
												fclose(plik8);
												free(temp);*/

	div_ComplexFloatD << < grid, block >> > (P_dev, K_dev);
	cudaDeviceSynchronize();
	err = cudaGetLastError();

	/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("phaseBefore4FFT.bin", "wb");
	fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	//4
	//ifftShift_1(D, D, n_proj, P_dev);
	//ifftShift_2(D, D, n_proj, P_dev);
	e = cufftExecC2C(plan1, P_dev, P_dev, CUFFT_INVERSE);
	cudaDeviceSynchronize();


	/*float* temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	cudaError_t cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	FILE* plik8 = fopen("psi.bin", "wb");
	fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	/*temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	cudaStatus1 = cudaMemcpy(temp, sPhase_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("sPhase_dev.bin", "wb");
	fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	D = D / 2;
	block = { (unsigned)D,(unsigned)1, (unsigned)1 };
	grid = { (unsigned)D, (unsigned)n_proj, (unsigned)1 };

	unwrapVolkovCut << < grid, block >> > (phaseExpanded_dev, sPhase_dev, P_dev, 2 * PI, 4 * D * D);
	/// !!!! temporary use of phaseExpanded_dev - here is unwrapped phase (Dminus x Dminus)
	cudaDeviceSynchronize();
	err = cudaGetLastError();

	D = Dorig;

	block = { (unsigned)D,(unsigned)1, (unsigned)1 };
	grid = { (unsigned)D, (unsigned)n_proj, (unsigned)1 };

	addCuttedEdges << < grid, block >> > (sinoPhase_dev, tempPh_dev, phaseExpanded_dev, (Dorig - Dminus) / 2);
	// result to sinoPhase_dev

	dim3 blockMed(D, 1, 1);
	dim3 gridMed(4, n_proj, 1);
	unwrapEdges << < gridMed, blockMed >> > (sinoPhase_dev, 4);
	unwrapEdges << < gridMed, blockMed >> > (sinoPhase_dev, 3);
	unwrapEdges << < gridMed, blockMed >> > (sinoPhase_dev, 2);
	unwrapEdges << < gridMed, blockMed >> > (sinoPhase_dev, 1);

	medianFiter << < gridMed, blockMed >> > (sinoPhase_dev, 3);
	medianFiter << < gridMed, blockMed >> > (sinoPhase_dev, 2);
	medianFiter << < gridMed, blockMed >> > (sinoPhase_dev, 1);

	cudaFree(phaseExpanded_dev); phaseExpanded_dev = NULL;
	cudaFree(tempAmp_dev); tempAmp_dev = NULL;
	cudaFree(tempPh_dev); tempPh_dev = NULL;
	cudaFree(sPhase_dev); sPhase_dev = NULL;
	cudaFree(P_dev); P_dev = NULL;
	cudaFree(K_dev); K_dev = NULL;
	cufftDestroy(plan1); plan1 = NULL;

	return 0;
}

int fastPhaseUnwrappingWithVolkov2(cufftComplex* complexFieldCSG_dev, float* sinoAmp_dev, float* sinoPhase_dev)
{
	cufftHandle plan1 = NULL, plan2 = NULL;
	int err;
	cudaError_t status;
	cufftResult_t e;

	cufftComplex* P_dev = NULL;

	float* K_dev;
	int e1, e2, e3;


	float* tempAmp_dev = NULL, * tempPh_dev = NULL;
	float* phaseExpanded_dev = NULL;

	e1 = cudaMalloc((void**)&tempAmp_dev, D * D * n_proj * sizeof(float));
	e2 = cudaMalloc((void**)&tempPh_dev, D * D * n_proj * sizeof(float));
	e3 = cudaMalloc((void**)&phaseExpanded_dev, 2 * D * 2 * D * n_proj * sizeof(float));
	if ((e1 != 0) || (e2 != 0) || (e3 != 0) || (phaseExpanded_dev == NULL) || (tempAmp_dev == NULL) || (tempPh_dev == NULL))
		return 84;


	err = complexFieldToComplexPhaseAndAmp(complexFieldCSG_dev, sinoPhComplex_dev, sinoAmp_dev);
	if (err != cudaSuccess)
		return err;

	if (referenceCalculated) {
		complexFieldReferenceSubtraction(sinoPhComplex_dev, sinoPhRefComplex_dev);
		if (err != cudaSuccess)
			return err;
	}

	err = complexPhaseToPhase(sinoPhComplex_dev, tempPh_dev);
	if (err != cudaSuccess)
		return err;

	//float* temp2 = (float*)malloc( 2*D*D*n_proj *sizeof(float));
	//cudaError_t cudaStatus101 = cudaMemcpy(temp2, sinoPhComplex_dev, 2*D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	//FILE* plik18;
	//plik18 = fopen("sinoPhComplex.bin", "wb");
	//fwrite(temp2, sizeof(float), 2*D*D*n_proj, plik18);
	//fclose(plik18);
	//free(temp2);

	//float* temp = (float*)malloc(D*D*n_proj * sizeof(float));
	//cudaError_t cudaStatus1 = cudaMemcpy(temp, sinoAmp_dev, D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	//FILE* plik8;
	//plik8 = fopen("amp.bin", "wb");
	//fwrite(temp, sizeof(float), D*D*n_proj, plik8);
	//fclose(plik8);
	//free(temp);

	//float* temp = (float*)malloc(D*D*n_proj * sizeof(float));
	//cudaError_t cudaStatus1 = cudaMemcpy(temp, tempPh_dev, D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	//FILE* plik8;
	//plik8 = fopen("fazaPrzedVolkov.bin", "wb");
	//fwrite(temp, sizeof(float), D*D*n_proj, plik8);
	//fclose(plik8);
	//free(temp);

	dim3 block(D, 1, 1);
	dim3 grid(D, n_proj, 1);
	volkovExpand << < grid, block >> > (phaseExpanded_dev, tempPh_dev, 0);

	/*float* temp = (float*)malloc(D*D*n_proj * sizeof(float));
	cudaError_t cudaStatus1 = cudaMemcpy(temp, sinoAmp_dev, D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	FILE* plik8;
	plik8 = fopen("amplituda.bin", "wb");
	fwrite(temp, sizeof(float), D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/
	/*float* temp = (float*)malloc(4 * D*D*n_proj *sizeof(float));
	cudaError_t cudaStatus1 = cudaMemcpy(temp, phaseExpanded_dev, 4*D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	FILE* plik8;
	plik8 = fopen("VolkovExpanded_phase.bin", "wb");
	fwrite(temp, sizeof(float), 4*D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	/*float* temp2 = (float*)malloc( D*D*n_proj *sizeof(float));
	cudaError_t cudaStatus101 = cudaMemcpy(temp2, complexFieldCSG_dev, D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	FILE* plik18;
	plik18 = fopen("complexFieldCSG.bin", "wb");
	fwrite(temp2, sizeof(float), D*D*n_proj, plik18);
	fclose(plik18);
	free(temp2);*/
	int Dorig = D;

	D = D * 2;

	cufftComplex* P2_dev = NULL;

	e1 = cudaMalloc((void**)&K_dev, D * D * sizeof(float));
	e3 = cudaMalloc((void**)&P_dev, D * D * n_proj * sizeof(cufftComplex));
	e2 = cudaMalloc((void**)&P2_dev, D * D * n_proj * sizeof(cufftComplex));

	if ((e1 != 0) || (e2 != 0) || (e3 != 0) || (P2_dev == NULL) || (P_dev == NULL) || (K_dev == NULL))
		return 82;

	/* Create a 3D FFT ComplexToComplex plan. */
	int NN_D[2] = { D,D };

	int idist_D = D * D;
	int odist_D = D * D;
	int inembed_D[] = { D, D };
	int onembed_D[] = { D, D };
	int istride_D = 1;
	int ostride_D = 1;

	e = cufftPlanMany(&plan1, 2, NN_D, inembed_D, istride_D, idist_D, onembed_D, ostride_D, odist_D, CUFFT_C2C, n_proj);
	if (e != cudaSuccess)
		return 81; // FFT plan creation error

	block = { (unsigned)D, (unsigned)1, (unsigned)1 };
	grid = { (unsigned)D, (unsigned)n_proj, (unsigned)1 };

	calc_matrix_K_2D(K_dev);
	//ifftShift_1_float(K_dev);
	//ifftShift_2_float(K_dev);
												//float * temp = (float*)malloc(D*D * sizeof(float));
												//cudaError_t cudaStatus1 = cudaMemcpy(temp, K_dev, D*D * sizeof(float), cudaMemcpyDeviceToHost);
												//FILE* plik8 = fopen("vecK.bin", "wb");
												//fwrite(temp, sizeof(float), D*D, plik8);
												//fclose(plik8);
												//free(temp);

												/*temp = (float*)malloc(D*D*n_proj * sizeof(float));
												cudaStatus1 = cudaMemcpy(temp, phaseExpanded_dev, D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
												plik8 = fopen("VolkovExpanded_phase2.bin", "wb");
												fwrite(temp, sizeof(float), D*D*n_proj, plik8);
												fclose(plik8);
												free(temp);*/

	P_vec_fromPhase2 << < grid, block >> > (P2_dev, phaseExpanded_dev);
	cudaDeviceSynchronize();
	err = cudaGetLastError();

	//temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	//cudaStatus1 = cudaMemcpy(temp, P2_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	//plik8 = fopen("phaseBefore1FFT.bin", "wb");
	//fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	//fclose(plik8);
	//free(temp);

	//1
	e = cufftExecC2C(plan1, P2_dev, P_dev, CUFFT_FORWARD);
	cudaDeviceSynchronize();
	ifftShift_1(D, D, n_proj, P_dev);
	ifftShift_2(D, D, n_proj, P_dev);

	//temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	//cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	//plik8 = fopen("phaseAfter1FFT.bin", "wb");
	//fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	//fclose(plik8);
	//free(temp);


	product_ComplexFloatD << < grid, block >> > (P_dev, K_dev);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	//temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	//cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	//plik8 = fopen("afterMult.bin", "wb");
	//fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	//fclose(plik8);
	//free(temp);

//2

	ifftShift_1(D, D, n_proj, P_dev);
	ifftShift_2(D, D, n_proj, P_dev);
	e = cufftExecC2C(plan1, P_dev, P_dev, CUFFT_INVERSE);
	cudaDeviceSynchronize();
	//signCorrection << < grid, block >> > (P_dev);


														//temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
														//cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
														//plik8 = fopen("phaseAfter2FFT.bin", "wb");
														//fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
														//fclose(plik8);
														//free(temp);

	Imag_z1byz2_v2 << < grid, block >> > (P_dev, P2_dev, D * D);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	//Imag_z1byz2 << < grid, block >> > (P_dev, sPhase_dev, D*D);


	//temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	//cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	//plik8 = fopen("phaseBefore3FFT.bin", "wb");
	//fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	//fclose(plik8);
	//free(temp);

	//3
	e = cufftExecC2C(plan1, P_dev, P_dev, CUFFT_FORWARD);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	ifftShift_1(D, D, n_proj, P_dev);
	ifftShift_2(D, D, n_proj, P_dev);

	//temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	//cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	//plik8 = fopen("phaseAfter3FFT.bin", "wb");
	//fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	//fclose(plik8);
	//free(temp);


	div_ComplexFloatD << < grid, block >> > (P_dev, K_dev);
	cudaDeviceSynchronize();
	err = cudaGetLastError();

	//temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	//cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	//plik8 = fopen("phaseBefore4FFT.bin", "wb");
	//fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	//fclose(plik8);
	//free(temp);

//4
	ifftShift_1(D, D, n_proj, P_dev);
	ifftShift_2(D, D, n_proj, P_dev);
	e = cufftExecC2C(plan1, P_dev, P_dev, CUFFT_INVERSE);
	cudaDeviceSynchronize();


	//temp = (float*)malloc(2 * D*D*n_proj * sizeof(float));
	//cudaStatus1 = cudaMemcpy(temp, P_dev, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	//plik8 = fopen("phaseAfter4FFT.bin", "wb");
	//fwrite(temp, sizeof(float), 2 * D*D*n_proj, plik8);
	//fclose(plik8);
	//free(temp);


	D = D / 2;
	block = { (unsigned)D,(unsigned)1, (unsigned)1 };
	grid = { (unsigned)D, (unsigned)n_proj, (unsigned)1 };

	unwrapVolkovCut_v2 << < grid, block >> > (sinoPhase_dev, P_dev, 2 * PI, 4 * D * D);
	cudaDeviceSynchronize();
	err = cudaGetLastError();

	//temp = (float*)malloc(D*D*n_proj * sizeof(float));
	//cudaStatus1 = cudaMemcpy(temp, sinoPhase_dev, D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	//plik8 = fopen("sinoPhaseUnw_dev.bin", "wb");
	//fwrite(temp, sizeof(float), D*D*n_proj, plik8);
	//fclose(plik8);
	//free(temp);


	cudaFree(phaseExpanded_dev); phaseExpanded_dev = NULL;
	cudaFree(tempAmp_dev); tempAmp_dev = NULL;
	cudaFree(tempPh_dev); tempPh_dev = NULL;
	cudaFree(P2_dev); P2_dev = NULL;
	cudaFree(P_dev); P_dev = NULL;
	cudaFree(K_dev); K_dev = NULL;
	cufftDestroy(plan1); plan1 = NULL;

	return 0;
}

#pragma endregion


void freeDeviceMemory()
{
	cudaFree(KO_dev); KO_dev = NULL;
	cudaFree(EW_dev); EW_dev = NULL;
	cudaFree(rf_dev); rf_dev = NULL;
	cudaFree(kn_dev); kn_dev = NULL;
	cudaFree(kxp_dev); kxp_dev = NULL;
	cudaFree(kyp_dev); kyp_dev = NULL;
	cudaFree(kxy_dev); kxy_dev = NULL;
	cudaFree(sinoAmp_dev); sinoAmp_dev = NULL;
	cudaFree(sinoPh_dev); sinoPh_dev = NULL;
	cudaFree(sinoPhComplex_dev); sinoPhComplex_dev = NULL;
	cudaFree(projections_dev); projections_dev = NULL;
	cudaFree(xShift_dev); xShift_dev = NULL;
	cudaFree(yShift_dev); yShift_dev = NULL;
	cufftDestroy(planFFTsino); planFFTsino = NULL;

	free(kxp_host); kxp_host = NULL;
	free(kyp_host); kyp_host = NULL;
}

void cleanDeviceMemoryBeforeDirectInversion()
{
	//cudaFree(KO_dev); KO_dev = NULL;
	//cudaFree(EW_dev); EW_dev = NULL;
	cudaFree(rf_dev); rf_dev = NULL;
	cudaFree(kn_dev); kn_dev = NULL;
	cudaFree(kxp_dev); kxp_dev = NULL;
	cudaFree(kyp_dev); kyp_dev = NULL;
	cudaFree(kxy_dev); kxy_dev = NULL;
	cudaFree(sinoAmp_dev); sinoAmp_dev = NULL;
	cudaFree(sinoPh_dev); sinoPh_dev = NULL;
	cudaFree(sinoPhComplex_dev); sinoPhComplex_dev = NULL;
	cudaFree(projections_dev); projections_dev = NULL;
	cudaFree(xShift_dev); xShift_dev = NULL;
	cudaFree(yShift_dev); yShift_dev = NULL;
	cufftDestroy(planFFTsino); planFFTsino = NULL;

	free(kxp_host); kxp_host = NULL;
	free(kyp_host); kyp_host = NULL;
}

int memoryAlloc()
{
	cudaMalloc((void**)&KO_dev, K_xy * K_xy * K_z * sizeof(cufftComplex));
	cudaMemset((void**)&KO_dev, 0, K_xy * K_xy * K_z * 2 * sizeof(float));
	cudaMalloc((void**)&EW_dev, K_xy * K_xy * K_z * sizeof(int));

	if (!preProcOnGPU) {
		cudaMalloc((void**)&sinoAmp_dev, Nx * Ny * n_proj * sizeof(float));
		cudaMalloc((void**)&sinoPh_dev, Nx * Ny * n_proj * sizeof(float));
		cudaMalloc((void**)&xShift_dev, n_proj * sizeof(int));
		cudaMalloc((void**)&yShift_dev, n_proj * sizeof(int));
		cudaMalloc((void**)&kxp_dev, n_proj * sizeof(float));
		cudaMalloc((void**)&kyp_dev, n_proj * sizeof(float));
	}
	cudaMalloc((void**)&projections_dev, Nx * Ny * n_proj * sizeof(cufftComplex));
	//cudaMalloc((void**)&projectionsInSizeK_dev, K_xy * K_xy * n_proj * sizeof(cufftComplex));
	//cudaMemset((void**)&projectionsInSizeK_dev, 0, K_xy * K_xy * n_proj * sizeof(cufftComplex));

	cudaMalloc((void**)&kn_dev, n_proj * sizeof(float));
	cudaMalloc((void**)&kz_dev, n_proj * sizeof(float));
	cudaMalloc((void**)&rf_dev, n_proj * sizeof(float));

	Nx_Fpmask = Nx / 32;
	if (Nx % 32 != 0)
		Nx_Fpmask++;
	if (!preProcOnGPU) {
		cudaMalloc((void**)&fpmask_dev, Nx_Fpmask * Ny * n_proj * sizeof(int));
		if (fpmask_dev == NULL)
			return 50; //memory allocation error
	}

	if ((KO_dev == NULL) || (EW_dev == NULL) || (sinoAmp_dev == NULL) || (sinoPh_dev == NULL) ||
		(projections_dev == NULL) || /*(projectionsInSizeK_dev == NULL) ||*/ (xShift_dev == NULL) ||
		(yShift_dev == NULL) || (kxp_dev == NULL) || (kyp_dev == NULL))
		return 50; //memory allocation error


	kxp_host = (float*)malloc(n_proj * sizeof(float));
	kyp_host = (float*)malloc(n_proj * sizeof(float));
	if ((kxp_host == NULL) || (kyp_host == NULL))
		return 50; //memory allocation error


	return 0;
}

// calculate kzz0
//__global__ void kernel_kzz0(cufftComplex* n_rec_dev, float nImm)
//{
//	register int index = blockIdx.x * blockDim.x + threadIdx.x;
//
//	if (n_rec_dev[index].x < nImm)
//		n_rec_dev[index].x = nImm;
//}
//
//int calculate_kzz0()
//{
//	int nThreads = n_rec_x;
//	dim3 nBlocks(n_rec_y*n_rec_z);
//	// execute the kernel
//	kernel_kzz0 << <nBlocks, nThreads >> > (n_rec_dev, n_imm);
//
//	cudaError_t cudaStatus1 = cudaGetLastError();
//
//	cudaError_t cudaStatus2 = cudaDeviceSynchronize();
//
//	if ((cudaStatus1 != cudaSuccess) && (cudaStatus2 != cudaSuccess)) {
//		return 8;
//	}
//	return 0;
//}

// approx Rytov 
__global__ void kernel_approxRytovKgen(float* sinoAmp_dev, float* sinoPh_dev, cufftComplex* projections_dev, int max)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;
	register float temp;

	if (index < max) {
		temp = sinoAmp_dev[index];
		temp = logf(temp);
		if (!isnan(temp))
			projections_dev[index].x = temp;
		else
			projections_dev[index].x = 0.0f; // logf(1.0f);

		if (!isnan(sinoPh_dev[index]))
			projections_dev[index].y = sinoPh_dev[index];
		else
			projections_dev[index].y = 0.0f;
	}

}

int approxRytovK_Gen()
{
	int nThreads = 512;

	int n = Nx * Ny * n_proj;
	int blocks = n / nThreads;
	int r = n % nThreads;
	if (r > 0) blocks++;
	dim3 nBlocks(blocks);
	// execute the kernel
	kernel_approxRytovKgen << <nBlocks, nThreads >> > (sinoAmp_dev, sinoPh_dev, projections_dev, n);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) && (cudaStatus2 != cudaSuccess)) {
		return 54;  //error in the function computing the Rytov approximation
	}
	return 0;
}

__device__
cuFloatComplex my_complex_exp(cuFloatComplex arg)
{
	cuFloatComplex res;
	float s, c;
	float e = expf(arg.x);
	sincosf(arg.y, &s, &c);
	res.x = c * e;
	res.y = s * e;
	return res;
}


__global__ void kernel_approxBornKgen(float* sinoAmp_dev, float* sinoPh_dev, cufftComplex* projections_dev, int max)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < max) {
		//thrust::complex<float> temp(0.0f, sinoPh_dev[index]);
		//register float re = thrust::exp(temp).real;
		//register float im = thrust::exp(temp).imag;
		register float re;
		register float im;

		cufftComplex temp;
		//float s, c;
		//float e = 1.0f; // expf(arg.x);
		__sincosf(sinoPh_dev[index], &im, &re); //sincosf(arg.y, &s, &c);
		//re = c * e;
		//im = s * e;


		register float amplitude = sinoAmp_dev[index];
		if (!isnan(amplitude)) amplitude = 1.0f;

		if ((!isnan(re)) && (!isnan(im))) {
			projections_dev[index].x = amplitude * re - 1.0f;
			projections_dev[index].y = amplitude * im;
		}
		else {
			projections_dev[index].x = amplitude /* *1 */ - 1.0f;
			projections_dev[index].y = /* amplitude * */ 0.0f;
		}
	}

}

int approxBornK_Gen()
{
	int nThreads = 512;

	int n = Nx * Ny * n_proj;
	int blocks = n / nThreads;
	int r = n % nThreads;
	if (r > 0) blocks++;
	dim3 nBlocks(blocks);
	// execute the kernel
	kernel_approxBornKgen << <nBlocks, nThreads >> > (sinoAmp_dev, sinoPh_dev, projections_dev, n);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) && (cudaStatus2 != cudaSuccess)) {
		return 54;  //error in the function computing the Rytov approximation
	}
	return 0;
}


__global__ void kernel_multComplexVector(cufftComplex* data_dev, float factor, int max)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < max) {
		data_dev[index].x *= factor;
		data_dev[index].y *= factor;
	}
}

int multComplexVector(int elements, cufftComplex* data_dev, float factor)
{
	int nThreads = 512;

	int blocks = elements / nThreads;
	int r = elements % nThreads;
	if (r > 0) blocks++;
	dim3 nBlocks(blocks);
	// execute the kernel
	kernel_multComplexVector << <nBlocks, nThreads >> > (data_dev, factor, elements);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) && (cudaStatus2 != cudaSuccess)) {
		return 55;  //error in the multComplexVector function
	}
	return 0;
}

//__global__ void kernel_approxRytovKgen2(float* sinoAmp_dev, float* sinoPh_dev, cufftComplex* projections_dev)
//{
//	register int index = blockIdx.x * blockDim.x + threadIdx.x;
//
//	projections_dev[index].x = logf(sinoAmp_dev[index]);
//	projections_dev[index].y = sinoPh_dev[index];
//}
//
//int approxRytovK_Gen2()
//{
//	int nThreads = Nx;
//	dim3 nBlocks(Ny * n_proj);
//	// execute the kernel
//	kernel_approxRytovKgen2 << <nBlocks, nThreads >> > (sinoAmp_dev, sinoPh_dev, projections_dev);
//
//	cudaError_t cudaStatus1 = cudaGetLastError();
//
//	cudaError_t cudaStatus2 = cudaDeviceSynchronize();
//
//	if ((cudaStatus1 != cudaSuccess) && (cudaStatus2 != cudaSuccess)) {
//		return 54;  //error in the function computing the Rytov approximation
//	}
//	return 0;
//}


__global__ void kernel_transferProjectionsToKbuffer(cufftComplex* projections_dev, cufftComplex* _KO_dev, int* _EW_dev,
	int* shiftX, int* shiftY, int Kxy, int Kz, float* rf_dev, unsigned int* fpmask_dev, int Nx_Fpmask,
	float* _kxp_dev, float* _kyp_dev, float* _kn_dev, float* _kz_dev, float _dkP, float _dkPz, int _preProcOnGPU)
{
	register int x_in = threadIdx.x;		// x within the projection
	register int y_in = blockIdx.x;			// y within the projection
	register int proj = blockIdx.y;			// index of the projection currently being processed
	register int in_proj_size = blockDim.x * gridDim.x;		//projection size
	register int out_proj_size = Kxy * Kxy;					// x*y size in the "horizontal" plane of K-space

	register int x_in_shifted = x_in - shiftY[proj];			// equivalent of shifting the input data in MATLAB
	if (x_in_shifted < 0) x_in_shifted += blockDim.x;			// implemented by computing the appropriate index
	if (x_in_shifted >= blockDim.x) x_in_shifted -= blockDim.x;	// if the index is below zero we take data "from the other side"
	register int y_in_shifted = y_in - shiftX[proj];			// and if it exceeds the frame size we wrap around by its size
	if (y_in_shifted < 0) y_in_shifted += gridDim.x;			// these operations correspond to the comcirc command applied in two dimensions
	if (y_in_shifted >= gridDim.x) y_in_shifted -= gridDim.x;

	register int in_index = in_proj_size * proj + blockDim.x * y_in_shifted + x_in_shifted;
	// computing the input image index; threads already read the image after the projection shifts, i.e. we have a circle

	register int deltaX = (Kxy - blockDim.x) / 2 - shiftY[proj];	// shifts of the Ewald spheres in the x and y plane
	register int deltaY = (Kxy - gridDim.x) / 2 - shiftX[proj];		// in K-space

	register int out_index;// = out_proj_size * proj + Kxy * (y_in + deltaY) + x_in + deltaX;
	// this index allows writing projection data into a buffer of size Kxy by Kxy, for now 2D
	// y_in + deltaY - the equivalent of the y index in the K-space plane
	// x_in + deltaX - the equivalent of the x index in the K-space plane

	register int x_r = threadIdx.x - blockDim.x / 2; // x and y coordinates of the Ewald spheres; they allow selecting only
	register int y_r = blockIdx.x - gridDim.x / 2;	// those threads holding data (inside the circle rf) that are actually to be placed into K

	register unsigned int mask = 0;
	if (!_preProcOnGPU)
		mask |= fpmask_dev[gridDim.x * Nx_Fpmask * proj + Nx_Fpmask * y_in_shifted + (x_in_shifted >> 5)] & (0x01 << (x_in_shifted & 0x1F));
	// index : y_size * x_size * projection number + x_size * y + "int number"
	else
		mask = 1;

	register float a, b;
	register int ai, bi;

	if (((x_r * x_r) + (y_r * y_r)) < (rf_dev[blockIdx.y] * rf_dev[blockIdx.y]) // select threads inside the circle of the given radius
		&& ((y_in + deltaY) < Kxy) && ((x_in + deltaX) < Kxy) && (mask > 0)) {
		ai = x_in + deltaX - (Kxy >> 1);
		a = (float)ai * _dkP + _kyp_dev[proj];
		a *= a;
		bi = y_in + deltaY - (Kxy >> 1);
		b = (float)bi * _dkP + _kxp_dev[proj];
		b *= b;

		a = sqrtf((_kn_dev[proj]) * (_kn_dev[proj]) - a - b); // teraz a = kzz
		if (a > 0) { // discard evanescent waves
			b = 4 * PI * a; // here b is only a temporary store for the next two lines

			out_index = roundf((a - _kz_dev[proj]) / _dkPz) + ceilf((Kz + 1) / 2);
			out_index = out_index * out_proj_size + Kxy * (y_in + deltaY) + x_in + deltaX;

			atomicAdd(&_KO_dev[out_index].x, -b * projections_dev[in_index].y);
			atomicAdd(&_KO_dev[out_index].y, b * projections_dev[in_index].x);
			atomicAdd(&_EW_dev[out_index], 1);
		}
	}
}

int transferProjectionsToKbuffer()
{
	int nThreads = Nx;
	dim3 nBlocks(Ny, n_proj, 1);
	// execute the kernel

						/*		float* temp = (float*)malloc(Nx*Ny * 2 * n_proj * sizeof(float));
								memset(temp, 0, Nx*Ny *n_proj * 2 * sizeof(float));
								cudaError_t cudaStatus3 = cudaMemcpy(temp, projections_dev, Nx*Ny*n_proj * 2 * sizeof(float), cudaMemcpyDeviceToHost);
								FILE* plik8;
								if (preProcOnGPU)
									plik8 = fopen("projectionsPreprocOnGPU.bin", "wb");
								else
									plik8 = fopen("projections.bin", "wb");
								fwrite(temp, sizeof(float), 2 * Nx*Ny*n_proj, plik8);
								fclose(plik8);
								free(temp);*/

								//temp = (float*)malloc(Nx*Ny * 2 * n_proj * sizeof(float));
								//memset(temp, 0, Nx*Ny *n_proj * 2 * sizeof(float));
								//cudaStatus3 = cudaMemcpy(temp, KO_dev, Nx*Ny*n_proj * 2 * sizeof(float), cudaMemcpyDeviceToHost);
								//if (preProcOnGPU)
								//	plik8 = fopen("KO_przed_PreprocOnGPU.bin", "wb");
								//else
								//	plik8 = fopen("KO_przed.bin", "wb");
								//fwrite(temp, sizeof(float), 2 * Nx*Ny*n_proj, plik8);
								//fclose(plik8);
								//free(temp);


	kernel_transferProjectionsToKbuffer << <nBlocks, nThreads >> > (projections_dev, KO_dev, EW_dev,
		xShift_dev, yShift_dev, K_xy, K_z, rf_dev, fpmask_dev, Nx_Fpmask, kxp_dev, kyp_dev, kn_dev, kz_dev, dkP, dkPz, preProcOnGPU);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) && (cudaStatus2 != cudaSuccess)) {
		return 57;  //error in the transferProjectionsToPreKbuffer function
	}

	/*int *temp = (int*)malloc( n_proj* sizeof(int));
	cudaError_t cudaStatus3 = cudaMemcpy(temp, xShift_dev, n_proj * sizeof(int), cudaMemcpyDeviceToHost);
	FILE *plik8;
	if (preProcOnGPU)
		plik8 = fopen("xShift_PreprocOnGPU.bin", "wb");
	else
		plik8 = fopen("xShift.bin", "wb");
	fwrite(temp, sizeof(int), n_proj, plik8);
	fclose(plik8);
	free(temp);

	temp = (int*)malloc(n_proj * sizeof(int));
	cudaStatus3 = cudaMemcpy(temp, yShift_dev, n_proj * sizeof(int), cudaMemcpyDeviceToHost);

	if (preProcOnGPU)
		plik8 = fopen("yShift_PreprocOnGPU.bin", "wb");
	else
		plik8 = fopen("yShift.bin", "wb");
	fwrite(temp, sizeof(int), n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	/*float* temp = (float*)malloc(n_proj * sizeof(float));
	memset(temp, 0, n_proj * sizeof(float));
	cudaError_t cudaStatus3 = cudaMemcpy(temp, kxp_dev, n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	FILE *plik8;
	if (preProcOnGPU)
		plik8 = fopen("kxp_po_PreprocOnGPU.bin", "wb");
	else
		plik8 = fopen("kxp_po.bin", "wb");
	fwrite(temp, sizeof(float), n_proj, plik8);
	fclose(plik8);
	free(temp);

	temp = (float*)malloc(n_proj * sizeof(float));
	memset(temp, 0, n_proj * sizeof(float));
	cudaStatus3 = cudaMemcpy(temp, kyp_dev, n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	if (preProcOnGPU)
		plik8 = fopen("kyp_po_PreprocOnGPU.bin", "wb");
	else
		plik8 = fopen("kyxp_po.bin", "wb");
	fwrite(temp, sizeof(float), n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	//temp = (float*)malloc(K_xy*K_xy * 2 * K_xy * sizeof(float));
	//memset(temp, 0, K_xy*K_xy * 2 * K_xy * sizeof(float));
	//cudaStatus3 = cudaMemcpy(temp, KO_dev, K_xy*K_xy * 2 * K_xy * sizeof(float), cudaMemcpyDeviceToHost);
	//if (preProcOnGPU)
	//	plik8 = fopen("KO_po_PreprocOnGPU.bin", "wb");
	//else
	//	plik8 = fopen("KO_po.bin", "wb");
	//fwrite(temp, sizeof(float), K_xy*K_xy * 2 * K_xy, plik8);
	//fclose(plik8);
	//free(temp);

	//temp = (float*)malloc(K_xy*K_xy * 2 * K_xy * sizeof(float));
	//plik8 = fopen("KO_po_PreprocOnGPU.bin", "rb");
	//fread(temp, sizeof(float), K_xy*K_xy * 2 * K_xy, plik8);
	//fclose(plik8);
	//cudaStatus3 = cudaMemcpy(KO_dev, temp, K_xy*K_xy * 2 * K_xy * sizeof(float), cudaMemcpyHostToDevice);
	//free(temp);


	return 0;
}

__global__ void kernel_divide_KO_by_EW(cufftComplex* data_dev, int* EW_dev, int max)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < max) {
		if (EW_dev[index] > 1) {
			data_dev[index].x /= EW_dev[index];
			data_dev[index].y /= EW_dev[index];
		}
	}
}

int divide_KO_by_EW(int elements, cufftComplex* data_dev, int* EW_dev)
{
	int nThreads = 512;

	int blocks = elements / nThreads;
	int r = elements % nThreads;
	if (r > 0) blocks++;
	dim3 nBlocks(blocks);
	// execute the kernel
	kernel_divide_KO_by_EW << <nBlocks, nThreads >> > (data_dev, EW_dev, elements);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) && (cudaStatus2 != cudaSuccess)) {
		return 60;  //error in the divide_KO_by_EW function
	}
	return 0;
}

__global__ void kernel_transferHologramsToCSGbuffer(cufftComplex* holograms_dev, short int* holograms_in_dev, float* mask)
{
	register int index = (blockIdx.y * gridDim.x + blockIdx.x) * 2 * blockDim.x + 2 * threadIdx.x;
	register int indexMask = blockIdx.x * 2 * blockDim.x + 2 * threadIdx.x;

	holograms_dev[index].x = (float)holograms_in_dev[index] * mask[indexMask];
	holograms_dev[index].y = 0;
	holograms_dev[index + 1].x = (float)holograms_in_dev[index + 1] * mask[indexMask + 1];
	holograms_dev[index + 1].y = 0;

}

int transferHologramsToCSGbufferWithMaskW()
{
	//int nThreads = 512;

	//int n = Nx * Ny * n_proj;
	//int blocks = n / nThreads;
	//int r = n % nThreads;
	//if (r > 0) blocks++;
	//dim3 nBlocks(blocks);
	int nThreads = Nx_holo / 2;
	dim3 nBlocks(Ny_holo, n_proj, 1);
	// execute the kernel
	kernel_transferHologramsToCSGbuffer << <nBlocks, nThreads >> > (hologramsCSG_device, holograms_device, tukeyWin2D_W_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) && (cudaStatus2 != cudaSuccess)) {
		return 65;
	}
	return 0;
}

/*__global__ void kernel_findMaxInAllHolograms(cufftComplex* holograms_dev, int omitedLines, int Nx, int Ny)
{
	register int offsetToFirstPixel = omitedLines * Nx + Nx / 2;

	register int index = (blockIdx.y * Nx * Ny + blockIdx.x) * Nx + threadIdx.x + offsetToFirstPixel;

	cufftComplex z1 = holograms_dev[index];

	__shared__ float buf[2 * 1024];

	register int i_loc = threadIdx.x;

	buf[2 * i_loc] = z1.x * z1.x + z1.y*z1.y; //no square root because we are looking for the maximum
	buf[2 * i_loc + 1] = i_loc;

	for (int krok = blockDim.x >> 1; krok > 0; krok >>= 1)
	{
		__syncthreads();
		if (i_loc < krok)
			if (buf[2 * i_loc] < buf[2 * (i_loc + krok)]) {
				buf[2 * i_loc] = buf[2 * (i_loc + krok)];
				buf[2 * i_loc + 1] = buf[2 * (i_loc + krok) + 1];
			}
	}

	if (i_loc == 0) {
		holograms_dev[(blockIdx.y * Nx * Ny + blockIdx.x)].x = buf[0]; // value
		holograms_dev[(blockIdx.y * Nx * Ny + blockIdx.x) + 1].y = buf[1];  // horizontal index
	}
}*/


__global__ void kernel_findMaxIn0incidence(cufftComplex* holograms_dev, int omitedLines, int Nx, int Ny)
{
	register int offsetToFirstPixel = omitedLines * Nx + Nx / 2;

	register int index = (blockIdx.x) * Nx + threadIdx.x + offsetToFirstPixel;

	__shared__ float buf[2 * 1024];
	cufftComplex z1;
	register int i_loc = threadIdx.x;

	if (threadIdx.x < Nx / 2) {
		z1 = holograms_dev[index];
		//holograms_dev[blockIdx.x * Nx + threadIdx.x] = z1;
		buf[2 * i_loc] = z1.x * z1.x + z1.y * z1.y; //no square root because we are looking for the maximum
	}
	else {
		buf[2 * i_loc] = 0;
	}

	buf[2 * i_loc + 1] = i_loc;


	for (int krok = blockDim.x >> 1; krok > 0; krok >>= 1)
	{
		__syncthreads();
		if (i_loc < krok)
			if (buf[2 * i_loc] < buf[2 * (i_loc + krok)]) {
				buf[2 * i_loc] = buf[2 * (i_loc + krok)];
				buf[2 * i_loc + 1] = buf[2 * (i_loc + krok) + 1];
			}
	}

	if (i_loc == 0) {
		holograms_dev[(blockIdx.x)].x = buf[0]; // value
		holograms_dev[(blockIdx.x) + 1].y = buf[1];  // horizontal index
	}
}


__global__ void kernel_findMaxIn0incidence_fin(cufftComplex* holograms_dev, int valuesToSearch, int cuttedX, int cuttedY)
{
	register int i_loc = threadIdx.x;

	float wartosc;
	//float indeksPoziomy = framesIn[2 * index_out + 1];

	__shared__ float buf[2 * 1024];

	if (i_loc < valuesToSearch) {
		wartosc = holograms_dev[i_loc].x;
	}
	else {
		wartosc = 0;
	}

	buf[2 * i_loc] = wartosc;   // values
	buf[2 * i_loc + 1] = i_loc; //vertical indices

	for (int krok = blockDim.x >> 1; krok > 0; krok >>= 1)
	{
		__syncthreads();
		if (i_loc < krok)
			if (buf[2 * i_loc] < buf[2 * (i_loc + krok)]) {
				buf[2 * i_loc] = buf[2 * (i_loc + krok)];
				buf[2 * i_loc + 1] = buf[2 * (i_loc + krok) + 1];
			}
	}
	__syncthreads();
	int y = (int)buf[1];
	if (i_loc == y) {
		buf[2] = holograms_dev[i_loc].y;  //horizontal index
	}

	__syncthreads();
	if (i_loc == 0) {
		holograms_dev[0].x = buf[2] + cuttedX;// horizontal index
		holograms_dev[0].y = buf[1] + cuttedY; //vertical index
		holograms_dev[1].x = buf[0]; //value
		//framesOut[2 * blockIdx.x + 1].x = buf[2]; //horizontal index
		//framesOut[2 * blockIdx.x+1].y = 0; //nic ;)
	}
}

int findMaxIn0incidence()
{

	int nThreads = 1024;
	int percent = 2; // percent to cut from half of the image
	int omitedLines = Ny_holo / 2 + (int)round(Ny_holo * (float)percent / 200); //additional division by 2 because the percentage refers to one half
	dim3 nBlocks(Ny_holo - omitedLines, 1, 1);
	// execute the kernel
	kernel_findMaxIn0incidence << <nBlocks, nThreads >> > (hologramsCSG_device, omitedLines, Nx_holo, Ny_holo);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) && (cudaStatus2 != cudaSuccess)) {
		return 68;
	}

	dim3 nThreads2(1024, 1, 1);  /// max hologram size 2048, less than half of the hologram - less than 1024
	dim3 nBlocks2(1, 1, 1);   //only one block - because we compute only from the first frame
	kernel_findMaxIn0incidence_fin << <nBlocks2, nThreads2 >> > (hologramsCSG_device, Ny_holo - omitedLines, Nx_holo / 2, omitedLines);

	cudaStatus1 = cudaDeviceSynchronize();

	if (cudaStatus1 != cudaSuccess) {
		return 69; // Error in the kernel function finalizing computations after FFT!
	}

	return 0;
}

__global__ void kernel_findMaxInAllFrames(cufftComplex* croppedPartsHoloCSG_device, int D, float* peaks1)
{
	register int index = blockIdx.y * D * D + blockIdx.x * D + threadIdx.x;

	__shared__ float buf[2 * 1024];
	cufftComplex z1;
	register int i_loc = threadIdx.x;

	z1 = croppedPartsHoloCSG_device[index];
	float wartosc = 0;

	if (threadIdx.x < D) {
		z1 = croppedPartsHoloCSG_device[index];
		wartosc = z1.x * z1.x + z1.y * z1.y; //no square root because we are looking for the maximum
	}
	else
		wartosc = 0;

	buf[2 * i_loc] = wartosc;   // values
	buf[2 * i_loc + 1] = i_loc; //horizontal indices

	//printf("index %3d, value %8f   \n", i_loc, wartosc);

	for (int krok = blockDim.x >> 1; krok > 0; krok >>= 1)
	{
		__syncthreads();
		if (i_loc < krok)
			if (buf[2 * i_loc] < buf[2 * (i_loc + krok)]) {
				buf[2 * i_loc] = buf[2 * (i_loc + krok)];
				buf[2 * i_loc + 1] = buf[2 * (i_loc + krok) + 1];
			}
	}

	if (i_loc == 0) {
		peaks1[2 * (blockIdx.y * D + blockIdx.x)] = buf[0]; // value
		//if (buf[1] > D) buf[1] = 0;
		peaks1[2 * (blockIdx.y * D + blockIdx.x) + 1] = buf[1];  // horizontal index
	}
}


__global__ void kernel_findMaxInAllFrames_fin(float* peaksStep1_dev, int* peaksStep2, int D)
{
	register int i_glob = blockIdx.x * D + threadIdx.x;

	register int i_loc = threadIdx.x;

	float wartosc;
	//float indeksPoziomy = framesIn[2 * index_out + 1];

	__shared__ float buf[2 * 1024];

	if (i_loc < D) {
		wartosc = peaksStep1_dev[2 * i_glob];
	}
	else {
		wartosc = 0;
	}

	buf[2 * i_loc] = wartosc;   // values
	buf[2 * i_loc + 1] = i_loc; //vertical indices

	for (int krok = blockDim.x >> 1; krok > 0; krok >>= 1)
	{
		__syncthreads();
		if (i_loc < krok)
			if (buf[2 * i_loc] < buf[2 * (i_loc + krok)]) {
				buf[2 * i_loc] = buf[2 * (i_loc + krok)];
				buf[2 * i_loc + 1] = buf[2 * (i_loc + krok) + 1];
			}
	}
	__syncthreads();
	int y = (int)buf[1];
	if (i_loc == y) {
		int x = peaksStep1_dev[2 * i_glob + 1];  //horizontal index
		peaksStep2[2 * blockIdx.x] = x;// horizontal index X
		peaksStep2[2 * blockIdx.x + 1] = y; //vertical index Y
	}
}

int findMaxInAllFrames()
{
	//
	// TO DO!!!!!!!!! Take care of bigger D than 256
	//
	int nThreads = 1024;
	dim3 nBlocks(D, n_proj, 1);
	//dim3 nBlocks(D, 1, 1); // temp
	// execute the kernel
	kernel_findMaxInAllFrames << <nBlocks, nThreads >> > (croppedPartsHoloCSG_device, D, peaksStep1_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 68;
	}

	dim3 nThreads2(nThreads, 1, 1);
	dim3 nBlocks2(n_proj, 1, 1);
	kernel_findMaxInAllFrames_fin << <nBlocks2, nThreads2 >> > (peaksStep1_dev, peaksStep2_dev, D);

	cudaStatus1 = cudaDeviceSynchronize();

	if (cudaStatus1 != cudaSuccess) {
		return 69;
	}

	return 0;
}
__global__ void kernel_circshift(cufftComplex* croppedPartsHoloCSG_device, int* peaksStep2, cufftComplex* complexFieldCSG_device)
{
	int frame = blockIdx.y;

	register int i_glob = (frame * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;

	int shiftX = blockDim.x / 2 - peaksStep2[2 * frame];
	int shiftY = gridDim.x / 2 - peaksStep2[2 * frame + 1];

	int newX = threadIdx.x + shiftX;
	if (newX >= blockDim.x) { newX -= blockDim.x; }
	if (newX < 0) { newX += blockDim.x; }

	int newY = blockIdx.x + shiftY;
	if (newY < 0) { newY += gridDim.x; }
	if (newY >= gridDim.x) { newY -= gridDim.x; }

	int address = frame * blockDim.x * gridDim.x + newY * blockDim.x + newX;

	if ((address >= 0) && (address < (blockDim.x * gridDim.x * gridDim.y))) {
		complexFieldCSG_device[frame * blockDim.x * gridDim.x + newY * blockDim.x + newX].x = shiftX;// newX; //croppedPartsHoloCSG_device[i_glob];
		complexFieldCSG_device[frame * blockDim.x * gridDim.x + newY * blockDim.x + newX].y = shiftY; // newY;
		//complexFieldCSG_device[frame*blockDim.x*gridDim.x + newY * blockDim.x + newX] = croppedPartsHoloCSG_device[i_glob];
	}
}


int circshift()
{

	int nThreads = D;
	dim3 nBlocks(D, n_proj, 1);
	//dim3 nBlocks(D, 1, 1); // temp
	// execute the kernel
	kernel_circshift << <nBlocks, nThreads >> > (croppedPartsHoloCSG_device, peaksStep2_dev, complexFieldCSG_device);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 75;
	}

	/*	dim3 nThreads2(nThreads, 1, 1);
		dim3 nBlocks2(n_proj, 1, 1);
		kernel_circshiftY << <nBlocks2, nThreads2 >> > (peaksStep1_dev, peaksStep2_dev, D);

		cudaStatus1 = cudaDeviceSynchronize();

		if (cudaStatus1 != cudaSuccess) {
			return 76;
		}*/

	return 0;
}

__global__ void kernel_circshift_v2(cufftComplex* croppedPartsHoloCSG_device, int* peaksStep2, cufftComplex* complexFieldCSG_device)
{
	int frame = blockIdx.y;

	register int i_glob = (frame * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	register int fr_offset = gridDim.x * blockDim.x * frame;

	int shiftX = blockDim.x / 2 - peaksStep2[2 * frame];
	int shiftY = gridDim.x / 2 - peaksStep2[2 * frame + 1];

	int srcX = threadIdx.x - shiftX;
	if (srcX >= blockDim.x) { srcX = srcX - blockDim.x; }
	if (srcX < 0) { srcX = srcX + blockDim.x; }
	if (srcX >= blockDim.x) { srcX = srcX - blockDim.x; }
	if (srcX < 0) { srcX = srcX + blockDim.x; }

	int srcY = blockIdx.x - shiftY;
	if (srcY < 0) { srcY = srcY + gridDim.x; }
	if (srcY >= gridDim.x) { srcY = srcY - gridDim.x; }
	if (srcY < 0) { srcY = srcY + gridDim.x; }
	if (srcY >= gridDim.x) { srcY = srcY - gridDim.x; }

	int srcIndex = fr_offset + blockDim.x * srcY + srcX;
	if (srcIndex < (gridDim.y * gridDim.x * blockDim.x)) {
		complexFieldCSG_device[i_glob] = croppedPartsHoloCSG_device[fr_offset + blockDim.x * srcY + srcX];
	}
	//else
	//{
	//	complexFieldCSG_device[i_glob].x = srcX;
	//	complexFieldCSG_device[i_glob].y = srcY;
	//}

}


int circshift_v2()
{

	int nThreads = D;
	dim3 nBlocks(D, n_proj, 1);
	//dim3 nBlocks(D, 1, 1); // temp
	// execute the kernel
	kernel_circshift_v2 << <nBlocks, nThreads >> > (croppedPartsHoloCSG_device, peaksStep2_dev, complexFieldCSG_device);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 75;
	}

	return 0;
}

__global__ void kernel_transferToCroppedHoloCSGbufferWith_sMask(cufftComplex* croppedHolo_dev, cufftComplex* holograms_in_dev, float* mask, int D, int Nx, int Ny)
{
	register int index_out = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	register int indexMask = blockIdx.x * blockDim.x + threadIdx.x;

	register int peakX = holograms_in_dev[0].x;
	register int peakY = holograms_in_dev[0].y;

	register int index_in = Nx * Ny * blockIdx.y + Nx * blockIdx.x + (peakY - D / 2) * Nx + (peakX - D / 2) + threadIdx.x;

	croppedHolo_dev[index_out].x = holograms_in_dev[index_in].x * mask[indexMask];
	croppedHolo_dev[index_out].y = holograms_in_dev[index_in].y * mask[indexMask];

}

int transferToCroppedHoloCSGbufferWith_sMask()
{
	//int nThreads = 512;

	//int n = Nx * Ny * n_proj;
	//int blocks = n / nThreads;
	//int r = n % nThreads;
	//if (r > 0) blocks++;
	//dim3 nBlocks(blocks);
	int nThreads = D;
	dim3 nBlocks(D, n_proj, 1);
	// execute the kernel
	kernel_transferToCroppedHoloCSGbufferWith_sMask << <nBlocks, nThreads >> > (croppedPartsHoloCSG_device, hologramsCSG_device, sMask_dev, D, Nx_holo, Ny_holo);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 72;
	}
	return 0;
}

__global__ void kernel_div2(cufftComplex* matrix3D, float factor)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	matrix3D[index].x /= factor;
	matrix3D[index].y /= factor;
}

int normalize_3Dmatrix(cufftComplex* matrix3D, float factor)
{
	int nThreads = D;
	dim3 nBlocks(D * n_proj);

	kernel_div2 << <nBlocks, nThreads >> > (matrix3D, factor);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 79;
	}
	return 0;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
/// EXPORTED FUNCTIONS
///////////////////////////////////////////////////////////////////////////

EXPORTED_FUNCTION int prepareKspaceGeneration(int _K_xy, int _K_z, float _dx, float _n_imm, int _n_proj, int _Nx, int _Ny,
	float _dkP, float _dkPz, float _NA, float* _lambda_all, float* _kxp, float* _kyp, int _approxBornNotRytov)
{
	cudaError_t cudaStatus1, cudaStatus2, cudaStatus3, cudaStatus4, cudaStatus5, cudaStatus6;
	dx = _dx;
	n_proj = _n_proj;
	dkP = _dkP;
	dkPz = _dkPz;
	K_xy = _K_xy;
	K_z = _K_z;
	n_imm = _n_imm;
	Nx = _Nx;
	Ny = _Ny;
	approxBornNotRytov = _approxBornNotRytov;

	int e = memoryAlloc();
	if (e != 0) return e;

	if (!preProcOnGPU) {
		memcpy(kxp_host, _kxp, n_proj * sizeof(float));
		memcpy(kyp_host, _kyp, n_proj * sizeof(float));

		// this can be moved later and performed while the slowest kernel function is executing

		int* shiftX = (int*)malloc(n_proj * sizeof(int));
		int* shiftY = (int*)malloc(n_proj * sizeof(int));

		float dxP = 1 / (dx * Nx);  // frequency sample(padded projection)

		for (int i = 0; i < n_proj; i++)
		{
			shiftX[i] = (int)roundf(kxp_host[i] / dxP);
			shiftY[i] = (int)roundf(kyp_host[i] / dxP);
		}

		cudaStatus3 = cudaMemcpy(xShift_dev, shiftX, n_proj * sizeof(int), cudaMemcpyHostToDevice);
		cudaStatus4 = cudaMemcpy(yShift_dev, shiftY, n_proj * sizeof(int), cudaMemcpyHostToDevice);

		free(shiftX); shiftX = NULL;
		free(shiftY); shiftY = NULL;


		cudaStatus5 = cudaMemcpy(kxp_dev, kxp_host, n_proj * sizeof(float), cudaMemcpyHostToDevice);
		cudaStatus6 = cudaMemcpy(kyp_dev, kyp_host, n_proj * sizeof(float), cudaMemcpyHostToDevice);

		if ((cudaStatus3 != cudaSuccess) || (cudaStatus4 != cudaSuccess) ||
			(cudaStatus5 != cudaSuccess) || (cudaStatus6 != cudaSuccess))
			return 59; // kxp and kyp or shiftX and shiftY vectors not transferred to GPU
	}
	else
	{
		cudaStatus5 = cudaMemcpy(kxp_host, kxp_dev, n_proj * sizeof(float), cudaMemcpyDeviceToHost);
		cudaStatus6 = cudaMemcpy(kyp_host, kyp_dev, n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	}


#pragma region radiusOfEwaldssSphere

	float* rf_host = (float*)malloc(n_proj * sizeof(float));
	if ((rf_dev == NULL) || (rf_host == NULL))
		return 51; //memory allocation error

	for (int i = 0; i < n_proj; i++)
		rf_host[i] = _NA / _lambda_all[i] / dkP;

	cudaStatus1 = cudaMemcpy(rf_dev, rf_host, n_proj * sizeof(float), cudaMemcpyHostToDevice);
	free(rf_host);

	if (cudaStatus1 != cudaSuccess)
		return 52; // data transfer to GPU error

#pragma endregion

#pragma region wavevectors

	float* kn_host = (float*)malloc(n_proj * sizeof(float));
	float* kz_host = (float*)malloc(n_proj * sizeof(float));
	if ((kn_dev == NULL) || (kn_host == NULL) || (kz_dev == NULL) || (kz_host == NULL))
		return 51; //memory allocation error

	float kn_m = 0;
	for (int i = 0; i < n_proj; i++) {
		kn_host[i] = n_imm / _lambda_all[i];
		kz_host[i] = sqrtf((kn_host[i] * kn_host[i]) - (kxp_host[i] * kxp_host[i]) - (kyp_host[i] * kyp_host[i]));
		kn_m += kn_host[i];
	}
	kn_mean = kn_m / n_proj;

	cudaStatus1 = cudaMemcpy(kn_dev, kn_host, n_proj * sizeof(float), cudaMemcpyHostToDevice);
	cudaStatus2 = cudaMemcpy(kz_dev, kz_host, n_proj * sizeof(float), cudaMemcpyHostToDevice);
	free(kn_host); kn_host = NULL;
	free(kz_host); kz_host = NULL;

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess))
		return 52; // data transfer to GPU error

#pragma endregion
#pragma region FFTplan
	int NN[2] = { Nx,Ny };

	int idist = Nx * Ny;
	int odist = Nx * Ny;
	int inembed[] = { Nx, Ny };
	int onembed[] = { Nx, Ny };
	int istride = 1;
	int ostride = 1;

	cufftResult ee = cufftPlanMany(&planFFTsino, 2, NN, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, n_proj);
	if (ee != cudaSuccess)
		return 55; // FFT plan creation error
#pragma endregion
	return 0;
}

EXPORTED_FUNCTION int sendSinograms(float* _sinoAmp, float* _sinoPh)
{
	cudaError_t cudaStatus1, cudaStatus2, cudaStatus3, cudaStatus4;
	if ((sinoAmp_dev != NULL) && (sinoPh_dev != NULL) && (_sinoAmp != NULL) && (_sinoPh != NULL)) {
		cudaStatus1 = cudaMemcpy(sinoAmp_dev, _sinoAmp, Nx * Ny * n_proj * sizeof(float), cudaMemcpyHostToDevice);
		cudaStatus2 = cudaMemcpy(sinoPh_dev, _sinoPh, Nx * Ny * n_proj * sizeof(float), cudaMemcpyHostToDevice);
	}
	else
		return 53; // sinograms not transferred to GPU

	if ((cudaStatus1 == cudaSuccess) && (cudaStatus2 == cudaSuccess))
		return 0;
	else
		return 53; // sinograms not transferred to GPU

}

EXPORTED_FUNCTION int sendFpmask_logical(unsigned char* Fpmask)
{
	unsigned int* fpmask_host = NULL;

	fpmask_host = (unsigned int*)malloc(Nx_Fpmask * Ny * n_proj * sizeof(int));
	int noOfInt = 0;
	int noOfBit = 0;

	for (int k = 0; k < n_proj; k++)
		for (int j = 0; j < Ny; j++) {
			for (int i = 0; i < Nx; i++) {
				noOfInt = i / 32;
				noOfBit = i % 32;
				if (Fpmask[Nx * Ny * k + Nx * j + i] == 0)
					fpmask_host[Nx_Fpmask * Ny * k + Nx_Fpmask * j + noOfInt] &= ~(1UL << noOfBit);
				else
					fpmask_host[Nx_Fpmask * Ny * k + Nx_Fpmask * j + noOfInt] |= (1UL << noOfBit);
			}
		}

	cudaError_t cudaStatus = cudaMemcpy(fpmask_dev, fpmask_host, Nx_Fpmask * Ny * n_proj * sizeof(int), cudaMemcpyHostToDevice);

	if ((cudaStatus != cudaSuccess))
		return 58; //error transferring the mask to the GPU

	free(fpmask_host);
	return 0;
}

EXPORTED_FUNCTION int generateKspace()
{
	LARGE_INTEGER countPerSec, tim1, tim2, tim3, tim4, tim5, tim6, tim7, tim8, tim9, tim10;

	int error;
	QueryPerformanceFrequency(&countPerSec);
	QueryPerformanceCounter(&tim1);
	// aproksymacja
	if (approxBornNotRytov == 1)
		error = approxBornK_Gen();
	else
		error = approxRytovK_Gen();
	QueryPerformanceCounter(&tim2);
	// ifftshift in 2D - done as in MATLAB - necessity to be verified
	if (error == 0)
		error = ifftShift_1(Nx, Ny, n_proj, projections_dev);
	QueryPerformanceCounter(&tim3);
	if (error == 0)
		error = ifftShift_2(Nx, Ny, n_proj, projections_dev);
	QueryPerformanceCounter(&tim4);
	// FFT
	if (error == 0)
		error = cufftExecC2C(planFFTsino, projections_dev, projections_dev, CUFFT_FORWARD);
	QueryPerformanceCounter(&tim5);
	// fftshift 2D because it is done as in MATLAB - necessity to be verified
	if (error == 0)
		error = ifftShift_1(Nx, Ny, n_proj, projections_dev);
	QueryPerformanceCounter(&tim6);
	if (error == 0)
		error = ifftShift_2(Nx, Ny, n_proj, projections_dev);
	QueryPerformanceCounter(&tim7);
	if (error == 0)
		error = multComplexVector(Nx * Ny * n_proj, projections_dev, dx * dx);

	QueryPerformanceCounter(&tim8);
	// transfer of data from the projection into a buffer of K-space size in 2D

								/*float* temp = (float*)malloc(Nx*Ny * 2 * n_proj * sizeof(float));
								memset(temp, 0, Nx*Ny *n_proj * 2 * sizeof(float));
								cudaError_t cudaStatus3 = cudaMemcpy(temp, projections_dev, Nx*Ny*n_proj * 2 * sizeof(float), cudaMemcpyDeviceToHost);
								FILE* plik8;
								if (preProcOnGPU)
									plik8 = fopen("projectionsPreprocOnGPU.bin", "wb");
								else
									plik8 = fopen("projections.bin", "wb");
								fwrite(temp, sizeof(float), 2 * Nx*Ny*n_proj, plik8);
								fclose(plik8);
								free(temp);*/
								// read
										/*	float* temp = (float*)malloc(Nx*Ny * 2 * n_proj * sizeof(float));
											FILE* plik8;
											plik8 = fopen("projections.bin", "rb");
											fread(temp, sizeof(float), 2 * Nx*Ny*n_proj, plik8);
											fclose(plik8);
											cudaError_t cudaStatus3 = cudaMemcpy( projections_dev, temp, Nx*Ny*n_proj * 2 * sizeof(float), cudaMemcpyHostToDevice);
											free(temp);*/



	if (error == 0)
		error = transferProjectionsToKbuffer();

	QueryPerformanceCounter(&tim9);
	if (error == 0)
		error = divide_KO_by_EW(K_xy * K_xy * K_z, KO_dev, EW_dev);

	//float * temp = (float*)malloc(K_xy*K_xy * 2 * K_xy * sizeof(float));
	//memset(temp, 0, K_xy*K_xy * 2 * K_xy * sizeof(float));
	//cudaError_t cudaStatus3 = cudaMemcpy(temp, KO_dev, K_xy*K_xy * 2 * K_xy * sizeof(float), cudaMemcpyDeviceToHost);
	//FILE *plik8;
	//if (preProcOnGPU)
	//	plik8 = fopen("KO_poDivByEW_PreprocOnGPU_nd.bin", "wb");
	//else
	//	plik8 = fopen("KO_poDivByEW.bin", "wb");
	//fwrite(temp, sizeof(float), K_xy*K_xy * 2 * K_xy, plik8);
	//fclose(plik8);
	//free(temp);

/*		int * tempI = (int*)malloc(K_xy*K_xy * K_xy * sizeof(int));
		memset(tempI, 0, K_xy*K_xy * K_xy * sizeof(int));
		cudaError_t cudaStatus13 = cudaMemcpy(tempI, EW_dev, K_xy*K_xy * K_xy * sizeof(int), cudaMemcpyDeviceToHost);
		FILE *plik9;
		if (preProcOnGPU)
			plik9 = fopen("EW_PreprocOnGPU.bin", "wb");
		else
			plik9 = fopen("EW.bin", "wb");
		fwrite(tempI, sizeof(int), K_xy*K_xy * K_xy, plik9);
		fclose(plik9);
		free(tempI);*/

	QueryPerformanceCounter(&tim10);

	double j_total = (double)(tim10.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;
	double* times = (double*)malloc(15 * sizeof(double));
	times[0] = j_total;
	times[1] = (double)(tim2.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;
	times[2] = (double)(tim3.QuadPart - tim2.QuadPart) / countPerSec.QuadPart * 1000;
	times[3] = (double)(tim4.QuadPart - tim3.QuadPart) / countPerSec.QuadPart * 1000;
	times[4] = (double)(tim5.QuadPart - tim4.QuadPart) / countPerSec.QuadPart * 1000;
	times[5] = (double)(tim6.QuadPart - tim5.QuadPart) / countPerSec.QuadPart * 1000;
	times[6] = (double)(tim7.QuadPart - tim6.QuadPart) / countPerSec.QuadPart * 1000;
	times[7] = (double)(tim8.QuadPart - tim7.QuadPart) / countPerSec.QuadPart * 1000;
	times[8] = (double)(tim9.QuadPart - tim8.QuadPart) / countPerSec.QuadPart * 1000;
	times[9] = (double)(tim10.QuadPart - tim9.QuadPart) / countPerSec.QuadPart * 1000;

	//size_t totalM, freeM;
	//cudaMemGetInfo(&freeM, &totalM);

	//freeM /= 1024; //w KB
	//freeM /= 1024; // w MB

	//totalM /= 1024; //w KB
	//totalM /= 1024; // w MB

	cudaFree(kz_dev);

	return 0;
}

EXPORTED_FUNCTION int takeDataHoloCSGdev(float* complexData_host)
{
	cudaError_t cudaStatus1 = cudaMemcpy(complexData_host, complexFieldCSG_device, 2 * D * D * n_proj * sizeof(float), cudaMemcpyDeviceToHost);

	if (cudaStatus1 != cudaSuccess) {
		return 56;
	}

	//freeDeviceMemory();

	return 0;
}


EXPORTED_FUNCTION int takeDataKspace(float* complexData_host)
{
	cudaError_t cudaStatus1 = cudaMemcpy(complexData_host, KO_dev, 2 * K_xy * K_xy * K_z * sizeof(float), cudaMemcpyDeviceToHost);

	if (cudaStatus1 != cudaSuccess) {
		return 56;
	}

	//freeDeviceMemory();

	return 0;
}

EXPORTED_FUNCTION int takePeaksInt(int* peaks_host)
{
	cudaError_t cudaStatus1 = cudaMemcpy(peaks_host, peaksStep2_dev, 2 * n_proj * sizeof(int), cudaMemcpyDeviceToHost);

	if (cudaStatus1 != cudaSuccess) {
		return 57;
	}

	//freeDeviceMemory();

	return 0;
}

EXPORTED_FUNCTION int takePeaksFloat(float* peaks_host)
{
	cudaError_t cudaStatus1 = cudaMemcpy(peaks_host, peaksStep1_dev, 2 * D * n_proj * sizeof(float), cudaMemcpyDeviceToHost);

	if (cudaStatus1 != cudaSuccess) {
		return 57;
	}

	//freeDeviceMemory();

	return 0;
}

EXPORTED_FUNCTION int takeDataEW(int* EW_host)
{
	cudaError_t cudaStatus1 = cudaMemcpy(EW_host, EW_dev, K_xy * K_xy * K_z * sizeof(int), cudaMemcpyDeviceToHost);

	if (cudaStatus1 != cudaSuccess) {
		return 57;
	}

	//freeDeviceMemory();

	return 0;
}


EXPORTED_FUNCTION int takeProjections(float* complexData_host)
{
	cudaError_t cudaStatus1 = cudaMemcpy(complexData_host, projections_dev, 2 * Nx * Ny * n_proj * sizeof(float), cudaMemcpyDeviceToHost);

	if (cudaStatus1 != cudaSuccess) {
		return 58;
	}

	//freeDeviceMemory();

	return 0;
}

EXPORTED_FUNCTION int freeGPUmemory()
{
	freeDeviceMemory();

	return 0;
}

int tukeyWindow(float* data, int L, int r)
{
	int i = 0;

	if ((data != NULL) && (L > r)) {
		for (i = 0; i < (r / 2); i++)
			data[i] = 0.5f * (1 + cos((2 * PI / r) * (i - (r / 2))));
		for (; i < (L - (r / 2)); i++)
			data[i] = 1;
		for (; i < L; i++)
			data[i] = 0.5f * (1 + cos((2 * PI / r) * (i + 1 - L + (r / 2))));
	}
	else return 100;

	return 0;
}

int tukeyWindow2D(float* tukey2D, int Nx, float rx, int Ny, float ry)
{
	float* tukey1D_x = (float*)malloc(Nx * sizeof(float));
	if (tukey1D_x == NULL) return 62;

	tukeyWindow(tukey1D_x, Nx, Nx * rx);

	float* tukey1D_y = (float*)malloc(Ny * sizeof(float));
	if (tukey1D_y == NULL) return 62;

	tukeyWindow(tukey1D_y, Ny, Ny * ry);

	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			tukey2D[j + i * Ny] = tukey1D_x[i] * tukey1D_y[j];

	free(tukey1D_x);
	free(tukey1D_y);

	return 0;
}


int circTukey2D(float* tukey2D, int D, float alfa)
{
	int interpolation = 10;

	float* tukey1D_x = (float*)malloc(interpolation * D * sizeof(float));
	if (tukey1D_x == NULL) return 62;

	int err = tukeyWindow(tukey1D_x, interpolation * D, interpolation * D * alfa);
	if (err) return 63;

	int mid = D / 2;
	float r = 0;
	int index = 0;

	for (int i = 0; i < D; i++)
		for (int j = 0; j < D; j++) {
			r = sqrtf(((i - mid) * (i - mid)) + ((j - mid) * (j - mid)));
			index = interpolation * mid - roundf(interpolation * r);
			if (index < 0) index = 0;
			tukey2D[j + i * D] = tukey1D_x[index];
		}
	free(tukey1D_x);

	return 0;
}

int allocBuffersForPreproc(short int* hologram)
{
	cudaError_t cudaStatus1, cudaStatus2, cudaStatus3, cudaStatus4, cudaStatus5, cudaStatus6, cudaStatus7,
		cudaStatus8, cudaStatus9, cudaStatus10, cudaStatus11, cudaStatus12, cudaStatus13, cudaStatus14, cudaStatus15,
		cudaStatus16, cudaStatus17;

	cudaStatus9 = cudaMalloc((void**)&sinoAmp_dev, D * D * n_proj * sizeof(float));
	cudaStatus10 = cudaMalloc((void**)&sinoPh_dev, D * D * n_proj * sizeof(float));
	cudaStatus17 = cudaMalloc((void**)&sinoPhComplex_dev, D * D * n_proj * sizeof(cufftComplex));

	cudaStatus12 = cudaMalloc((void**)&kxp_dev, n_proj * sizeof(float));
	cudaStatus13 = cudaMalloc((void**)&kyp_dev, n_proj * sizeof(float));
	cudaStatus16 = cudaMalloc((void**)&kxy_dev, n_proj * sizeof(float));

	cudaStatus14 = cudaMalloc((void**)&xShift_dev, n_proj * sizeof(int));
	cudaStatus15 = cudaMalloc((void**)&yShift_dev, n_proj * sizeof(int));

	cudaStatus1 = cudaMalloc((void**)&holograms_device, X * Y * n_proj * sizeof(short int));
	cudaStatus2 = cudaMalloc((void**)&hologramsCSG_device, X * Y * n_proj * sizeof(cufftComplex));
	cudaStatus3 = cudaMalloc((void**)&tukeyWin2D_W_dev, Nx_holo * Ny_holo * sizeof(float));
	cudaStatus4 = cudaMalloc((void**)&croppedPartsHoloCSG_device, D * D * n_proj * sizeof(cufftComplex));
	cudaStatus5 = cudaMalloc((void**)&sMask_dev, D * D * sizeof(float));
	cudaStatus11 = cudaMalloc((void**)&mask_dev, D * D * sizeof(float));
	cudaStatus6 = cudaMalloc((void**)&peaksStep1_dev, 2 * D * n_proj * sizeof(int));
	cudaStatus7 = cudaMalloc((void**)&peaksStep2_dev, 2 * n_proj * sizeof(int));
	cudaStatus8 = cudaMalloc((void**)&complexFieldCSG_device, D * D * n_proj * sizeof(cufftComplex));

	if ((holograms_device != NULL) && (cudaStatus1 == cudaSuccess)) {
		cudaStatus1 = cudaMemcpy(holograms_device, hologram, X * Y * n_proj * sizeof(short int), cudaMemcpyHostToDevice);
		if (cudaStatus1 != cudaSuccess)
			return 61; // holograms not transferred to GPU
	}
	else
		return 61; // holograms not transferred to GPU

	if ((hologramsCSG_device == NULL) || (cudaStatus2 != cudaSuccess))
		return 64;

	if ((tukeyWin2D_W_dev == NULL) || (cudaStatus3 != cudaSuccess))
		return 67;

	if ((mask_dev == NULL) || (cudaStatus11 != cudaSuccess))
		return 101;

	if ((croppedPartsHoloCSG_device == NULL) || (cudaStatus4 != cudaSuccess))
		return 70;

	if ((peaksStep1_dev == NULL) || (cudaStatus6 != cudaSuccess))
		return 73;

	if ((peaksStep2_dev == NULL) || (cudaStatus7 != cudaSuccess))
		return 74;

	if ((complexFieldCSG_device == NULL) || (cudaStatus8 != cudaSuccess))
		return 77;

	if ((sinoAmp_dev == NULL) || (cudaStatus9 != cudaSuccess))
		return 80;

	if ((sinoPh_dev == NULL) || (cudaStatus10 != cudaSuccess))
		return 81;

	if ((sinoPhComplex_dev == NULL) || (cudaStatus17 != cudaSuccess))
		return 107;

	if ((kxp_dev == NULL) || (cudaStatus12 != cudaSuccess))
		return 106;

	if ((kyp_dev == NULL) || (cudaStatus13 != cudaSuccess))
		return 113;

	if ((xShift_dev == NULL) || (cudaStatus14 != cudaSuccess))
		return 109;

	if ((yShift_dev == NULL) || (cudaStatus15 != cudaSuccess))
		return 110;

	if ((kxy_dev == NULL) || (cudaStatus16 != cudaSuccess))
		return 111;

	return 0;
}

int allocBuffersForPreproc_Reference(short int* hologram)
{
	cudaError_t cudaStatus1, cudaStatus2, cudaStatus3, cudaStatus4, cudaStatus5, cudaStatus6, cudaStatus7,
		cudaStatus8, cudaStatus9, cudaStatus10, cudaStatus11, cudaStatus12, cudaStatus13, cudaStatus14, cudaStatus15,
		cudaStatus16, cudaStatus17;

	cudaStatus9 = cudaMalloc((void**)&sinoAmpRef_dev, D * D * n_proj * sizeof(float));
	cudaStatus10 = cudaMalloc((void**)&sinoPhRef_dev, D * D * n_proj * sizeof(float));
	cudaStatus17 = cudaMalloc((void**)&sinoPhRefComplex_dev, D * D * n_proj * sizeof(cufftComplex));

	/*cudaStatus12 = cudaMalloc((void**)&kxp_dev, n_proj * sizeof(float));
	cudaStatus13 = cudaMalloc((void**)&kyp_dev, n_proj * sizeof(float));
	cudaStatus16 = cudaMalloc((void**)&kxy_dev, n_proj * sizeof(float));

	cudaStatus14 = cudaMalloc((void**)&xShift_dev, n_proj * sizeof(int));
	cudaStatus15 = cudaMalloc((void**)&yShift_dev, n_proj * sizeof(int));*/

	cudaStatus1 = cudaMalloc((void**)&holograms_device, X * Y * n_proj * sizeof(short int));
	cudaStatus2 = cudaMalloc((void**)&hologramsCSG_device, X * Y * n_proj * sizeof(cufftComplex));
	cudaStatus3 = cudaMalloc((void**)&tukeyWin2D_W_dev, Nx_holo * Ny_holo * sizeof(float));
	cudaStatus4 = cudaMalloc((void**)&croppedPartsHoloCSG_device, D * D * n_proj * sizeof(cufftComplex));
	cudaStatus5 = cudaMalloc((void**)&sMask_dev, D * D * sizeof(float));
	cudaStatus11 = cudaMalloc((void**)&mask_dev, D * D * sizeof(float));
	cudaStatus6 = cudaMalloc((void**)&peaksStep1_dev, 2 * D * n_proj * sizeof(int));
	cudaStatus7 = cudaMalloc((void**)&peaksStep2_dev, 2 * n_proj * sizeof(int));
	cudaStatus8 = cudaMalloc((void**)&complexFieldCSG_device, D * D * n_proj * sizeof(cufftComplex));

	if ((holograms_device != NULL) && (cudaStatus1 == cudaSuccess)) {
		cudaStatus1 = cudaMemcpy(holograms_device, hologram, X * Y * n_proj * sizeof(short int), cudaMemcpyHostToDevice);
		if (cudaStatus1 != cudaSuccess)
			return 61; // holograms not transferred to GPU
	}
	else
		return 61; // holograms not transferred to GPU

	if ((hologramsCSG_device == NULL) || (cudaStatus2 != cudaSuccess))
		return 64;

	if ((tukeyWin2D_W_dev == NULL) || (cudaStatus3 != cudaSuccess))
		return 67;

	if ((mask_dev == NULL) || (cudaStatus11 != cudaSuccess))
		return 101;

	if ((croppedPartsHoloCSG_device == NULL) || (cudaStatus4 != cudaSuccess))
		return 70;

	if ((peaksStep1_dev == NULL) || (cudaStatus6 != cudaSuccess))
		return 73;

	if ((peaksStep2_dev == NULL) || (cudaStatus7 != cudaSuccess))
		return 74;

	if ((complexFieldCSG_device == NULL) || (cudaStatus8 != cudaSuccess))
		return 77;

	if ((sinoAmpRef_dev == NULL) || (cudaStatus9 != cudaSuccess))
		return 80;

	if ((sinoPhRef_dev == NULL) || (cudaStatus10 != cudaSuccess))
		return 81;

	if ((sinoPhRefComplex_dev == NULL) || (cudaStatus17 != cudaSuccess))
		return 105;

	/*if ((kxp_dev == NULL) || (cudaStatus12 != cudaSuccess))
		return 106;

	if ((kyp_dev == NULL) || (cudaStatus13 != cudaSuccess))
		return 107;

	if ((xShift_dev == NULL) || (cudaStatus14 != cudaSuccess))
		return 109;

	if ((yShift_dev == NULL) || (cudaStatus15 != cudaSuccess))
		return 110;

	if ((kxy_dev == NULL) || (cudaStatus16 != cudaSuccess))
		return 111;*/

	return 0;
}

__global__ void kernel_sum_fin(float* resultsMid, float* results, int D)
{
	register int index_out = blockIdx.x * D + threadIdx.x;

	__shared__ float buf[1024];

	register int i_loc = threadIdx.x;

	if (i_loc < D)
		buf[i_loc] = resultsMid[index_out];
	else
		buf[i_loc] = 0;

	for (int krok = blockDim.x >> 1; krok > 0; krok >>= 1)
	{
		__syncthreads();
		if (i_loc < krok)
			buf[i_loc] += buf[i_loc + krok];
	}
	__syncthreads();
	if (i_loc == 0)
		results[blockIdx.x] = buf[0]; // value
}

__global__ void kernel_sum_fin_b2(float* resultsMid, float* results, int D)
{
	register int index_out = blockIdx.x * D + threadIdx.x;

	__shared__ float buf[1024];

	register int i_loc = threadIdx.x;

	if (i_loc < D)
		buf[i_loc] = resultsMid[index_out] * threadIdx.x;
	else
		buf[i_loc] = 0;

	for (int krok = blockDim.x >> 1; krok > 0; krok >>= 1)
	{
		__syncthreads();
		if (i_loc < krok)
			buf[i_loc] += buf[i_loc + krok];
	}
	__syncthreads();
	if (i_loc == 0)
		results[blockIdx.x] = buf[0]; // value
}


__global__ void kernel_sum(float* framesIn, float* resultsMid, int D)
{
	register int index = blockIdx.y * D * D + blockIdx.x * D + threadIdx.x;

	__shared__ float buf[1024];

	register int i_loc = threadIdx.x;

	if (i_loc < D)
		buf[i_loc] = framesIn[index];
	else
		buf[i_loc] = 0;


	for (int krok = blockDim.x >> 1; krok > 0; krok >>= 1)
	{
		__syncthreads();
		if (i_loc < krok)
			buf[i_loc] += buf[i_loc + krok];
	}

	if (i_loc == 0)
		resultsMid[blockIdx.y * gridDim.x + blockIdx.x] = buf[0]; // value
}

__global__ void kernel_sum_b1(float* framesIn, float* resultsMid, int D)
{
	register int index = blockIdx.y * D * D + blockIdx.x * D + threadIdx.x;

	__shared__ float buf[1024];

	register int i_loc = threadIdx.x;

	if (i_loc < D)
		buf[i_loc] = framesIn[index] * threadIdx.x;
	else
		buf[i_loc] = 0;


	for (int krok = blockDim.x >> 1; krok > 0; krok >>= 1)
	{
		__syncthreads();
		if (i_loc < krok)
			buf[i_loc] += buf[i_loc + krok];
	}

	if (i_loc == 0)
		resultsMid[blockIdx.y * gridDim.x + blockIdx.x] = buf[0]; // value
}

__global__ void fillPhase_testOnly(float* phaseOut_dev)
{
	register int index_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int index_glob = blockIdx.y * blockDim.x * gridDim.x + index_loc;

	phaseOut_dev[index_glob] = threadIdx.x;
}

__global__ void move2frameTo1_testOnly(float* phaseOut_dev)
{
	register int index_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int index_glob = blockIdx.y * blockDim.x * gridDim.x + index_loc;

	phaseOut_dev[index_glob] = phaseOut_dev[index_glob + blockDim.x * gridDim.x];
}

__global__ void kernel_determine_fj_fi_f0(float* fj_fi_f0, float* a, float* b1, float* b2, float* b3)
{
	//a[0] a[1] a[2] a[3] a[4] a[5]
	//a11_  12_  13_  22_  23_  33
	float delta = a[0] * (a[3] * a[5] - a[4] * a[4]) + a[1] * (a[4] * a[2] - a[1] * a[5]) + a[2] * (a[1] * a[4] - a[3] * a[2]);

	int i = threadIdx.x;

	// fj
	fj_fi_f0[3 * i] = (b1[i] * (a[3] * a[5] - a[4] * a[4]) + a[1] * (a[4] * b3[i] - b2[i] * a[5]) + a[2] * (b2[i] * a[4] - a[3] * b3[i])) / delta;
	// fi 
	fj_fi_f0[3 * i + 1] = (a[0] * (b2[i] * a[5] - a[4] * b3[i]) + b1[i] * (a[4] * a[2] - a[1] * a[5]) + a[2] * (a[1] * b3[i] - b2[i] * a[2])) / delta;
	// f0
	fj_fi_f0[3 * i + 2] = (a[0] * (a[3] * b3[i] - b2[i] * a[4]) + a[1] * (b2[i] * a[2] - a[1] * b3[i]) + b1[i] * (a[1] * a[4] - a[3] * a[2])) / delta;
}

__global__ void phaseCorrection(float* phaseInOut_dev, float* fj_fi_f0)
{
	register int index_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int index_glob = blockIdx.y * blockDim.x * gridDim.x + index_loc;

	int i = blockIdx.y; // n_proj

	phaseInOut_dev[index_glob] = phaseInOut_dev[index_glob] - (fj_fi_f0[3 * i] * threadIdx.x + fj_fi_f0[3 * i + 1] * blockIdx.x + fj_fi_f0[3 * i + 2]);
}

int detrend2D(float* sPh_dev)
{
	long long int an1 = D;
	long long int an2 = D;
	long long int a11 = an2 * an1 * (an1 + 1) * (2 * an1 + 1) / 6;
	long long int a12 = an1 * (an1 + 1) * an2 * (an2 + 1) / 4;
	long long int a13 = an2 * an1 * (an1 + 1) / 2;
	long long int a22 = an1 * an2 * (an2 + 1) * (2 * an2 + 1) / 6;
	long long int a23 = an1 * an2 * (an2 + 1) / 2;
	long long int a33 = an1 * an2;
	int size = sizeof(long long int);

	float* a_params_host;
	a_params_host = (float*)malloc(6 * sizeof(float));
	a_params_host[0] = a11;
	a_params_host[1] = a12;
	a_params_host[2] = a13;
	a_params_host[3] = a22;
	a_params_host[4] = a23;
	a_params_host[5] = a33;


	float* a11_12_13_22_23_33 = NULL;
	cudaError_t cudaStatus5 = cudaMalloc((void**)&a11_12_13_22_23_33, 6 * sizeof(float));
	if ((a11_12_13_22_23_33 == NULL) || (cudaStatus5 != cudaSuccess)) return 93;

	cudaError_t cudaStatus9 = cudaMemcpy(a11_12_13_22_23_33, a_params_host, 6 * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus9 != cudaSuccess) return 94;

	float* fj_fi_f0_params = NULL;
	cudaError_t cudaStatus6 = cudaMalloc((void**)&fj_fi_f0_params, 3 * n_proj * sizeof(float));
	if ((fj_fi_f0_params == NULL) || (cudaStatus6 != cudaSuccess)) return 95;

	float* temp_sum = NULL;
	cudaError_t cudaStatus1 = cudaMalloc((void**)&temp_sum, D * n_proj * sizeof(float));
	if ((temp_sum == NULL) || (cudaStatus1 != cudaSuccess)) return 85;

	float* results_sum_b3 = NULL;
	cudaStatus1 = cudaMalloc((void**)&results_sum_b3, n_proj * sizeof(float));
	if (results_sum_b3 == NULL) return 86;

	float* results_sum_b1 = NULL;
	cudaStatus1 = cudaMalloc((void**)&results_sum_b1, n_proj * sizeof(float));
	if (results_sum_b1 == NULL) return 88;

	float* results_sum_b2 = NULL;
	cudaStatus1 = cudaMalloc((void**)&results_sum_b2, n_proj * sizeof(float));
	if (results_sum_b2 == NULL) return 89;
	// B3
	int numThreads = 1024;
	if (D < 513) numThreads = 512;
	if (D < 257) numThreads = 256;
	dim3 nThreads(numThreads, 1, 1);  /// resultant tomogram height
	dim3 nBlocks(D, n_proj, 1);

	kernel_sum << <nBlocks, nThreads >> > (sPh_dev, temp_sum, D);

	cudaError_t cudaStatus = cudaDeviceSynchronize();

	if (cudaStatus != cudaSuccess) {
		return 87;
	}

	dim3 nThreads2(numThreads, 1, 1);  /// resultant tomogram height
	dim3 nBlocks2(n_proj, 1, 1);
	kernel_sum_fin << <nBlocks2, nThreads2 >> > (temp_sum, results_sum_b3, D);

	cudaStatus = cudaDeviceSynchronize();

	if (cudaStatus != cudaSuccess) {
		return 90;
	}
	// end B3

	// B1
	kernel_sum_b1 << <nBlocks, nThreads >> > (sPh_dev, temp_sum, D);

	cudaStatus = cudaDeviceSynchronize();

	if (cudaStatus != cudaSuccess) {
		return 89;
	}

	kernel_sum_fin << <nBlocks2, nThreads2 >> > (temp_sum, results_sum_b1, D);

	cudaStatus = cudaDeviceSynchronize();

	if (cudaStatus != cudaSuccess) {
		return 91;
	}
	// end B1

	// B2
	kernel_sum << <nBlocks, nThreads >> > (sPh_dev, temp_sum, D);

	cudaStatus = cudaDeviceSynchronize();

	if (cudaStatus != cudaSuccess) {
		return 92;
	}

	kernel_sum_fin_b2 << <nBlocks2, nThreads2 >> > (temp_sum, results_sum_b2, D);

	cudaStatus = cudaDeviceSynchronize();

	if (cudaStatus != cudaSuccess) {
		return 90;
	}
	// end B2

	/// determine params: fj, fi, f0 for each frame
	dim3 nThreads3(n_proj, 1, 1);  /// resultant tomogram height
	dim3 nBlocks3(1, 1, 1);
	kernel_determine_fj_fi_f0 << <nBlocks3, nThreads3 >> > (fj_fi_f0_params, a11_12_13_22_23_33, results_sum_b1,
		results_sum_b2, results_sum_b3);

	/// end of determine params: fj, fi, f0
	dim3 block5(D, 1, 1);
	dim3 grid5(D, n_proj, 1);

	phaseCorrection << <grid5, block5 >> > (sPh_dev, fj_fi_f0_params);



	//!!!!!!!!!!!!!!!!!!!! cudaFree !!!!!!!!!!!!!!!!!!!!
	cudaFree(a11_12_13_22_23_33);
	cudaFree(fj_fi_f0_params);
	cudaFree(temp_sum);
	cudaFree(results_sum_b3);
	cudaFree(results_sum_b1);
	cudaFree(results_sum_b2);
	return 0;
}

__global__ void amplitudeProc(float* amp_dev, float* mask, float* sums)
{
	register int index_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int index_glob = blockIdx.y * blockDim.x * gridDim.x + index_loc;

	float amplitude = amp_dev[index_glob];

	float meanValue = sums[blockIdx.y] / (blockDim.x * gridDim.x);

	amplitude = amplitude / meanValue;

	amp_dev[index_glob] = (amplitude - 1) * mask[index_loc] + 1;
}

__global__ void phaseProc(float* ph_dev, float* mask, float* correctors)
{
	register int index_loc = blockIdx.x * blockDim.x + threadIdx.x;
	register int index_glob = blockIdx.y * blockDim.x * gridDim.x + index_loc;

	float phase = ph_dev[index_glob] - correctors[blockIdx.y];

	ph_dev[index_glob] = mask[index_loc] * phase;
}

__global__ void kernel_sort(float* phaseIn, float* maskIn, float* results_sort, int D, int do_NNC)
{
	register int index1 = blockIdx.x * D * D + threadIdx.x + D * (D >> 1); //horizontal
	register int index2 = blockIdx.x * D * D + threadIdx.x * D + (D >> 1); //vertical
	register int index3 = blockIdx.x * D * D + threadIdx.x * (D + 1); //diagonal1
	register int index4 = blockIdx.x * D * D + threadIdx.x * (D - 1) + (D - 1);  //diagonal2

	extern __shared__ float buf[];

	register int i_loc = threadIdx.x;
	register double temp = 0;

	register int i_loc_x4 = 4 * i_loc;

	buf[i_loc_x4] = phaseIn[index1] * maskIn[threadIdx.x + D * (D >> 1)];
	buf[i_loc_x4 + 1] = phaseIn[index2] * maskIn[threadIdx.x * D + (D >> 1)];
	buf[i_loc_x4 + 2] = phaseIn[index3] * maskIn[threadIdx.x * (D + 1)];
	buf[i_loc_x4 + 3] = phaseIn[index4] * maskIn[threadIdx.x * (D - 1) + (D - 1)];

	__syncthreads();


	for (int i = 0; i < 4 * D; i++) {

		if (buf[i_loc_x4] > buf[i_loc_x4 + 1]) {
			temp = buf[i_loc_x4];
			buf[i_loc_x4] = buf[i_loc_x4 + 1];
			buf[i_loc_x4 + 1] = temp;
		}

		if (buf[i_loc_x4 + 1] > buf[i_loc_x4 + 2]) {
			temp = buf[i_loc_x4 + 1];
			buf[i_loc_x4 + 1] = buf[i_loc_x4 + 2];
			buf[i_loc_x4 + 2] = temp;
		}

		if (buf[i_loc_x4 + 2] > buf[i_loc_x4 + 3]) {
			temp = buf[i_loc_x4 + 2];
			buf[i_loc_x4 + 2] = buf[i_loc_x4 + 3];
			buf[i_loc_x4 + 3] = temp;
		}
		__syncthreads();  //nessesary because "+4" sample is from the next warp
		if ((i_loc_x4 + 4) < 4 * D)
			if (buf[i_loc_x4 + 3] > buf[i_loc_x4 + 4]) {
				temp = buf[i_loc_x4 + 3];
				buf[i_loc_x4 + 3] = buf[i_loc_x4 + 4];
				buf[i_loc_x4 + 4] = temp;
			}
		__syncthreads();   //nessesary because "+4" sample is from the next warp
	}

	if (threadIdx.x == 0)
		if (do_NNC == -1)
			results_sort[blockIdx.x] = buf[(int)roundf(0.95f * 4 * D)];
		else
			results_sort[blockIdx.x] = buf[(int)roundf(0.05f * 4 * D)];
}

int amplAndPhaseProc(float* sAmp_dev, float* sPh_dev, float* mask)
{
	// sum to: mean(amp(:))
	float* temp_sum = NULL;
	cudaError_t cudaStatus1 = cudaMalloc((void**)&temp_sum, D * n_proj * sizeof(float));
	if ((temp_sum == NULL) || (cudaStatus1 != cudaSuccess)) return 96;

	float* results_sum_or_sort = NULL;
	cudaStatus1 = cudaMalloc((void**)&results_sum_or_sort, n_proj * sizeof(float));
	if (results_sum_or_sort == NULL) return 97;

	int numThreads = 1024;
	if (D < 513) numThreads = 512;
	if (D < 257) numThreads = 256;
	dim3 nThreads(numThreads, 1, 1);  /// resultant tomogram height
	dim3 nBlocks(D, n_proj, 1);

	kernel_sum << <nBlocks, nThreads >> > (sAmp_dev, temp_sum, D);

	cudaError_t cudaStatus = cudaDeviceSynchronize();

	if (cudaStatus != cudaSuccess) {
		return 98;
	}

	dim3 nThreads2(numThreads, 1, 1);  /// resultant tomogram height
	dim3 nBlocks2(n_proj, 1, 1);
	kernel_sum_fin << <nBlocks2, nThreads2 >> > (temp_sum, results_sum_or_sort, D);

	cudaStatus = cudaDeviceSynchronize();

	if (cudaStatus != cudaSuccess) {
		return 99;
	}

	/*	float* temp = (float*)malloc( n_proj *sizeof(float));
		cudaError_t cudaStatus11 = cudaMemcpy(temp, results_sum, n_proj * sizeof(float), cudaMemcpyDeviceToHost);
		//FILE* plik8;
		//plik8 = fopen("testyMediany.bin", "wb");
		//fwrite(temp, sizeof(float), n_proj, plik8);
		//fclose(plik8);
		free(temp);*/

		//amp = amp. / mean(amp(:));  % normalized amplitude
		//amp = (amp - 1).*mask + 1; % apply Tukey mask

	dim3 block2(D, 1, 1);
	dim3 grid2(D, n_proj, 1);
	amplitudeProc << <grid2, block2 >> > (sAmp_dev, mask, results_sum_or_sort);

	cudaStatus = cudaDeviceSynchronize();

	if (cudaStatus != cudaSuccess) {
		return 102;
	}

	//Phase:


	/// debug olny
	/*float* debugOnly;
	cudaStatus1 = cudaMalloc((void**)&debugOnly, n_proj * D * 4 * sizeof(float));

	float* debugOnly2;
	cudaStatus1 = cudaMalloc((void**)&debugOnly2, n_proj * D * 4 * sizeof(float));*/
	///


	dim3 block3(D, 1, 1);
	dim3 grid3(n_proj, 1, 1);
	kernel_sort << <grid3, block3, 4 * D * sizeof(float) >> > (sPh_dev, mask, results_sum_or_sort, D, do_NNC);

	cudaStatus = cudaDeviceSynchronize();

	cudaStatus = cudaGetLastError();

	if (cudaStatus != cudaSuccess) {
		return 103;
	}

	/*float* temp = (float*)malloc( n_proj *D *4*sizeof(float));
	cudaError_t cudaStatus11 = cudaMemcpy(temp, debugOnly, n_proj * D*4* sizeof(float), cudaMemcpyDeviceToHost);
	FILE* plik8;
	plik8 = fopen("poSortowaniu.bin", "wb");
	fwrite(temp, sizeof(float), n_proj*D*4, plik8);
	fclose(plik8);
	free(temp);

	float* temp2 = (float*)malloc(n_proj *D * 4 * sizeof(float));
	cudaError_t cudaStatus13 = cudaMemcpy(temp2, debugOnly2, n_proj * D * 4 * sizeof(float), cudaMemcpyDeviceToHost);
	FILE* plik9;
	plik9 = fopen("przedSortowaniem.bin", "wb");
	fwrite(temp, sizeof(float), n_proj*D * 4, plik9);
	fclose(plik9);
	free(temp2);*/

	/*
	temp = (float*)malloc(n_proj * sizeof(float));
	cudaStatus11 = cudaMemcpy(temp, results_sum_or_sort, n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("wynikiMediany.bin", "wb");
	fwrite(temp, sizeof(float), n_proj, plik8);
	fclose(plik8);
	free(temp);*/


	dim3 block4(D, 1, 1);
	dim3 grid4(D, n_proj, 1);
	phaseProc << <grid4, block4 >> > (sPh_dev, mask, results_sum_or_sort);

	cudaStatus = cudaDeviceSynchronize();

	cudaStatus = cudaGetLastError();

	if (cudaStatus != cudaSuccess) {
		return 104;
	}

	cudaFree(temp_sum);
	cudaFree(results_sum_or_sort);

	return 0;
}

__global__ void kernel_ray_XY(int nproj, float* kxp_dev, float* kyp_dev, float* kxy_dev, int* xShift_dev, int* yShift_dev, float multiplier, int* peaks)
{
	__shared__ float buf[1024];

	if (threadIdx.x < nproj) {
		int dy = peaks[2 * threadIdx.x] - peaks[0];
		int dx = peaks[2 * threadIdx.x + 1] - peaks[1];

		xShift_dev[threadIdx.x] = dx;
		yShift_dev[threadIdx.x] = dy;

		float kxp = dx * multiplier;
		float kyp = dy * multiplier;

		kyp_dev[threadIdx.x] = kyp;
		kxp_dev[threadIdx.x] = kxp;

		buf[threadIdx.x] = sqrt((kxp * kxp) + (kyp * kyp));
		kxy_dev[threadIdx.x] = buf[threadIdx.x];
	}
	else
	{
		buf[threadIdx.x] = 0;
	}


	for (int krok = blockDim.x >> 1; krok > 0; krok >>= 1)
	{
		__syncthreads();
		if (threadIdx.x < krok)
			buf[threadIdx.x] += buf[threadIdx.x + krok];
	}
	__syncthreads();
	if (threadIdx.x == 0)
		kxy_dev[0] = buf[0] / (nproj - 1); // value
}

int ray_XY(float lambda, float cam_pix, float M)
{
	float system_scaling_factor = cam_pix * Nx_holo / M;
	float multiplier = 1 / (system_scaling_factor);//lambda / (system_scaling_factor * n_imm);


	int nThr = n_proj;
	if (nThr <= 128)
		nThr = 128;
	else if (nThr <= 256)
		nThr = 256;
	else if (nThr <= 512)
		nThr = 512;
	else if (nThr <= 1024)
		nThr = 1024;
	dim3 nThreads(nThr, 1, 1);  /// resultant tomogram height
	dim3 nBlocks(1, 1, 1);

	kernel_ray_XY << <nBlocks, nThreads >> > (n_proj, kxp_dev, kyp_dev, kxy_dev, xShift_dev, yShift_dev, multiplier, peaksStep2_dev);

	cudaError_t cudaStatus = cudaGetLastError();

	if (cudaStatus != cudaSuccess) {
		return 108;
	}
	return 0;
}

__global__ void kernel_takeReference(float* sAmp_dev, float* sPh_dev, float* sAmpR_dev, float* sPhR_dev)
{
	int index = blockIdx.y * blockDim.x * gridDim.x + blockIdx.x * blockDim.x + threadIdx.x;

	//sPh_dev[index] = sPh_dev[index] - sPhR_dev[index];
	sAmp_dev[index] = sAmp_dev[index] / sAmpR_dev[index];
}

int takeReference(float* sAmp_dev, float* sPh_dev, float* sAmpR_dev, float* sPhR_dev)
{
	dim3 nThreads(D, 1, 1);
	dim3 nBlocks(D, n_proj, 1);

	kernel_takeReference << <nBlocks, nThreads >> > (sAmp_dev, sPh_dev, sAmpR_dev, sPhR_dev);

	cudaError_t cudaStatus = cudaGetLastError();

	if (cudaStatus != cudaSuccess) {
		return 112;
	}
	return 0;
}


// High level functions
EXPORTED_FUNCTION int preproc_sendHologram(short int* hologram, int _X, int _Y, int _nproj, float _NA, float lambda,
	float cam_pix, float M, float _n_imm, int _do_NNC, float _cosFactor, float _fftWindowsScale)
{
	X = _X;
	Y = _Y;
	n_proj = _nproj;
	NA = _NA;
	n_imm = _n_imm;
	do_NNC = _do_NNC;
	preProcOnGPU = true;

	/// crop holograms - option N/A yet - dimension after crop: Nx_holo i Ny_holo
	Nx_holo = X;
	Ny_holo = Y;

	///several calculations
	//int fftWindowsScale = 1;
	D = 2 * round(1 / (1 - 0.1) * _fftWindowsScale * NA / lambda * Nx_holo * cam_pix / M);
	//D = D - 8;
	Dminus = D - 8;
	//int Dmax = 2 * round(1 / (1 - 0.1)*fftWindowsScale*NA / 0.8 * 2048*cam_pix / M);
	float downsampling = ((float)Nx_holo) / D;
	dx = cam_pix * downsampling / M;
	//#define timingsPreproc
#ifdef timingsPreproc
	cudaDeviceProp prop;

	cudaGetDeviceProperties(&prop, 0);
	char nazwaGPU[80];
	strcpy(nazwaGPU, prop.name);

	LARGE_INTEGER countPerSec, tim1, tim2, tim3, tim4, tim5, tim6, tim7, tim8, tim9, tim10, tim11, tim12, tim13;
	QueryPerformanceFrequency(&countPerSec);
	QueryPerformanceCounter(&tim1);
#endif 


	/// memory alloc and FFT config
	int err = allocBuffersForPreproc(hologram);
	if (err != 0)
		return err;

#ifdef timingsPreproc
	QueryPerformanceCounter(&tim2);
#endif

#pragma region FFTplanHolo
	int NN[2] = { Nx_holo,Ny_holo };

	int idist = Nx_holo * Ny_holo;
	int odist = Nx_holo * Ny_holo;
	int inembed[] = { Nx_holo, Ny_holo };
	int onembed[] = { Nx_holo, Ny_holo };
	int istride = 1;
	int ostride = 1;

	cufftResult ee = cufftPlanMany(&planFFTholo, 2, NN, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, n_proj);
	if (ee != cudaSuccess)
		return 66; // FFT plan creation error
#pragma endregion

#pragma region FFTplanComplexField
	int NN_D[2] = { D,D };

	int idist_D = D * D;
	int odist_D = D * D;
	int inembed_D[] = { D, D };
	int onembed_D[] = { D, D };
	int istride_D = 1;
	int ostride_D = 1;

	ee = cufftPlanMany(&planFFTcomplexField, 2, NN_D, inembed_D, istride_D, idist_D, onembed_D, ostride_D, odist_D, CUFFT_C2C, n_proj);
	if (ee != cudaSuccess)
		return 78; // FFT plan creation error
#pragma endregion

#ifdef timingsPreproc
	QueryPerformanceCounter(&tim3);
#endif

	// mask W
	float* tukey2D_W = (float*)malloc(Nx_holo * Ny_holo * sizeof(float));
	tukeyWindow2D(tukey2D_W, Nx_holo, 0.25, Ny_holo, 0.25); // it was 0.25 in MATLAB
	cudaError_t cudaStatus1 = cudaMemcpy(tukeyWin2D_W_dev, tukey2D_W, Nx_holo * Ny_holo * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus1 != cudaSuccess)
		return 71;

	//transfer from short int input buffer to CSG buffer for FFT
	err = transferHologramsToCSGbufferWithMaskW();
	if (err != 0)
		return err;

#ifdef timingsPreproc
	QueryPerformanceCounter(&tim4);
#endif

	ifftShift_1(Nx_holo, Ny_holo, n_proj, hologramsCSG_device);
	ifftShift_2(Nx_holo, Ny_holo, n_proj, hologramsCSG_device);

	err = cufftExecC2C(planFFTholo, hologramsCSG_device, hologramsCSG_device, CUFFT_FORWARD);
	if (err != 0)
		return err;

	ifftShift_1(Nx_holo, Ny_holo, n_proj, hologramsCSG_device);
	ifftShift_2(Nx_holo, Ny_holo, n_proj, hologramsCSG_device);

	err = findMaxIn0incidence();
	if (err != 0)
		return err;

#ifdef timingsPreproc
	QueryPerformanceCounter(&tim5);
#endif

	float* sMask = (float*)malloc(D * D * sizeof(float));
	circTukey2D(sMask, D, 0.1f);
	cudaStatus1 = cudaMemcpy(sMask_dev, sMask, D * D * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus1 != cudaSuccess)
		return 71;

	cudaStatus1 = cudaMemcpy(sMask_dev, sMask, D * D * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus1 != cudaSuccess)
		return 71;

	err = transferToCroppedHoloCSGbufferWith_sMask();
	if (err != 0)
		return err;

#ifdef timingsPreproc
	QueryPerformanceCounter(&tim6);
#endif

	err = findMaxInAllFrames();
	if (err != 0)
		return err;

	//int* temp = (int*)malloc(2*n_proj *sizeof(int));
	//cudaError_t cudaStatus19 = cudaMemcpy(temp, peaksStep2_dev, 2*n_proj * sizeof(int), cudaMemcpyDeviceToHost);
	//FILE* plik8;
	//plik8 = fopen("peaksStep2_ND.bin", "wb");
	//fwrite(temp, sizeof(int), 2*n_proj, plik8);
	//fclose(plik8);
	//free(temp);


#ifdef timingsPreproc
	QueryPerformanceCounter(&tim7);
#endif

	err = circshift_v2();
	if (err != 0)
		return err;

	//float* temp2 = (float*)malloc(2*D*D*n_proj *sizeof(float));
	//cudaError_t cudaStatus11 = cudaMemcpy(temp2, complexFieldCSG_device, 2 * D*D*n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	//FILE* plik9;
	//plik9 = fopen("afterComcirc_ND.bin", "wb");
	//fwrite(temp2, sizeof(float), 2 * D*D*n_proj , plik9);
	//fclose(plik9);
	//free(temp2);

#ifdef timingsPreproc
	QueryPerformanceCounter(&tim8);
#endif

	ifftShift_1(D, D, n_proj, complexFieldCSG_device);
	ifftShift_2(D, D, n_proj, complexFieldCSG_device);

	err = cufftExecC2C(planFFTcomplexField, complexFieldCSG_device, complexFieldCSG_device, CUFFT_INVERSE);
	if (err != 0)
		return err;

	err = normalize_3Dmatrix(complexFieldCSG_device, (float)D * D);
	if (err != 0)
		return err;

	ifftShift_1(D, D, n_proj, complexFieldCSG_device);
	ifftShift_2(D, D, n_proj, complexFieldCSG_device);

#ifdef timingsPreproc
	QueryPerformanceCounter(&tim9);
#endif
	err = fastPhaseUnwrappingWithVolkov2(complexFieldCSG_device, sinoAmp_dev, sinoPh_dev);
	if (err != 0)
		return err;

#ifdef timingsPreproc
	QueryPerformanceCounter(&tim10);
#endif
	if (referenceCalculated) {
		err = takeReference(sinoAmp_dev, sinoPh_dev, sinoAmpRef_dev, sinoPhRef_dev);
		if (err != 0)
			return err;
	}


	/*float* temp = (float*)malloc(n_proj *D * D * sizeof(float));
	cudaError_t cudaStatus19 = cudaMemcpy(temp, sinoAmp_dev, n_proj * D * D * sizeof(float), cudaMemcpyDeviceToHost);
	FILE* plik8;
	plik8 = fopen("SinoAmpPrzedDetrend.bin", "wb");
	fwrite(temp, sizeof(float), D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);

	temp = (float*)malloc(n_proj *D * D * sizeof(float));
	cudaStatus19 = cudaMemcpy(temp, sinoPh_dev, n_proj * D * D * sizeof(float), cudaMemcpyDeviceToHost);
	plik8 = fopen("SinoPhPrzedDetrend.bin", "wb");
	fwrite(temp, sizeof(float), D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

	err = detrend2D(sinoPh_dev);
	if (err != 0)
		return err;



#ifdef timingsPreproc
	QueryPerformanceCounter(&tim11);
#endif

	//mask
	float* mask = (float*)malloc(D * D * sizeof(float));
	tukeyWindow2D(mask, D, _cosFactor, D, _cosFactor);
	cudaStatus1 = cudaMemcpy(mask_dev, mask, D * D * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus1 != cudaSuccess)
		return 102;

	//FILE* plik8;
	//plik8 = fopen("mask.bin", "wb");
	//fwrite(mask, sizeof(float), D*D, plik8);
	//fclose(plik8);

	err = amplAndPhaseProc(sinoAmp_dev, sinoPh_dev, mask_dev);
	if (err != 0)
		return err;

#ifdef timingsPreproc
	QueryPerformanceCounter(&tim12);
#endif

	err = ray_XY(lambda, cam_pix, M);
	if (err != 0)
		return err;

#ifdef timingsPreproc
	QueryPerformanceCounter(&tim13);
#endif

	/*float* temp = (float*)malloc(n_proj *sizeof(float));
	cudaError_t cudaStatus19 = cudaMemcpy(temp, kyp_dev, n_proj * sizeof(float), cudaMemcpyDeviceToHost);
	FILE* plik8;
	plik8 = fopen("VolkovExpanded_phase.bin", "wb");
	fwrite(temp, sizeof(float), 4*D*D*n_proj, plik8);
	fclose(plik8);
	free(temp);*/

#ifdef timingsPreproc
	struct tm newtime;
	time_t now = time(0);
	errno_t eeee = localtime_s(&newtime, &now);
	char buffer[256];
	strftime(buffer, 80, "%c", &newtime); //"%x - %I:%M%p"
	FILE* plik = NULL;
	fopen_s(&plik, "timingsGPUpreproc.txt", "at");

	if (plik != NULL) {

		char* hostname;
		hostname = getenv("COMPUTERNAME");

		double time1 = (double)(tim2.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;
		double time2 = (double)(tim3.QuadPart - tim2.QuadPart) / countPerSec.QuadPart * 1000;
		double time3 = (double)(tim4.QuadPart - tim3.QuadPart) / countPerSec.QuadPart * 1000;
		double time4 = (double)(tim5.QuadPart - tim4.QuadPart) / countPerSec.QuadPart * 1000;
		double time5 = (double)(tim6.QuadPart - tim5.QuadPart) / countPerSec.QuadPart * 1000;
		double time6 = (double)(tim7.QuadPart - tim6.QuadPart) / countPerSec.QuadPart * 1000;
		double time7 = (double)(tim8.QuadPart - tim7.QuadPart) / countPerSec.QuadPart * 1000;
		double time8 = (double)(tim9.QuadPart - tim8.QuadPart) / countPerSec.QuadPart * 1000;
		double time9 = (double)(tim10.QuadPart - tim9.QuadPart) / countPerSec.QuadPart * 1000;
		double time10 = (double)(tim11.QuadPart - tim10.QuadPart) / countPerSec.QuadPart * 1000;
		double time11 = (double)(tim12.QuadPart - tim11.QuadPart) / countPerSec.QuadPart * 1000;
		double time12 = (double)(tim13.QuadPart - tim12.QuadPart) / countPerSec.QuadPart * 1000;

		double timeTot = (double)(tim13.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;

		fprintf(plik, "\n\nAnaliza wykonana w dniu:\n %s\n", buffer);
		fprintf(plik, "na karcie graficznej: %s\n", nazwaGPU);
		fprintf(plik, "komputer: %s\n\n", hostname);

		fprintf(plik, "\ntotal preprocessing time: \t\t%10.2f", timeTot);
		fprintf(plik, "\n  - time 1:    \t%10.2f", time1);
		fprintf(plik, "\n  - time 2:	\t%10.2f", time2);
		fprintf(plik, "\n  - time 3:    \t%10.2f", time3);
		fprintf(plik, "\n  - time 4:    \t%10.2f", time4);
		fprintf(plik, "\n  - time 5:    \t%10.2f", time5);
		fprintf(plik, "\n  - time 6:	\t%10.2f", time6);
		fprintf(plik, "\n  - time 7:    \t%10.2f", time7);
		fprintf(plik, "\n  - time 8:    \t%10.2f", time8);
		fprintf(plik, "\n  - time 9:    \t%10.2f", time9);
		fprintf(plik, "\n  - time 10:	\t%10.2f", time10);
		fprintf(plik, "\n  - time 11:    \t%10.2f", time11);
		fprintf(plik, "\n  - time 12:    \t%10.2f", time12);

		fclose(plik);
	}
#endif												
	//free
	cudaFree(holograms_device);     holograms_device = NULL;
	cudaFree(hologramsCSG_device); hologramsCSG_device = NULL;
	cudaFree(tukeyWin2D_W_dev); tukeyWin2D_W_dev = NULL;
	cudaFree(croppedPartsHoloCSG_device); croppedPartsHoloCSG_device = NULL;
	cudaFree(sMask_dev); sMask_dev = NULL;
	cudaFree(mask_dev); mask_dev = NULL;
	cudaFree(peaksStep1_dev); peaksStep1_dev = NULL;
	cudaFree(peaksStep2_dev); peaksStep2_dev = NULL;
	cudaFree(complexFieldCSG_device); complexFieldCSG_device = NULL;

	free(tukey2D_W);
	free(sMask);
	free(mask);
	cufftDestroy(planFFTholo); planFFTholo = NULL;
	cufftDestroy(planFFTcomplexField); planFFTcomplexField = NULL;

	return 0;
}

EXPORTED_FUNCTION int HL_removeReference()
{
	referenceCalculated = false;
	cudaFree(sinoAmpRef_dev); sinoAmpRef_dev = NULL;
	cudaFree(sinoPhRef_dev); sinoPhRef_dev = NULL;
	cudaFree(sinoPhRefComplex_dev); sinoPhRefComplex_dev = NULL;
	return 0;
}

EXPORTED_FUNCTION int HL00_removeReference()
{
	return HL_removeReference();
}

EXPORTED_FUNCTION int HL_addReference(short int* hologram, int _X, int _Y, int _nproj, float _NA, float lambda,
	float cam_pix, float M, float _n_imm, int _do_NNC, float _fftWindowsScale)
{
	X = _X;
	Y = _Y;
	n_proj = _nproj;
	NA = _NA;
	n_imm = _n_imm;
	do_NNC = _do_NNC;
	referenceCalculated = true;

	//FILE* plik8;
	//plik8 = fopen("referencja_ODT_2021-9-9_12-57-22.bin", "wb");
	//fwrite(hologram, sizeof(short int), X * Y * n_proj, plik8);
	//fclose(plik8);
	cudaInitDev();

	/// crop holograms - option N/A yet - dimension after crop: Nx_holo i Ny_holo
	Nx_holo = X;
	Ny_holo = Y;

	///several calculations
	//int fftWindowsScale = 1;
	D = 2 * round(1 / (1 - 0.1) * _fftWindowsScale * NA / lambda * Nx_holo * cam_pix / M);

	Dminus = D - 8;
	//int Dmax = 2 * round(1 / (1 - 0.1)*fftWindowsScale*NA / 0.8 * 2048*cam_pix / M);
	float downsampling = ((float)Nx_holo) / D;
	dx = cam_pix * downsampling / M;



	/// memory alloc and FFT config
	int err = allocBuffersForPreproc_Reference(hologram);
	if (err != 0)
		return err;

#pragma region FFTplanHolo
	int NN[2] = { Nx_holo,Ny_holo };

	int idist = Nx_holo * Ny_holo;
	int odist = Nx_holo * Ny_holo;
	int inembed[] = { Nx_holo, Ny_holo };
	int onembed[] = { Nx_holo, Ny_holo };
	int istride = 1;
	int ostride = 1;

	cufftResult ee = cufftPlanMany(&planFFTholo, 2, NN, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, n_proj);
	if (ee != cudaSuccess)
		return 66; // FFT plan creation error
#pragma endregion

#pragma region FFTplanComplexField
	int NN_D[2] = { D,D };

	int idist_D = D * D;
	int odist_D = D * D;
	int inembed_D[] = { D, D };
	int onembed_D[] = { D, D };
	int istride_D = 1;
	int ostride_D = 1;

	ee = cufftPlanMany(&planFFTcomplexField, 2, NN_D, inembed_D, istride_D, idist_D, onembed_D, ostride_D, odist_D, CUFFT_C2C, n_proj);
	if (ee != cudaSuccess)
		return 78; // FFT plan creation error
#pragma endregion


	// mask W
	float* tukey2D_W = (float*)malloc(Nx_holo * Ny_holo * sizeof(float));
	tukeyWindow2D(tukey2D_W, Nx_holo, 0.25, Ny_holo, 0.25); // it was 0.25 in MATLAB
	cudaError_t cudaStatus1 = cudaMemcpy(tukeyWin2D_W_dev, tukey2D_W, Nx_holo * Ny_holo * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus1 != cudaSuccess)
		return 71;

	//transfer from short int input buffer to CSG buffer for FFT
	err = transferHologramsToCSGbufferWithMaskW();
	if (err != 0)
		return err;


	ifftShift_1(Nx_holo, Ny_holo, n_proj, hologramsCSG_device);
	ifftShift_2(Nx_holo, Ny_holo, n_proj, hologramsCSG_device);

	err = cufftExecC2C(planFFTholo, hologramsCSG_device, hologramsCSG_device, CUFFT_FORWARD);
	if (err != 0)
		return err;

	ifftShift_1(Nx_holo, Ny_holo, n_proj, hologramsCSG_device);
	ifftShift_2(Nx_holo, Ny_holo, n_proj, hologramsCSG_device);

	err = findMaxIn0incidence();
	if (err != 0)
		return err;


	float* sMask = (float*)malloc(D * D * sizeof(float));
	circTukey2D(sMask, D, 0.1f);
	cudaStatus1 = cudaMemcpy(sMask_dev, sMask, D * D * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus1 != cudaSuccess)
		return 71;

	cudaStatus1 = cudaMemcpy(sMask_dev, sMask, D * D * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus1 != cudaSuccess)
		return 71;

	err = transferToCroppedHoloCSGbufferWith_sMask();
	if (err != 0)
		return err;


	err = findMaxInAllFrames();
	if (err != 0)
		return err;


	err = circshift_v2();
	if (err != 0)
		return err;


	ifftShift_1(D, D, n_proj, complexFieldCSG_device);
	ifftShift_2(D, D, n_proj, complexFieldCSG_device);

	err = cufftExecC2C(planFFTcomplexField, complexFieldCSG_device, complexFieldCSG_device, CUFFT_INVERSE);
	if (err != 0)
		return err;

	err = normalize_3Dmatrix(complexFieldCSG_device, (float)D * D);
	if (err != 0)
		return err;

	ifftShift_1(D, D, n_proj, complexFieldCSG_device);
	ifftShift_2(D, D, n_proj, complexFieldCSG_device);

	err = complexFieldToComplexPhaseAndAmp(complexFieldCSG_device, sinoPhRefComplex_dev, sinoAmpRef_dev);
	if (err != 0)
		return err;

	/*err = fastPhaseUnwrappingWithVolkov(complexFieldCSG_device, sinoAmpRef_dev, sinoPhRef_dev);
	if (err != 0)
		return err;*/

		//free
	cudaFree(holograms_device);     holograms_device = NULL;
	cudaFree(hologramsCSG_device); hologramsCSG_device = NULL;
	cudaFree(tukeyWin2D_W_dev); tukeyWin2D_W_dev = NULL;
	cudaFree(croppedPartsHoloCSG_device); croppedPartsHoloCSG_device = NULL;
	cudaFree(sMask_dev); sMask_dev = NULL;
	cudaFree(mask_dev); mask_dev = NULL;
	cudaFree(peaksStep1_dev); peaksStep1_dev = NULL;
	cudaFree(peaksStep2_dev); peaksStep2_dev = NULL;
	cudaFree(complexFieldCSG_device); complexFieldCSG_device = NULL;

	free(tukey2D_W);
	free(sMask);
	cufftDestroy(planFFTholo); planFFTholo = NULL;
	cufftDestroy(planFFTcomplexField); planFFTcomplexField = NULL;

	//xMallocRaport();

	return 0;
}

// Backward-compatibility wrapper for legacy name HL00_addReference
EXPORTED_FUNCTION int HL00_addReference(short int* hologram, int _X, int _Y, int _nproj, float _NA, float lambda,
	float cam_pix, float M, float _n_imm, int _do_NNC, float _fftWindowsScale)
{
	return HL_addReference(hologram, _X, _Y, _nproj, _NA, lambda,
		cam_pix, M, _n_imm, _do_NNC, _fftWindowsScale);
}

EXPORTED_FUNCTION int HL00to02_FromPreprocToGenKO(short int* hologram, int X, int Y, int _nproj, float NA, float lambda,
	float cam_pix, float M, float _n_imm, int _do_NNC, int* _K_xy, int Kspace_oversampling_z, float _cosFactor, float _fftWindowsScale, int _approxBornNotRytov)
{
	LARGE_INTEGER countPerSec, tim1, tim2, tim3, tim4, tim5;
#ifdef Save_timings
	cudaDeviceProp prop;

	cudaGetDeviceProperties(&prop, 0);
	char nazwaGPU[80];
	strcpy(nazwaGPU, prop.name);

	QueryPerformanceFrequency(&countPerSec);
	QueryPerformanceCounter(&tim1);
#endif //  Save_timings
	if (!referenceCalculated)
		cudaInitDev();

	int error = preproc_sendHologram(hologram, X, Y, _nproj, NA, lambda, cam_pix, M, _n_imm, _do_NNC, _cosFactor, _fftWindowsScale);
	if (error) return error;


	//FILE* plik8;
	//plik8 = fopen("hologram_ODT_2021-9-9_12-57-22.bin", "wb");
	//fwrite(hologram, sizeof(short int), X * Y * n_proj, plik8);
	//fclose(plik8);
#ifdef Save_timings
	QueryPerformanceCounter(&tim2);
#endif //  Save_timings

	float kxy_mean;
	cudaError_t cudaStatus5 = cudaMemcpy(&kxy_mean, kxy_dev, 1 * sizeof(float), cudaMemcpyDeviceToHost);

	float OTF_diameter = 2 * (NA / lambda + kxy_mean);

	//dkP = 1 / (dx_projection*N_projection_padded);
	dkP = 1 / (dx * D);
	dkPz = dkP / Kspace_oversampling_z; // simplified for preproc on GPU

	int Kxy = (int)round(OTF_diameter / dkP / 2) * 2;
	*_K_xy = Kxy;

	int Kz;
	Kz = Kspace_oversampling_z * Kxy;

	float* lambda_all = (float*)malloc(n_proj * sizeof(float));
	for (int i = 0; i < n_proj; i++)
		lambda_all[i] = lambda;
#ifdef Save_timings
	QueryPerformanceCounter(&tim3);
#endif //  Save_timings
	error = prepareKspaceGeneration(Kxy, Kz, dx, n_imm, n_proj, D, D,
		dkP, dkPz, NA, lambda_all, kxp_dev, kyp_dev, _approxBornNotRytov);
	if (error) return error;
#ifdef Save_timings
	QueryPerformanceCounter(&tim4);
#endif //  Save_timings
	error = generateKspace();
	if (error) return error;
#ifdef Save_timings
	QueryPerformanceCounter(&tim5);

	struct tm newtime;
	time_t now = time(0);
	errno_t eeee = localtime_s(&newtime, &now);
	char buffer[256];
	strftime(buffer, 80, "%c", &newtime); //"%x - %I:%M%p"
	FILE* plik = NULL;
	fopen_s(&plik, "timingsGPU.txt", "at");

	if (plik != NULL) {

		char* hostname;
		hostname = getenv("COMPUTERNAME");

		double time1 = (double)(tim2.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;
		double time2 = (double)(tim3.QuadPart - tim2.QuadPart) / countPerSec.QuadPart * 1000;
		double time3 = (double)(tim4.QuadPart - tim3.QuadPart) / countPerSec.QuadPart * 1000;
		double time4 = (double)(tim5.QuadPart - tim4.QuadPart) / countPerSec.QuadPart * 1000;
		double timeTot = (double)(tim5.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;

		fprintf(plik, "\n\nAnaliza wykonana w dniu:\n %s\n", buffer);
		fprintf(plik, "na karcie graficznej: %s\n", nazwaGPU);
		fprintf(plik, "komputer: %s\n\n", hostname);

		fprintf(plik, "\nTotal time HL00to02_FromPreprocToGenKO: \t\t%10.2f", timeTot);
		fprintf(plik, "\n  - preprocessing time:     \t%10.2f", time1);
		fprintf(plik, "\n  - minor CPU computations: \t%10.2f", time2);
		fprintf(plik, "\n  - time prepareKspaceGen:  \t%10.2f", time3);
		fprintf(plik, "\n  - time generateKspace:    \t%10.2f", time4);

		fclose(plik);
	}
#endif //  Save_timings
	free(lambda_all);


	return 0;
}

// after preprocessing the following are available: dx, D, Nx, Ny
EXPORTED_FUNCTION int HL01_setParams(int _K_xy, int _K_z, float _dx, float _n_imm, int _n_proj, int _Nx, int _Ny,
	float _dkP, float _dkPz, float _NA, float* _lambda_all, float* _kxp, float* _kyp, int _approxBornNotRytov)
{
	preProcOnGPU = false;
	LARGE_INTEGER countPerSec, tim1, tim2, tim3;
#ifdef Save_timings
	cudaDeviceProp prop;

	cudaGetDeviceProperties(&prop, 0);
	char nazwaGPU[80];
	strcpy(nazwaGPU, prop.name);

	QueryPerformanceFrequency(&countPerSec);
	QueryPerformanceCounter(&tim1);
#endif //  Save_timings
	int error = cudaInitDev();
	if (error) return error;
#ifdef Save_timings
	QueryPerformanceCounter(&tim2);
#endif //  Save_timings

	error = prepareKspaceGeneration(_K_xy, _K_z, _dx, _n_imm, _n_proj, _Nx, _Ny,
		_dkP, _dkPz, _NA, _lambda_all, _kxp, _kyp, _approxBornNotRytov);
	if (error) return error;
#ifdef Save_timings
	QueryPerformanceCounter(&tim3);

	struct tm newtime;
	time_t now = time(0);
	errno_t eeee = localtime_s(&newtime, &now);
	char buffer[256];
	strftime(buffer, 80, "%c", &newtime); //"%x - %I:%M%p"
	FILE* plik = NULL;
	fopen_s(&plik, "timingsGPU.txt", "at");

	if (plik != NULL) {

		char* hostname;
		hostname = getenv("COMPUTERNAME");

		double time1 = (double)(tim2.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;
		double time2 = (double)(tim3.QuadPart - tim2.QuadPart) / countPerSec.QuadPart * 1000;
		double timeTot = (double)(tim3.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;

		fprintf(plik, "\n\nAnaliza wykonana w dniu:\n %s\n", buffer);
		fprintf(plik, "na karcie graficznej: %s\n", nazwaGPU);
		fprintf(plik, "komputer: %s\n\n", hostname);

		fprintf(plik, "\nTotal time HL_01_setParams: \t\t%10.2f", timeTot);
		fprintf(plik, "\n  - GPU initialization time: \t%10.2f", time1);
		fprintf(plik, "\n  - time prepareKspaceGen:   \t%10.2f", time2);

		fclose(plik);
	}
#endif //  Save_timings

	return 0;
}

EXPORTED_FUNCTION int HL02_sendDataAndGenerateKO(float* _sinoAmp, float* _sinoPh, unsigned char* _FpmaskLogical)
{
	LARGE_INTEGER countPerSec, tim1, tim2, tim3, tim4;
#ifdef Save_timings
	QueryPerformanceFrequency(&countPerSec);
	QueryPerformanceCounter(&tim1);
#endif //  Save_timings


	int error = sendSinograms(_sinoAmp, _sinoPh);
	if (error) return error;
#ifdef Save_timings
	QueryPerformanceCounter(&tim2);
#endif //  Save_timings
	error = sendFpmask_logical(_FpmaskLogical);
	if (error) return error;
#ifdef Save_timings
	QueryPerformanceCounter(&tim3);
#endif //  Save_timings
	error = generateKspace();
	if (error) return error;
#ifdef Save_timings
	QueryPerformanceCounter(&tim4);


	FILE* plik = NULL;
	fopen_s(&plik, "timingsGPU.txt", "at");

	if (plik != NULL) {

		double time1 = (double)(tim2.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;
		double time2 = (double)(tim3.QuadPart - tim2.QuadPart) / countPerSec.QuadPart * 1000;
		double time3 = (double)(tim4.QuadPart - tim3.QuadPart) / countPerSec.QuadPart * 1000;
		double timeTot = (double)(tim4.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;

		fprintf(plik, "\n\nTotal time HL02_sendDataAndGenKO: \t%10.2f", timeTot);
		fprintf(plik, "\n  - time sendSinograms:                \t%10.2f", time1);
		fprintf(plik, "\n  - time sendFpmask_logical:           \t%10.2f", time2);
		fprintf(plik, "\n  - time generateKspace:               \t%10.2f", time3);

		fprintf(plik, "\n\n_______________________________________");
		fprintf(plik, "\nInfo o pliku:");
		fprintf(plik, "\n - number of projections:    %d", n_proj);
		fprintf(plik, "\n - projection size:      %d x %d", Nx, Ny);
		fprintf(plik, "\n_______________________________________");

		fclose(plik);
	}
#endif //  Save_timings
	return 0;
}
EXPORTED_FUNCTION int HL02_B_optionTakeKO(float* _complexData_host)
{
	int error = takeDataKspace(_complexData_host);
	if (error) return error;

	return 0;
}

EXPORTED_FUNCTION int HL02_B_optionTakeEW(int* _ew_host)
{
	int error = takeDataEW(_ew_host);
	if (error) return error;

	return 0;
}

EXPORTED_FUNCTION int HL02_C_optionFreeGPUmemory()
{
	freeGPUmemory();

	return 0;
}

EXPORTED_FUNCTION int takeUnwrappedFrame(float* unwrapped_host)
{
	cudaError_t cudaStatus1 = cudaMemcpy(unwrapped_host, sinoPh_dev, D * D * n_proj * sizeof(float), cudaMemcpyDeviceToHost);

	if (cudaStatus1 != cudaSuccess) {
		return 56;
	}

	return 0;
}

EXPORTED_FUNCTION int HL00_B_optionTakeSinoAmpRef(float* sinAmpRef_host)
{
	cudaError_t cudaStatus1 = cudaMemcpy(sinAmpRef_host, sinoAmpRef_dev, D * D * n_proj * sizeof(float), cudaMemcpyDeviceToHost);

	if (cudaStatus1 != cudaSuccess) {
		return 500;
	}

	return 0;
}

EXPORTED_FUNCTION int HL00_B_optionTakeSinoPhRef(float* sinPhRef_host)
{
	cudaError_t cudaStatus1 = cudaMemcpy(sinPhRef_host, sinoPhRef_dev, D * D * n_proj * sizeof(float), cudaMemcpyDeviceToHost);

	if (cudaStatus1 != cudaSuccess) {
		return 500;
	}

	return 0;
}