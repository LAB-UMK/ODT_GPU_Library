/**
 * @file Reconstruction.cu
 * @brief GPU kernels and exported functions for tomographic reconstruction.
 *
 * Implements the reconstruction stage of the LaODT processing pipeline:
 * Direct Inversion (DI) and the iterative Gerchberg-Papoulis (GP) algorithm,
 * together with the associated CUDA kernels, GPU memory management and
 * 3D cuFFT transforms.
 *
 * Part of the ODT_GPU library.
 * Distributed under the GNU General Public License v3.0 (GPL-3.0).
 * See the LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html for details.
 *
 * @see https://github.com/LAB-UMK/ODT_GPU_Library
 */

// Reconstruction.cu : Defines the exported functions for the DLL application.
//
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
#include <stdlib.h>
#include <string.h> 

#define NRANK 3

extern int K_xy;
extern int K_z;
extern float dx;
extern int Nx;
extern bool preProcOnGPU;

int do_TC = 1;
int do_NNC = 1;
float n_imm = 1.335f;
float relaxM = 1;
float betaM = 1;
float kn_mean = 0;
int objShift = 0;
float dxo = 0;
int n_EW8bit = 0;
int relaxGP = 1;
int betaGP = 1;
int nGPi = 30;
bool objSupportAv = false;  /// object support available - parameter changed in "sendData_objectSupport" function,
/// if user is sending empty matrix, objSupportAv=false and object suport procedure is not called in iterations

cufftHandle plan;
int n[NRANK];

cufftComplex* n_rec_dev = NULL;
cufftComplex* KOi_dev = NULL;
extern cufftComplex* KO_dev;
extern int* EW_dev;
float* objSupport_dev = NULL;
cufftComplex* EKO_dev = NULL;
int* EKOind_dev = NULL;

extern void cleanDeviceMemoryBeforeDirectInversion();

/////////////////////
/// Kernel Functions
/////////////////////

// Transparency Constraint
__global__ void kernel_TC(cufftComplex* n_rec_dev)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	n_rec_dev[index].y = 0;
}

int transparencyConstraint()
{
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);
	// execute the kernel
	kernel_TC << <nBlocks, nThreads >> > (n_rec_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 7;
	}
	return 0;
}


// Non-negativity constraint
__global__ void kernel_NNC(cufftComplex* n_rec_dev, float nImm)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (n_rec_dev[index].x < nImm)
		n_rec_dev[index].x = nImm;
}

int nonNegativityConstraint()
{
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);
	// execute the kernel
	kernel_NNC << <nBlocks, nThreads >> > (n_rec_dev, n_imm);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 8;
	}
	return 0;
}


// Non-Positivity constraint
__global__ void kernel_NPC(cufftComplex* n_rec_dev, float nImm)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (n_rec_dev[index].x > nImm)
		n_rec_dev[index].x = nImm;
}

int nonPositivityConstraint()
{
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);
	// execute the kernel
	kernel_NPC << <nBlocks, nThreads >> > (n_rec_dev, n_imm);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 9;
	}
	return 0;

}
// transparency and Non-negativity constraint

__global__ void kernel_TCandNNC(cufftComplex* n_rec_dev, float nImm)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (n_rec_dev[index].x < nImm)
		n_rec_dev[index].x = nImm;

	n_rec_dev[index].y = 0;
}

int transparencyAndNonNegativityConstraint()
{
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);
	// execute the kernel
	kernel_TCandNNC << <nBlocks, nThreads >> > (n_rec_dev, n_imm);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 8;
	}
	return 0;
}


// transparency and Non-positivity constraint

__global__ void kernel_TCandNPC(cufftComplex* n_rec_dev, float nImm)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (n_rec_dev[index].x > nImm)
		n_rec_dev[index].x = nImm;

	n_rec_dev[index].y = 0;
}

int transparencyAndNonPositivityConstraint()
{
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);
	// execute the kernel
	kernel_TCandNPC << <nBlocks, nThreads >> > (n_rec_dev, n_imm);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 8;
	}
	return 0;
}


// transparency and Non-negativity constraint and Rytov

__global__ void kernel_TCandNNCandRytov(cufftComplex* n_rec_dev, float nImm, float param1, float param2)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	float re = n_rec_dev[index].x;
	float im = n_rec_dev[index].y;

	if (re < nImm)
		re = nImm;

	im = 0;

	float re2 = re * re - im * im;
	float im2 = 2 * re * im;

	n_rec_dev[index].x = param1 * (1 - re2 / param2);
	n_rec_dev[index].y = param1 * (-im2 / param2);
}

int transparencyAndNonNegativityConstraintandRytov()
{
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);

	float param1 = (2 * PI * kn_mean) * (2 * PI * kn_mean);
	float param2 = n_imm * n_imm;

	// execute the kernel
	kernel_TCandNNCandRytov << <nBlocks, nThreads >> > (n_rec_dev, n_imm, param1, param2);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 8;
	}
	return 0;
}


// transparency and Non-positivity constraint and Rytov

__global__ void kernel_TCandNPCandRytov(cufftComplex* n_rec_dev, float nImm, float param1, float param2)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	float re = n_rec_dev[index].x;
	float im = n_rec_dev[index].y;

	if (re > nImm)
		re = nImm;

	im = 0;

	float re2 = re * re - im * im;
	float im2 = 2 * re * im;

	n_rec_dev[index].x = param1 * (1 - re2 / param2);
	n_rec_dev[index].y = param1 * (-im2 / param2);
}

int transparencyAndNonPositivityConstraintAndRytov()
{
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);

	float param1 = (2 * PI * kn_mean) * (2 * PI * kn_mean);
	float param2 = n_imm * n_imm;

	// execute the kernel
	kernel_TCandNPCandRytov << <nBlocks, nThreads >> > (n_rec_dev, n_imm, param1, param2);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 8;
	}
	return 0;
}


// spatial support
__global__ void kernel_spatialSupport(cufftComplex* n_rec_dev, float* objSupport_dev, float n_imm, float betaM)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	n_rec_dev[index].x = (n_rec_dev[index].x - n_imm) * (objSupport_dev[index] * betaM + (1 - betaM)) + n_imm;
	n_rec_dev[index].y = n_rec_dev[index].y * (objSupport_dev[index] * betaM + (1 - betaM));
}

__global__ void kernel_spatialSupport_forBetaM_1(cufftComplex* n_rec_dev, float* objSupport_dev, float n_imm, float betaM)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	n_rec_dev[index].x = (n_rec_dev[index].x - n_imm) * objSupport_dev[index] + n_imm;
	n_rec_dev[index].y = n_rec_dev[index].y * objSupport_dev[index];
}

int spatialSupport()
{
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);
	// execute the kernel
	if ((betaM > 0.99999) && (betaM < 1.00001))
		kernel_spatialSupport_forBetaM_1 << <nBlocks, nThreads >> > (n_rec_dev, objSupport_dev, n_imm, betaM);
	else
		kernel_spatialSupport << <nBlocks, nThreads >> > (n_rec_dev, objSupport_dev, n_imm, betaM);

	betaM = relaxM * betaM;

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 9;
	}
	return 0;
}


__global__ void kernel_approxRytov(cufftComplex* n_rec_dev, float param1, float param2)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	float re = n_rec_dev[index].x;
	float im = n_rec_dev[index].y;

	float re2 = re * re - im * im;
	float im2 = 2 * re * im;

	n_rec_dev[index].x = param1 * (1 - re2 / param2);
	n_rec_dev[index].y = param1 * (-im2 / param2);
}

int approxRytov()
{
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);

	float param1 = (2 * PI * kn_mean) * (2 * PI * kn_mean);
	float param2 = n_imm * n_imm;

	kernel_approxRytov << <nBlocks, nThreads >> > (n_rec_dev, param1, param2);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 10;
	}
	return 0;
}

__global__ void kernel_ifftShift_3(cufftComplex* data_copmlex_dev)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	register int offset = blockDim.x * gridDim.x;

	float re = data_copmlex_dev[index].x;
	float im = data_copmlex_dev[index].y;

	data_copmlex_dev[index].x = data_copmlex_dev[index + offset].x;
	data_copmlex_dev[index].y = data_copmlex_dev[index + offset].y;

	data_copmlex_dev[index + offset].x = re;
	data_copmlex_dev[index + offset].y = im;
}

int ifftShift_3(int x, int y, int z, cufftComplex* data_copmlex_dev)
{
	int nThreads = x;
	dim3 nBlocks(y * z / 2);

	kernel_ifftShift_3 << <nBlocks, nThreads >> > (data_copmlex_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 11;
	}
	return 0;
}

__global__ void kernel_ifftShift_1(cufftComplex* data_complex_dev)
{
	register int index = blockIdx.y * 2 * blockDim.x * gridDim.x + blockIdx.x * 2 * blockDim.x + threadIdx.x;

	register int offset = blockDim.x;

	float re = data_complex_dev[index].x;
	float im = data_complex_dev[index].y;

	data_complex_dev[index].x = data_complex_dev[index + offset].x;
	data_complex_dev[index].y = data_complex_dev[index + offset].y;

	data_complex_dev[index + offset].x = re;
	data_complex_dev[index + offset].y = im;
}

int ifftShift_1(int x, int y, int z, cufftComplex* data_complex_dev)
{
	int nThreads = x / 2;
	dim3 nBlocks(x, z, 1);

	kernel_ifftShift_1 << <nBlocks, nThreads >> > (data_complex_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 11;
	}
	return 0;
}

__global__ void kernel_ifftShift_2(cufftComplex* data_complex_dev)
{
	register int index = blockIdx.y * blockDim.x * 2 * gridDim.x + blockIdx.x * blockDim.x + threadIdx.x;

	register int offset = blockDim.x * gridDim.x;

	float re = data_complex_dev[index].x;
	float im = data_complex_dev[index].y;

	data_complex_dev[index].x = data_complex_dev[index + offset].x;
	data_complex_dev[index].y = data_complex_dev[index + offset].y;

	data_complex_dev[index + offset].x = re;
	data_complex_dev[index + offset].y = im;
}

__global__ void kernel_ifftShift_2_big(cufftComplex* data_complex_dev)
{
	register int index = blockIdx.y * blockDim.x * 2 * gridDim.x + blockIdx.x * blockDim.x + threadIdx.x;

	register int offset = 2 * blockDim.x * gridDim.x;

	float re = data_complex_dev[2 * index].x;
	float im = data_complex_dev[2 * index].y;
	float re2 = data_complex_dev[2 * index + 1].x;
	float im2 = data_complex_dev[2 * index + 1].y;

	data_complex_dev[2 * index].x = data_complex_dev[2 * index + offset].x;
	data_complex_dev[2 * index].y = data_complex_dev[2 * index + offset].y;
	data_complex_dev[2 * index + 1].x = data_complex_dev[2 * index + 1 + offset].x;
	data_complex_dev[2 * index + 1].y = data_complex_dev[2 * index + 1 + offset].y;

	data_complex_dev[2 * index + offset].x = re;
	data_complex_dev[2 * index + offset].y = im;
	data_complex_dev[2 * index + 1 + offset].x = re2;
	data_complex_dev[2 * index + 1 + offset].y = im2;
}

int ifftShift_2(int x, int y, int z, cufftComplex* data_complex_dev)
{
	int nThreads = x;
	dim3 nBlocks(x / 2, z, 1);

	if (nThreads < 1025)
		kernel_ifftShift_2 << <nBlocks, nThreads >> > (data_complex_dev);
	else {
		nThreads = nThreads / 2;
		kernel_ifftShift_2_big << <nBlocks, nThreads >> > (data_complex_dev);
	}


	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 11;
	}
	return 0;
}

#pragma region fftshiftInt
__global__ void kernel_ifftShift_3_int(int* data_dev)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	register int offset = blockDim.x * gridDim.x;

	int data = data_dev[index];

	data_dev[index] = data_dev[index + offset];

	data_dev[index + offset] = data;
}

int ifftShift_3_int(int x, int y, int z, int* data_dev)
{
	int nThreads = x;
	dim3 nBlocks(y * z / 2);

	kernel_ifftShift_3_int << <nBlocks, nThreads >> > (data_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 11;
	}
	return 0;
}

__global__ void kernel_ifftShift_1_int(int* data_dev)
{
	register int index = blockIdx.y * 2 * blockDim.x * gridDim.x + blockIdx.x * 2 * blockDim.x + threadIdx.x;

	register int offset = blockDim.x;

	int data = data_dev[index];

	data_dev[index] = data_dev[index + offset];

	data_dev[index + offset] = data;
}

int ifftShift_1_int(int x, int y, int z, int* data_dev)
{
	int nThreads = x / 2;
	dim3 nBlocks(x, z, 1);

	kernel_ifftShift_1_int << <nBlocks, nThreads >> > (data_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 11;
	}
	return 0;
}

__global__ void kernel_ifftShift_2_int(int* data_dev)
{
	register int index = blockIdx.y * blockDim.x * 2 * gridDim.x + blockIdx.x * blockDim.x + threadIdx.x;

	register int offset = blockDim.x * gridDim.x;

	float data = data_dev[index];

	data_dev[index] = data_dev[index + offset];

	data_dev[index + offset] = data;
}

int ifftShift_2_int(int x, int y, int z, int* data_dev)
{
	int nThreads = x;
	dim3 nBlocks(x / 2, z, 1);

	kernel_ifftShift_2_int << <nBlocks, nThreads >> > (data_dev);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 11;
	}
	return 0;
}
#pragma endregion fftshiftInt

__global__ void kernel_mult(cufftComplex* matrix3D, float factor)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	matrix3D[index].x *= factor;
	matrix3D[index].y *= factor;
}

int mult_3Dmatrix(cufftComplex* matrix3D, float factor)
{
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);

	kernel_mult << <nBlocks, nThreads >> > (matrix3D, factor);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 15;
	}
	return 0;
}

__global__ void kernel_div(cufftComplex* matrix3D, float factor)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	matrix3D[index].x /= factor;
	matrix3D[index].y /= factor;
}

int div_3Dmatrix(cufftComplex* matrix3D, float factor)
{
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);

	kernel_div << <nBlocks, nThreads >> > (matrix3D, factor);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 15;
	}
	return 0;
}

//__global__ void kernel_temp(cufftComplex* temp, int n_EW8bit)
//{
//	register int index = blockIdx.x * blockDim.x + threadIdx.x;
//
//	if (index < n_EW8bit) {
//		temp[index].x = 100;
//		temp[index].y = 100;
//	}
//}

__global__ void kernel_insertProjData(cufftComplex* _KO_dev, cufftComplex* _EKO_dev, int* EKOind_dev, int n_EW8bit)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	register int indexKO;

	if (index < n_EW8bit) {
		indexKO = EKOind_dev[index];

		_KO_dev[indexKO].x = _EKO_dev[index].x;
		_KO_dev[indexKO].y = _EKO_dev[index].y;
	}
}

int insertProjData(cufftComplex* _KOi, cufftComplex* _EKO_dev, int* _EKOind_dev, int _n_EW8bit)
{

	int nThreads = 1024;
	int nBlocks = n_EW8bit / nThreads;
	nBlocks++;

	kernel_insertProjData << <nBlocks, nThreads >> > (_KOi, _EKO_dev, _EKOind_dev, _n_EW8bit);

	betaGP = relaxGP * betaGP;

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 16;
	}
	return 0;
}

__global__ void kernel_approxRytovProj2Rec(cufftComplex* n_rec_dev, float param1, float param2)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	float re = n_rec_dev[index].x;
	float im = n_rec_dev[index].y;

	cufftComplex val;
	val.x = 1 - re / param2;
	val.y = -im / param2;

	//sqrtf(complex val)
	float radius = cuCabsf(val);
	float cosA = val.x / radius;
	cufftComplex out;
	out.x = sqrtf(radius * (cosA + 1.0) / 2.0);
	out.y = sqrtf(radius * (1.0 - cosA) / 2.0);
	// signbit should be false if x.y is negative
	if (signbit(val.y))
		out.y *= -1.0;

	n_rec_dev[index].x = param1 * out.x;
	n_rec_dev[index].y = param1 * out.y;
}

int approxRytovProj2Rec()
{
	//n_rec = n_imm * sqrt(1 - n_rec / (2 * pi*mean(kn)) ^ 2); % o(x, y, z)->n(x, y, z)
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);

	float param1 = n_imm;
	float param2 = (2 * PI * kn_mean) * (2 * PI * kn_mean);

	kernel_approxRytovProj2Rec << <nBlocks, nThreads >> > (n_rec_dev, param1, param2);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 10;
	}
	return 0;
}

__global__ void kernel_approxRytovProj2RecWithDiv(cufftComplex* n_rec_dev, float param1, float param2, float factor)
{
	register int index = blockIdx.x * blockDim.x + threadIdx.x;

	float re = n_rec_dev[index].x / factor;
	float im = n_rec_dev[index].y / factor;

	cufftComplex val;
	val.x = 1 - re / param2;
	val.y = -im / param2;

	//sqrtf(complex val)
	float radius = cuCabsf(val);
	float cosA = val.x / radius;
	cufftComplex out;
	out.x = sqrtf(radius * (cosA + 1.0) / 2.0);
	out.y = sqrtf(radius * (1.0 - cosA) / 2.0);
	// signbit should be false if x.y is negative
	if (signbit(val.y))
		out.y *= -1.0;

	n_rec_dev[index].x = param1 * out.x;
	n_rec_dev[index].y = param1 * out.y;
}

int approxRytovProj2RecWithDiv(float factor)
{
	//n_rec = n_imm * sqrt(1 - n_rec / (2 * pi*mean(kn)) ^ 2); % o(x, y, z)->n(x, y, z)
	int nThreads = K_xy;
	dim3 nBlocks(K_xy * K_z);

	float param1 = n_imm;
	float param2 = (2 * PI * kn_mean) * (2 * PI * kn_mean);

	kernel_approxRytovProj2RecWithDiv << <nBlocks, nThreads >> > (n_rec_dev, param1, param2, factor);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 10;
	}
	return 0;
}



void freeMemory()
{

	cudaFree(KO_dev);
	cudaFree(EW_dev);
	if (objSupport_dev != NULL) cudaFree(objSupport_dev);
	cudaFree(EKO_dev);
	cudaFree(EKOind_dev);

	//cudaFree((n_rec_dev);
	cudaFree(KOi_dev);
	cufftDestroy(plan);
	n_EW8bit = 0;
}

int transfer_KO_and_EW(cufftComplex* _KO_dev, int* _EWint_dev)
{
	LARGE_INTEGER countPerSec, tim1, tim2;

	QueryPerformanceFrequency(&countPerSec);
	QueryPerformanceCounter(&tim1);

	cudaError_t cudaStatus;

	int* EWint_host;
	cufftComplex* KO_host;
	int* EKOind_host = NULL; // Sample index in K-space, used for zero-removal compression

	cudaStatus = cudaHostAlloc((void**)&KO_host, K_xy * K_xy * K_z * sizeof(cufftComplex), cudaHostAllocDefault);
	if (cudaStatus != cudaSuccess) return 401;

	cudaStatus = cudaHostAlloc((void**)&EWint_host, K_xy * K_xy * K_z * sizeof(int), cudaHostAllocDefault);
	if (cudaStatus != cudaSuccess) return 402;

	cudaStatus = cudaHostAlloc((void**)&EKOind_host, K_xy * K_xy * K_z * sizeof(int), cudaHostAllocDefault);
	if (cudaStatus != cudaSuccess) return 403;

	cudaStatus = cudaMemcpyAsync(KO_host, _KO_dev, 2 * K_xy * K_xy * K_z * sizeof(float), cudaMemcpyDeviceToHost, 0);
	if (cudaStatus != cudaSuccess) return 404;

	cudaStatus = cudaMemcpyAsync(EWint_host, _EWint_dev, K_xy * K_xy * K_z * sizeof(int), cudaMemcpyDeviceToHost, 0);
	if (cudaStatus != cudaSuccess) return 405;

	cudaDeviceSynchronize();

	float factor = pow(dxo, 3);


	int j = 0;
	for (int i = 0; i < K_xy * K_xy * K_z; i++) {
		if (EWint_host[i]) {
			KO_host[j].x = KO_host[i].x / factor;
			KO_host[j].y = KO_host[i].y / factor;
			//EKO_host[j].y = 100;
			EKOind_host[j] = i;
			j++;
		}
	}
	n_EW8bit = j;

#ifdef Save_data2
	FILE* file1 = NULL;
	fopen_s(&file1, "KO_data_afterFFTshift_transfer.bin", "wb");
	if (file1) {
		int size = 2 * j;
		fwrite(&size, sizeof(int), 1, file1);
		fwrite(KO_host, sizeof(float), 2 * j, file1);
		fclose(file1);
	}
	fopen_s(&file1, "EW8bit_data_transfer.bin", "wb");
	if (file1) {
		int size = j;
		fwrite(&size, sizeof(int), 1, file1);
		fwrite(EKOind_host, sizeof(int), j, file1);
		fclose(file1);
	}
#endif
	cudaStatus = cudaMalloc((void**)&EKO_dev, j * sizeof(cufftComplex));
	if (cudaStatus != 0) return 409;

	cudaStatus = cudaMemcpyAsync(EKO_dev, KO_host, j * sizeof(cufftComplex), cudaMemcpyHostToDevice, 0);
	if (cudaStatus != 0) return 406;

	cudaStatus = cudaMalloc((void**)&EKOind_dev, j * sizeof(int));
	if (cudaStatus != 0) return 407;

	cudaStatus = cudaMemcpyAsync(EKOind_dev, EKOind_host, j * sizeof(int), cudaMemcpyHostToDevice, 0);
	if (cudaStatus != 0) return 408;

	cudaFreeHost(EWint_host);
	cudaFreeHost(KO_host);
	cudaFreeHost(EKOind_host);

	QueryPerformanceCounter(&tim2);
	double j_total = (double)(tim2.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;

	return 0;
}

/////////////////////
/// Exported Functions
/////////////////////
EXPORTED_FUNCTION int cudaInitDev()
{
	int count = 0;

	cudaGetDeviceCount(&count);
	if (count == 0) return 1;

	if (count > 0) {
		cudaSetDevice(0);
		//OutputDebugStringA("CUDA initialized.\n");
		return 0;
	}
	else
		return 2;
}

EXPORTED_FUNCTION int setDims_n_rec(int x, int y, int z)
{
	K_xy = x;
	K_xy = y;
	K_z = z;
	int BATCH = 1;

	//size_t free;
	//size_t total;
	//cudaMemGetInfo(&free, &total);

	int n[NRANK] = { z, y, x };
	cudaMalloc((void**)&n_rec_dev, x * y * z * sizeof(cufftComplex));
	KOi_dev = n_rec_dev;
	//cudaMalloc((void**)&KOi_dev, x * y * z * sizeof(cufftComplex));
	cudaMalloc((void**)&objSupport_dev, x * y * z * sizeof(float));

	cufftResult_t error = cufftPlanMany(&plan, NRANK, n,
		NULL, 1, x * y * z, // *inembed, istride, idist 
		NULL, 1, x * y * x, // *onembed, ostride, odist
		CUFFT_C2C, BATCH);

	int freeMB;
	int totalMB;

	memoryInfo(&freeMB, &totalMB);

	//cufftResult_t error = cufftPlan3d(&plan, z, y, x, CUFFT_C2C);

	if (error != CUFFT_SUCCESS)
		return error;

	if ((n_rec_dev == NULL) || (KOi_dev == NULL) || (objSupport_dev == NULL))
		return 3;
	else
		return 0;
}

EXPORTED_FUNCTION int setParams(int _do_TC, int _do_NNC, float _n_imm, float _relaxM,
	float _betaM, float _kn_mean, int _objshift, float _dxo, int _nGPi)
{
	do_TC = _do_TC;
	do_NNC = _do_NNC;
	n_imm = _n_imm;
	relaxM = _relaxM;
	betaM = _betaM;
	kn_mean = _kn_mean;
	objShift = _objshift;
	dxo = _dxo;
	nGPi = _nGPi;

	return 0;
}

EXPORTED_FUNCTION int sendData_n_rec(float* n_rec_host)
{
	cudaError_t cudaStatus;

	if (n_rec_dev != NULL) {
		cudaStatus = cudaMemcpy(n_rec_dev, n_rec_host, 2 * K_xy * K_xy * K_z * sizeof(float), cudaMemcpyHostToDevice);

#ifdef Save_data	
		FILE* file1 = NULL;
		fopen_s(&file1, "n_rec_data.bin", "wb");
		if (file1) {
			fwrite(n_rec_host, sizeof(float), 2 * n_rec_x * n_rec_y * n_rec_z, file1);
			fclose(file1);
		}
#endif
	}
	else
		return 4;

	if (cudaStatus != cudaSuccess) {
		return 5;
	}
	return 0;
}


EXPORTED_FUNCTION int sendData_objectSupport(float* _objSupport)
{
	cudaError_t cudaStatus;

	if (_objSupport == NULL) {
		objSupportAv = false;
		return 0;
	}

	if (n_rec_dev != NULL) {
		cudaStatus = cudaMemcpy(objSupport_dev, _objSupport, K_xy * K_xy * K_z * sizeof(float), cudaMemcpyHostToDevice);
		objSupportAv = true;

#ifdef Save_data
		FILE* file1 = NULL;
		fopen_s(&file1, "obj_support_data.bin", "wb");
		if (file1) {
			fwrite(_objSupport, sizeof(float), n_rec_x * n_rec_y * n_rec_z, file1);
			fclose(file1);
		}
#endif
	}
	else
		return 4;

	if (cudaStatus != cudaSuccess) {
		return 5;
	}
	return 0;
}

EXPORTED_FUNCTION int sendData_KO_Re_Im_and_EW8bit(float* _KO, unsigned char* EW8bit, int n)
{
	cudaError_t cudaStatus;
	cufftComplex* EKO_host = NULL;
	cufftComplex* KO = (cufftComplex*)_KO;
	int* EKOind_host = NULL;

#ifdef Save_data
	FILE* file1 = NULL;
	fopen_s(&file1, "KO_data_afterFFTshift.bin", "wb");
	if (file1) {
		fwrite(_KO, sizeof(float), 2 * n_rec_x * n_rec_y * n_rec_z, file1);
		fclose(file1);
	}
	fopen_s(&file1, "EW8bit_data.bin", "wb");
	if (file1) {
		fwrite(EW8bit, sizeof(float), n_rec_x * n_rec_y * n_rec_z, file1);
		fclose(file1);
	}
#endif


	EKO_host = (cufftComplex*)malloc(n * sizeof(cufftComplex));
	EKOind_host = (int*)malloc(n * sizeof(int));

	cudaStatus = cudaMalloc((void**)&EKO_dev, n * sizeof(cufftComplex));
	cudaStatus = cudaMalloc((void**)&EKOind_dev, n * sizeof(int));

	float factor = pow(dxo, 3);
	if ((EKO_host != NULL) && (EKOind_host != NULL) && (EKO_dev != NULL) && (EKO_dev != NULL)) {
		int j = 0;
		for (int i = 0; i < K_xy * K_xy * K_z; i++) {
			if (EW8bit[i]) {
				EKO_host[j].x = KO[i].x / factor;
				EKO_host[j].y = KO[i].y / factor;
				//EKO_host[j].y = 100;
				EKOind_host[j] = i;
				j++;
			}
		}

#ifdef Save_data2
		FILE* file1 = NULL;
		fopen_s(&file1, "KO_data_afterFFTshift_send.bin", "wb");
		if (file1) {
			int size = 2 * j;
			fwrite(&size, sizeof(int), 1, file1);
			fwrite(EKO_host, sizeof(float), 2 * j, file1);
			fclose(file1);
		}
		fopen_s(&file1, "EW8bit_data_send.bin", "wb");
		if (file1) {
			int size = j;
			fwrite(&size, sizeof(int), 1, file1);
			fwrite(EKOind_host, sizeof(int), j, file1);
			fclose(file1);
		}
#endif
		if (j == n) {
			n_EW8bit = n;
			cudaStatus = cudaMemcpy(EKO_dev, EKO_host, n * sizeof(cufftComplex), cudaMemcpyHostToDevice);
			if (cudaStatus == 0) {
				cudaStatus = cudaMemcpy(EKOind_dev, EKOind_host, n * sizeof(int), cudaMemcpyHostToDevice);
			}
		}
		else {
			free(EKO_host);
			free(EKOind_host);
			return 16;
		}
	}
	else {
		free(EKO_host);
		free(EKOind_host);
		return 4;
	}

	free(EKO_host);
	free(EKOind_host);
	if (cudaStatus != cudaSuccess) {
		return 5;
	}
	return 0;
}

EXPORTED_FUNCTION int takeComplexData(float* complexData_host)
{
	cudaError_t cudaStatus1 = cudaMemcpy(complexData_host, n_rec_dev, 2 * K_xy * K_xy * K_z * sizeof(float), cudaMemcpyDeviceToHost);

	if (cudaStatus1 != cudaSuccess) {
		return 6;
	}

#ifdef Save_data
	FILE* file1 = NULL;
	fopen_s(&file1, "results_data_2iter.bin", "wb");
	if (file1) {
		fwrite(complexData_host, sizeof(float), 2 * n_rec_x * n_rec_y * n_rec_z, file1);
		fclose(file1);
	}
#endif

	freeMemory();

	return 0;
}

__global__ void kernel_absSingleFrame(float* output_dev, cufftComplex* input_n_rec_dev, unsigned long int offset)
{
	register unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	float re = input_n_rec_dev[index + offset].x;
	float im = input_n_rec_dev[index + offset].y;

	cufftComplex val;
	val.x = re;
	val.y = im;

	output_dev[index] = cuCabsf(val);
}

int absSingleFrame(float* output_dev, cufftComplex* input_n_rec_dev, unsigned long int offset)
{
	int nThreads = K_xy;
	dim3 nBlocks = K_xy;

	kernel_absSingleFrame << <nBlocks, nThreads >> > (output_dev, input_n_rec_dev, offset);

	cudaError_t cudaStatus1 = cudaGetLastError();

	cudaError_t cudaStatus2 = cudaDeviceSynchronize();

	if ((cudaStatus1 != cudaSuccess) || (cudaStatus2 != cudaSuccess)) {
		return 117;
	}
	return 0;
}

int takeDataSingleSlice(float* complexData_host, int noOfSlice)
{
	if (noOfSlice < 1)
		noOfSlice++;

	float* frame_temp_dev;
	cudaError_t cudaStatus1 = cudaMalloc((void**)&frame_temp_dev, K_xy * K_xy * sizeof(float));
	if (cudaStatus1 != cudaSuccess) {
		freeMemory();
		return 3;
	}

	unsigned long int offset = (K_xy * K_xy * (noOfSlice - 1));

	absSingleFrame(frame_temp_dev, n_rec_dev, offset);

	cudaStatus1 = cudaMemcpy(complexData_host, frame_temp_dev, K_xy * K_xy * sizeof(float), cudaMemcpyDeviceToHost);

	cudaFree(frame_temp_dev);

	if (cudaStatus1 != cudaSuccess) {
		freeMemory();
		return 6;
	}

	freeMemory();

	return 0;
}

int performDI()
{
	cufftResult_t err;

	// fftshift 3D for DI
	int error = ifftShift_3(K_xy, K_xy, K_z, KOi_dev);

	if (error == 0)
		error = ifftShift_2(K_xy, K_xy, K_z, KOi_dev);
	if (error == 0)
		error = ifftShift_1(K_xy, K_xy, K_z, KOi_dev);

	// FFT shift EW - for transfer to spare matrix
	if (nGPi > 0) {
		if (error == 0)
			error = ifftShift_3_int(K_xy, K_xy, K_z, EW_dev);

		if (error == 0)
			error = ifftShift_2_int(K_xy, K_xy, K_z, EW_dev);
		if (error == 0)
			error = ifftShift_1_int(K_xy, K_xy, K_z, EW_dev);

		if ((error == 0))
			error = transfer_KO_and_EW(KOi_dev, EW_dev);
	}

	if (error == 0) {
		err = cufftExecC2C(plan, KOi_dev, n_rec_dev, CUFFT_INVERSE);
		if (err != CUFFT_SUCCESS) error = 14;
		cudaDeviceSynchronize();
	}

	//if (error == 0) {
	//	float factor = (float)K_xy * K_xy * K_z;
	//	error = div_3Dmatrix(n_rec_dev, factor);
	//}


	if (error == 0) {
		float factor = pow(dxo, 3.0f);
		factor = factor * K_xy * K_xy * K_z;
		error = div_3Dmatrix(n_rec_dev, factor);
	}

	if ((!objShift) && (error == 0))
		error = ifftShift_3(K_xy, K_xy, K_z, n_rec_dev);

	if (error == 0)
		error = approxRytovProj2Rec();

	//float * temp = (float*)malloc(K_xy*K_xy * 2 * K_xy * sizeof(float));
	//memset(temp, 0, K_xy*K_xy * 2 * K_xy * sizeof(float));
	//cudaError_t cudaStatus3 = cudaMemcpy(temp, KOi_dev, K_xy*K_xy * 2 * K_xy * sizeof(float), cudaMemcpyDeviceToHost);
	//FILE *plik8;
	//if (preProcOnGPU)
	//	plik8 = fopen("KO_wHL03_PreprocOnGPU.bin", "wb");
	//else
	//	plik8 = fopen("KO_wHL03.bin", "wb");
	//fwrite(temp, sizeof(float), K_xy*K_xy * 2 * K_xy, plik8);
	//fclose(plik8);
	//free(temp);

	return error;
}

EXPORTED_FUNCTION int cudaProcessing()
{
	int error = 0;
	cufftResult_t err;
	int fasterMode = 0;
	//#define Save_timings_GP

	LARGE_INTEGER countPerSec, tim1, tim2, tim3, tim4, tim5, tim6, tim7, tim8, tim9, tim10, tim11, tim12, tim13, tim14, tim15, tim16, tim17, tim18;

	QueryPerformanceFrequency(&countPerSec);
	QueryPerformanceCounter(&tim1);
	if ((!objSupportAv) && do_TC && ((do_NNC * do_NNC) > 0)) {// faster option
		fasterMode = 1;
		for (int i = 0; i < nGPi; i++) { // nGPi loop
			QueryPerformanceCounter(&tim3);
			if ((do_NNC > 0) && (error == 0)) {
				error = transparencyAndNonNegativityConstraintandRytov();
			}
			else if ((do_NNC < 0) && (error == 0)) {
				error = transparencyAndNonPositivityConstraintAndRytov();
			}
			QueryPerformanceCounter(&tim4);
			if ((!objShift) && (error == 0))
				error = ifftShift_3(K_xy, K_xy, K_z, n_rec_dev);
			QueryPerformanceCounter(&tim5);
			if (error == 0) {
				err = cufftExecC2C(plan, n_rec_dev, KOi_dev, CUFFT_FORWARD);
				if (err != CUFFT_SUCCESS) error = 14;
				cudaDeviceSynchronize();
			}
			QueryPerformanceCounter(&tim6);
			if (error == 0)
				error = insertProjData(KOi_dev, EKO_dev, EKOind_dev, n_EW8bit);
			QueryPerformanceCounter(&tim7);
			if (error == 0) {
				err = cufftExecC2C(plan, KOi_dev, n_rec_dev, CUFFT_INVERSE);
				if (err != CUFFT_SUCCESS) error = 14;
				cudaDeviceSynchronize();
			}
			QueryPerformanceCounter(&tim8);
			if ((!objShift) && (error == 0))
				error = ifftShift_3(K_xy, K_xy, K_z, n_rec_dev);
			QueryPerformanceCounter(&tim9);
			float factor = (float)K_xy * K_xy * K_z;
			if (error == 0)
				error = approxRytovProj2RecWithDiv(factor);
			QueryPerformanceCounter(&tim10);
		}
	}
	else {
		fasterMode = 0;
		for (int i = 0; i < nGPi; i++) { // nGPi loop
			QueryPerformanceCounter(&tim3);
			if (do_TC && (do_NNC > 0) && (error == 0)) {
				error = transparencyAndNonNegativityConstraint();
			}
			else if (do_TC && (do_NNC < 0) && (error == 0)) {
				error = transparencyAndNonPositivityConstraint();
			}
			else {

				if (do_TC && (error == 0))		//transparency constraint
					error = transparencyConstraint();

				if ((do_NNC > 0) && (error == 0))		//nonNegativity constraint
					error = nonNegativityConstraint();
				else if ((do_NNC < 0) && (error == 0))
					error = nonPositivityConstraint();
			}
			QueryPerformanceCounter(&tim4);
			if ((error == 0) && (objSupportAv))
				error = spatialSupport();
			QueryPerformanceCounter(&tim6);
			if (error == 0)
				error = approxRytov();
			QueryPerformanceCounter(&tim7);
			if ((!objShift) && (error == 0))
				error = ifftShift_3(K_xy, K_xy, K_z, n_rec_dev);
			QueryPerformanceCounter(&tim8);
			if (error == 0) {
				err = cufftExecC2C(plan, n_rec_dev, KOi_dev, CUFFT_FORWARD);
				if (err != CUFFT_SUCCESS) error = 14;
				cudaDeviceSynchronize();
			}
			QueryPerformanceCounter(&tim9);
			// FFT shift 3D:
			/*if (error == 0)
				error = ifftShift_3(KOi_dev);
			if (error == 0)
				error = ifftShift_2(KOi_dev);
			if (error == 0)
				error = ifftShift_1(KOi_dev);*/
			QueryPerformanceCounter(&tim10);
			//if (error == 0) {
			//	float factor = pow(dxo,3);
			//	error = mult_3Dmatrix(KOi_dev, factor);
			//}
			QueryPerformanceCounter(&tim11);
			if (error == 0)
				error = insertProjData(KOi_dev, EKO_dev, EKOind_dev, n_EW8bit);
			QueryPerformanceCounter(&tim12);
			// FFT shift 3D:
			/*if (error == 0)
				error = ifftShift_3(KOi_dev);
			if (error == 0)
				error = ifftShift_2(KOi_dev);
			if (error == 0)
				error = ifftShift_1(KOi_dev);*/
			QueryPerformanceCounter(&tim13);
			if (error == 0) {
				err = cufftExecC2C(plan, KOi_dev, n_rec_dev, CUFFT_INVERSE);
				if (err != CUFFT_SUCCESS) error = 14;
				cudaDeviceSynchronize();
			}
			QueryPerformanceCounter(&tim14);
			if (error == 0) {
				float factor = (float)K_xy * K_xy * K_z;
				error = div_3Dmatrix(n_rec_dev, factor);
			}
			QueryPerformanceCounter(&tim15);
			//if (error == 0) {
			//	float factor = pow(dxo, 3.0f);
			//	error = div_3Dmatrix(n_rec_dev, factor);
			//}
			QueryPerformanceCounter(&tim16);
			if ((!objShift) && (error == 0))
				error = ifftShift_3(K_xy, K_xy, K_z, n_rec_dev);
			QueryPerformanceCounter(&tim17);
			if (error == 0)
				error = approxRytovProj2Rec();
			QueryPerformanceCounter(&tim18);
		}
	}



	QueryPerformanceCounter(&tim2);
	//store the number of clock ticks in tim2

	double j_total = (double)(tim2.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;
	//convert to milliseconds
	double timings[16];

	timings[0] = (double)(tim4.QuadPart - tim3.QuadPart) / countPerSec.QuadPart * 1000;
	timings[1] = (double)(tim5.QuadPart - tim4.QuadPart) / countPerSec.QuadPart * 1000;
	timings[2] = (double)(tim6.QuadPart - tim5.QuadPart) / countPerSec.QuadPart * 1000;
	timings[3] = (double)(tim7.QuadPart - tim6.QuadPart) / countPerSec.QuadPart * 1000;
	timings[4] = (double)(tim8.QuadPart - tim7.QuadPart) / countPerSec.QuadPart * 1000;
	timings[5] = (double)(tim9.QuadPart - tim8.QuadPart) / countPerSec.QuadPart * 1000;
	timings[6] = (double)(tim10.QuadPart - tim9.QuadPart) / countPerSec.QuadPart * 1000;
	if (fasterMode == 0) {
		timings[7] = (double)(tim11.QuadPart - tim10.QuadPart) / countPerSec.QuadPart * 1000;
		timings[8] = (double)(tim12.QuadPart - tim11.QuadPart) / countPerSec.QuadPart * 1000;
		timings[9] = (double)(tim13.QuadPart - tim12.QuadPart) / countPerSec.QuadPart * 1000;
		timings[10] = (double)(tim14.QuadPart - tim13.QuadPart) / countPerSec.QuadPart * 1000;
		timings[11] = (double)(tim15.QuadPart - tim14.QuadPart) / countPerSec.QuadPart * 1000;
		timings[12] = (double)(tim16.QuadPart - tim15.QuadPart) / countPerSec.QuadPart * 1000;
		timings[13] = (double)(tim17.QuadPart - tim16.QuadPart) / countPerSec.QuadPart * 1000;
		timings[14] = (double)(tim18.QuadPart - tim17.QuadPart) / countPerSec.QuadPart * 1000;
	}
	timings[15] = j_total;

	float t = (float)j_total;

#ifdef Save_timings_GP

	struct tm newtime;
	time_t now = time(0);
	errno_t eeee = localtime_s(&newtime, &now);
	char buffer[256];
	strftime(buffer, 80, "%c", &newtime); //"%x - %I:%M%p"
	FILE* plik = NULL;
	fopen_s(&plik, "timingsGP_loop.txt", "at");

	if (plik != NULL) {
		if (fasterMode == 0) {
			fprintf(plik, "\n\nTime of all iterations in the GP loop:  %8.2f [ms]\n %d iterations performed\n\n", t, nGPi);

			fprintf(plik, "\n 1. transparencyConstraint      %8.3f [ms]", timings[0]);
			fprintf(plik, "\n 2. nonNegativityConstraint     %8.3f [ms]", timings[1]);
			fprintf(plik, "\n 3. spatialSupport              %8.3f [ms]", timings[2]);
			fprintf(plik, "\n 4. approxRytov                 %8.3f [ms]", timings[3]);
			fprintf(plik, "\n 5. ifftShift 3rdDim            %8.3f [ms]", timings[4]);
			fprintf(plik, "\n 6. FFT                         %8.3f [ms]", timings[5]);
			fprintf(plik, "\n 7. fftShift 3D                 %8.3f [ms]", timings[6]);
			fprintf(plik, "\n 8. mult 3D matrix              %8.3f [ms]", timings[7]);
			fprintf(plik, "\n 9. insert proj. data           %8.3f [ms]", timings[8]);
			fprintf(plik, "\n10. fftShift 3D                 %8.3f [ms]", timings[9]);
			fprintf(plik, "\n11. IFFT                        %8.3f [ms]", timings[10]);
			fprintf(plik, "\n12. div 3D matrix               %8.3f [ms]", timings[11]);
			fprintf(plik, "\n13. div 3D matrix 2             %8.3f [ms]", timings[12]);
			fprintf(plik, "\n14. ifftShift 3rdDim            %8.3f [ms]", timings[13]);
			fprintf(plik, "\n15. approxRytovProj2Rec         %8.3f [ms]", timings[14]);
		}
		else {
			fprintf(plik, "\n\nTime of all iterations in the GP loop:  %8.2f [ms]\n %d iterations performed\n\n", t, nGPi);

			fprintf(plik, "\n 1. TC, NNC/NPC, Rytov          %8.3f [ms]", timings[0]);
			fprintf(plik, "\n 2. fft shift                   %8.3f [ms]", timings[1]);
			fprintf(plik, "\n 3. FFT                         %8.3f [ms]", timings[2]);
			fprintf(plik, "\n 4. insert proj. data           %8.3f [ms]", timings[3]);
			fprintf(plik, "\n 5. IFFT                        %8.3f [ms]", timings[4]);
			fprintf(plik, "\n 6. fft shift                   %8.3f [ms]", timings[5]);
			fprintf(plik, "\n 7. div and Rytov               %8.3f [ms]", timings[6]);
		}

		fclose(plik);
	}
#endif //  Save_timings 

	return error;
}

EXPORTED_FUNCTION int cudaFFT_n_rec_2_KOi()
{
	if (cufftExecC2C(plan, n_rec_dev, KOi_dev, CUFFT_FORWARD) == CUFFT_SUCCESS)
		return 0;
	else return 13;
}

EXPORTED_FUNCTION int memoryInfo(int* freeMB, int* totalMB)
{
	size_t total;
	size_t free;

	cudaMemGetInfo(&free, &total);

	*totalMB = (int)(total / 1048576);
	*freeMB = (int)(free / 1048576);

	return 0;
}

int prepareDirectInversion(float* _objSupport)
{
	int BATCH = 1;

	int n[NRANK] = { K_z, K_xy, K_xy };
	KOi_dev = KO_dev;
	n_rec_dev = KOi_dev;

	cufftResult_t error = cufftPlanMany(&plan, NRANK, n,
		NULL, 1, K_xy * K_xy * K_z, // *inembed, istride, idist 
		NULL, 1, K_xy * K_xy * K_z, // *onembed, ostride, odist
		CUFFT_C2C, BATCH);

	int freeMB;
	int totalMB;

	memoryInfo(&freeMB, &totalMB);
	//cufftResult_t error = cufftPlan3d(&plan, z, y, x, CUFFT_C2C);

	if (error != CUFFT_SUCCESS)
		return error;

	//cudaMalloc((void**)&KOi_dev, x * y * z * sizeof(cufftComplex));

	cudaError_t cudaStatus;

	if (_objSupport == NULL) {
		objSupportAv = false;
		objSupport_dev = NULL;
	}
	else {
		cudaMalloc((void**)&objSupport_dev, K_xy * K_xy * K_z * sizeof(float));
		if (objSupport_dev != NULL) {
			cudaStatus = cudaMemcpy(objSupport_dev, _objSupport, K_xy * K_xy * K_z * sizeof(float), cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess)
				return 5;
			objSupportAv = true;
		}
		else
			return 3;
	}

	return 0;
}

EXPORTED_FUNCTION int HL03_setParamsAndStartDIandGP(int _nGPi, int _do_TC, int _do_NNC, float _relaxM,
	float _betaM, float _kn_mean, int _objShift, float* _objSupport)
{

#ifdef Save_timings	
	LARGE_INTEGER countPerSec, tim1, tim2, tim3, tim4, tim5, tim6, tim7, tim8, tim9, tim10, tim11, tim12, tim13, tim14, tim15, tim16, tim17, tim18;

	QueryPerformanceFrequency(&countPerSec);
	QueryPerformanceCounter(&tim1);
#endif //  Save_timings
	cleanDeviceMemoryBeforeDirectInversion();
#ifdef Save_timings
	QueryPerformanceCounter(&tim2);
#endif //  Save_timings
	do_TC = _do_TC;
	do_NNC = _do_NNC;
	relaxM = _relaxM;
	betaM = _betaM;
	if (_kn_mean > -0.001f) kn_mean = _kn_mean;
	objShift = _objShift;
	dxo = dx * (float)Nx / K_xy;
	nGPi = _nGPi;
#ifdef Save_timings
	QueryPerformanceCounter(&tim3);
#endif //  Save_timings
	int error = prepareDirectInversion(_objSupport);
#ifdef Save_timings
	QueryPerformanceCounter(&tim4);
#endif //  Save_timings
	if (error == 0)
		error = performDI();
#ifdef Save_timings	
	QueryPerformanceCounter(&tim5);
#endif //  Save_timings
	if ((error == 0) && (nGPi > 0))
		error = cudaProcessing();

	if ((objShift) && (error == 0)) {
		if (error == 0)
			error = ifftShift_3(K_xy, K_xy, K_z, n_rec_dev);
		if (error == 0)
			error = ifftShift_2(K_xy, K_xy, K_z, n_rec_dev);
		if (error == 0)
			error = ifftShift_1(K_xy, K_xy, K_z, n_rec_dev);
	}
#ifdef Save_timings	
	QueryPerformanceCounter(&tim6);

	double time[5];

	time[0] = (double)(tim2.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;
	time[1] = (double)(tim3.QuadPart - tim2.QuadPart) / countPerSec.QuadPart * 1000;
	time[2] = (double)(tim4.QuadPart - tim3.QuadPart) / countPerSec.QuadPart * 1000;
	time[3] = (double)(tim5.QuadPart - tim4.QuadPart) / countPerSec.QuadPart * 1000;
	time[4] = (double)(tim6.QuadPart - tim5.QuadPart) / countPerSec.QuadPart * 1000;

	FILE* plik = NULL;
	fopen_s(&plik, "timingsGPU.txt", "at");

	if (plik != NULL) {

		double timeTot = (double)(tim6.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;

		fprintf(plik, "\n\nTotal time 03_setParamsAndStartDIandGP:        \t%10.2f", timeTot);
		fprintf(plik, "\n  - time cleanDeviceMemoryBeforeDirectInversion: \t%10.2f", time[0]);
		fprintf(plik, "\n  - time of variable reassignment:               \t%10.2f", time[1]);
		fprintf(plik, "\n  - time prepareDirectInversion:                 \t%10.2f", time[2]);
		fprintf(plik, "\n  - time performDI:                              \t%10.2f", time[3]);
		fprintf(plik, "\n  - time GP:                                     \t%10.2f", time[4]);

		fclose(plik);
	}
#endif //  Save_timings
	return error;
}

EXPORTED_FUNCTION int HL04_takeReconstructionAndFreeMemory(float* complexData_host)
{
	LARGE_INTEGER countPerSec, tim1, tim2, tim3;
#ifdef Save_timings

	QueryPerformanceFrequency(&countPerSec);
	QueryPerformanceCounter(&tim1);
#endif //  Save_timings
	int error = takeComplexData(complexData_host);
#ifdef Save_timings
	QueryPerformanceCounter(&tim2);
#endif //  Save_timings
	freeMemory();
#ifdef Save_timings
	QueryPerformanceCounter(&tim3);

	FILE* plik = NULL;
	fopen_s(&plik, "timingsGPU.txt", "at");

	if (plik != NULL) {

		double time1 = (double)(tim2.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;
		double time2 = (double)(tim3.QuadPart - tim2.QuadPart) / countPerSec.QuadPart * 1000;
		double timeTot = (double)(tim3.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;

		fprintf(plik, "\n\nTotal time HL04_takeReconstructionAndFreeMemory: \t%10.2f", timeTot);
		fprintf(plik, "\n  - time of data transfer:                         \t%10.2f", time1);
		fprintf(plik, "\n  - time of memory deallocation:                   \t%10.2f", time2);

		fprintf(plik, "\n\n_______________________________________");
		fprintf(plik, "\nInfo o pliku:");
		fprintf(plik, "\n - K-space:          %d x %d x %d", K_xy, K_xy, K_z);
		fprintf(plik, "\n - K-space [MB]:     %d", K_xy * K_xy * K_z * 8 / (1024 * 1024));
		fprintf(plik, "\n - number of iterations:         %d", nGPi);
		fprintf(plik, "\n_______________________________________");

		fclose(plik);
	}
#endif //  Save_timings
	return error;
}


EXPORTED_FUNCTION int HL04a_takeRecSingleSliceAndFreeMemory(float* complexData_host, int noOfSlice)
{
	LARGE_INTEGER countPerSec, tim1, tim2, tim3;
#ifdef Save_timings

	QueryPerformanceFrequency(&countPerSec);
	QueryPerformanceCounter(&tim1);
#endif //  Save_timings
	int error = takeDataSingleSlice(complexData_host, noOfSlice);
#ifdef Save_timings
	QueryPerformanceCounter(&tim2);
#endif //  Save_timings
	freeMemory();
#ifdef Save_timings
	QueryPerformanceCounter(&tim3);

	FILE* plik = NULL;
	fopen_s(&plik, "timingsGPU.txt", "at");

	if (plik != NULL) {

		double time1 = (double)(tim2.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;
		double time2 = (double)(tim3.QuadPart - tim2.QuadPart) / countPerSec.QuadPart * 1000;
		double timeTot = (double)(tim3.QuadPart - tim1.QuadPart) / countPerSec.QuadPart * 1000;

		fprintf(plik, "\n\nTotal time HL04_takeReconstructionAndFreeMemory: \t%10.2f", timeTot);
		fprintf(plik, "\n  - time of data transfer:                         \t%10.2f", time1);
		fprintf(plik, "\n  - time of memory deallocation:                   \t%10.2f", time2);

		fprintf(plik, "\n\n_______________________________________");
		fprintf(plik, "\nInfo o pliku:");
		fprintf(plik, "\n - K-space:          %d x %d x %d", K_xy, K_xy, K_z);
		fprintf(plik, "\n - K-space [MB]:     %d", K_xy * K_xy * K_z * 8 / (1024 * 1024));
		fprintf(plik, "\n - number of iterations:         %d", nGPi);
		fprintf(plik, "\n_______________________________________");

		fclose(plik);
	}
#endif //  Save_timings
	return error;
}