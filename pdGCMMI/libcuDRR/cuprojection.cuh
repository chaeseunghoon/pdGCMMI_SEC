#pragma warning(disable:4244)

#include <cuda_runtime.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda_texture_types.h>//texture<> depend head file
#include <texture_fetch_functions.h>

#include "libcuDRRdll.h"

__global__ void kernelEtriDRRplaneX(float* pfPrjs, float* pVol);
__global__ void kernelEtriDRRplaneY(float* pfPrjs, float* pVol);
__global__ void kernelEtriDRRplaneZ(float* pfPrjs, float* pVol);
__global__ void kernelNormalization2D32f(float* img);
__global__ void kernelEtriShrink(float* dst, float* src);

bool etri_set_common_params(SPARAMS_DRR* params);
int etri_set_individual_params(SPARAMS_DRR* params, float* trans);
bool cuda_foward_projection(float* d_prj, float* d_vol, int w, int h, int dir);
bool normalization_digital_projeciton_image(float* d_prj, int w, int h, float norm_factor);
bool cuda_volume_shrink(float* dst, float* src, int old_x, int old_y, int old_z, int new_x, int new_y, int new_z, float pitch_x, float pitch_y, float pitch_z, int scale_x, int scale_y, int scale_z);

__device__ __constant__ int C_PRJ_U, C_PRJ_V;
__device__ __constant__ int C_VOL_X, C_VOL_Y, C_VOL_Z;
__device__ __constant__ float C_ORIGIN_X, C_ORIGIN_Y, C_ORIGIN_Z;
__device__ __constant__ int C_OLD_VOL_X, C_OLD_VOL_Y, C_OLD_VOL_Z;
__device__ __constant__ int C_SCALE_X, C_SCALE_Y, C_SCALE_Z;
__device__ __constant__ float C_PITCH_X, C_PITCH_Y, C_PITCH_Z;
__device__ __constant__ int C_ROI_X, C_ROI_Y;
__device__ __constant__ float C_SOURCE_POS_X, C_SOURCE_POS_Y, C_SOURCE_POS_Z;
__device__ __constant__ float C_DETECTOR_ORI_X, C_DETECTOR_ORI_Y, C_DETECTOR_ORI_Z;	//detector origin
__device__ __constant__ float C_DETECTOR_U_X, C_DETECTOR_U_Y, C_DETECTOR_U_Z;	//detector u vector
__device__ __constant__ float C_DETECTOR_V_X, C_DETECTOR_V_Y, C_DETECTOR_V_Z;	//detector v vector
__device__ __constant__ float C_N_FACTOR;

