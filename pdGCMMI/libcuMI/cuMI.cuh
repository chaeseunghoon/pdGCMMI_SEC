#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "libcuMIdll.h"

__device__ __constant__ float C_TRANS_X, C_TRANS_Y, C_COS_VAL, C_SIN_VAL;
__device__ __constant__ float C_N_FACTOR;
__device__ __constant__ int C_ROI_X, C_ROI_Y, C_IMG_W, C_IMG_H;
__device__ __constant__ int C_HISTBIN;
__device__ __constant__ int C_SUM;

__global__ void kernelEtriRotationShift(unsigned short* trans, unsigned short* ori);
__global__ void kernelNormalization2D16u(unsigned short* img);
__global__ void kernelJointHistogram(float* hist, float* img, unsigned short* ref, unsigned char* mask);
__global__ void kernelNormJointHistogram(float* d_hist);
__global__ void kernelSumX(float* d_hist, float* pSum);
__global__ void kernelSumY(float* d_hist, float* pSum);
__global__ void kernelGetMI(float* d_hist, float* prob_x, float* prob_y, float* mi_val);
__global__ void kernelGetMI_X(float* d_hist, float* prob_x, float* prob_y, float* mi_val);
__global__ void kernelGetMI_Y(float* d_hist, float* prob_x, float* prob_y, float* mi_val);

bool img_rotation_shift(unsigned short* d_trans_ref, unsigned short* d_ref, int w, int h, int roi_w, int roi_h, float trans_x, float trans_y, float rot);
bool etri_set_common_params(SPARAMS_MI* params);
bool normalization_projeciton_image(unsigned short* d_img, int roi_w, int roi_h, int w, int h, float norm_factor);
bool joint_histogram(float* hist, float* img, unsigned short* ref, int roi_w, int roi_h, unsigned char* mask, int bin);
bool normalization_joint_histogram(float *d_hist, int bin);
float get_MI_value(float* d_hist, float* prob_x, float* prob_y, float *d_mi, int bin);
