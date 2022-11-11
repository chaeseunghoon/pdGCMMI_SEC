#define DLL_EXPORT

#include "pch.h"

#include <stdio.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "libcuDRRdll.h"
#include "cuprojection.cuh"

extern "C" DLL_DESC etriDRR* create_etriDRR(SPARAMS_DRR* params, float *vol_h, float *vol_d)
{
	etriDRR* class_drr = new etriDRR();
	
	class_drr->set_etriDRR(params, vol_h, vol_d);

	return class_drr;
}

extern "C" DLL_DESC etriDRR * create_etriDRR_class(SPARAMS_DRR * params)
{
	etriDRR* class_drr = new etriDRR();

	class_drr->init_etriDRR(params);

	return class_drr;
}

// etriDRR class에 volume 데이터 입력
// DEVICE_TYPE == HOST: DEVICE 메모리 할당 후 입력
// DEVICE_TYPE == DEVICE: 배열 위치만 연결
// SHRINK_ON_OFF == 0: volume의 W, H, D가 512이하가 되도록 조절
extern "C" DLL_DESC bool set_vol_etriDRR(etriDRR * etri_drr, float* vol, DEVICE_TYPE type)
{
	etri_drr->set_volume(vol, type);

	return true;
}

extern "C" DLL_DESC float* get_device_DRR_image(etriDRR* etri_drr, float* geo)
{
	etri_drr->digital_foward_projection(geo);

	return etri_drr->get_device_image();
}

extern "C" DLL_DESC float* get_image_device_to_host(etriDRR * etri_drr)
{
	return etri_drr->get_host_image();
}


etriDRR::etriDRR(void)
{
	d_vol = NULL;
	d_prj = NULL;
	params = NULL;
	max = NULL;
	min = NULL;
	nBufferSize = 0;
	nppMinMaxBuffer = NULL;
}


float* etriDRR::get_device_image(void)
{
	return d_prj;
}

float* etriDRR::get_host_image(void)
{
	float* img = NULL;
	img = (float*)malloc(sizeof(float) * static_cast<size_t>(params->prj_u) * params->prj_v);

	cudaMemcpy(img, d_prj, sizeof(float) * static_cast<size_t>(params->prj_u) * params->prj_v, cudaMemcpyDeviceToHost);

	return img;
}

bool etriDRR::set_volume(float* vol, DEVICE_TYPE d_type)
{
	bool bRtn = true;

	int t_xd = params->vol_x, t_yd = params->vol_y, t_zd = params->vol_z;
	int scale_z = 1, scale_y = 1, scale_x = 1;
	size_t vol_size;
	float* tmp_vol;

	// SHRINK ON일경우 축소비율 확인
	if (params->shink == SHRINK_ON_OFF::ON)
	{
		//image dimension의 축소 by 1/2
		if (t_zd > params->shrink_size) {
			do {
				t_zd /= 2;
				scale_z *= 2;
			} while (t_zd > params->shrink_size);
		}
		if (t_yd > params->shrink_size) {
			do {
				t_yd /= 2;
				scale_y *= 2;
			} while (t_yd > params->shrink_size);
		}
		if (t_xd > params->shrink_size) {
			do {
				t_xd /= 2;
				scale_x *= 2;
			} while (t_xd > params->shrink_size);
		}
	}

	vol_size = static_cast<size_t>(t_xd) * t_yd * t_zd;
	
	if (d_vol != NULL)
	{
		// 기존 vol과 크기가 달려질 경우 메모리 해제 후 재할당, 기존 크기와 동일하면 그냥 재사용
		if (!(params->vol_x == t_xd && params->vol_y == t_yd && params->vol_z == t_zd))
		{
			cudaFree(d_vol);
			etri_cuda_device_alloc(&d_vol, sizeof(float) * vol_size);
		}
	}
	else
		etri_cuda_device_alloc(&d_vol, sizeof(float) * vol_size);

	// scale이 모두 1인 경우
	// shrink 과정 없이 내부 gpu 할당
	if (scale_x == 1 && scale_y == 1 && scale_z == 1)
	{
		// 입력된 vol 데이터 위치에 따라 복사 방법 변경
		if (d_type == DEVICE_TYPE::DEVICE)
			cudaMemcpy(d_vol, vol, sizeof(float) * vol_size, cudaMemcpyDeviceToDevice);
		else
			cudaMemcpy(d_vol, vol, sizeof(float) * vol_size, cudaMemcpyHostToDevice);
	}
	else
	{
		//shrink가 발생될 경우 
		//volume 데이터를 Device에 저장후 Cuda Shrink 수행
		if (d_type == DEVICE_TYPE::HOST)
		{
			etri_cuda_device_alloc(&tmp_vol, sizeof(float) * static_cast<size_t>(params->vol_x) * params->vol_y * params->vol_z);
			cudaMemcpy(tmp_vol, vol, sizeof(float) * static_cast<size_t>(params->vol_x) * params->vol_y * params->vol_z, cudaMemcpyHostToDevice);
			cuda_volume_shrink(d_vol, tmp_vol, params->vol_x, params->vol_y, params->vol_z, t_xd, t_yd, t_zd, params->pitch_x, params->pitch_y, params->pitch_z, scale_x, scale_y, scale_z);
			cudaFree(tmp_vol);
		}
		else // volume 데이터가 이미 Device일 경우 바로 cuda shrink 수행
			cuda_volume_shrink(d_vol, vol, params->vol_x, params->vol_y, params->vol_z, t_xd, t_yd, t_zd, params->pitch_x, params->pitch_y, params->pitch_z, scale_x, scale_y, scale_z);
	}

	// 실제 내부에서 사용할 volume size update
	params->vol_x = t_xd;
	params->vol_y = t_yd;
	params->vol_z = t_zd;
	params->pitch_z *= scale_z;
	params->pitch_y *= scale_y;
	params->pitch_x *= scale_x;

	return bRtn;
}

bool etriDRR::digital_foward_projection(float* trans)
{
	// DRR for CT gantry rotation
	// only gantry rotation beta is allowed.
	// identical with z-axis rotation of getDRR2
	bool bRtn = true;
	int dir;
	float norm_factor;
	float max_pixel;

	dir = etri_set_individual_params(params, trans);
	cudaMemset(d_prj, 0, sizeof(float) * static_cast<size_t>(params->prj_u) * params->prj_v);

	cuda_foward_projection(d_prj, d_vol, params->roi_w, params->roi_h, dir);
	
	nppsMinMax_32f(d_prj, params->prj_u * params->prj_v, min, max, nppMinMaxBuffer);
	cudaMemcpy(&max_pixel, max, sizeof(float), cudaMemcpyDeviceToHost);

	norm_factor = ((2 << (params->prj_bits_count - 1)) - 1) / max_pixel;
	//norm_factor = max_pixel;
	normalization_digital_projeciton_image(d_prj, params->roi_w, params->roi_h, norm_factor);

	return bRtn;

}

bool etriDRR::set_etriDRR(SPARAMS_DRR* input_params, float* vol_h, float* vol_d)
{
	size_t vol_size;
	size_t prj_size;

	d_vol = NULL;
	params = input_params;

	vol_size = static_cast<size_t>(params->vol_x) * params->vol_y * params->vol_z;
	prj_size = static_cast<size_t>(params->prj_u) * params->prj_v;

	if (vol_d != NULL)
	{
		d_vol = vol_d;
	}
	else
	{
		if (vol_h != NULL)
		{
			etri_cuda_device_alloc(&d_vol, sizeof(float) * vol_size);
			cudaMemcpy(d_vol, vol_h, sizeof(float) * vol_size, cudaMemcpyHostToDevice);
		}
	}

	etri_cuda_device_alloc(&d_prj, sizeof(float) * prj_size);

	params->roi_x = 0;
	params->roi_y = 0;
	params->roi_w = params->prj_u;
	params->roi_h = params->prj_v;

	etri_set_common_params(params);

	nppsMinMaxGetBufferSize_32f(params->prj_u * params->prj_v, &nBufferSize);
	cudaMalloc(&nppMinMaxBuffer, nBufferSize);
	cudaMalloc(&max, sizeof(float));
	cudaMalloc(&min, sizeof(float));

	return true;
}


bool etriDRR::init_etriDRR(SPARAMS_DRR* input_params)
{
	size_t prj_size;

	params = input_params;
	d_vol = NULL;
	d_prj = NULL;

	prj_size = static_cast<size_t>(params->prj_u) * params->prj_v;

	etri_set_common_params(params);
	

	etri_cuda_device_alloc(&d_prj, sizeof(float) * prj_size);

	nppsMinMaxGetBufferSize_32f(params->prj_u * params->prj_v, &nBufferSize);
	cudaMalloc(&nppMinMaxBuffer, nBufferSize);
	cudaMalloc(&max, sizeof(float));
	cudaMalloc(&min, sizeof(float));

	return true;
}


bool etriDRR::etri_cuda_device_alloc(float** d_mem, size_t mem_size)
{
	cudaError_t error = cudaMalloc(d_mem, mem_size);
	if (error != cudaSuccess)
	{
		printf("Error in cudaMalloc us projection memory with %s\n", cudaGetErrorName(error));
		return false;
	}

	cudaMemset(*d_mem, 0, mem_size);

	return true;
}


bool etriDRR::etri_cuda_device_free(void)
{
	bool bRtn = true;
	if (d_vol != NULL) cudaFree(d_vol);
	if (d_prj != NULL) cudaFree(d_prj);

	cudaFree(nppMinMaxBuffer);
	cudaFree(max);
	cudaFree(min);
	return bRtn;
}
