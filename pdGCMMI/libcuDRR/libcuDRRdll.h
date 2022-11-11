#pragma once
#include "pch.h"
#include <npp.h> //max

enum class DEVICE_TYPE
{
	HOST = 0,
	DEVICE = 1
};

enum class SHRINK_ON_OFF
{
	OFF = 0,
	ON = 1
};

typedef struct __struct_params_DRR__
{
	//infomation for initial geometry
//	GEO_METHOD method;
//	GEO_TYPE type;
	SHRINK_ON_OFF shink;
	float sid, sod;
	float depth_offset;
	float elevation;
	int reverse_gantry_angle;

	//infomation for volume
	int vol_x, vol_y, vol_z;
	int shrink_size;
	float pitch_x, pitch_y, pitch_z;
	float origin_x, origin_y, origin_z;
	
	//infomation for projection image
	int prj_u, prj_v;
	float pitch_u, pitch_v;
	float offset_u, offset_v;
	int prj_bits_count;

	int roi_x, roi_y, roi_w, roi_h;

}SPARAMS_DRR;


class etriDRR {
public:
	etriDRR(void);

	bool init_etriDRR(SPARAMS_DRR* params);
	bool set_etriDRR(SPARAMS_DRR* params, float* vol_h, float* vol_d);
	bool set_volume(float* vol, DEVICE_TYPE d_type);
	bool digital_foward_projection(float* geo);
	bool etri_cuda_device_alloc(float** d_mem, size_t mem_size); 
	bool etri_cuda_device_free(void);
	float* get_host_image(void);
	float* get_device_image(void);

private:
	SPARAMS_DRR* params;
	float* d_vol;
	float* d_prj;

	Npp8u* nppMinMaxBuffer;
	int nBufferSize;
	float* max;
	float* min;
};

extern "C" {
	DLL_DESC etriDRR* create_etriDRR(SPARAMS_DRR* params, float* vol_h, float* vol_d);
	DLL_DESC etriDRR* create_etriDRR_class(SPARAMS_DRR* params);
	DLL_DESC bool set_vol_etriDRR(etriDRR* params, float* vol, DEVICE_TYPE type);
	DLL_DESC float *get_device_DRR_image(etriDRR *etri_drr, float *geo);
	DLL_DESC float* get_image_device_to_host(etriDRR* etri_drr);
}
