#pragma once
#include "pch.h"
#include <npp.h> //max
#pragma warning(disable:4244)

typedef struct __struct_params_MI__
{
	float trans_max, trans_min;
	int hist_bins;
	int w, h;
	float pitch_x, pitch_y;
	int roi_x, roi_y, roi_w, roi_h;
	unsigned char* mask;
	int mask_sum;
	int prj_bits_count;
	int iter_bound;

}SPARAMS_MI;


class etriMI {
public:
	etriMI(void);

	bool set_etriMI(SPARAMS_MI* params);
	bool set_device_images(float* d_prj, unsigned short* h_ref);
	bool matching_MI(float *trs_u, float *trs_v, float *rot);
	void write_device_image_u16(char* fname, unsigned short* d_img, int w, int h);
	void write_device_image_f32(char* fname, float* d_img, int w, int h);

	private:
	SPARAMS_MI* params;
	float* d_img;
	unsigned short* d_ref;
	unsigned short* d_trans_ref;
	unsigned char* d_mask;
	float* d_hist;
	float* prob_y;
	float* prob_x;

// projection image의 최대값 최소값을 찾기 위한 변수
	unsigned short* max_v;
	unsigned short* min_v;
	float* d_mi;

	Npp8u* nppMinMaxBuffer;
	int nBufferSize;
};

extern "C" {

	DLL_DESC etriMI* create_etriMI(SPARAMS_MI* params);
	DLL_DESC void set_device_MI_images(etriMI* class_mi, float* d_prj, unsigned short* h_ref);
	DLL_DESC void get_infomation_MI(etriMI* class_mi, float *trans_u, float *trans_v, float *rot);
}

