#define DLL_EXPORT

#include "pch.h"
#include <stdio.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math.h>
#include "libcuMIdll.h"
#include "cuMI.cuh"

extern "C" DLL_DESC etriMI * create_etriMI(SPARAMS_MI * params)
{
	etriMI* class_mi = new etriMI();

	class_mi->set_etriMI(params);

	return class_mi;
}


extern "C" DLL_DESC void set_device_MI_images(etriMI * class_mi, float* d_prj, unsigned short* h_ref)
{
	class_mi->set_device_images(d_prj, h_ref);
}

extern "C" DLL_DESC void get_infomation_MI(etriMI * class_mi, float *trans_u, float *trans_v, float *rot )
{
	class_mi->matching_MI(trans_u, trans_v, rot);
}

etriMI::etriMI(void)
{
	d_img = NULL;
	d_ref = NULL;
	d_mask = NULL;
	d_trans_ref = NULL;
	d_hist = NULL;
	prob_x = NULL;
	prob_y = NULL;
	params = NULL;
	max_v = NULL;
	min_v = NULL;
	nBufferSize = 0;
	nppMinMaxBuffer = NULL;
	
}
bool etriMI::set_etriMI(SPARAMS_MI* params_mi)
{
	bool bRtn = true;

	params = params_mi;

	if (d_trans_ref == NULL)
		cudaMalloc(&d_trans_ref, sizeof(unsigned short) * static_cast<size_t>(params->w) * params->h);
	if (d_hist == NULL)
		cudaMalloc(&d_hist, sizeof(float) * params->hist_bins * params->hist_bins);

	if (prob_x == NULL) cudaMalloc(&prob_x, sizeof(float) * params->hist_bins);

	if (prob_y == NULL) cudaMalloc(&prob_y, sizeof(float) * params->hist_bins);
	
	if (d_mask == NULL) cudaMalloc(&d_mask, sizeof(unsigned char) * static_cast<size_t>(params->w) * params->h);
	
	cudaMemcpy(d_mask, params->mask, sizeof(unsigned char) * static_cast<size_t>(params->w) * params->h, cudaMemcpyHostToDevice);

	etri_set_common_params(params);

	nppsMinMaxGetBufferSize_16u(params->w * params->h, &nBufferSize);
	cudaMalloc(&nppMinMaxBuffer, nBufferSize);
	cudaMalloc(&max_v, sizeof(unsigned short));
	cudaMalloc(&min_v, sizeof(unsigned short));

	cudaMalloc(&d_mi, sizeof(float));


	return bRtn;
}

bool etriMI::set_device_images(float* d_prj, unsigned short* h_ref)
{
	bool bRtn = true;
	size_t img_size = static_cast<size_t>(params->w) * params->h;

	if (d_ref == NULL)
	{
		cudaMalloc(&d_ref, sizeof(unsigned short) * img_size);
		cudaMalloc(&d_trans_ref, sizeof(unsigned short) * img_size);
	}

	cudaMemcpy(d_ref, h_ref, sizeof(unsigned short) * img_size, cudaMemcpyHostToDevice);
	
	d_img = d_prj;

	return bRtn;
}

bool etriMI::matching_MI(float *trs_u, float *trs_v, float *rot)
{
	bool bRtn = true;
	
	
	float MAX_TRS = params->trans_max / params->pitch_x; //mm단위를 pixel단위로 변환
	float MIN_TRS = params->trans_min / params->pitch_x; //mm단위를 pixel단위로 변환
	float MAX_ANG = (asin(params->trans_max * 2 * 2 / params->pitch_x / params->w)); //  ~ MAX_TRS
	float MIN_ANG = (asin(params->trans_min * 2 * 2 / params->pitch_x / params->w)); //  ~ MIN_TRS

	//자신+12-nearest neighbor의 6D 좌표값
	float trs_step_x[7] = { 0,				MAX_TRS,		0,				-MAX_TRS,		0,				0,				0 };
	float trs_step_y[7] = { 0,				0,				MAX_TRS,		0,				-MAX_TRS,		0,				0 };
	float ang_step[7] = { 0,				0,				0,				0,				0,				MAX_ANG,		-MAX_ANG };

	float trs_x_sum = 0.0;
	float trs_y_sum = 0.0;
	float rot_sum = 0.0;

	int i;
	unsigned short max_pixel;
	float norm_factor;
	float mi;
	float mi_max;
	int mi_max_index;
	int iter = 0;
	char buf[256];

	cudaMemset(d_hist, 0, sizeof(float) * params->hist_bins * params->hist_bins);
	
	mi_max = 0;
	do {
		mi_max_index = 0;

		for (i = 0; i < 7; i++)
		{
			if (iter != 0 && i == 0) continue;

			//img_rotation_shift(d_trans_ref, d_ref, params->w, params->h, params->roi_w, params->roi_h, trs_x_sum + trs_step_x[i], trs_y_sum + trs_step_y[i], rot_sum + ang_step[i]);
			img_rotation_shift(d_trans_ref, d_ref, params->w, params->h, params->roi_w, params->roi_h, 50, 50, 45*(3.14/180));

			nppsMinMax_16u(d_trans_ref, params->w * params->h, min_v, max_v, nppMinMaxBuffer);
			cudaMemcpy(&max_pixel, max_v, sizeof(unsigned short), cudaMemcpyDeviceToHost);
		
			norm_factor = ((2 << (params->prj_bits_count - 1)) - 1) / (float)max_pixel;
			normalization_projeciton_image(d_trans_ref, params->roi_w, params->roi_h, params->w, params->h, norm_factor);

			sprintf_s(buf, "E:\\sec\\13\\ref\\ref_%d.raw", i);
			write_device_image_u16(buf, d_trans_ref, params->w, params->h);

			joint_histogram(d_hist, d_img, d_trans_ref, params->roi_w, params->roi_h, d_mask, params->hist_bins);

			cudaMemset(prob_x, 0, sizeof(float) * params->hist_bins);
			cudaMemset(prob_y, 0, sizeof(float) * params->hist_bins);

			normalization_joint_histogram(d_hist, params->hist_bins);

			cudaMemset(d_mi, 0, sizeof(float));

			mi = get_MI_value(d_hist, prob_x, prob_y, d_mi, params->hist_bins);

			if (mi > mi_max)
			{
				mi_max_index = i;
				mi_max = mi;
			}
		}
		trs_x_sum += trs_step_x[mi_max_index];
		trs_y_sum += trs_step_y[mi_max_index];
		rot_sum += ang_step[mi_max_index];

		iter++;
	} while (((trs_step_x[1] > MIN_TRS || (ang_step[3] > MIN_ANG)) && (iter < params->iter_bound)));

	*trs_u = trs_x_sum * params->pitch_x;
	*trs_v = trs_y_sum * params->pitch_y;
	*rot = rot_sum;

	return bRtn;
}

void etriMI::write_device_image_u16(char* fname, unsigned short* d_img, int w, int h)
{
	FILE* fp = NULL;

	unsigned short *h_img = NULL;
	size_t img_size = static_cast<size_t>(w) * h;
	
	h_img = (unsigned short*)malloc(sizeof(unsigned short) * img_size);

	cudaMemcpy(h_img, d_img, sizeof(unsigned short) * img_size, cudaMemcpyDeviceToHost);

	fopen_s(&fp, fname, "wb");
	fwrite(h_img, sizeof(unsigned short), img_size, fp);
	fclose(fp);

	free(h_img);
	//printf("u16 image save: %s (%d, %d)\n", fname, w, h);
}


void etriMI::write_device_image_f32(char* fname, float* d_img, int w, int h)
{
	FILE* fp = NULL;

	float* h_img = NULL;
	size_t img_size = static_cast<size_t>(w) * h;

	h_img = (float*)malloc(sizeof(float) * img_size);

	cudaMemcpy(h_img, d_img, sizeof(float) * img_size, cudaMemcpyDeviceToHost);

	fopen_s(&fp, fname, "wb");
	fwrite(h_img, sizeof(float), img_size, fp);
	fclose(fp);

	free(h_img);
	printf("f32 image save: %s (%d, %d)\n", fname, w, h);
}