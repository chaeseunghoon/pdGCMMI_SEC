#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <io.h>
#include <math.h>
#include <time.h>
#include <cuda_runtime.h>
#include "utils_params.h"
#include "utils_PM.h"
#include "utils_proj_image.h"
#include "libcuDRRdll.h"
#include "libcuMIdll.h"

#define __DEBUG__
#pragma warning(disable: 4819)

enum class VOL_TYPE
{
	UINT16 = 0,
	FLOAT32 = 1,
};

bool set_parameter(char* fname, SPARAMS* params);
void print_error(const char* meg);
void create_params_DRR(SPARAMS* src, SPARAMS_DRR* dst);
float *load_volume(SPARAMS* params, VOL_TYPE type);
float **create_drr_info(SPARAMS* params);
void create_params_MI(SPARAMS* src, SPARAMS_MI* dst);
bool load_projection_image(char* fname, unsigned short *img, SPARAMS* params);
int count_number(int val);
void init_mask(SPARAMS* params, float limit);

int main(int argc, char* argv[])
{
	char param_file[] = "E:\\SEC\\13\\gcmparams_gcm.ini"; 
	char proj_path[] = "E:\\SEC\\13\\projections";

	SPARAMS params;
	SPARAMS_DRR params_DRR;
	SPARAMS_MI params_MI;
	
	etriDRR* etri_prj;
	etriMI* etri_mi;
	float* vol = NULL;
	float* d_prj = NULL;
	float** geo = NULL;
	unsigned short* img;
	char buf[MAX_LENGTH];
	char temp_number[100];
	int i, j;
	int cnt, idx;
	int number_len;
	float trans_u, trans_v, rot;
	FILE* fp = NULL;
	float* gcm = NULL;
	int makePMmode = 0;
	cudaError_t error;

	clock_t start1, end1;
	
	/* Demo를 위한 parameter 로딩*/
	if (!set_parameter(param_file, &params))
	{
		print_error("set parameter");
	}

	if (makePMmode == 1)
	{
		float* pm;
		char buf[256];
		printf("make_pm\n");
		sprintf_s(buf, "E:\\SEC\\13\\PrjAngle_00.txt");	
		pm = makeProjectionMatrix(&params);
		write_PM(buf, pm, params.N_prj, params.reverse_gantry_angle);

		return 0;
	}
	/* Demo를 위한 Volume 데이터 로딩*/
	vol = load_volume(&params, VOL_TYPE::UINT16);
	if (vol == NULL)
	{
		printf("does not found file: %s\n", params.buf_vol_file);
		exit(0);
	}
	/* MI를 위한 ROI Mask 생성 */
	init_mask(&params, 1.0f);

	/* libDRR을 위한 파라미터 구조체 생성*/
	create_params_DRR(&params, &params_DRR);

	/* libMI를 위한 파라미터 구조체 생성 */
	create_params_MI(&params, &params_MI);

	/* params 정보를 갖는 DRR용 class 생성 */
	etri_prj = create_etriDRR_class(&params_DRR);

	/* etri_prj에 volume 데이터 입력 */
	/* vol의 저장공간이 HOST일 경우와 DEVICE(GPU)일 경우를 구분 */
	set_vol_etriDRR(etri_prj, vol, DEVICE_TYPE::HOST);

	//volume데이터가 Device에 올린후 메모리 해제
	//DEVICE_TYPE이 HOST일 경우만 수행해야 함
	cudaFreeHost(vol);

	etri_mi = create_etriMI(&params_MI);

	/* planar CT 에 맞게 Geo 변환 */
	geo = create_drr_info(&params);
	if (geo == NULL)
	{
		return 0;
	}
	//write_PM("E:\\xavis\\H160\\case3\\prj_00_i.txt", params.pm, params.N_prj, params.reverse_gantry_angle);
	//exit(0);

	gcm = (float*)malloc(sizeof(float) * 3 * params.N_prj);

	i = 0;
	start1 = clock();

	// MI를 위한 projection영상 저장 공간
	error = cudaHostAlloc((void**)&img, sizeof(unsigned short) * static_cast<size_t>(params.prj_u) * params.prj_v, 0);
	float* tmp;
	char save_buf[256];
	FILE* fp1;

	for (i = 0; i < params.N_prj; i++)
	{
		printf("%d\n", i);
		// DRR 영상의 디바이스 메모리 주소 리턴
		d_prj = get_device_DRR_image(etri_prj, geo[i]);
		tmp = get_image_device_to_host(etri_prj);
		sprintf_s(save_buf,"E:\\SEC\\13\\drr\\%d.raw", i);
		fopen_s(&fp1, save_buf, "wb");
		fwrite(tmp, sizeof(float), static_cast<size_t>(params.prj_u) * params.prj_v, fp1);
		fclose(fp1);

		// projeciton 파일명 로딩용
		idx = i;
		cnt = count_number(idx);
		number_len = params.prj_file_digit;

		if (number_len < cnt) number_len = cnt;


		for (j = 0; j < number_len; j++) {
			if (j < cnt) {
				temp_number[(number_len - 1) - j] = (idx % 10) + '0';
				idx /= 10;
			}
			else
				temp_number[(number_len - 1) - j] = '0';
		}

		if (number_len >= 0 && number_len < 100)
			temp_number[number_len] = '\0';
		else
			temp_number[0] = '\0';

		sprintf_s(buf, "%s\\%s%s.%s", proj_path, params.prj_file_prefix, temp_number, params.prj_file_ext);
		load_projection_image(buf, img, &params);

		// host image to device
		set_device_MI_images(etri_mi, d_prj, img);

		// mi 계산
		get_infomation_MI(etri_mi, &trans_u, &trans_v, &rot);
		gcm[i * 3 + 0] = trans_u;
		gcm[(i * 3) + 1] = trans_v;
		gcm[(i * 3) + 2] = rot;
	}

	// projection image loading에 사용된 변수 메모리 해제
	cudaFreeHost(img);

	convertedtrans2PM(gcm, params.pm, params.N_prj, params.prj_u, params.prj_v, params.pitch_u);
	write_PM(params.output_pm, params.pm, params.N_prj, params.reverse_gantry_angle);
	write_trans("E:\\SEC\\13\\gcm01.txt", gcm, params.N_prj);

	end1 = clock();
	printf("run-time : %.3f s \n", (float)(end1 - start1) / CLOCKS_PER_SEC);

	free(gcm);

	for (i = 0; i < params.N_prj; i++)
		free(geo[i]);
	free(geo);

	return 0;
}

int count_number(int val) {
	int cnt = 0;
	while (val > 0) {
		cnt++;
		val /= 10;
	}
	return cnt;
}


void init_mask(SPARAMS* params, float limit )
{
	float cy, cx;
	float r, r_y, r_x;
	unsigned char* mask;
	int cnt = 0;
	int x, y;
	float dis;

	cx = params->prj_u / 2.f;
	cy = params->prj_v / 2.f;

	/* 영상에 포함되는 최대 반지름*/
	r_x = (cx < (params->prj_u - cx)) ? cx : (params->prj_u - cx);
	r_y = (cy < (params->prj_v - cy)) ? cy : (params->prj_v - cy);

	/* r 의 반지름을 갖는 원이 mask가 됨*/
	r = (r_x < r_y) ? r_x : r_y;
	r *= limit;

	params->roi_x = int(cx - r);
	params->roi_y = int(cy - r);
	params->roi_w = int(2 * r);
	params->roi_h = int(2 * r);

	mask = (unsigned char*)malloc(sizeof(unsigned char) * params->prj_u * params->prj_v);
	memset(mask, 0, sizeof(unsigned char) * params->prj_u * params->prj_v);

	//mask 생성
	for (y = params->roi_y; y < params->roi_y + params->roi_h; y++)
	{
		for (x = params->roi_x; x < params->roi_x + params->roi_w; x++)
		{
			dis = sqrt(((cx - x) * (cx - x)) + ((cy - y) * (cy - y)));

			if (dis < r) {
				mask[(y * params->prj_u) + x] = 1;
				cnt++;
			}
		}
	}

	params->mask_sum = cnt;
	params->mask = mask;
}


bool load_projection_image(char* fname, unsigned short *img, SPARAMS* params)
{
	//unsigned short* img;
	size_t img_size;
	FILE* fp = NULL;
	bool bRtn = true;

	img_size = static_cast<size_t>(params->prj_u) * params->prj_v;

	fopen_s(&fp, fname, "rb");
	if (fp == NULL)
	{
		printf("file open error: %s\n", fname);
		exit(0);
	}
	//img = (unsigned short*)malloc(sizeof(unsigned short) * img_size);

	fread(img, sizeof(unsigned short), img_size, fp);

	fclose(fp);

	return bRtn;

}
float **create_drr_info(SPARAMS* params)
{
	float** geo = NULL;
	int direction = -1;
	float sid, sod_00, sod_01, angle;
	float ppd = 0.0;
	int i;

	if (params->method == GEO_METHOD::PM)
	{
		printf("PM mode\n");

		params->pm = loadProjectionMatrix(params->buf_pm_file);
		geo = convertedPM2Geo(params->pm, params->N_prj, params->pitch_u);
		params->elevation = (geo[0][3] * 180.) / PI;
		params->sod = geo[0][1];
		params->sid = geo[0][6];
	}
	else
	{
		printf("Numeric mode\n");
		params->pm = makeProjectionMatrix(params);
		
		angle = params->elevation * (PI / 180.);

		sod_00 = params->sod;
		sod_01 = 0.;
		sid = params->sid;
		ppd = 0.0;

		if (params->type == GEO_TYPE::PLANARCT)
		{
			sod_00 = params->sod * cos(angle);
			sod_01 = params->sod * sin(angle);
			sid = params->sid * cos(angle);
			ppd = params->sid * sin(angle);
			angle = -90 * (PI / 180.);
		}

		geo = (float**)malloc(sizeof(float*) * params->N_prj);
		if (geo == NULL)return NULL;

		for (i = 0; i < params->N_prj; i++)
		{
#pragma warning(disable:6386 6385)
			geo[i] = (float*)malloc(sizeof(float) * 9);
			if (geo[i] == NULL)
			{
				free(geo);
				return NULL;
			}
			//trans
			geo[i][0] = 0.;
			geo[i][1] = sod_00;
			geo[i][2] = sod_01;

			//rot
			geo[i][3] = angle;
			geo[i][4] = 0.;

 			if (params->reverse_gantry_angle != 1) direction = 1;
			geo[i][5] = float(i * (PI / 180.) * 360. / params->N_prj * direction);

			// sid mm
			geo[i][6] = sid;
			//u0 mm
			geo[i][7] = ((params->prj_u / 2.) + params->offset_u) * params->pitch_u;;
			//v0 mm 
			geo[i][8] = ((params->prj_v / 2.) + params->offset_v) * params->pitch_v + ppd;
			
		}
	}
	return geo;
}

void print_error(const char* meg)
{

#ifdef __DEBUG__
	printf("[ERROR] %s\n", meg);
#endif

}
float* load_volume(SPARAMS* params, VOL_TYPE type)
{
	unsigned short* vol_uint16 = NULL;
	float* vol_float32 = NULL;
	size_t vol_size;
	int i, j, t;
	FILE* fp = NULL;
	cudaError_t error;

	vol_size = static_cast<size_t>(params->vol_x) * params->vol_y * params->vol_z;

	if(_access(params->buf_vol_file, 4))
	{
		print_error("volume image loading error");
		return NULL;
	}

	error = cudaHostAlloc((void **)&vol_float32, sizeof(float) * vol_size, 0);

	if (vol_float32 == NULL)
		return NULL;
	for (size_t idx = 0; idx < vol_size; idx++)
		vol_float32[idx] = 0.f;

	fopen_s(&fp, params->buf_vol_file, "rb");
	if (fp == NULL)
	{
		print_error("volume file open error");
		return NULL;
	}
	if (type == VOL_TYPE::UINT16)
	{
		vol_uint16 = (unsigned short*)malloc(sizeof(unsigned short) * vol_size);
		if (vol_uint16 == NULL)
		{
			fclose(fp);
			return NULL;
		}
		fread(vol_uint16, sizeof(unsigned short), vol_size, fp);
		for (i = 0; i < params->vol_z; i++)
			//if (i > 124)// && i < 205)
			{
				for (j = 0; j < params->vol_y; j++)
					for (t = 0; t < params->vol_x; t++)
#pragma warning(disable:6385 6386)
						vol_float32[(i * params->vol_y * params->vol_x) + (j * params->vol_x) + t] = (float)(vol_uint16[(i * params->vol_y * params->vol_x) + (j * params->vol_x) + t]);
			}
		free(vol_uint16);
		
	}
	else
		fread(vol_float32, sizeof(float), vol_size, fp);

	fclose(fp);

	return vol_float32;
}

void create_params_MI(SPARAMS* src, SPARAMS_MI* dst)
{
	dst->w = src->prj_u;
	dst->h = src->prj_v;
	dst->trans_max = src->trans_max;
	dst->trans_min = src->trans_min;
	dst->hist_bins = src->hist_bin;
	dst->iter_bound = src->iter_bound;

	dst->pitch_x = src->pitch_u;
	dst->pitch_y = src->pitch_v;

	dst->prj_bits_count = src->prj_bits_count;
	dst->roi_x = src->roi_x;
	dst->roi_y = src->roi_y;
	dst->roi_w = src->roi_w;
	dst->roi_h = src->roi_h;
	dst->mask = src->mask;
	dst->mask_sum = src->mask_sum;
}
void create_params_DRR(SPARAMS* src, SPARAMS_DRR* dst)
{
	dst->sid = src->sid;
	dst->sod = src->sod;
	dst->depth_offset = src->depth_offset;
	dst->elevation = src->elevation;
	dst->reverse_gantry_angle = src->reverse_gantry_angle;

	//infomation for volume
	dst->vol_x = src->vol_x;
	dst->vol_y = src->vol_y;
	dst->vol_z = src->vol_z;

	dst->pitch_x = src->pitch_x;
	dst->pitch_y = src->pitch_y;
	dst->pitch_z = src->pitch_z;
	dst->origin_x = src->origin_x;
	dst->origin_y = src->origin_y;
	dst->origin_z = src->origin_z;

	dst->shink = SHRINK_ON_OFF::ON;
	dst->shrink_size = 512;

	//infomation for projection image
	dst->prj_u = src->prj_u;
	dst->prj_v = src->prj_v;
	dst->pitch_u = src->pitch_u;
	dst->pitch_v = src->pitch_v;

	dst->roi_x = src->roi_x;
	dst->roi_y = src->roi_y;
	dst->roi_w = src->roi_w;
	dst->roi_h = src->roi_h;

	dst->offset_u = src->offset_u;
	dst->offset_v = src->offset_v;
	dst->prj_bits_count = src->prj_bits_count;
	
}
bool set_parameter(char* fname, SPARAMS* params)
{
	FILE* fp = NULL;
	char buf[128];

	if (_access(fname, 4))
	{
		print_error("parameter file access error");

		return false;
	}

	fopen_s(&fp, fname, "r");
	if (fp == NULL) {
		print_error("file open error");

		return false;
	}

	if (getProfileString(fp, "Geometry", "method", buf, 128))
		return false;
	if (strncmp(buf, "NUMERIC", strlen(buf)) == 0)
		params->method = GEO_METHOD::NUMERIC;
	else
		params->method = GEO_METHOD::PM;

	if (getProfileString(fp, "Geometry", "type", buf, 128))
		params->type = GEO_TYPE::CT;
	else
	{
		if (strncmp(buf, "OBLIQUECT", strlen(buf)) == 0)
			params->type = GEO_TYPE::OBLIQUECT;
		else if (strncmp(buf, "PLANARCT", strlen(buf)) == 0)
			params->type = GEO_TYPE::PLANARCT;
		else
			params->type = GEO_TYPE::CT;
	}

	if (params->method == GEO_METHOD::NUMERIC)
	{
		if (getIniFloat(fp, "Geometry", "elevation", &params->elevation))
			params->elevation = 0.0;

		if (getIniFloat(fp, "Geometry", "SOD", &params->sod))
		{
			print_error("the SOD value does not exist: ex) SOD = 620.0");
			return false;
		}
		if (getIniFloat(fp, "Geometry", "SID", &params->sid)) {
			print_error("the SID value does not exist: ex) SID = 10.0");
			return false;
		}
	}
	else
	{
		if (getProfileString(fp, "Geometry", "pm_file_name", params->buf_pm_file, MAX_LENGTH))
		{
			print_error("the pm_file_name value does not exist: ex) pm_file_name = PrjAngle.txt");
			return false;
		}
	}


	if (getIniFloat(fp, "Geometry", "depthOffset", &params->depth_offset))
		params->depth_offset = 0.0;
	
	int nReverse;
	if (!getIniInt(fp, "Geometry", "reverseAngle", &nReverse))
		params->reverse_gantry_angle = nReverse == 1 ? 1 : 0;
	else
		params->reverse_gantry_angle = 0;


	// Volume
	// volume size
	if (getProfileString(fp, "Volume", "vol_file_name", params->buf_vol_file, MAX_LENGTH)) {
		print_error("the vol_file_name value does not exist: ex) vol_file_name = vol.raw");
		return false;
	}
	if (getIniInt(fp, "Volume", "N_x", &params->vol_x))
	{
		print_error("the volume size X does not exist: ex) N_x = 256");
		return false;
	}
	if (getIniInt(fp, "Volume", "N_y", &params->vol_y))
	{
		print_error("the volume size Y does not exist: ex) N_y = 256");
		return false;
	}
	if (getIniInt(fp, "Volume", "N_z", &params->vol_z))
	{
		print_error("the volume size Z does not exist: ex) N_z = 256");
		return false;
	}
	
	//voxel space
	if (getIniFloat(fp, "Volume", "dx", &params->pitch_x))
	{
		print_error("the volume pitch X does not exist: ex) dx = 0.5");
		return false;
	}
	if (getIniFloat(fp, "Volume", "dy", &params->pitch_y))
	{
		print_error("the volume pitch Y does not exist: ex) dy = 0.5");
		return false;
	}
	if (getIniFloat(fp, "Volume", "dz", &params->pitch_z))
	{
		print_error("the volume pitch Z does not exist: ex) dz = 0.5");
		return false;
	}
	
	//location of zero index
	if (getIniFloat(fp, "Volume", "X_origin", &params->origin_x))
		params->origin_x = -(float(params->vol_x) * float(params->pitch_x)) / 2.f;
	if (getIniFloat(fp, "Volume", "Y_origin", &params->origin_y))
		params->origin_y = -(float(params->vol_y) * float(params->pitch_y)) / 2.f;
	if (getIniFloat(fp, "Volume", "Z_origin", &params->origin_z))
		params->origin_z = -(float(params->vol_z) * float(params->pitch_z)) / 2.f;


	if (getIniInt(fp, "Detector", "N_u", &params->prj_u))
	{
		print_error("the detector size u does not exist: ex) N_u = 256");
		return false;
	}
	if (getIniInt(fp, "Detector", "N_v", &params->prj_v))
	{
		print_error("the detector size v does not exist: ex) N_v = 256");
		return false;
	}

	//voxel space
	if (getIniFloat(fp, "Detector", "du", &params->pitch_u))
	{
		print_error("the detector pitch u does not exist: ex) dx = 0.5");
		return false;
	}
	if (getIniFloat(fp, "Detector", "dv", &params->pitch_v))
	{
		print_error("the detector pitch v does not exist: ex) dy = 0.5");
		return false;
	}

	if (params->method == GEO_METHOD::NUMERIC)
	{
		if (getIniFloat(fp, "Detector", "offset_u", &params->offset_u))
			params->offset_u = 0.0f;
		if (getIniFloat(fp, "Detector", "offset_v", &params->offset_v))
			params->offset_v = 0.0f;
	}
	if (getIniInt(fp, "Detector", "detector_image_bits", &params->prj_bits_count))
		params->prj_bits_count = 16;


	if (getIniInt(fp, "Gcm", "iter_bound", &params->iter_bound))
		params->iter_bound = 50;
	if (getIniInt(fp, "Gcm", "hist_bin", &params->hist_bin))
		params->hist_bin = 128;

	if (getIniFloat(fp, "Gcm", "trans_max", &params->trans_max))
		params->trans_max = 0.2f;
	if (getIniFloat(fp, "Gcm", "trans_min", &params->trans_min))
		params->trans_min = 0.025f;
	if (getProfileString(fp, "Gcm", "output_PM", params->output_pm, MAX_LENGTH))
	{
		print_error("the file name for new PM  does not exist: ex) output_PM = E:\\test.txt");
		return false;
	}

	//Projection
	//file path
	if (getProfileString(fp, "Projection", "File_path", params->prj_file_path, MAX_LENGTH))
	{
		print_error("the projection file path does not exist: ex) File_path = c:\\projections\\");
		return false;
	}
	if (getIniInt(fp, "Projection", "N_proj", &params->N_prj))
	{
		print_error("the number of prjection image does not exist: ex) N_prj = 360");
		return false;
	}
	//file명
	if (getProfileString(fp, "Projection", "File_prefix", params->prj_file_prefix, 128))
	{
		print_error("the file name prefix for projection images does not exist: ex) File_prefix = ViewImage");
		return false;
	}
	if (getProfileString(fp, "Projection", "File_ext", params->prj_file_ext, 10))
	{

		print_error("the file name extenstion for projection images does not exist: ex) File_ext = raw");
		return false;
	}
	if (getIniInt(fp, "Projection", "File_digit", &params->prj_file_digit))
	{

		print_error("the file digit for projection images does not exist: ex) File_digit = 4");
		return false;
	}

	fclose(fp);

	return true;
}
