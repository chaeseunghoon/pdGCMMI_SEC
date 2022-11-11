#pragma once
#pragma warning(disable:4244)

#define MAX_LENGTH 1024

enum class GEO_METHOD
{
	PM = 0,
	NUMERIC = 1
};

enum class GEO_TYPE
{
	CT = 0,
	OBLIQUECT = 1,
	PLANARCT = 2
};

typedef struct __struct_params__
{
	//infomation for initial geometry
	GEO_METHOD method;
	GEO_TYPE type;
	float sid, sod;
	float depth_offset;
	float elevation;
	int reverse_gantry_angle;
	char buf_pm_file[MAX_LENGTH];

	//infomation for volume
	char buf_vol_file[MAX_LENGTH];
	int vol_x, vol_y, vol_z;
	float pitch_x, pitch_y, pitch_z;
	float origin_x, origin_y, origin_z;

	//information for GCM
	float trans_max;
	float trans_min;
	int hist_bin;
	int iter_bound;
	char output_pm[MAX_LENGTH];

	//infomation for projection image
	int prj_u, prj_v;
	float pitch_u, pitch_v;
	float offset_u, offset_v;
	int N_prj;
	int prj_bits_count;
	char prj_file_path[MAX_LENGTH];
	char prj_file_prefix[128];
	char prj_file_ext[10];
	int prj_file_digit;

	int roi_x, roi_y, roi_w, roi_h;
	unsigned char* mask;
	int mask_sum;

	float* pm;

}SPARAMS;


int getProfileString(FILE* fp, const char* section, const char* key, char* ret, int size);
int seekSection(FILE* fp, const char* section);
char* trimString(char* str);
bool isBlank(char ch);
bool isComment(char ch);
int getNextKey(FILE* fp, char* retKey, char* retValue);
int getIniInt(FILE* fp, const char* section, const char* key, int* ret);
int getIniFloat(FILE* fp, const char* section, const char* key, float* ret);
