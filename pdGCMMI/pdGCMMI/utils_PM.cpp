#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <cstdlib> //abort()
#include "utils_PM.h"

float* loadProjectionMatrix(char* fName)
{
	FILE* fp = NULL;
	char ch;
	int ready = 0;
	int cnt;
	char buf[256];

	float* pm;

	double tmp;
	int dump_int;
	int num = 0;
	int stop = 0;
	int i;

	char* pStr;

	fopen_s(&fp, fName, "rt");
	if (fp == NULL)
	{
		printf("do not find projection file : %s\n", fName);
		return NULL;
	}


	while (1) {
		ch = fgetc(fp);

		if (ch == EOF || stop == 1) break;
		if (ch == '#')
		{
			while (1) {
				ch = fgetc(fp);
				if (ch == EOF)
				{
					stop = 1;
					break;
				}
				if (ch == '\n') break;
			}

			pStr = fgets(buf, 256, fp);
			if (pStr == NULL)
				return NULL;

			num = atoi(buf);
			break;
		}

	}
	if (num == 0)
	{
		printf("projection matrix error\n");
		fclose(fp);
		return NULL;
	}

	pm = (float*)malloc(sizeof(float) * 12 * num);

	cnt = 0;

	while (1) {
		if (ready != 2)
		{
			ch = fgetc(fp);
			if (ch == EOF) break;
			if (ch == '@')
				ready = 1;
			if (ready == 1 && ch == ' ')
				ready = 2;
		}
		else
		{
			fscanf_s(fp, "%d", &dump_int);
			fscanf_s(fp, "%lf", &tmp);
			fscanf_s(fp, "%lf", &tmp);

			for (i = 0; i < 3; i++)
			{
				fscanf_s(fp, "%f %f %f %f", &pm[(cnt * 12) + (i * 4)], &pm[(cnt * 12) + (i * 4) + 1], &pm[(cnt * 12) + (i * 4) + 2], &pm[(cnt * 12) + (i * 4) + 3]);
			}

			ready = 0;
			cnt++;
			if (cnt >= num)
				break;
		}

	}
	fclose(fp);

	return pm;
}


float* makeProjectionMatrix(SPARAMS* params)
{
	int idx;
	float view_angle;
	float angle_step;
	float* pm = NULL;
	float sc_off[3][3] = {
		{params->sid / params->pitch_u, -(params->prj_u / 2.f + params->offset_u), 0},
		{0, -(params->prj_v / 2.f + params->offset_v), params->sid / params->pitch_v},
		{0, -1, 0}
	};
	float Rx[3][3] = {
		{1,0,0},
		{0,1,0},
		{0,0,1}
	};
	float Rz[3][3] = {
		{1,0,0},
		{0,1,0},
		{0,0,1}
	};
	float R[3][4];
	float tilt = 0.0;
	float ppd;
	float v0;

	float T[3] = { 0, -params->sod, 0 };
	int direction = 1;
	if (params->reverse_gantry_angle == 1) direction = -1;

	angle_step = (360.f / params->N_prj) / 180.f * PI * direction;

	pm = (float*)malloc(sizeof(float) * 12 * params->N_prj);
	if (pm == NULL)
		return NULL;

	float rad = params->elevation * PI / 180.f;

	if (params->type == GEO_TYPE::PLANARCT)
	{
		tilt = rad;
		rad = -90.f * PI / 180.f;
		T[1] = -params->sod * cos(tilt);
		T[2] = -params->sod * sin(tilt);

		ppd = params->sid * sin(tilt) / params->pitch_v;
		v0 = params->prj_v / 2.f + params->offset_v + ppd;

		sc_off[0][0] = params->sid * cos(tilt) / params->pitch_u;
		sc_off[1][1] = -v0;
		sc_off[1][2] = params->sid * cos(tilt) / params->pitch_u;
	}
	Rx[1][1] = cos(rad);
	Rx[1][2] = sin(rad);
	Rx[2][1] = -sin(rad);
	Rx[2][2] = cos(rad);


	for (idx = 0; idx < params->N_prj; idx++)
	{
		view_angle = idx * angle_step;
		Rz[0][0] = cos(view_angle);
		Rz[0][1] = sin(view_angle);
		Rz[1][0] = -sin(view_angle);
		Rz[1][1] = cos(view_angle);
		matrixMultiplication(Rx, Rz, R);

		R[0][3] = T[0];
		R[1][3] = T[1];
		R[2][3] = T[2];

		genPM(sc_off, R, &pm[idx * 12]);
	}

	return pm;
}


/**
*@fn void matrixMultiplication(float a[][3], float b[][3], float c[][4])
* @brief PM을 만들기 위한 object회전에 관련된 행렬 곱
* @details b를 먼저 회전하고 a를 회전 진행, 결과는 3x4 행렬 c에 저장되며, c의 4열은 object의 이동값을 입력하는 공간으로 행렬곱에서는 사용하지 않음
* @param a x 축에 대한 회전 행렬(3x3)
* @param b z 축에 대한 회전 행렬(3x3)
* @param c 행렬곱 저장 변수 3x4 배열로 구성
*/
void matrixMultiplication(float a[][3], float b[][3], float c[][4])
{
	int i, j, k;
	float sum;

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			sum = 0.0;
			for (k = 0; k < 3; k++)
				sum += a[i][k] * b[k][j];
			c[i][j] = sum;
		}
	}

}

/**
*@fn void genPM(float a[][3], float b[][4], float* pm)
* @brief 최종 3x4크기의 PM을 만드는 함수
* @param a PM을 Geometry정보가 저장된 행렬
* @param b object의 회전, 이동정보과 저장된 행렬
* @param pm 한 view angle의 PM이 저장되는 변수
*/
void genPM(float a[][3], float b[][4], float* pm)
{
	int i, j, k;
	float sum;
	int cnt = 0;

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 4; j++)
		{
			sum = 0.0;
			for (k = 0; k < 3; k++)
				sum += a[i][k] * b[k][j];
			pm[cnt++] = sum;
		}
	}
}


/**
*@fn void write_PM(const char* fName, float* pm, int num, int bReverseGantryAngle)
* @brief PM을 파일로 저장하는 함수
* @param fName PM을 저장할 파일명(경로 포함)
* @param pm PM 변수(12* projections 수)
* @param num projections 수
* @param bReverseGantryAngle 회전 방향
*/
void write_PM(const char* fName, float* pm, int num, int bReverseGantryAngle)
{
	FILE* fp = NULL;
	int t, m, n;
	struct tm _t;
	char datetime[256];
	time_t tnow;
	int direction = 1;
	float angle_step;
	float azimuth;

	tnow = time(NULL);
	localtime_s(&_t, &tnow);

	fopen_s(&fp, fName, "wt");
	if (fp == NULL)
	{
		printf("projection matrix save faile: %s\r\n", fName);
		return;
	}

	sprintf_s(datetime, "%04d. %02d. %02d. %02d:%02d", _t.tm_year + 1900, _t.tm_mon + 1, _t.tm_mday, _t.tm_hour, _t.tm_min);

	printf("ETRI projtable version iterative \n%s \n\n # format : angle / entries of projection matrices\n%d\n\n", datetime, num);

	fprintf_s(fp, "ETRI projtable version iterative \n%s \n\n # format : angle / entries of projection matrices\n%d\n\n", datetime, num);


	if (bReverseGantryAngle == 1) direction = -1;

	angle_step = (360. / num) / 180. * PI * direction;

	for (t = 0; t < num; t++)
	{
		azimuth = t * angle_step;

		fprintf(fp, "@ %d\n", t);
		fprintf(fp, "%lf %lf\n", azimuth, 0.0);

		for (n = 0; n < 3; n++) {
			for (m = 0; m < 4; m++) fprintf(fp, "%+11.10e ", pm[(12 * t) + (4 * n) + m]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
	printf(" Projection matrix save : %s\n", fName);

}


void write_Geo(const char* fName, float** geo, int num)
{
	FILE* fp;
	int i, j;

	fopen_s(&fp, fName, "wt");

	fprintf(fp, "trans_x\t trans_y\t trans_z\t rot_x\t rot_y\t rot_z\t SID\t offset_u\t offset_v\n");

	for (i = 0; i < num; i++)
	{
		for (j = 0; j < 8; j++)
			fprintf(fp, "%.3f\t ", geo[i][j]);
		fprintf(fp, "%.3f\n ", geo[i][8]);
	}
	fclose(fp);
}

void normalizationmatrix(float* ori, float pm[][4], int row, int col)
{
	double norm2;
	int i, j;

	norm2 = sqrt((ori[(2 * 4) + 0] * ori[(2 * 4) + 0]) + (ori[(2 * 4) + 1] * ori[(2 * 4) + 1]) + (ori[(2 * 4) + 2] * ori[(2 * 4) + 2]));

	for (i = 0; i < row; i++)
		for (j = 0; j < col; j++)
			pm[i][j] = ori[(i * 4) + j] / norm2;

}

float abscheckverticality(float a[][1], float b[][1], int row, int col)
{
	//a = 3x1 -> 1x3으로 계산
	int i;
	float sum;

	sum = 0.0;

	//a'를 구해야 하지만 단순화 해서 계산
	for (i = 0; i < row; i++)
	{
		sum += (a[i][0] * b[i][0]);
	}

	if (sum < 0)
		sum += -1.;

	return sum;
}


void matrix_element_mul_3x1(float data[][1], float val)
{
	int i, j;

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 1; j++)
		{
			data[i][j] *= val;
		}
	}
}


void matrix_element_mul_3x3(float data[][3], float val)
{
	int i, j;

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			data[i][j] *= val;
		}
	}
}

float** convertedPM2Geo(float* matrix, int num, float pitch)
{
	float pm[3][4], R[3][3], T[3][1], basis1[3][1], basis2[3][1], basis3[3][1];
	float** geo = NULL;

	double radian;
	int i, n;
	double norm2;
	double u0, v0, FaCosParity, FaSin, FbCos, FbSinParity;
	double parity, theta, fa, fb;
	double det_R;
	double cosA, cosB, cosC, sinA, sinB, sinC, argB, argA, argC;

	geo = (float**)malloc(sizeof(float*) * num);
	if (geo == NULL)
		return NULL;

	for (i = 0; i < num; i++)
	{
		geo[i] = (float*)malloc(sizeof(float) * 9);
		if (geo[i] == NULL)
			return NULL;

		normalizationmatrix(&matrix[i * 12], pm, 3, 4);

		for (n = 0; n < 3; n++)
			R[1][n] = -pm[2][n];

		T[1][0] = -pm[2][3];

		for (n = 0; n < 3; n++)
			basis2[n][0] = R[1][n];

		basis3[0][0] = 1;
		basis3[1][0] = 0;
		basis3[2][0] = 0;

		if (abscheckverticality(basis2, basis3, 3, 1) < 0.58)
		{
			basis1[0][0] = (basis3[2][0] * basis2[1][0]) - (basis3[1][0] * basis2[2][0]);
			basis1[1][0] = (basis3[0][0] * basis2[2][0]) - (basis3[2][0] * basis2[0][0]);
			basis1[2][0] = (basis3[1][0] * basis2[0][0]) - (basis3[0][0] * basis2[1][0]);
		}

		basis3[0][0] = 0;
		basis3[1][0] = 1;
		basis3[2][0] = 0;

		if (abscheckverticality(basis2, basis3, 3, 1) < 0.58)
		{
			basis1[0][0] = (basis3[2][0] * basis2[1][0]) - (basis3[1][0] * basis2[2][0]);
			basis1[1][0] = (basis3[0][0] * basis2[2][0]) - (basis3[2][0] * basis2[0][0]);
			basis1[2][0] = (basis3[1][0] * basis2[0][0]) - (basis3[0][0] * basis2[1][0]);
		}

		basis3[0][0] = 0;
		basis3[1][0] = 0;
		basis3[2][0] = 1;

		if (abscheckverticality(basis2, basis3, 3, 1) < 0.58)
		{
			basis1[0][0] = (basis3[2][0] * basis2[1][0]) - (basis3[1][0] * basis2[2][0]);
			basis1[1][0] = (basis3[0][0] * basis2[2][0]) - (basis3[2][0] * basis2[0][0]);
			basis1[2][0] = (basis3[1][0] * basis2[0][0]) - (basis3[0][0] * basis2[1][0]);
		}

		norm2 = sqrt((basis1[0][0] * basis1[0][0])
			+ (basis1[1][0] * basis1[1][0])
			+ (basis1[2][0] * basis1[2][0]));
		basis1[0][0] /= norm2;
		basis1[1][0] /= norm2;
		basis1[2][0] /= norm2;

		basis3[0][0] = (basis2[2][0] * basis1[1][0]) - (basis2[1][0] * basis1[2][0]);
		basis3[1][0] = (basis2[0][0] * basis1[2][0]) - (basis2[2][0] * basis1[0][0]);
		basis3[2][0] = (basis2[1][0] * basis1[0][0]) - (basis2[0][0] * basis1[1][0]);

		u0 = -(pm[0][0] * basis2[0][0] + pm[0][1] * basis2[1][0] + pm[0][2] * basis2[2][0]);
		v0 = -((pm[1][0] * basis2[0][0]) + (pm[1][1] * basis2[1][0]) + (pm[1][2] * basis2[2][0]));
		FaCosParity = pm[0][0] * basis1[0][0] + pm[0][1] * basis1[1][0] + pm[0][2] * basis1[2][0];
		FaSin = -(pm[0][0] * basis3[0][0] + pm[0][1] * basis3[1][0] + pm[0][2] * basis3[2][0]);

		FbCos = pm[1][0] * basis3[0][0] + pm[1][1] * basis3[1][0] + pm[1][2] * basis3[2][0];
		FbSinParity = pm[1][0] * basis1[0][0] + pm[1][1] * basis1[1][0] + pm[1][2] * basis1[2][0];

		if (fabs(FbCos) > 1.0E-5)
		{
			parity = SIGN(FaCosParity / FbCos);
			if (FbCos > 0)
			{
				theta = atan(parity * FbSinParity / FbCos);
				fa = FaCosParity * parity / cos(theta);
				fb = FbCos / cos(theta);
			}
			else
			{
				theta = atan(parity * FbSinParity / FbCos) + PI;
				fa = FaCosParity * parity / cos(theta);
				fb = FbCos / cos(theta);
			}
		}
		else
		{
			parity = SIGN(FbSinParity / FaSin);
			if (FbCos > 0)
			{
				radian = parity * FbCos / FbSinParity;
				//theta = acot(parity*FbCos/FbSinParity);
				theta = atan(1. / radian);

				fa = FaSin / sin(theta);
				fb = FbSinParity * parity / sin(theta);
			}
			else
			{
				radian = parity * FbCos / FbSinParity;
				//theta = acot(parity*FbCos/FbSinParity)+PI;
				theta = atan(1. / radian) + PI;
				fa = FaSin / sin(theta);
				fb = FbSinParity * parity / sin(theta);
			}
		}

		// T(1), T(3) 결정
		T[0][0] = (pm[0][3] + u0 * T[1][0]) / fa;
		T[2][0] = (pm[1][3] + v0 * T[1][0]) / fb;

		// [R] 결정
		R[0][0] = cos(theta) * parity * basis1[0][0] - sin(theta) * basis3[0][0];
		R[0][1] = cos(theta) * parity * basis1[1][0] - sin(theta) * basis3[1][0];
		R[0][2] = cos(theta) * parity * basis1[2][0] - sin(theta) * basis3[2][0];
		R[2][0] = cos(theta) * basis3[0][0] + sin(theta) * parity * basis1[0][0];
		R[2][1] = cos(theta) * basis3[1][0] + sin(theta) * parity * basis1[1][0];
		R[2][2] = cos(theta) * basis3[2][0] + sin(theta) * parity * basis1[2][0];

		// Parity를 항상 1로 만듬
		// 여기부터 theta는 의미없을 듯
		det_R = (R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]))
			- (R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]))
			+ (R[0][2] * (R[1][1] * R[2][1] - R[1][1] * R[2][1]));

		if (det_R < 0)
		{
			matrix_element_mul_3x3(R, -1);
			matrix_element_mul_3x1(T, -1);
		}

		sinB = -R[0][2];
		cosB = sqrt(1 - sinB * sinB);

		sinA = R[1][2] / cosB;
		cosA = R[2][2] / cosB;

		cosC = R[0][0] / cosB;
		sinC = R[0][1] / cosB;

		argB = atan(sinB / cosB);
		if (cosA > 0)
		{
			argA = atan(sinA / cosA);
		}

		else if (cosA < 0)
		{
			if (sinA > 0)
				argA = atan(sinA / cosA) + PI;
			else
				argA = atan(sinA / cosA) - PI;
		}
		else
		{
			if (sinA > 0)
				argA = PI / 2.;
			else
				argA = -PI / 2.;
		}

		if (cosC > 0)
		{
			argC = atan(sinC / cosC);
		}
		else if (cosC < 0)
		{
			if (sinC > 0)
				argC = atan(sinC / cosC) + PI;
			else
				argC = atan(sinC / cosC) - PI;
		}
		else
		{
			if (sinC > 0)
				argC = PI / 2.;
			else
				argC = -PI / 2.;
		}
/*
		geo[i][0] = -T[0][0];
		geo[i][1] = -T[1][0];
		geo[i][2] = -T[2][0];
		geo[i][3] = -argA;
		geo[i][4] = -argB;
		geo[i][5] = -argC;
		geo[i][6] = (fa + fb) / 2. * pitch;
		geo[i][7] = u0 * pitch;
		geo[i][8] = v0 * pitch;
*/
		geo[i][0] = T[0][0];
		geo[i][1] = T[1][0];
		geo[i][2] = T[2][0];
		geo[i][3] = argA;
		geo[i][4] = argB;
		geo[i][5] = argC;
		geo[i][6] = (fa + fb) / 2. * pitch;
		geo[i][7] = u0 * pitch;
		geo[i][8] = v0 * pitch;
		
	}

	return geo;
}


// trans[0] : trans_u, trans[1] : trans_v, trans[2] : angle
void convertedtrans2PM(float* trans, float* old, int num, int ud, int vd, double pitch)
{
	int i;
	float R[2][2], T[2][1], CV[3][3];
	//double** R, ** T, ** CV;
	//double** sub_pm;
	//float sub_pm[3][4];
	float trans_v, trans_u, angle;
	float primary_angle, secondary_angle;
	float uc, vc;

	if (old == NULL || trans == NULL)
	{
		printf("data is NULL\n");
		return;
	}

	uc = ud / 2.;
	vc = vd / 2.;

	for (i = 0; i < num; i++)
	{
		primary_angle = 0.0;
		secondary_angle = 0.0;

		trans_u = -trans[(i * 3) + 0] / pitch;
		trans_v = -trans[(i * 3) + 1] / pitch;
		angle = trans[(i * 3) + 2];


		// Counter clockwise object 회전 (v, u) -> (v', u') :  R (u,v)' = (u',v')'
		R[0][0] = cos(angle);
		R[0][1] = -sin(angle);
		R[1][0] = sin(angle);
		R[1][1] = cos(angle);

		// T = [uc,vc]'-R*([uc,vc]'+[trans_u, trans_v]') 계산
		T[0][0] = uc - (R[0][0] * (uc + trans_u) + R[0][1] * (vc + trans_v));
		T[1][0] = vc - (R[1][0] * (uc + trans_u) + R[1][1] * (vc + trans_v));

		CV[0][0] = R[0][0]; CV[0][1] = R[0][1]; CV[0][2] = T[0][0];
		CV[1][0] = R[1][0]; CV[1][1] = R[1][1]; CV[1][2] = T[1][0];
		CV[2][0] = 0; CV[2][1] = 0; CV[2][2] = 1;

		updateMatrix(CV, &old[i * 12]);
	}
}

void updateMatrix(float a[3][3], float* pm)
{
	int row = 3;
	int col = 4;
	int i, j, k;
	float temp[12];
	float sum;

	memcpy(temp, pm, sizeof(float) * 12);

	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
		{
			sum = 0.0;
			for (k = 0; k < row; k++)
				sum += a[i][k] * temp[(k * 4) + j];
			pm[(i * 4) + j] = sum;
		}
	}
}


void write_trans(const char* fName, float* data, int num)
{
	FILE* fp;
	int i;
	auto errorno = fopen_s(&fp, fName, "wt");

	if (fp != NULL)
	{
		fprintf(fp, "#idx, offset u(mm), offset v(mm), rotation(radian)\n");
		for (i = 0; i < num; i++)
		{
			fprintf(fp, "%d, %lf, %lf, %lf\n", i, data[i * 3], data[(i * 3) + 1], data[(i * 3) + 2]);
		}
		fclose(fp);
	}
	else
	{
		std::abort();
	}

}
