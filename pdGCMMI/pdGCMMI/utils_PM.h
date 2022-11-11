#pragma once
#pragma warning(disable:4244)

#include "utils_params.h"

#define PI 3.14159265359f

#define SIGN(A) (A==0)?0:(A<0)?-1:1

float* loadProjectionMatrix(char* fName);
float* makeProjectionMatrix(SPARAMS* params);

void write_PM(const char* fName, float* pm, int num, int bReverseGantryAngle);
void genPM(float a[][3], float b[][4], float* pm);
void matrixMultiplication(float a[][3], float b[][3], float c[][4]);
void normalizationmatrix(float* ori, float pm[][4], int row, int col);
float abscheckverticality(float a[][1], float b[][1], int row, int col);
void matrix_element_mul_3x1(float data[][1], float val);
void matrix_element_mul_3x3(float data[][3], float val);
float** convertedPM2Geo(float* matrix, int num, float pitch);


void updateMatrix(float a[3][3], float* pm);
// trans[0] : trans_u, trans[1] : trans_v, trans[2] : angle
void convertedtrans2PM(float* trans, float* old, int num, int ud, int vd, double pitch);
void write_trans(const char* fName, float* data, int num);
void write_Geo(const char* fName, float** geo, int num);