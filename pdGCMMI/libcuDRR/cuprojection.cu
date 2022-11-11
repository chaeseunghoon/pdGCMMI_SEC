#include "cuprojection.cuh"
#include <math.h>

bool etri_set_common_params(SPARAMS_DRR *params)
{
    bool bRtn = true;
    
	cudaMemcpyToSymbol(C_PRJ_U, &(params->prj_u), sizeof(int));
	cudaMemcpyToSymbol(C_PRJ_V, &(params->prj_v), sizeof(int));

	cudaMemcpyToSymbol(C_VOL_X, &(params->vol_x), sizeof(int));
	cudaMemcpyToSymbol(C_VOL_Y, &(params->vol_y), sizeof(int));
	cudaMemcpyToSymbol(C_VOL_Z, &(params->vol_z), sizeof(int));

	cudaMemcpyToSymbol(C_PITCH_X, &(params->pitch_x), sizeof(float));
	cudaMemcpyToSymbol(C_PITCH_Y, &(params->pitch_y), sizeof(float));
	cudaMemcpyToSymbol(C_PITCH_Z, &(params->pitch_z), sizeof(float));

    cudaMemcpyToSymbol(C_ORIGIN_X, &(params->origin_x), sizeof(float));
    cudaMemcpyToSymbol(C_ORIGIN_Y, &(params->origin_y), sizeof(float));
    cudaMemcpyToSymbol(C_ORIGIN_Z, &(params->origin_z), sizeof(float));

	cudaMemcpyToSymbol(C_ROI_X, &(params->roi_x), sizeof(int));
	cudaMemcpyToSymbol(C_ROI_Y, &(params->roi_y), sizeof(int));

	return bRtn;
}

int etri_set_individual_params(SPARAMS_DRR* params, float* trans)
{
    float Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz;
    float t_xray_x, t_xray_y, t_xray_z; //임시
    float t_detector_origin_x, t_detector_origin_y, t_detector_origin_z; // 임시
    //디텍터의 각 방향 픽셀 pitch의 벡터값
    float detector_du_x, detector_du_y, detector_du_z, detector_dv_x, detector_dv_y, detector_dv_z; 
    // 소스의 회전/이동 위치
    float xray_x, xray_y, xray_z;
    // 디텍터 중심의 회전/이동 위치
    float detector_origin_x, detector_origin_y, detector_origin_z;
    
    // ray path 방향 결정
    float x_m, y_m, z_m;
    int dir = 0;
    float diff_x, diff_y, diff_z;

    t_xray_x = t_xray_y = t_xray_z = 0.0; //디텍터 (0,0)좌표의 위치
    t_detector_origin_x = -trans[7];// *pParam->du;
    t_detector_origin_y = -trans[6]; // SID mm단위로 입력
    t_detector_origin_z = -trans[8];// *pParam->dv;

    // 소스-디텍터 평행이동
    //t_xray_x += trans[0]; t_xray_y += trans[1]; t_xray_z += trans[2];
    //t_detector_origin_x += trans[0]; t_detector_origin_y += trans[1]; t_detector_origin_z += trans[2];

    t_xray_x -= trans[0]; t_xray_y -= trans[1]; t_xray_z -= trans[2];
    t_detector_origin_x -= trans[0]; t_detector_origin_y -= trans[1]; t_detector_origin_z -= trans[2];

    // 소스-디텍터 회전행렬(Z<-Y<-X 순서)
    Rxx = cos(trans[5]) * cos(trans[4]);
    Rxy = cos(trans[5]) * sin(trans[4]) * sin(trans[3]) - sin(trans[5]) * cos(trans[3]);
    Rxz = cos(trans[5]) * sin(trans[4]) * cos(trans[3]) + sin(trans[5]) * sin(trans[3]);

    Ryx = sin(trans[5]) * cos(trans[4]);
    Ryy = sin(trans[5]) * sin(trans[4]) * sin(trans[3]) + cos(trans[5]) * cos(trans[3]);
    Ryz = sin(trans[5]) * sin(trans[4]) * cos(trans[3]) - cos(trans[5]) * sin(trans[3]);

    Rzx = -sin(trans[4]);
    Rzy = cos(trans[4]) * sin(trans[3]);
    Rzz = cos(trans[4]) * cos(trans[3]);

    //소스의 변화된 위치
    xray_x = Rxx * t_xray_x + Rxy * t_xray_y + Rxz * t_xray_z;
    xray_y = Ryx * t_xray_x + Ryy * t_xray_y + Ryz * t_xray_z;
    xray_z = Rzx * t_xray_x + Rzy * t_xray_y + Rzz * t_xray_z;

    //디텍터 origin의 변화된 위치
    detector_origin_x = Rxx * t_detector_origin_x + Rxy * t_detector_origin_y + Rxz * t_detector_origin_z;
    detector_origin_y = Ryx * t_detector_origin_x + Ryy * t_detector_origin_y + Ryz * t_detector_origin_z;
    detector_origin_z = Rzx * t_detector_origin_x + Rzy * t_detector_origin_y + Rzz * t_detector_origin_z;

    //디턱터상의 u-방향 단위벡터
    detector_du_x = Rxx * params->pitch_u;
    detector_du_y = Ryx * params->pitch_u;
    detector_du_z = Rzx * params->pitch_u;

    //디턱터상의 v-방향 단위벡터
    detector_dv_x = Rxz * params->pitch_v;
    detector_dv_y = Ryz * params->pitch_v;
    detector_dv_z = Rzz * params->pitch_v;

    cudaMemcpyToSymbol(C_SOURCE_POS_X, &(xray_x), sizeof(float));
    cudaMemcpyToSymbol(C_SOURCE_POS_Y, &(xray_y), sizeof(float));
    cudaMemcpyToSymbol(C_SOURCE_POS_Z, &(xray_z), sizeof(float));

    cudaMemcpyToSymbol(C_DETECTOR_ORI_X, &(detector_origin_x), sizeof(float));
    cudaMemcpyToSymbol(C_DETECTOR_ORI_Y, &(detector_origin_y), sizeof(float));
    cudaMemcpyToSymbol(C_DETECTOR_ORI_Z, &(detector_origin_z), sizeof(float));

    cudaMemcpyToSymbol(C_DETECTOR_U_X, &(detector_du_x), sizeof(float));
    cudaMemcpyToSymbol(C_DETECTOR_U_Y, &(detector_du_y), sizeof(float));
    cudaMemcpyToSymbol(C_DETECTOR_U_Z, &(detector_du_z), sizeof(float));

    cudaMemcpyToSymbol(C_DETECTOR_V_X, &(detector_dv_x), sizeof(float));
    cudaMemcpyToSymbol(C_DETECTOR_V_Y, &(detector_dv_y), sizeof(float));
    cudaMemcpyToSymbol(C_DETECTOR_V_Z, &(detector_dv_z), sizeof(float));



    x_m = detector_origin_x + (params->prj_v / 2.0) * detector_dv_x + (params->prj_u / 2.0) * detector_du_x;
    y_m = detector_origin_y + (params->prj_v / 2.0) * detector_dv_y + (params->prj_u / 2.0) * detector_du_y;
    z_m = detector_origin_z + (params->prj_v / 2.0) * detector_dv_z + (params->prj_u / 2.0) * detector_du_z;


    //소스-디텍터 방향(ray)과 최소각을 갖는 축 찾기용
    diff_x = fabs(xray_x - x_m);
    diff_y = fabs(xray_y - y_m);
    diff_z = fabs(xray_z - z_m);

    //raysum 방향 결정
    if (diff_x > diff_y) {
        if (diff_x > diff_z) dir = 1; else dir = 3;
    }
    else {
        if (diff_y > diff_z) dir = 2; else dir = 3;
    }

    return dir;
}


__global__ void kernelEtriDRRplaneX(float* pfPrjs, float* pVol)
{
    int u = blockDim.x * blockIdx.x + threadIdx.x + C_ROI_X;
    int v = blockDim.y * blockIdx.y + threadIdx.y + C_ROI_Y;

    float ray_start_x, ray_end_x, ray_start_y, ray_end_y, ray_start_z, ray_end_z, ray_temp;
    float x_m, y_m, z_m;
    int plane_start, plane_end;
    float raysum = 0.0f;
    float ratio_k;
    int x;
    float ray_y, ray_z;
    float diff_x, diff_y, diff_z;

    if ((u < C_PRJ_U) && (v < C_PRJ_U))
    {
        //디텍터 픽셀의 위치 계산
        x_m = C_DETECTOR_ORI_X + v * C_DETECTOR_V_X + u * C_DETECTOR_U_X;
        y_m = C_DETECTOR_ORI_Y + v * C_DETECTOR_V_Y + u * C_DETECTOR_U_Y;
        z_m = C_DETECTOR_ORI_Z + v * C_DETECTOR_V_Z + u * C_DETECTOR_U_Z;

        //소스-디텍터 픽셀 사이의 raysum 계산
        //소스와 디텍터 위치를 object 인덱스로 변경 (double 값)
        /*
        ray_start_x = x_m / C_PITCH_X + (int)(C_VOL_X / 2.f);
        ray_end_x = C_SOURCE_POS_X / C_PITCH_X + (int)(C_VOL_X / 2.f);
        ray_start_y = y_m / C_PITCH_Y + (int)(C_VOL_Y / 2.f);
        ray_end_y = C_SOURCE_POS_Y / C_PITCH_Y + (int)(C_VOL_Y / 2.f);
        ray_start_z = z_m / C_PITCH_Z + (int)(C_VOL_Z / 2.f);
        ray_end_z = C_SOURCE_POS_Z/ C_PITCH_Z + (int)(C_VOL_Z / 2.f);
        */
        ray_start_x = (int)((x_m - C_ORIGIN_X) / C_PITCH_X);
        ray_end_x = (int)((C_SOURCE_POS_X - C_ORIGIN_X)/ C_PITCH_X);
        ray_start_y = (int)((y_m - C_ORIGIN_Y) / C_PITCH_Y);
        ray_end_y = (int)((C_SOURCE_POS_Y - C_ORIGIN_Y) / C_PITCH_Y);
        ray_start_z = (int)((z_m - C_ORIGIN_Z) / C_PITCH_Z);
        ray_end_z = (int)((C_SOURCE_POS_Z - C_ORIGIN_Z) / C_PITCH_Z);

        diff_x = fabs(C_SOURCE_POS_X - x_m);
        diff_y = fabs(C_SOURCE_POS_Y - y_m);
        diff_z = fabs(C_SOURCE_POS_Z - z_m);

        //시작과 끝 재조정
        if (ray_start_x > ray_end_x) {

            ray_temp = ray_start_x;
            ray_start_x = ray_end_x;
            ray_end_x = ray_temp;

            ray_temp = ray_start_y;
            ray_start_y = ray_end_y;
            ray_end_y = ray_temp;

            ray_temp = ray_start_z;
            ray_start_z = ray_end_z;
            ray_end_z = ray_temp;
        }

        //object 시작점 및 끝점 계산
        if (ray_start_x < 0.0) plane_start = 0;
        else plane_start = (int)ray_start_x + 1; // 강제로 +1 ; 어떤 경우에도 division by zero 해소
        if (ray_end_x > (C_VOL_X - 1)) plane_end = (C_VOL_X - 1);
        else plane_end = (int)ray_end_x - 1; //강제로 -1

        //object 시작에서 끝까지 raysum
        for (x = plane_start; x <= plane_end; x++) {

            //ray와 x-번째 plane의 교차점을 비례식으로 계산
            ratio_k = (ray_end_x - x) / (x - ray_start_x); // division by zero는 일어나지 않음
            ray_y = (ray_end_y + ratio_k * ray_start_y) / (1 + ratio_k);
            ray_z = (ray_end_z + ratio_k * ray_start_z) / (1 + ratio_k);

            //y-z 평면에서 interpolation을 위한 주변 4픽셀 위치 및 weight 계산
            if (((x - C_VOL_X / 2.f) * (x - C_VOL_X / 2.f) + (ray_y - C_VOL_Y / 2.f) * (ray_y - C_VOL_Y / 2.f)) < (C_VOL_X / 2.f * C_VOL_Y / 2.f)) {
                if ((ray_y >= 0) && (ray_y < C_VOL_Y - 1) && (ray_z >= 0) && (ray_z < C_VOL_Z - 1)) {
                   raysum += pVol[int(ray_z) * C_VOL_X * C_VOL_Y + int(ray_y) * C_VOL_X + x];
                    
                }
            }
        }

        raysum *= sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z) / fabs(diff_x);

        raysum *= C_PITCH_X;

        pfPrjs[v * C_PRJ_U + u] = raysum;
    }

}

__global__ void kernelEtriDRRplaneY(float* pfPrjs, float* pVol)
{
    int u = blockDim.x * blockIdx.x + threadIdx.x + C_ROI_X;
    int v = blockDim.y * blockIdx.y + threadIdx.y + C_ROI_Y;

    float ray_start_x, ray_end_x, ray_start_y, ray_end_y, ray_start_z, ray_end_z, ray_temp;
    float x_m, y_m, z_m;
    int plane_start, plane_end;
    float raysum = 0.0;
    float ratio_k;
    int y;
    float ray_x, ray_z;
    float diff_x, diff_y, diff_z;

    if ((u < C_PRJ_U) && (v < C_PRJ_V))
    {        
        //디텍터 픽셀의 위치 계산
        x_m = C_DETECTOR_ORI_X + v * C_DETECTOR_V_X + u * C_DETECTOR_U_X;
        y_m = C_DETECTOR_ORI_Y + v * C_DETECTOR_V_Y + u * C_DETECTOR_U_Y;
        z_m = C_DETECTOR_ORI_Z + v * C_DETECTOR_V_Z + u * C_DETECTOR_U_Z;

        //소스-디텍터 픽셀 사이의 raysum 계산
        //소스와 디텍터 위치를 object 인덱스로 변경 (double 값)
        /*
        ray_start_x = x_m / C_PITCH_X + (int)(C_VOL_X / 2.f);
        ray_end_x = C_SOURCE_POS_X / C_PITCH_X + (int)(C_VOL_X / 2.f);
        ray_start_y = y_m / C_PITCH_Y + (int)(C_VOL_Y / 2.f);
        ray_end_y = C_SOURCE_POS_Y / C_PITCH_Y + (int)(C_VOL_Y / 2.f);
        ray_start_z = z_m / C_PITCH_Z + (int)(C_VOL_Z / 2.f);
        ray_end_z = C_SOURCE_POS_Z / C_PITCH_Z + (int)(C_VOL_Z / 2.f);
        */
        ray_start_x = (int)((x_m - C_ORIGIN_X) / C_PITCH_X);
        ray_end_x = (int)((C_SOURCE_POS_X - C_ORIGIN_X) / C_PITCH_X);
        ray_start_y = (int)((y_m - C_ORIGIN_Y) / C_PITCH_Y);
        ray_end_y = (int)((C_SOURCE_POS_Y - C_ORIGIN_Y) / C_PITCH_Y);
        ray_start_z = (int)((z_m - C_ORIGIN_Z) / C_PITCH_Z);
        ray_end_z = (int)((C_SOURCE_POS_Z - C_ORIGIN_Z) / C_PITCH_Z);

        diff_x = fabs(C_SOURCE_POS_X - x_m);
        diff_y = fabs(C_SOURCE_POS_Y - y_m);
        diff_z = fabs(C_SOURCE_POS_Z - z_m);

        //시작과 끝 재조정
        if (ray_start_y > ray_end_y) {

            ray_temp = ray_start_x;
            ray_start_x = ray_end_x;
            ray_end_x = ray_temp;

            ray_temp = ray_start_y;
            ray_start_y = ray_end_y;
            ray_end_y = ray_temp;

            ray_temp = ray_start_z;
            ray_start_z = ray_end_z;
            ray_end_z = ray_temp;
        }

        //object 시작점 및 끝점 계산
        if (ray_start_y < 0.0) plane_start = 0;
        else plane_start = (int)ray_start_y + 1; // 강제로 +1 ; 어떤 경우에도 division by zero 해소
        if (ray_end_y > (C_VOL_Y - 1)) plane_end = (C_VOL_Y - 1);
        else plane_end = (int)ray_end_y - 1; //강제로 -1

        for (y = plane_start; y <= plane_end; y++) {

            //ray와 y-번째 plane의 교차점을 비례식으로 계산
            ratio_k = (ray_end_y - y) / (y - ray_start_y); // division by zero는 일어나지 않음

            ray_x = (ray_end_x + ratio_k * ray_start_x) / (1 + ratio_k);
            ray_z = (ray_end_z + ratio_k * ray_start_z) / (1 + ratio_k);

            //x-z 평면에서 interpolation을 위한 주변 4픽셀 위치 및 weight 계산
            //x-방향 계산
            if (((ray_x - C_VOL_X / 2.f) * (ray_x - C_VOL_X / 2.f) + (y - C_VOL_Y / 2.f) * (y - C_VOL_Y / 2.f)) < (C_VOL_X / 2.f * C_VOL_Y / 2.f)) {
                if ((ray_x >= 0) && (ray_x < C_VOL_X - 1) && (ray_z >= 0) && (ray_z < C_VOL_Z - 1)) {
                   raysum += pVol[int(ray_z) * C_VOL_X * C_VOL_Y + y * C_VOL_X + int(ray_x)];
                   
                }
            }
        }

        raysum *= sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z) / fabs(diff_y);
        raysum *= C_PITCH_Y;

        pfPrjs[v * C_PRJ_U + u] = raysum;

    }
}


__global__ void kernelEtriDRRplaneZ(float* pfPrjs, float * pVol)
{
    int u = blockDim.x * blockIdx.x + threadIdx.x + C_ROI_X;
    int v = blockDim.y * blockIdx.y + threadIdx.y + C_ROI_Y;

    float ray_start_x, ray_end_x, ray_start_y, ray_end_y, ray_start_z, ray_end_z, ray_temp;
    float x_m, y_m, z_m;
    int plane_start, plane_end;
    float raysum = 0.0;
    float ratio_k;
    int z;
    float ray_x, ray_y;
    float diff_x, diff_y, diff_z;

    if ((u < C_PRJ_U) && (v < C_PRJ_V))
    {
        //디텍터 픽셀의 위치 계산
        x_m = C_DETECTOR_ORI_X + v * C_DETECTOR_V_X + u * C_DETECTOR_U_X;
        y_m = C_DETECTOR_ORI_Y + v * C_DETECTOR_V_Y + u * C_DETECTOR_U_Y;
        z_m = C_DETECTOR_ORI_Z + v * C_DETECTOR_V_Z + u * C_DETECTOR_U_Z;

        //소스-디텍터 픽셀 사이의 raysum 계산
        //소스와 디텍터 위치를 object 인덱스로 변경 (double 값)
        /*
        ray_start_x = x_m / C_PITCH_X + (int)(C_VOL_X / 2.f);
        ray_end_x = C_SOURCE_POS_X / C_PITCH_X + (int)(C_VOL_X / 2.f);
        ray_start_y = y_m / C_PITCH_Y + (int)(C_VOL_Y / 2.f);
        ray_end_y = C_SOURCE_POS_Y / C_PITCH_Y + (int)(C_VOL_Y / 2.f);
        ray_start_z = z_m / C_PITCH_Z + (int)(C_VOL_Z / 2.f);
        ray_end_z = C_SOURCE_POS_Z / C_PITCH_Z + (int)(C_VOL_Z / 2.f);
        */
        ray_start_x = (int)((x_m - C_ORIGIN_X) / C_PITCH_X);
        ray_end_x = (int)((C_SOURCE_POS_X - C_ORIGIN_X) / C_PITCH_X);
        ray_start_y = (int)((y_m - C_ORIGIN_Y) / C_PITCH_Y);
        ray_end_y = (int)((C_SOURCE_POS_Y - C_ORIGIN_Y) / C_PITCH_Y);
        ray_start_z = (int)((z_m - C_ORIGIN_Z) / C_PITCH_Z);
        ray_end_z = (int)((C_SOURCE_POS_Z - C_ORIGIN_Z) / C_PITCH_Z);

        diff_x = fabs(C_SOURCE_POS_X - x_m);
        diff_y = fabs(C_SOURCE_POS_Y - y_m);
        diff_z = fabs(C_SOURCE_POS_Z - z_m);

        //시작과 끝 재조정
        if (ray_start_z > ray_end_z) {

            ray_temp = ray_start_x;
            ray_start_x = ray_end_x;
            ray_end_x = ray_temp;

            ray_temp = ray_start_y;
            ray_start_y = ray_end_y;
            ray_end_y = ray_temp;

            ray_temp = ray_start_z;
            ray_start_z = ray_end_z;
            ray_end_z = ray_temp;
        }

        //object 시작점 및 끝점 계산
        if (ray_start_z < 0.0) plane_start = 0;
        else plane_start = (int)ray_start_z + 1; // 강제로 +1 ; 아래의 division by zero 해소
        if (ray_end_z > (C_VOL_Z - 1)) plane_end = (C_VOL_Z - 1);
        else plane_end = (int)ray_end_z - 1; //강제로 -1

        for (z = plane_start; z <= plane_end; z++) {

            //ray와 z-번째 plane의 교차점을 비례식으로 계산
            ratio_k = (ray_end_z - z) / (z - ray_start_z); // division by zero problem 없음
            ray_y = (ray_end_y + ratio_k * ray_start_y) / (1 + ratio_k);
            ray_x = (ray_end_x + ratio_k * ray_start_x) / (1 + ratio_k);

            //x-y 평면에서 interpolation을 위한 주변 4픽셀 위치 및 weight 계산
            if (((ray_x - C_VOL_X / 2.f) * (ray_x - C_VOL_X / 2.f) + (ray_y - C_VOL_Y / 2.f) * (ray_y - C_VOL_Y / 2.f)) < (C_VOL_X / 2.f * C_VOL_Y / 2.f)) {
                if ((ray_x >= 0) && (ray_x < C_VOL_X - 1) && (ray_y >= 0) && (ray_y < C_VOL_Y - 1)) {
                   raysum += pVol[z * C_VOL_X * C_VOL_Y + int(ray_y) * C_VOL_X + int(ray_x)];
                    
                }
            }
        }

        raysum *= sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z) / fabs(diff_z);
        raysum *= C_PITCH_Z;

        pfPrjs[v * C_PRJ_U + u] = raysum;
    }
}


bool cuda_volume_shrink(float* dst, float *src, int old_x, int old_y, int old_z, int new_x, int new_y, int new_z, float p_x, float p_y, float p_z, int s_x, int s_y, int s_z)
{
    bool bRtn = true;

    float new_pitch_x = p_x * s_x;
    float new_pitch_y = p_y * s_y;
    float new_pitch_z = p_z * s_z;

    int TILE_DIM_X;
    int TILE_DIM_Y;
    int TILE_DIM_Z;

    TILE_DIM_X = 16;
    TILE_DIM_Y = 16;
    TILE_DIM_Z = 4;

    dim3 dimGrid(int(ceilf(new_x / TILE_DIM_X)), int(ceilf(new_y / TILE_DIM_Y)), int(ceilf(new_z / TILE_DIM_Z)));
    dim3 dimBlock(TILE_DIM_X, TILE_DIM_Y, TILE_DIM_Z);


    cudaMemcpyToSymbol(C_VOL_X, &(new_x), sizeof(int));
    cudaMemcpyToSymbol(C_VOL_Y, &(new_y), sizeof(int));
    cudaMemcpyToSymbol(C_VOL_Z, &(new_z), sizeof(int));

    cudaMemcpyToSymbol(C_SCALE_X, &(s_x), sizeof(int));
    cudaMemcpyToSymbol(C_SCALE_Y, &(s_y), sizeof(int));
    cudaMemcpyToSymbol(C_SCALE_Z, &(s_z), sizeof(int));


    cudaMemcpyToSymbol(C_OLD_VOL_X, &(old_x), sizeof(int));
    cudaMemcpyToSymbol(C_OLD_VOL_Y, &(old_y), sizeof(int));
    cudaMemcpyToSymbol(C_OLD_VOL_Z, &(old_z), sizeof(int));

    cudaMemcpyToSymbol(C_PITCH_X, &(new_pitch_x), sizeof(float));
    cudaMemcpyToSymbol(C_PITCH_Y, &(new_pitch_y), sizeof(float));
    cudaMemcpyToSymbol(C_PITCH_Z, &(new_pitch_z), sizeof(float));

    kernelEtriShrink << <dimGrid, dimBlock >> > (dst, src);

    return bRtn;
}

bool cuda_foward_projection(float* d_prj, float* d_vol, int w, int h, int dir)
{
    bool bRtn = true;
    int TILE_DIM_X;
    int TILE_DIM_Y;

    TILE_DIM_X = 16;
    TILE_DIM_Y = 16;
    dim3 dimGrid(int(ceilf(w / TILE_DIM_X)), int(ceilf(h / TILE_DIM_Y)));
    dim3 dimBlock(TILE_DIM_X, TILE_DIM_Y);

    if (dir == 1)
        kernelEtriDRRplaneX << <dimGrid, dimBlock >> > (d_prj, d_vol);
    else if (dir == 2)
        kernelEtriDRRplaneY << <dimGrid, dimBlock >> > (d_prj, d_vol);
    else
        kernelEtriDRRplaneZ << <dimGrid, dimBlock >> > (d_prj, d_vol);

    cudaDeviceSynchronize();

    return bRtn;
}

bool normalization_digital_projeciton_image(float* d_prj, int w, int h, float norm_factor)
{
    bool bRtn = true;
    int TILE_DIM_X;
    int TILE_DIM_Y;

    TILE_DIM_X = 16;
    TILE_DIM_Y = 16;
    dim3 dimGrid(int(ceilf(w / TILE_DIM_X)), int(ceilf(h / TILE_DIM_Y)));
    dim3 dimBlock(TILE_DIM_X, TILE_DIM_Y);

    cudaMemcpyToSymbol(C_N_FACTOR, &norm_factor, sizeof(float));

    kernelNormalization2D32f<<< dimGrid, dimBlock>>>(d_prj);

    return bRtn;
}


__global__ void kernelEtriShrink(float* dst, float* src)
{
    int idx_x = blockIdx.x * blockDim.x + threadIdx.x;
    int idx_y = blockIdx.y * blockDim.y + threadIdx.y;
    int idx_z = blockIdx.z * blockDim.z + threadIdx.z;

    int src_x = idx_x * C_SCALE_X;
    int src_y = idx_y * C_SCALE_Y;
    int src_z = idx_z * C_SCALE_Z;

    if ((idx_x < C_VOL_X) && (idx_y < C_VOL_Y) && (idx_z < C_VOL_Z))
    {
        if ((src_x < C_OLD_VOL_X) && (src_y < C_OLD_VOL_Y) && (src_z < C_OLD_VOL_Z))
            dst[(idx_z * C_VOL_X * C_VOL_Y) + (idx_y * C_VOL_X) + idx_x] = src[(src_z * C_OLD_VOL_X * C_OLD_VOL_Y) + (src_y * C_OLD_VOL_X) + src_x];
        else
            dst[(idx_z * C_VOL_X * C_VOL_Y) + (idx_y * C_VOL_X) + idx_x] = 0.f;
    }
}
__global__ void kernelNormalization2D32f(float* img)
{
    int u = blockDim.x * blockIdx.x + threadIdx.x + C_ROI_X;
    int v = blockDim.y * blockIdx.y + threadIdx.y + C_ROI_Y;

    int max_value = (2 << 15) - 1;

    if ((u < C_PRJ_U) && (v < C_PRJ_V))
    {
        //img[v * C_PRJ_U + u] = max_value * exp(-(img[v * C_PRJ_U + u]) / C_N_FACTOR);
        img[v * C_PRJ_U + u] = img[v * C_PRJ_U + u] * C_N_FACTOR;
    }
}

