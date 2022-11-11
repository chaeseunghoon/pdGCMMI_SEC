#include "cuMI.cuh"

bool img_rotation_shift(unsigned short* d_trans_ref, unsigned short* d_ref, int w, int h, int roi_w, int roi_h, float trs_x, float trs_y, float rot)
{
	bool bRtn = true;

    float cc, ss;
    int TILE_DIM_X = 16;
    int TILE_DIM_Y = 16;
    dim3 dimGrid(int(ceilf(roi_w / TILE_DIM_X)), int(ceilf(roi_h / TILE_DIM_Y)));
    dim3 dimBlock(TILE_DIM_X, TILE_DIM_Y);

    cc = cos(rot);
    ss = sin(rot);

    cudaMemcpyToSymbol(C_TRANS_X, &(trs_x), sizeof(float));
    cudaMemcpyToSymbol(C_TRANS_Y, &(trs_y), sizeof(float));
    cudaMemcpyToSymbol(C_COS_VAL, &(cc), sizeof(float));
    cudaMemcpyToSymbol(C_SIN_VAL, &(ss), sizeof(float));
    cudaMemset(d_trans_ref, 0, sizeof(unsigned short) * static_cast<size_t>(w) * h);

    kernelEtriRotationShift << <dimGrid, dimBlock >> > (d_trans_ref, d_ref);

    return bRtn;
}

bool etri_set_common_params(SPARAMS_MI *params)
{
    bool bRtn = true;

    cudaMemcpyToSymbol(C_IMG_W, &(params->w), sizeof(int));
    cudaMemcpyToSymbol(C_IMG_H, &(params->h), sizeof(int));
    cudaMemcpyToSymbol(C_ROI_X, &(params->roi_x), sizeof(int));
    cudaMemcpyToSymbol(C_ROI_Y, &(params->roi_y), sizeof(int));
    cudaMemcpyToSymbol(C_HISTBIN, &(params->hist_bins), sizeof(int));
    cudaMemcpyToSymbol(C_SUM, &(params->mask_sum), sizeof(int));

    return bRtn;
}


bool normalization_projeciton_image(unsigned short* d_img, int roi_w, int roi_h, int w, int h, float norm_factor)
{
    bool bRtn = true;
    int TILE_DIM_X;
    int TILE_DIM_Y;

    TILE_DIM_X = 16;
    TILE_DIM_Y = 16;
    dim3 dimGrid(int(ceilf(roi_w / TILE_DIM_X)), int(ceilf(roi_h / TILE_DIM_Y)));
    dim3 dimBlock(TILE_DIM_X, TILE_DIM_Y);

    cudaMemcpyToSymbol(C_N_FACTOR, &norm_factor, sizeof(float));

    kernelNormalization2D16u << < dimGrid, dimBlock >> > (d_img);
    return bRtn;
}

bool joint_histogram(float* hist, float* img, unsigned short* ref, int roi_w, int roi_h, unsigned char* mask, int bin)
{
    bool bRtn = true;

    int TILE_DIM_X = 16;
    int TILE_DIM_Y = 16;
    dim3 dimGrid(int(ceilf(roi_w / TILE_DIM_X)), int(ceilf(roi_h / TILE_DIM_Y)));
    dim3 dimBlock(TILE_DIM_X, TILE_DIM_Y);

    cudaMemset(hist, 0, static_cast<size_t>(bin) * bin * sizeof(float));
    kernelJointHistogram << < dimGrid, dimBlock >> > (hist, img, ref, mask);

    return bRtn;

}

bool normalization_joint_histogram(float* d_hist, int bin)
{
    bool bRtn = true;

    dim3 dimBlock(256, 1, 1);
    dim3 dimGrid;
    dimGrid.x = ((bin * bin) + dimBlock.x - 1) / dimBlock.x;

    kernelNormJointHistogram << <dimGrid, dimBlock >> > (d_hist);

    return bRtn;
}

float get_MI_value(float* d_hist, float *prob_x, float *prob_y, float *d_mi, int bin)
{
    float mi, mi_x, mi_y, norm_mi;

    int TILE_DIM_X = 16;
    int TILE_DIM_Y = 16;
    dim3 dimGrid2(int(ceilf(bin / TILE_DIM_X)), int(ceilf(bin / TILE_DIM_Y)));
    dim3 dimBlock2(TILE_DIM_X, TILE_DIM_Y);
    kernelSumX << <dimGrid2, dimBlock2 >> > (d_hist, prob_x);
    kernelSumY << <dimGrid2, dimBlock2 >> > (d_hist, prob_y);
    kernelGetMI << <dimGrid2, dimBlock2 >> > (d_hist, prob_x, prob_y, d_mi);
    
    cudaMemcpy(&mi, d_mi, sizeof(float), cudaMemcpyDeviceToHost);
    kernelGetMI_X << <dimGrid2, dimBlock2 >> > (d_hist, prob_x, prob_y, d_mi);
    cudaMemcpy(&mi_x, d_mi, sizeof(float), cudaMemcpyDeviceToHost);

    kernelGetMI_Y << <dimGrid2, dimBlock2 >> > (d_hist, prob_x, prob_y, d_mi);
    cudaMemcpy(&mi_y, d_mi, sizeof(float), cudaMemcpyDeviceToHost);

    norm_mi = mi / min(mi_x, mi_y);
    
    return norm_mi;
}

__global__ void kernelGetMI(float* d_hist, float* prob_x, float* prob_y, float* mi_val)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    float val;

    if (x < C_HISTBIN && y < C_HISTBIN)
    {
        if (prob_x[x] * prob_y[y] > 1.0e-10)
        {
            if (d_hist[(y * C_HISTBIN) + x] > 1.0e-10)
            {
                val = d_hist[(y * C_HISTBIN) + x] * log(d_hist[(y * C_HISTBIN) + x] / (prob_x[x] * prob_y[y]));
                atomicAdd(mi_val, val);

            }
        }
    }
}


__global__ void kernelGetMI_X(float* d_hist, float* prob_x, float* prob_y, float* mi_val)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    float val;

    if (x < C_HISTBIN && y < C_HISTBIN)
    {
        if (prob_x[x] * prob_y[y] > 1.0e-10)
        {
            if (d_hist[(y * C_HISTBIN) + x] > 1.0e-10)
            {
                val = -d_hist[(y * C_HISTBIN) + x] * log(prob_x[x]);
                atomicAdd(mi_val, val);
            }
        }
    }
}

__global__ void kernelGetMI_Y(float* d_hist, float* prob_x, float* prob_y, float* mi_val)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    float val;
    if (x < C_HISTBIN && y < C_HISTBIN)
    {
        if (prob_x[x] * prob_y[y] > 1.0e-10)
        {
            if (d_hist[(y * C_HISTBIN) + x] > 1.0e-10)
            {
                val = -d_hist[(y * C_HISTBIN) + x] * log(prob_y[y]);
                atomicAdd(mi_val, val);
            }
        }
    }
}

__global__ void kernelSumX(float* d_hist, float* pSum)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x < C_HISTBIN && y < C_HISTBIN)
    {
        atomicAdd(&pSum[x], d_hist[(y * C_HISTBIN) + x]);
    }
}

__global__ void kernelSumY(float* d_hist, float* pSum)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x < C_HISTBIN && y < C_HISTBIN)
    {
        atomicAdd(&pSum[y], d_hist[(y * C_HISTBIN) + x]);
    }
}

__global__ void kernelNormJointHistogram(float* d_hist)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i > C_HISTBIN * C_HISTBIN)	return;

    d_hist[i] = d_hist[i] / C_SUM;
}


__global__ void kernelJointHistogram(float* hist, float* img, unsigned short* ref, unsigned char* mask)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x + C_ROI_X;
    int y = blockDim.y * blockIdx.y + threadIdx.y + C_ROI_Y;
    
    int p0, p1;

    if ((x < C_IMG_W) && (y < C_IMG_H) && mask[(y * C_IMG_W) + x] != 0)
    {
        p0 = (int)(img[y* C_IMG_W+x] / (unsigned short(65536. / float(C_HISTBIN))));
        p1 = (int)(ref[y* C_IMG_W+x] / (unsigned short(65536. / float(C_HISTBIN))));

        atomicAdd(&hist[p0 * C_HISTBIN + p1], 1.0);
    }

}

__global__ void kernelNormalization2D16u(unsigned short* img)
{
    int u = blockDim.x * blockIdx.x + threadIdx.x + C_ROI_X;
    int v = blockDim.y * blockIdx.y + threadIdx.y + C_ROI_Y;
    int pixel;
    int max_value = (2 << 15) - 1;
    float logair = __log10f(float(max_value));

    if ((u < C_IMG_W) && (v < C_IMG_H))
    {
        pixel = int(img[v * C_IMG_W + u] * C_N_FACTOR);

        if (pixel > max_value) pixel = max_value;

        if (pixel == 0)
            pixel = max_value;
        else
        {
            pixel = int((__log10f(max_value / float(pixel)) / logair) * max_value);
        }

        img[v * C_IMG_W + u] = unsigned short(pixel);
    }
}

__global__ void kernelEtriRotationShift(unsigned short* trans, unsigned short* ori)
{
    int u = blockDim.x * blockIdx.x + threadIdx.x + C_ROI_X;
    int v = blockDim.y * blockIdx.y + threadIdx.y + C_ROI_Y;
    float ori_x, ori_y;
    int x_1, y_1, x_2, y_2;
    float xr, yr;

    float center_x, center_y;
    float pixel;
    unsigned short max_value;
    //   unsigned short max_d;
    //   float max_rate;
    //   max_d = ((2 << 15) - 1);
    max_value = ((2 << 15) - 1);
    //   max_rate = max_value / max_d;

    center_x = C_IMG_W / 2.0;
    center_y = C_IMG_H/ 2.0;


    if ((u < C_IMG_W) && (v < C_IMG_H))
    {
        // rotation -> translation
        ori_x = (center_x + (v - C_TRANS_Y  - center_y) * C_SIN_VAL + (u - C_TRANS_X - center_x) * C_COS_VAL);
        ori_y = (center_y + (v - C_TRANS_Y - center_y) * C_COS_VAL - (u - C_TRANS_X - center_x) * C_SIN_VAL);

        if ((ori_y >= 0 && ori_y <= C_IMG_H - 1) && (ori_x >= 0 && ori_x <= C_IMG_W - 1))
        {

            x_1 = (int)ori_x;
            y_1 = (int)ori_y;

            xr = ori_x - x_1;
            yr = ori_y - y_1;

            x_2 = x_1 + 1.;
            if (x_2 >= C_IMG_W)
                x_2 = C_IMG_W - 1;
            y_2 = y_1 + 1.;
            if (y_2 >= C_IMG_H)
                y_2 = C_IMG_H - 1;

            pixel = ori[(y_1 * C_IMG_W) + x_1] * (1. - xr) * (1. - yr)
                + ori[(y_1 * C_IMG_W) + x_2] * xr * (1. - yr)
                + ori[(y_2 * C_IMG_W) + x_1] * (1. - xr) * yr
                + ori[(y_2 * C_IMG_W) + x_2] * xr * yr;

            if (pixel > max_value) pixel = max_value;

            trans[v * C_IMG_W + u] = unsigned short(pixel); // (max_d - pixel)* max_rate);
        }
        else
        {
            trans[v * C_IMG_W + u] = unsigned short(0);
        }
    }
}
