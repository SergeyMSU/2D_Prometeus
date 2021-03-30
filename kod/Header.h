#pragma once
#include "Cell.h"
struct Cell;

extern __device__ double get_square(const double& x1, const double& y1, const double& dx1, const double& dy1, const double& x2, const double& y2, //
    const double& dx2, const double& dy2, double& n1, double& n2, double& dist);

extern __device__ double HLLC_2d_Korolkov_b_s(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L,//
    const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, //
    double* P, double& PQ, const double& n1, const double& n2, const double& rad, int metod);

extern __global__ void Cuda_main_HLLDQ(int* NN, double* X, double* Y, int* Size,//
    double* RO1, double* RO2, double* Q1, double* Q2, double* P1, double* P2, double* U1, double* U2, double* V1, double* V2,//
    int* SOSED, int* L, int* R, double* T, double* T_do, int step_, double DX, double DY, int metod = 0);