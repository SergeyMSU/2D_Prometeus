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

__global__ void Cuda_main_5_komponent(int* NN, double* X, double* Y, int* Size,//
    double* RO1, double* RO2, double* Q1, double* Q2, double* P1, double* P2, double* U1, double* U2, double* V1, double* V2,//
    double* RO1_H1, double* RO2_H1, double* P1_H1, double* P2_H1, double* U1_H1, double* U2_H1, double* V1_H1, double* V2_H1,//
    double* RO1_H2, double* RO2_H2, double* P1_H2, double* P2_H2, double* U1_H2, double* U2_H2, double* V1_H2, double* V2_H2,//
    double* RO1_H3, double* RO2_H3, double* P1_H3, double* P2_H3, double* U1_H3, double* U2_H3, double* V1_H3, double* V2_H3,//
    double* RO1_H4, double* RO2_H4, double* P1_H4, double* P2_H4, double* U1_H4, double* U2_H4, double* V1_H4, double* V2_H4,//
    int* SOSED, int* L, int* R, double* T, double* T_do, int step_, double DX, double DY, int metod);

__global__ void Cuda_main_HLLDQ_M_K(int* NN, double* X, double* Y, int* Size,//
    double* RO1, double* RO2, double* Q1, double* Q2, double* P1, double* P2, double* U1, double* U2, double* V1, double* V2,//
    double* I_u, double* I_v, double* I_T, int* SOSED, int* L, int* R, double* T, double* T_do, int step_, double DX, double DY, int metod);