#include "Header.h"
#include "math.h"


/// ¬се функции описанные в этой файле должны иметь прототипы в файле Header.h
__device__ double polar_angle(double x, double y)
{
    if (x < 0)
    {
        return atan(y / x) + 1.0 * PI;
    }
    else if (x > 0 && y >= 0)
    {
        return atan(y / x);
    }
    else if (x > 0 && y < 0)
    {
        return atan(y / x) + 2.0 * PI;
    }
    else if (y > 0 && x >= 0 && x <= 0)
    {
        return PI / 2.0;
    }
    else if (y < 0 && x >= 0 && x <= 0)
    {
        return  3.0 * PI / 2.0;
    }
    return 0.0;
}

__device__ double polar_perenos(const double& x1, const double& y1, const double& x2, const double& y2, double& u, double& v)
{
    double phi1 = polar_angle(x1, y1);
    double phi2 = polar_angle(x2, y2);
    double fr = u * cos(phi1) + v * sin(phi1);
    double ff = -u * sin(phi1) + v * cos(phi1);
    u = fr * cos(phi2) - ff * sin(phi2);
    v = fr * sin(phi2) + ff * cos(phi2);
}

__device__ double get_square(const double& x1, const double& y1, const double& dx1, const double& dy1, const double& x2, const double& y2, //
    const double& dx2, const double& dy2,  double& n1, double& n2, double& dist)
{
    if (fabs(fabs(x1 - x2) - dx1 - dx2) < geo)
    {
        n1 = (x2 - x1) / fabs(x1 - x2);
        n2 = 0.0;
        dist = min(dx1, dx2);
        return 2.0 * min(dy1, dy2);
    }
    else if (fabs(fabs(y1 - y2) - dy1 - dy2) < geo)
    {
        n2 = (y2 - y1) / fabs(y1 - y2);
        n1 = 0.0;
        dist = min(dy1, dy2);
        return 2.0 * min(dx1, dx2);
    }
    else
    {
        printf("Error:  get_square: %lf, %lf, %lf, %lf, %lf,  %lf,  %lf,  %lf\n", //
            x1, y1, x2, y2, dx1, dy1, dx2, dy2);
    }
    return 0.0;
}

__device__ double HLLC_2d_Korolkov_b_s(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L,//
    const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, //
    double* P, double& PQ, const double& n1, const double& n2, const double& rad, int metod)
    // BestSeries
    // Ћучший работающий 2д распадник
    //
    //  ¬ывод:
    // P[1]       // —корости
    // P[2]
    // P[0]       // ћасса
    // P[3]       // Ёнерги€
{
    double t1 = -n2;
    double t2 = n1;

    double u1, v1, u2, v2;
    u1 = v1_L * n1 + v2_L * n2;
    v1 = v1_L * t1 + v2_L * t2;
    u2 = v1_R * n1 + v2_R * n2;
    v2 = v1_R * t1 + v2_R * t2;

    double sqrtroL = sqrt(ro_L);
    double sqrtroR = sqrt(ro_R);
    double cL = sqrt(ggg * p_L / ro_L);
    double cR = sqrt(ggg * p_R / ro_R);


    double uu_L = (kv(v1_L) + kv(v2_L)) / 2.0;
    double uu_R = (kv(v1_R) + kv(v2_R)) / 2.0;



    double SL = min(u1, u2) - max(cL, cR);
    double SR = max(u1, u2) + max(cL, cR);

    double suR = (SR - u2);
    double suL = (SL - u1);

    double SM = (suR * ro_R * u2 - suL * ro_L * u1 - p_R + p_L) //
        / (suR * ro_R - suL * ro_L);

    double PTT = (suR * ro_R * p_L - suL * ro_L * p_R + ro_L * ro_R * suR * suL * (u2 - u1)) / (suR * ro_R - suL * ro_L);

    double UU = max(fabs(SL), fabs(SR));
    double time = kurant * rad / UU;

    double FL[5], FR[5], UL[5], UR[5];

    double e1 = p_L / g1 + ro_L * uu_L;
    double e2 = p_R / g1 + ro_R * uu_R;


    FL[0] = ro_L * u1;
    FL[1] = ro_L * u1 * u1 + p_L;
    FL[2] = ro_L * u1 * v1;
    FL[3] = (e1 + p_L) * u1;
    FL[4] = Q_L * u1;

    if (SL >= 0.0)
    {
        P[1] = n1 * FL[1] + t1 * FL[2];     // —корости
        P[2] = n2 * FL[1] + t2 * FL[2];
        P[0] = FL[0];                       // ћасса
        P[3] = FL[3];                       // Ёнерги€
        PQ = FL[4];
        return time;
    }

    FR[0] = ro_R * u2;
    FR[1] = ro_R * u2 * u2 + p_R;
    FR[2] = ro_R * u2 * v2;
    FR[3] = (e2 + p_R) * u2;
    FR[4] = Q_R * u2;

    if (SR <= 0.0)
    {
        P[1] = n1 * FR[1] + t1 * FR[2];     // —корости
        P[2] = n2 * FR[1] + t2 * FR[2];
        P[0] = FR[0];                       // ћасса
        P[3] = FR[3];                       // Ёнерги€
        PQ = FR[4];
        return time;
    }

    UL[0] = ro_L;
    UL[1] = ro_L * u1;
    UL[2] = ro_L * v1;
    UL[3] = e1;
    UL[4] = Q_L;

    UR[0] = ro_R;
    UR[1] = ro_R * u2;
    UR[2] = ro_R * v2;
    UR[3] = e2;
    UR[4] = Q_R;

    if (metod == 0)
    {
        double  PO[5];
        for (int i = 0; i < 5; i++)
        {
            PO[i] = (SR * FL[i] - SL * FR[i] + SR * SL * (UR[i] - UL[i])) / (SR - SL);
        }

        P[1] = n1 * PO[1] + t1 * PO[2];     // —корости
        P[2] = n2 * PO[1] + t2 * PO[2];
        P[0] = PO[0];                       // ћасса
        P[3] = PO[3];                       // Ёнерги€
        PQ = PO[4];
        return time;
    }


    double ro_LL = ro_L * (SL - u1) / (SL - SM);
    double ro_RR = ro_R * (SR - u2) / (SR - SM);
    double Q_LL = Q_L * (SL - u1) / (SL - SM);
    double Q_RR = Q_R * (SR - u2) / (SR - SM);


    double UZ0 = (SR * UR[0] - SL * UL[0] + FL[0] - FR[0]) / (SR - SL);
    double UZ1 = (SR * UR[1] - SL * UL[1] + FL[1] - FR[1]) / (SR - SL);
    double UZ2 = (SR * UR[2] - SL * UL[2] + FL[2] - FR[2]) / (SR - SL);
    double UZ3 = (SR * UR[3] - SL * UL[3] + FL[3] - FR[3]) / (SR - SL);
    double UZ4 = (SR * UR[4] - SL * UL[4] + FL[4] - FR[4]) / (SR - SL);
    double vzL, vzR, vLL, vRR, ppLR, ee1, ee2;


    double suRm = suR / (SR - SM);
    double suLm = suL / (SL - SM);
    double rzR = ro_R * suRm;
    double rzL = ro_L * suLm;

    double ptzR = p_R + ro_R * suR * (SM - u2);
    double ptzL = p_L + ro_L * suL * (SM - u1);
    double ptz = (ptzR + ptzL) / 2.0;


    if( fabs(v1 - v2) > 0.1)
    {
        vLL = v1;
        vRR = v2;
    }
    else
    {
        vRR = UZ2 / UZ0;
        vLL = vRR;
    }


    ee2 = e2 * suRm + (ptz * SM - p_R * u2) / (SR - SM);
    ee1 = e1 * suLm + (ptz * SM - p_L * u1) / (SL - SM);


    double  ULL[5], URR[5], PO[5];
    ULL[0] = ro_LL;
    ULL[1] = ro_LL * SM;
    ULL[2] = ro_LL * vLL;
    ULL[3] = ee1;
    ULL[4] = Q_LL;

    URR[0] = ro_RR;
    URR[1] = ro_RR * SM;
    URR[2] = ro_RR * vRR;
    URR[3] = ee2;
    URR[4] = Q_RR;

    if (SL < 0.0 && SM >= 0.0)
    {
        for (int i = 0; i < 5; i++)
        {
            PO[i] = FL[i] + SL * ULL[i] - SL * UL[i];
        }
    }
    else if (SR > 0.0 && SM < 0.0)
    {
        for (int i = 0; i < 5; i++)
        {
            PO[i] = FR[i] + SR * URR[i] - SR * UR[i];
        }
    }

    P[1] = n1 * PO[1] + t1 * PO[2];     // —корости
    P[2] = n2 * PO[1] + t2 * PO[2];
    P[0] = PO[0];                       // ћасса
    P[3] = PO[3];                       // Ёнерги€
    PQ = PO[4];

    return time;
}

__global__ void Cuda_main_HLLDQ(int* NN, double* X, double* Y, int* Size,//
    double* RO1, double* RO2, double* Q1, double* Q2, double* P1, double* P2, double* U1, double* U2, double* V1, double* V2,//
    int* SOSED, int* L, int* R, double* T, double* T_do, int step_, double DX, double DY, int metod)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x; // √лобальный индекс текущего потока
    if (index > * NN - 1)
    {
        return;
    }
    double x, y, ro, p, u, v, Q;
    int size = Size[index];
    double dx = (DX / pow(2, size - 1)) / 2.0;   // ѕоловина длины €чейки
    double dy = (DY / pow(2, size - 1)) / 2.0;   // ѕоловина ширины €чейки
    int l = L[index];
    int r = R[index];
    x = X[index];
    y = Y[index];
    ro = RO1[index];
    p = P1[index];
    u = U1[index];
    v = V1[index];
    Q = Q1[index];
    double radius = sqrt(kv(x) + kv(y));


    if (radius <= Distant) // || (ddd <= 4.0 && x > -5 && x < 0) ) //(ddd < 5.76 || ddd2 <= 2.0) //1.5
    {
        RO2[index] = ro;
        P2[index] = p;
        U2[index] = u;
        V2[index] = v;
        Q2[index] = Q;
    }
    else
    {
        double PQ = 0.0;
        double n1 = 0.0;
        double n2 = 0.0;
        double dist = 0.0;
        double P[4] = { 0.0 };
        P[0] = P[1] = P[2] = P[3] = 0.0;
        double Potok[5] = { 0.0 };
        Potok[0] = Potok[1] = Potok[2] = Potok[3] = Potok[4] = 0.0;
        double tmin = 10000000;
        double Volume = dx * dy * 4.0;
        int ii = 0;
        double x2, y2, dx2, dy2, ro2, p2, u2, v2, Q_2, size2;
        double roC = 1.0; 
        double pC = 1.0; 
        double uC = Velosity_inf;
        double vC = 0.0;
        double QC = 100.0;
        double u1_polar, v1_polar;

        for (int i = l; i <= r; i++)
        {
            ii = SOSED[i];
            if (ii >= 0)
            {
                x2 = X[ii];
                y2 = Y[ii];
                size2 = Size[ii];
                dx2 = (DX / pow(2, size2 - 1)) / 2.0;   // ѕоловина длины €чейки
                dy2 = (DY / pow(2, size2 - 1)) / 2.0;   // ѕоловина ширины €чейки
                ro2 = RO1[ii];
                p2 = P1[ii];
                u2 = U1[ii];
                v2 = V1[ii];
                Q_2 = Q1[ii];
                double S = get_square(x, y, dx, dy, x2, y2, dx2, dy2, n1, n2, dist);

                u1_polar = u;
                v1_polar = v;

                if (radius < 100)
                {
                    polar_perenos(x, y, x + n1 * dx, y + n2 * dy, u1_polar, v1_polar);
                    polar_perenos(x2, y2, x2 - n1 * dx2, y2 - n2 * dy2, u2, v2);
                }
                

                tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro, Q, p, u1_polar, v1_polar, ro2, Q_2, p2, u2, v2, P, PQ, n1, n2, dist, metod));

                for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[4] = Potok[4] + PQ * S;
            }
            else if (ii == -1)
            {
                double S = dy * 2.0;
                n1 = 1.0;
                n2 = 0.0;
                dist = dx;

                tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro, Q, p, u, v, roC, QC, pC, uC, vC, P, PQ, n1, n2, dist, metod));

                for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[4] = Potok[4] + PQ * S;
            }
            else if (ii == -2)
            {
                double S = dy * 2.0;
                n1 = -1.0;
                n2 = 0.0;
                dist = dx;

                double uu = u;
                if (uu > Velosity_inf && y < 300)
                {
                    uu = Velosity_inf;
                }

                tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro, Q, p, u, v, ro, Q, p, uu, v, P, PQ, n1, n2, dist, metod));

                for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[4] = Potok[4] + PQ * S;
            }
            else if (ii == -3)
            {
                double S = dx * 2.0;
                n1 = 0.0;
                n2 = 1.0;
                dist = dy;

                tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro, Q, p, u, v, ro, Q, p, u, v, P, PQ, n1, n2, dist, metod));

                for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[4] = Potok[4] + PQ * S;
            }
            else if (ii == -4)
            {
                double S = dx * 2.0;
                n1 = 0.0;
                n2 = -1.0;
                dist = dy;

                u1_polar = u;
                v1_polar = v;

                if (radius < 100)
                {
                    u1_polar = u;
                    v1_polar = v;
                    polar_perenos(x, y, x + n1 * dx, y + n2 * dy, u1_polar, v1_polar);
                }


                tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro, Q, p, u1_polar, v1_polar, ro, Q, p, //
                                                        u1_polar, -v1_polar, P, PQ, n1, n2, dist, metod));

                for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[4] = Potok[4] + PQ * S;
            }
            else
            {
                printf("Error 12438wedew4353jdyu. Ne doljni suda popadat = %d \n", ii);
            }
        }

        double ro3, p3, u3, v3, Q33;

        ro3 = ro - *T_do * (Potok[0] / Volume + ro * v / y);
        Q33 = Q - (*T_do / Volume) * Potok[4] - *T_do * Q * v / y;
        if (ro3 <= 0)
        {
            printf("Problemsssss  ro < 0! %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", x, y, dx, dy, ro, p, u, v, Q);
            ro3 = 0.00001;
        }
        u3 = (ro * u - *T_do * (Potok[1] / Volume + ro * v * u / y)) / ro3;
        v3 = (ro * v - *T_do * (Potok[2] / Volume + ro * v * v / y)) / ro3;
        p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - *T_do * (Potok[3] / Volume + //
            + v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y)) - //
            0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
        if (p3 <= 0)
        {
            p3 = 0.000001;
        }

        Q2[index] = Q33;
        RO2[index] = ro3;
        P2[index] = p3;
        U2[index] = u3;
        V2[index] = v3;

        if (*T > tmin)
        {
            *T = tmin;
            __threadfence();
        }
    }

}

__global__ void Cuda_main_5_komponent(int* NN, double* X, double* Y, int* Size,//
    double* RO1, double* RO2, double* Q1, double* Q2, double* P1, double* P2, double* U1, double* U2, double* V1, double* V2,//
    double* RO1_H1, double* RO2_H1, double* P1_H1, double* P2_H1, double* U1_H1, double* U2_H1, double* V1_H1, double* V2_H1,//
    double* RO1_H2, double* RO2_H2, double* P1_H2, double* P2_H2, double* U1_H2, double* U2_H2, double* V1_H2, double* V2_H2,//
    double* RO1_H3, double* RO2_H3, double* P1_H3, double* P2_H3, double* U1_H3, double* U2_H3, double* V1_H3, double* V2_H3,//
    double* RO1_H4, double* RO2_H4, double* P1_H4, double* P2_H4, double* U1_H4, double* U2_H4, double* V1_H4, double* V2_H4,//
    int* SOSED, int* L, int* R, double* T, double* T_do, int step_, double DX, double DY, int metod)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x; // √лобальный индекс текущего потока
    if (index > * NN - 1)
    {
        return;
    }
    double x, y, ro, p, u, v, Q;
    double ro_H1, p_H1, u_H1, v_H1;
    double ro_H2, p_H2, u_H2, v_H2;
    double ro_H3, p_H3, u_H3, v_H3;
    double ro_H4, p_H4, u_H4, v_H4;
    int size = Size[index];
    double dx = (DX / pow(2, size - 1)) / 2.0;   // ѕоловина длины €чейки
    double dy = (DY / pow(2, size - 1)) / 2.0;   // ѕоловина ширины €чейки
    int l = L[index];
    int r = R[index];
    x = X[index];
    y = Y[index];
    ro = RO1[index];
    p = P1[index];
    u = U1[index];
    v = V1[index];
    ro_H1 = RO1_H1[index];
    p_H1 = P1_H1[index];
    u_H1 = U1_H1[index];
    v_H1 = V1_H1[index];
    ro_H2 = RO1_H2[index];
    p_H2 = P1_H2[index];
    u_H2 = U1_H2[index];
    v_H2 = V1_H2[index];
    ro_H3 = RO1_H3[index];
    p_H3 = P1_H3[index];
    u_H3 = U1_H3[index];
    v_H3 = V1_H3[index];
    ro_H4 = RO1_H4[index];
    p_H4 = P1_H4[index];
    u_H4 = U1_H4[index];
    v_H4 = V1_H4[index];
    Q = Q1[index];

    double radius = sqrt(kv(x) + kv(y));


    double PQ = 0.0;
    double n1 = 0.0;
    double n2 = 0.0;
    double dist = 0.0;
    double P[4] = { 0.0 };
    P[0] = P[1] = P[2] = P[3] = 0.0;
    double Potok[5] = { 0.0 };
    Potok[0] = Potok[1] = Potok[2] = Potok[3] = Potok[4] = 0.0;
    double Potok_H1[4] = { 0.0 };
    Potok_H1[0] = Potok_H1[1] = Potok_H1[2] = Potok_H1[3] = 0.0;
    double Potok_H2[4] = { 0.0 };
    Potok_H2[0] = Potok_H2[1] = Potok_H2[2] = Potok_H2[3] = 0.0;
    double Potok_H3[4] = { 0.0 };
    Potok_H3[0] = Potok_H3[1] = Potok_H3[2] = Potok_H3[3] = 0.0;
    double Potok_H4[4] = { 0.0 };
    Potok_H4[0] = Potok_H4[1] = Potok_H4[2] = Potok_H4[3] = 0.0;
    double tmin = 10000000;
    double Volume = dx * dy * 4.0;
    int ii = 0;
    double x2, y2, dx2, dy2, ro2, p2, u2, v2, Q_2, size2;
    double ro2_H1, p2_H1, u2_H1, v2_H1;
    double ro2_H2, p2_H2, u2_H2, v2_H2;
    double ro2_H3, p2_H3, u2_H3, v2_H3;
    double ro2_H4, p2_H4, u2_H4, v2_H4;
    double roC = 1.0;
    double pC = 1.0;
    double uC = Velosity_inf;
    double vC = 0.0;
    double QC = 100.0;
    double u1_polar, v1_polar;

    for (int i = l; i <= r; i++)
    {
        ii = SOSED[i];
        if (ii >= 0)
        {
            x2 = X[ii];
            y2 = Y[ii];
            size2 = Size[ii];
            dx2 = (DX / pow(2, size2 - 1)) / 2.0;   // ѕоловина длины €чейки
            dy2 = (DY / pow(2, size2 - 1)) / 2.0;   // ѕоловина ширины €чейки
            ro2 = RO1[ii];
            p2 = P1[ii];
            u2 = U1[ii];
            v2 = V1[ii];
            ro2_H1 = RO1_H1[ii];
            p2_H1 = P1_H1[ii];
            u2_H1 = U1_H1[ii];
            v2_H1 = V1_H1[ii];
            ro2_H2 = RO1_H2[ii];
            p2_H2 = P1_H2[ii];
            u2_H2 = U1_H2[ii];
            v2_H2 = V1_H2[ii];
            ro2_H3 = RO1_H3[ii];
            p2_H3 = P1_H3[ii];
            u2_H3 = U1_H3[ii];
            v2_H3 = V1_H3[ii];
            ro2_H4 = RO1_H4[ii];
            p2_H4 = P1_H4[ii];
            u2_H4 = U1_H4[ii];
            v2_H4 = V1_H4[ii];
            Q_2 = Q1[ii];
            double S = get_square(x, y, dx, dy, x2, y2, dx2, dy2, n1, n2, dist);

            u1_polar = u;
            v1_polar = v;

            if (radius < 100)
            {
                polar_perenos(x, y, x + n1 * dx, y + n2 * dy, u1_polar, v1_polar);
                polar_perenos(x2, y2, x2 - n1 * dx2, y2 - n2 * dy2, u2, v2);
            }

            if (radius > Distant)
            {
                tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H1, 1.0, p_H1, u_H1, v_H1, ro2_H1, 1.0, p2_H1, u2_H1, v2_H1, P, PQ, n1, n2, dist, metod));
                for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
                {
                    Potok_H1[k] = Potok_H1[k] + P[k] * S;
                }
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H2, 1.0, p_H2, u_H2, v_H2, ro2_H2, 1.0, p2_H2, u2_H2, v2_H2, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H2[k] = Potok_H2[k] + P[k] * S;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H3, 1.0, p_H3, u_H3, v_H3, ro2_H3, 1.0, p2_H3, u2_H3, v2_H3, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H3[k] = Potok_H3[k] + P[k] * S;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H4, 1.0, p_H4, u_H4, v_H4, ro2_H4, 1.0, p2_H4, u2_H4, v2_H4, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H4[k] = Potok_H4[k] + P[k] * S;
            }

            if (radius > Distant)
            {
                tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro, Q, p, u1_polar, v1_polar, ro2, Q_2, p2, u2, v2, P, PQ, n1, n2, dist, metod));

                for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[4] = Potok[4] + PQ * S;
            }
        }
        else if (ii == -1)
        {
            double S = dy * 2.0;
            n1 = 1.0;
            n2 = 0.0;
            dist = dx;

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H1, 1.0, p_H1, u_H1, v_H1, ro_H1, 1.0, p_H1, u_H1, v_H1, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H1[k] = Potok_H1[k] + P[k] * S;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H2, 1.0, p_H2, u_H2, v_H2, ro_H2, 1.0, p_H2, u_H2, v_H2, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H2[k] = Potok_H2[k] + P[k] * S;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H3, 1.0, p_H3, u_H3, v_H3, ro_H3, 1.0, p_H3, u_H3, v_H3, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H3[k] = Potok_H3[k] + P[k] * S;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H4, 1.0, p_H4, u_H4, v_H4, roC, 1.0, 0.5 * pC, Velosity_inf, 0.0, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H4[k] = Potok_H4[k] + P[k] * S;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro, Q, p, u, v, roC, QC, pC, uC, vC, P, PQ, n1, n2, dist, metod));

            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok[k] = Potok[k] + P[k] * S;
            }
            Potok[4] = Potok[4] + PQ * S;
        }
        else if (ii == -2)
        {
            double S = dy * 2.0;
            n1 = -1.0;
            n2 = 0.0;
            dist = dx;

            double uu = u;
            if (uu > Velosity_inf && y < 300)
            {
                uu = Velosity_inf;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H1, 1.0, p_H1, u_H1, v_H1, ro_H1, 1.0, p_H1, u_H1, v_H1, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H1[k] = Potok_H1[k] + P[k] * S;
            }
            

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H2, 1.0, p_H2, u_H2, v_H2, ro_H2, 1.0, p_H2, u_H2, v_H2, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H2[k] = Potok_H2[k] + P[k] * S;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H3, 1.0, p_H3, u_H3, v_H3, ro_H3, 1.0, p_H3, u_H3, v_H3, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H3[k] = Potok_H3[k] + P[k] * S;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H4, 1.0, p_H4, u_H4, v_H4, ro_H4, 1.0, p_H4, u_H4, v_H4, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H4[k] = Potok_H4[k] + P[k] * S;
            }


            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro, Q, p, u, v, ro, Q, p, uu, v, P, PQ, n1, n2, dist, metod));

            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok[k] = Potok[k] + P[k] * S;
            }
            Potok[4] = Potok[4] + PQ * S;
            
        }
        else if (ii == -3)
        {
            double S = dx * 2.0;
            n1 = 0.0;
            n2 = 1.0;
            dist = dy;

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H1, 1.0, p_H1, u_H1, v_H1, ro_H1, 1.0, p_H1, u_H1, v_H1, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H1[k] = Potok_H1[k] + P[k] * S;
            }
            

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H2, 1.0, p_H2, u_H2, v_H2, ro_H2, 1.0, p_H2, u_H2, v_H2, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H2[k] = Potok_H2[k] + P[k] * S;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H3, 1.0, p_H3, u_H3, v_H3, ro_H3, 1.0, p_H3, u_H3, v_H3, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H3[k] = Potok_H3[k] + P[k] * S;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H4, 1.0, p_H4, u_H4, v_H4, ro_H4, 1.0, p_H4, u_H4, v_H4, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H4[k] = Potok_H4[k] + P[k] * S;
            }


            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro, Q, p, u, v, ro, Q, p, u, v, P, PQ, n1, n2, dist, metod));

            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok[k] = Potok[k] + P[k] * S;
            }
            Potok[4] = Potok[4] + PQ * S;
            
        }
        else if (ii == -4)
        {
            double S = dx * 2.0;
            n1 = 0.0;
            n2 = -1.0;
            dist = dy;

            u1_polar = u;
            v1_polar = v;

            if (radius < 100)
            {
                u1_polar = u;
                v1_polar = v;
                polar_perenos(x, y, x + n1 * dx, y + n2 * dy, u1_polar, v1_polar);
            }

            if (radius > Distant)
            {
                tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H1, 1.0, p_H1, u_H1, v_H1, ro_H1, 1.0, p_H1, u_H1, -v_H1, P, PQ, n1, n2, dist, metod));
                for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
                {
                    Potok_H1[k] = Potok_H1[k] + P[k] * S;
                }
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H2, 1.0, p_H2, u_H2, v_H2, ro_H2, 1.0, p_H2, u_H2, -v_H2, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H2[k] = Potok_H2[k] + P[k] * S;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H3, 1.0, p_H3, u_H3, v_H3, ro_H3, 1.0, p_H3, u_H3, -v_H3, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H3[k] = Potok_H3[k] + P[k] * S;
            }

            tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro_H4, 1.0, p_H4, u_H4, v_H4, ro_H4, 1.0, p_H4, u_H4, -v_H4, P, PQ, n1, n2, dist, metod));
            for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
            {
                Potok_H4[k] = Potok_H4[k] + P[k] * S;
            }

            if (radius > Distant)
            {
                tmin = min(tmin, HLLC_2d_Korolkov_b_s(ro, Q, p, u1_polar, v1_polar, ro, Q, p, //
                    u1_polar, -v1_polar, P, PQ, n1, n2, dist, metod));

                for (int k = 0; k < 4; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[4] = Potok[4] + PQ * S;
            }
        }
        else
        {
            printf("Error 12438wedew4353jdyu. Ne doljni suda popadat = %d \n", ii);
        }
    }


    double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
    double U_H1, U_H2, U_H3, U_H4;
    double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
    double nu_H1, nu_H2, nu_H3, nu_H4;
    double q2_1, q2_2, q3;

    /// «десь остановились вчера!!!

    double ro3, p3, u3, v3, Q33;

    ro3 = ro - *T_do * (Potok[0] / Volume + ro * v / y);
    Q33 = Q - (*T_do / Volume) * Potok[4] - *T_do * Q * v / y;
    if (ro3 <= 0)
    {
        printf("Problemsssss  ro < 0! %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", x, y, dx, dy, ro, p, u, v, Q);
        ro3 = 0.00001;
    }
    u3 = (ro * u - *T_do * (Potok[1] / Volume + ro * v * u / y)) / ro3;
    v3 = (ro * v - *T_do * (Potok[2] / Volume + ro * v * v / y)) / ro3;
    p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - *T_do * (Potok[3] / Volume + //
        +v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y)) - //
        0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
    if (p3 <= 0)
    {
        p3 = 0.000001;
    }

    Q2[index] = Q33;
    RO2[index] = ro3;
    P2[index] = p3;
    U2[index] = u3;
    V2[index] = v3;

    if (*T > tmin)
    {
        *T = tmin;
        __threadfence();
    }
    

}

