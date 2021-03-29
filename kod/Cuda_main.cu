#include "Header.h"
#include "math.h"


/// ¬се функции описанные в этой файле должны иметь прототипы в файле Header.h

__global__ void Cuda_main_HLLDQ(int* NN, double* X, double* Y, int* Size,//
    double* RO1, double* RO2, double* Q1, double* Q2, double* P1, double* P2, double* U1, double* U2, double* V1, double* V2,//
    int* SOSED, int* L, int* R, double* T, double* T_do, int step_, double M_inf_, const double& DX, const double& DY, int metod = 0)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x; // √лобальный индекс текущего потока
    if (index > * NN - 1)
    {
        return;
    }
    double x, y, z, dx, dy, ro, p, u, v, Q;
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

    //double ddd = kv(y) + kv(z);
    //double ddd2 = kv(x + 0.8) + kv(y) + kv(z);
    //double dist3 = kv(x + 1.0) / kv(1.6) + kv(y) / kv(1.6) + kv(z) / kv(1.6);
    double dist3 = kv(x + 1.08) / kv(2.4) + kv(y) / kv(2.0) + kv(z) / kv(2.0);

    if (dist3 < 1.0) // || (ddd <= 4.0 && x > -5 && x < 0) ) //(ddd < 5.76 || ddd2 <= 2.0) //1.5
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
        double P[8] = { 0.0 };
        P[0] = P[1] = P[2] = P[3] = P[4] = P[5] = P[6] = P[7] = 0.0;
        double Potok[10] = { 0.0 };
        Potok[0] = Potok[1] = Potok[2] = Potok[3] = Potok[4] = Potok[5] = Potok[6] = Potok[7] = Potok[8] = Potok[9] = 0.0;
        double tmin = 1000;
        double Volume = dx * dy * 4.0;
        int ii = 0;
        double x2, y2, dx2, dy2, ro2, p2, u2, v2, sks, Q_2, size2;
        double roC = 1.0; // 8.2598; //  1.0;
        double pC = 1.0 / (ggg); // 1.0 / (ggg * M_inf * M_inf);
        double uC = -M_inf_; // -1.0;
        double vC = 0.0;
        double QC = 100.0;


        for (int i = l; i <= r; i++)
        {
            ii = SOSED[i];
            if (ii >= 0)
            {
                x2 = X[ii];
                y2 = Y[ii];
                size2 = Size[ii];
                double dx2 = (DX / pow(2, size - 1)) / 2.0;   // ѕоловина длины €чейки
                double dy2 = (DY / pow(2, size - 1)) / 2.0;   // ѕоловина ширины €чейки
                ro2 = RO1[ii];
                p2 = P1[ii];
                u2 = U1[ii];
                v2 = V1[ii];
                Q_2 = Q1[ii];
                double S = get_square(x, y, z, dx, dy, dz, x2, y2, z2, dx2, dy2, dz2, n1, n2, n3, dist);
                if (diver == true)
                {
                    sks = n1 * (bx + bx2) / 2.0 + n2 * (by + by2) / 2.0 + n3 * (bz + bz2) / 2.0;
                }
                else
                {
                    sks = 0.0;
                }
                Potok[8] = Potok[8] + sks * S;

                if (!kor_Sol || metod == 1)//(y * y + z * z < 225 && y2 * y2 + z2 * z2 < 225 && x > -15 && x2 > -15 && x < 8 && x2 < 8  && step_ > 10000)
                {
                    tmin = min(tmin, HLLDQ_Alexashov(ro, Q, p, u, v, w, bx, by, bz, ro2, Q_2, p2, u2, v2, w2, bx2, by2, bz2, P, PQ, n1, n2, n3, dist, metod));
                }
                else
                {
                    tmin = min(tmin, HLLDQ_Korolkov(ro, Q, p, u, v, w, bx, by, bz, ro2, Q_2, p2, u2, v2, w2, bx2, by2, bz2, P, PQ, n1, n2, n3, dist, metod));
                }
                for (int k = 0; k < 8; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[9] = Potok[9] + PQ * S;
            }
            else if (ii == -1)
            {
                double S = dy * dz * 4.0;
                n1 = 1.0;
                n2 = 0.0;
                n3 = 0.0;
                dist = dx;
                if (diver == true)
                {
                    sks = n1 * (bx + bxC) / 2.0 + n2 * (by + byC) / 2.0 + n3 * (bz + bzC) / 2.0;
                }
                else
                {
                    sks = 0.0;
                }
                /*double uu = u;
                if (uu < 0.0)
                {
                    uu = 0.0;
                }*/
                Potok[8] = Potok[8] + sks * S;
                if (!kor_Sol || metod == 1 || metod == 3)
                {
                    tmin = min(tmin, HLLDQ_Alexashov(ro, Q, p, u, v, w, bx, by, bz, roC, QC, pC, uC, vC, wC, bxC, byC, bzC, P, PQ, n1, n2, n3, dist, metod));
                }
                else
                {
                    tmin = min(tmin, HLLDQ_Korolkov(ro, Q, p, u, v, w, bx, by, bz, roC, QC, pC, uC, vC, wC, bxC, byC, bzC, P, PQ, n1, n2, n3, dist, metod));
                }
                //  ћожно вручную выписать потоки дл€ ускорени€ времени
                /*double b2R = kv(bxC) + kv(byC) + kv(bzC);
                double ptR = pC + b2R / 2.0;
                double upt2 = (kv(uC) + kv(vC) + kv(wC)) / 2.0;
                double sbv2 = uC * bxC + vC * byC + wC * bzC;
                double e2 = pC / g1 + roC * upt2 + b2R / 2.0;

                P[0] = roC * uC;
                P[1] = roC * uC * uC + ptR - kv(bxC);
                P[2] = roC * uC * vC - bxC * byC;
                P[3] = roC * uC * wC - bxC * bzC;
                P[7] = (e2 + ptR) * uC - bxC * sbv2;
                P[4] = 0.0;
                P[5] = uC * byC - vC * bxC;
                P[6] = uC * bzC - wC * bxC;*/

                for (int k = 0; k < 8; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[9] = Potok[9] + PQ * S;
            }
            else if (ii == -2)
            {
                double S = dy * dz * 4.0;
                n1 = -1.0;
                n2 = 0.0;
                n3 = 0.0;
                dist = dx;
                if (diver == true)
                {
                    sks = n1 * bx + n2 * by + n3 * bz;
                }
                else
                {
                    sks = 0.0;
                }
                Potok[8] = Potok[8] + sks * S;
                double uu = u;
                if (uu > -M_inf_ && step_ < 5000)
                {
                    uu = -M_inf_;
                }
                else if (uu > -0.01)
                {
                    uu = -0.01;
                }

                if (!kor_Sol || metod == 1 || metod == 3)
                {
                    tmin = min(tmin, HLLDQ_Alexashov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, p, uu, v, w, bx, by, bz, P, PQ, n1, n2, n3, dist, metod));
                }
                else
                {
                    tmin = min(tmin, HLLDQ_Korolkov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, p, uu, v, w, bx, by, bz, P, PQ, n1, n2, n3, dist, metod));
                }

                /*double t1, t2, t3, m1, m2, m3;
                double bx_L = bx / spi4;
                double by_L = by / spi4;
                double bz_L = bz / spi4;
                t1 = 0.0;
                t2 = 0.0;
                t3 = 1.0;
                m1 = 0.0;
                m2 = 1.0;
                m3 = 0.0;
                double u1 = uu * n1 + v * n2 + w * n3;
                double v1 = uu * t1 + v * t2 + w * t3;
                double w1 = uu * m1 + v * m2 + w * m3;
                double bn1, bt1, bm1;
                bn1 = bx_L * n1 + by_L * n2 + bz_L * n3;
                bt1 = bx_L * t1 + by_L * t2 + bz_L * t3;
                bm1 = bx_L * m1 + by_L * m2 + bz_L * m3;
                double uu_L = (kv(uu) + kv(v) + kv(w)) / 2.0;
                double bb_L = kv(bx_L) + kv(by_L) + kv(bz_L);
                double e1 = p / g1 + ro * uu_L + bb_L / 2.0;
                double pTL = p + bb_L / 2.0;

                double PO[9];

                PO[0] = ro * u1;
                PO[1] = ro * u1 * u1 + pTL - kv(bn1);
                PO[2] = ro * u1 * v1 - bn1 * bt1;
                PO[3] = ro * u1 * w1 - bn1 * bm1;
                PO[4] = (e1 + pTL) * u1 - bn1 * (u1 * bn1 + v1 * bt1 + w1 * bm1);
                PO[5] = 0.0;
                PO[6] = u1 * bt1 - v1 * bn1;
                PO[7] = u1 * bm1 - w1 * bn1;
                PO[8] = Q * u1;


                P[1] = n1 * PO[1] + t1 * PO[2] + m1 * PO[3];
                P[2] = n2 * PO[1] + t2 * PO[2] + m2 * PO[3];
                P[3] = n3 * PO[1] + t3 * PO[2] + m3 * PO[3];
                P[5] = spi4 * (n1 * PO[5] + t1 * PO[6] + m1 * PO[7]);
                P[6] = spi4 * (n2 * PO[5] + t2 * PO[6] + m2 * PO[7]);
                P[7] = spi4 * (n3 * PO[5] + t3 * PO[6] + m3 * PO[7]);
                P[0] = PO[0];
                P[4] = PO[4];
                PQ = PO[8];

                double SWAP = P[4];
                P[4] = P[5];
                P[5] = P[6];
                P[6] = P[7];
                P[7] = SWAP;*/

                for (int k = 0; k < 8; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[9] = Potok[9] + PQ * S;
            }
            else if (ii == -3)
            {
                double S = dx * dz * 4.0;
                n1 = 0.0;
                n2 = 1.0;
                n3 = 0.0;
                dist = dy;
                double uu = v;
                if (uu < 0.0)
                {
                    uu = 0.0;
                }
                if (diver == true)
                {
                    sks = n1 * (bx + bx) / 2.0 + n2 * (by + by) / 2.0 + n3 * (bz + bz) / 2.0;
                }
                else
                {
                    sks = 0.0;
                }
                Potok[8] = Potok[8] + sks * S;
                if (!kor_Sol || metod == 1 || metod == 3)
                {
                    tmin = min(tmin, HLLDQ_Alexashov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, pC, u, v, w, bx, by, bz, P, PQ, n1, n2, n3, dist, metod));
                }
                else
                {
                    tmin = min(tmin, HLLDQ_Korolkov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, pC, u, v, w, bx, by, bz, P, PQ, n1, n2, n3, dist, metod));
                }
                for (int k = 0; k < 8; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[9] = Potok[9] + PQ * S;
            }
            else if (ii == -4)
            {
                double S = dx * dz * 4.0;
                n1 = 0.0;
                n2 = -1.0;
                n3 = 0.0;
                dist = dy;
                if (diver == true)
                {
                    sks = n1 * (bx + bx) / 2.0 + n2 * (by + by) / 2.0 + n3 * (bz + bz) / 2.0;
                    //sks =  n2 * (by + by) / 2.0;
                }
                else
                {
                    sks = 0.0;
                }

                Potok[8] = Potok[8] + sks * S;
                if (!kor_Sol || metod == 1 || metod == 3)
                {
                    //tmin = min(tmin, HLLDQ_Alexashov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, p, u, -v, w, -bx, by, -bz, P, PQ, n1, n2, n3, dist, metod));
                    tmin = min(tmin, HLLDQ_Alexashov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, p, u, v, w, bx, by, bz, P, PQ, n1, n2, n3, dist, metod));
                }
                else
                {
                    //tmin = min(tmin, HLLDQ_Korolkov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, p, u, -v, w, -bx, by, -bz, P, PQ, n1, n2, n3, dist, metod));
                    tmin = min(tmin, HLLDQ_Korolkov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, p, u, v, w, bx, by, bz, P, PQ, n1, n2, n3, dist, metod));
                }

                for (int k = 0; k < 8; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[9] = Potok[9] + PQ * S;
            }
            else if (ii == -5)
            {
                double S = dy * dx * 4.0;
                n1 = 0.0;
                n2 = 0.0;
                n3 = 1.0;
                dist = dz;
                if (diver == true)
                {
                    sks = n1 * (bx + bx) / 2.0 + n2 * (by + by) / 2.0 + n3 * (bz + bz) / 2.0;
                }
                else
                {
                    sks = 0.0;
                }


                Potok[8] = Potok[8] + sks * S;
                if (!kor_Sol || metod == 1 || metod == 3)
                {
                    tmin = min(tmin, HLLDQ_Alexashov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, pC, u, v, w, bx, by, bz, P, PQ, n1, n2, n3, dist, metod));
                }
                else
                {
                    tmin = min(tmin, HLLDQ_Korolkov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, pC, u, v, w, bx, by, bz, P, PQ, n1, n2, n3, dist, metod));
                }

                for (int k = 0; k < 8; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[9] = Potok[9] + PQ * S;
            }
            else if (ii == -6)
            {
                double S = dy * dx * 4.0;
                n1 = 0.0;
                n2 = 0.0;
                n3 = -1.0;
                dist = dz;
                if (diver == true)
                {
                    //sks = n1 * (bx + bx) / 2.0 + n2 * (by + by) / 2.0 + n3 * (bz + bz) / 2.0;
                    sks = n1 * (bx + bx) / 2.0 + n2 * (by + by) / 2.0;
                }
                else
                {
                    sks = 0.0;
                }


                Potok[8] = Potok[8] + sks * S;
                if (!kor_Sol || metod == 1 || metod == 3)
                {
                    tmin = min(tmin, HLLDQ_Alexashov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, p, u, v, -w, bx, by, -bz, P, PQ, n1, n2, n3, dist, metod));
                    //tmin = min(tmin, HLLDQ_Alexashov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, pC, u, v, w, bx, by, bz, P, PQ, n1, n2, n3, dist, metod));
                }
                else
                {
                    tmin = min(tmin, HLLDQ_Korolkov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, p, u, v, -w, bx, by, -bz, P, PQ, n1, n2, n3, dist, metod));
                    //tmin = min(tmin, HLLDQ_Korolkov(ro, Q, p, u, v, w, bx, by, bz, ro, Q, pC, u, v, w, bx, by, bz, P, PQ, n1, n2, n3, dist, metod));
                }
                for (int k = 0; k < 8; k++)  // —уммируем все потоки в €чейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }
                Potok[9] = Potok[9] + PQ * S;
            }
            else
            {
                printf("Error 12438jdyu. Ne doljni suda popadat = %d \n", ii);
            }
        }

        double ro3, p3, u3, v3, w3, bx3, by3, bz3, Q33;

        Q33 = Q - *T_do * Potok[9] / Volume;
        ro3 = ro - *T_do * Potok[0] / Volume;
        if (ro3 <= 0.0)
        {
            printf("ERROR -  dssdbfhfshjskfutytqqazz\n");
            printf("%lf, %lf, %lf, %lf\n", x, y, z, ro3);
            ro3 = ro;
        }
        u3 = (ro * u - *T_do * (Potok[1] + (bx / cpi4) * Potok[8]) / Volume) / ro3;
        v3 = (ro * v - *T_do * (Potok[2] + (by / cpi4) * Potok[8]) / Volume) / ro3;
        w3 = (ro * w - *T_do * (Potok[3] + (bz / cpi4) * Potok[8]) / Volume) / ro3;
        bx3 = (bx - *T_do * (Potok[4] + u * Potok[8]) / Volume);
        by3 = (by - *T_do * (Potok[5] + v * Potok[8]) / Volume);
        bz3 = (bz - *T_do * (Potok[6] + w * Potok[8]) / Volume);
        p3 = ((U8(ro, p, u, v, w, bx, by, bz) - *T_do * (Potok[7] + (skk(u, v, w, bx, by, bz) / cpi4) * Potok[8])//
            / Volume) - 0.5 * ro3 * kvv(u3, v3, w3) - kvv(bx3, by3, bz3) / cpi8) * (ggg - 1.0);
        //u3 = (ro * u - *T_do * (Potok[1] + (bx) * Potok[8]) / Volume) / ro3;
        //v3 = (ro * v - *T_do * (Potok[2] + (by) * Potok[8]) / Volume) / ro3;
        //w3 = (ro * w - *T_do * (Potok[3] + (bz) * Potok[8]) / Volume) / ro3;
        //bx3 = (bx - *T_do * (Potok[4] + u * Potok[8]) / Volume);
        //by3 = (by - *T_do * (Potok[5] + v * Potok[8]) / Volume);
        //bz3 = (bz - *T_do * (Potok[6] + w * Potok[8]) / Volume);
        //p3 = ((U8(ro, p, u, v, w, bx, by, bz) - *T_do * (Potok[7] + (skk(u, v, w, bx, by, bz)) * Potok[8])//
        //    / Volume) - 0.5 * ro3 * kvv(u3, v3, w3) - kvv(bx3, by3, bz3)) * (ggg - 1.0);
        if (p3 <= 0)
        {
            p3 = 0.000001;
        }

        Q2[index] = Q33;
        RO2[index] = ro3;
        P2[index] = p3;
        U2[index] = u3;
        V2[index] = v3;
        W2[index] = w3;
        /*if (Q33 / ro3 > 50)
        {
            BX2[index] = 0.0;
            BY2[index] = 0.0;
            BZ2[index] = 0.0;
        }
        else
        {
            BX2[index] = bx3;
            BY2[index] = by3;
            BZ2[index] = bz3;
        }*/
        BX2[index] = bx3;
        BY2[index] = by3;
        BZ2[index] = bz3;

        if (*T > tmin)
        {
            *T = tmin;
            __threadfence();
        }
    }

}

