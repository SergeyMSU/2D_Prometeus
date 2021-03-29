#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "Cell.h"
#include "Header.h"
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include "Konstruktor.h"
#include "Cuda_main.cu"
#include "Kyb.h"

#define ER_S(x) printf("Standart error in kernel.cu: kod - x\n")


using namespace std;

cudaError_t addWithCuda(void);
__global__ void funk_time(double* T, double* T_do, double* TT, int* i);



__global__ void funk_time(double* T, double* T_do, double* TT, int* i)
{
    *T_do = *T;
    *TT = *TT + *T_do;
    *T = 10000000;
    *i = *i + 1;
    if (*i % 1000 == 0)
    {
        printf("i = %d,  TT = %lf \n", *i, *TT);
    }
    return;
}


int main()
{
    if (true)
    {
        // Add vectors in parallel.
        cudaError_t cudaStatus = addWithCuda();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "addWithCuda failed!");
            return 1;
        }

        // cudaDeviceReset must be called before exiting in order for profiling and
        // tracing tools such as Nsight and Visual Profiler to show complete traces
        // Эта штука должна быть вызвана
        cudaStatus = cudaDeviceReset();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceReset failed!\n");
            return 1;
        }
    }


    return 0;
}


cudaError_t addWithCuda()
{
    cudaError_t cudaStatus = cudaSuccess;

    Konstruktor K(30, 30, -100, 100, 100);
    K.Drobim(0, 0, 0, 50);
    K.Drobim(0, 0, 30, 80);

    int N = K.all_Kyb.size();          // Число ячеек
    cout << "All size = " << N << endl;
    int nn = K.get_size_conektiv();    // Число связей (размер массива связей)
    cout << "Connect = " << nn << endl;


    K.initial_condition();   // Заполнение начальными условиями


    int* host_sosed;
    int* dev_sosed;
    double* host_T, * host_T_do, * host_TT;
    double* dev_T, * dev_T_do, * dev_TT;
    int* host_i;
    int* dev_i;

    // Создаём массивы переменных\ячеек

    int* dev_l, * dev_r, * Nu, * dev_N;
    double* dev_x, * dev_y, * dev_z;
    double* dev_dx, * dev_dy, * dev_dz;
    double* dev_Q1, * host_Q1, * dev_Q2;
    double* dev_ro1, * dev_p1, * dev_u1, * dev_v1,  * dev_ro2, * dev_p2, * dev_u2, * dev_v2;
    double* host_x, * host_y;
    double* host_ro1, * host_p1, * host_u1, * host_v1;
    int* host_l, * host_r;

    host_T = (double*)malloc(sizeof(double));
    host_T_do = (double*)malloc(sizeof(double));
    host_TT = (double*)malloc(sizeof(double));
    host_i = (int*)malloc(sizeof(int));
    Nu = (int*)malloc(sizeof(int));

    int NNN;
    if (N % 256 == 0)
    {
        NNN = N / 256;
    }
    else
    {
        NNN = (int)(N / 256) + 1;
    }

    host_x = new double[N];
    host_y = new double[N];
    host_ro1 = new double[N];
    host_Q1 = new double[N];
    
    host_p1 = new double[N];
    host_u1 = new double[N];
    host_v1 = new double[N];
    host_l = new int[N];
    host_r = new int[N];

    *host_T = 10000000.0;
    *host_T_do = 0.00000001;
    *host_TT = 0.0;
    *host_i = 0;
    *Nu = N;

    host_sosed = new int[nn];

    // Заполнение массивов
    int c = 0;
    for (auto& i : K.all_Kyb)
    {
        for (auto& j : i->sosed)
        {
            host_sosed[c] = j->number;
            c++;
        }
    }

    int ll = 0;
    for (int i = 0; i < K.all_Kyb.size(); i++)
    {
        host_x[i] = K.all_Kyb[i]->x;
        host_y[i] = K.all_Kyb[i]->y;
        host_ro1[i] = K.all_Kyb[i]->ro;
        host_Q1[i] = K.all_Kyb[i]->Q;
        host_p1[i] = K.all_Kyb[i]->p;
        host_u1[i] = K.all_Kyb[i]->u;
        host_v1[i] = K.all_Kyb[i]->v;
        host_l[i] = ll;
        host_r[i] = ll + K.all_Kyb[i]->sosed.size() - 1;
        ll = ll + K.all_Kyb[i]->sosed.size();
    }

    // Выделение памяти на девайсе
    if (true)
    {
        cudaStatus = cudaMalloc((void**)&dev_x, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_y, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 2!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_z, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 3!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_dx, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_dy, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_dz, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_ro1, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_Q1, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_Q2, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }
        

        cudaStatus = cudaMalloc((void**)&dev_ro2, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_p1, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_p2, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_u1, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_u2, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_v1, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_v2, N * sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_l, N * sizeof(int));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_r, N * sizeof(int));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 1!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_sosed, nn * sizeof(int));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 2!");
            goto Error;
        }


        cudaStatus = cudaMalloc((void**)&dev_T, sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 3!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_T_do, sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 4!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_TT, sizeof(double));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 5!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_i, sizeof(int));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 6!");
            goto Error;
        }

        cudaStatus = cudaMalloc((void**)&dev_N, sizeof(int));
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed 6!");
            goto Error;
        }
    }


    // Копируем массивы с хоста на девайс
    if (true)
    {
        cudaStatus = cudaMemcpy(dev_sosed, host_sosed, nn * sizeof(int), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed -1 !");
            goto Error;
        }

        cudaStatus = cudaMemcpy(dev_x, host_x, N * sizeof(double), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 0!");
            goto Error;
        }

        cudaStatus = cudaMemcpy(dev_x, host_x, N * sizeof(double), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 0!");
            goto Error;
        }
        cudaStatus = cudaMemcpy(dev_y, host_y, N * sizeof(double), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 0!");
            goto Error;
        }
        cudaStatus = cudaMemcpy(dev_ro1, host_ro1, N * sizeof(double), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 0!");
            goto Error;
        }
        cudaStatus = cudaMemcpy(dev_Q1, host_Q1, N * sizeof(double), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 012!");
            goto Error;
        }
        
        cudaStatus = cudaMemcpy(dev_p1, host_p1, N * sizeof(double), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 0!");
            goto Error;
        }
        cudaStatus = cudaMemcpy(dev_u1, host_u1, N * sizeof(double), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 0!");
            goto Error;
        }
        cudaStatus = cudaMemcpy(dev_v1, host_v1, N * sizeof(double), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 0!");
            goto Error;
        }
        cudaStatus = cudaMemcpy(dev_l, host_l, N * sizeof(int), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 0!");
            goto Error;
        }
        cudaStatus = cudaMemcpy(dev_r, host_r, N * sizeof(int), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 0!");
            goto Error;
        }

        cudaStatus = cudaMemcpy(dev_T, host_T, sizeof(double), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 1!");
            goto Error;
        }

        cudaStatus = cudaMemcpy(dev_TT, host_TT, sizeof(double), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 2!");
            goto Error;
        }

        cudaStatus = cudaMemcpy(dev_T_do, host_T_do, sizeof(double), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 3!");
            goto Error;
        }

        cudaStatus = cudaMemcpy(dev_i, host_i, sizeof(int), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 4!");
            goto Error;
        }

        cudaStatus = cudaMemcpy(dev_N, Nu, sizeof(int), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed 4!");
            goto Error;
        }
    }

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }


    for (int i = 0; i < 100000; i = i + 2)  // Сколько шагов по времени делаем?
    {
        // запускаем add() kernel на GPU, передавая параметры
        Cuda_main_HLLDQ << <NNN, 256 >> > (dev_N, dev_x, dev_y, dev_z, dev_dx, dev_dy, dev_dz,//
            dev_ro1, dev_ro2, dev_Q1, dev_Q2, dev_p1, dev_p2, dev_u1, dev_u2, dev_v1, dev_v2,//
            dev_w1, dev_w2, dev_bx1, dev_by1, dev_bz1, dev_bx2, dev_by2, dev_bz2,//
            dev_sosed, dev_l, dev_r, dev_T, dev_T_do, i, MMM, true, true, 1);

        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel! 11111\n", cudaStatus);
            goto Error;
        }

        funk_time << <1, 1 >> > (dev_T, dev_T_do, dev_TT, dev_i);
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel! 222222\n", cudaStatus);
            goto Error;
        }

        Cuda_main_HLLDQ << <NNN, 256 >> > (dev_N, dev_x, dev_y, dev_z, dev_dx, dev_dy, dev_dz,//
            dev_ro2, dev_ro1, dev_Q2, dev_Q1, dev_p2, dev_p1, dev_u2, dev_u1, dev_v2, dev_v1,//
            dev_w2, dev_w1, dev_bx2, dev_by2, dev_bz2, dev_bx1, dev_by1, dev_bz1,//
            dev_sosed, dev_l, dev_r, dev_T, dev_T_do, i, MMM, true, true, 1);

        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel! 333333\n", cudaStatus);
            goto Error;
        }

        funk_time << <1, 1 >> > (dev_T, dev_T_do, dev_TT, dev_i);
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel! 4444444\n", cudaStatus);
            goto Error;
        }



        if ((i % 15000000000 == 0 && i > 1))
        {
            cout << "HLL "  << endl;
            if (true)
            {
                cudaStatus = cudaMemcpy(host_ro1, dev_ro1, N * sizeof(double), cudaMemcpyDeviceToHost);
                if (cudaStatus != cudaSuccess) {
                    fprintf(stderr, "cudaMemcpy failed!  3452\n");
                    goto Error;
                }
                cudaStatus = cudaMemcpy(host_p1, dev_p1, N * sizeof(double), cudaMemcpyDeviceToHost);
                if (cudaStatus != cudaSuccess) {
                    fprintf(stderr, "cudaMemcpy failed!  3452\n");
                    goto Error;
                }
                cudaStatus = cudaMemcpy(host_u1, dev_u1, N * sizeof(double), cudaMemcpyDeviceToHost);
                if (cudaStatus != cudaSuccess) {
                    fprintf(stderr, "cudaMemcpy failed!  3452\n");
                    goto Error;
                }
                cudaStatus = cudaMemcpy(host_v1, dev_v1, N * sizeof(double), cudaMemcpyDeviceToHost);
                if (cudaStatus != cudaSuccess) {
                    fprintf(stderr, "cudaMemcpy failed!  3452\n");
                    goto Error;
                }
                cudaStatus = cudaMemcpy(host_Q1, dev_Q1, N * sizeof(double), cudaMemcpyDeviceToHost);
                if (cudaStatus != cudaSuccess) {
                    fprintf(stderr, "cudaMemcpy failed!  3452\n");
                    goto Error;
                }
            }
            K.read_Cuda_massiv(host_ro1, host_p1, host_u1, host_v1, host_Q1);
            K.print_Tecplot();
        }
    }


    // Копирование массивов обратно
    if (true)
    {
        cudaStatus = cudaMemcpy(host_ro1, dev_ro1, N * sizeof(double), cudaMemcpyDeviceToHost);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed!  3452\n");
            goto Error;
        }
        cudaStatus = cudaMemcpy(host_p1, dev_p1, N * sizeof(double), cudaMemcpyDeviceToHost);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed!  3452\n");
            goto Error;
        }
        cudaStatus = cudaMemcpy(host_u1, dev_u1, N * sizeof(double), cudaMemcpyDeviceToHost);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed!  3452\n");
            goto Error;
        }
        cudaStatus = cudaMemcpy(host_v1, dev_v1, N * sizeof(double), cudaMemcpyDeviceToHost);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed!  3452\n");
            goto Error;
        }
        cudaStatus = cudaMemcpy(host_Q1, dev_Q1, N * sizeof(double), cudaMemcpyDeviceToHost);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed!  3452\n");
            goto Error;
        }
    }  

    K.read_Cuda_massiv(host_ro1, host_p1, host_u1, host_v1,  host_Q1);
    K.print_Tecplot();



Error:
    cudaFree(dev_sosed);
    cudaFree(dev_ro1);
    cudaFree(dev_ro2);
    cudaFree(dev_p1);
    cudaFree(dev_p2);
    cudaFree(dev_u1);
    cudaFree(dev_u2);
    cudaFree(dev_v1);
    cudaFree(dev_v2);
    cudaFree(dev_T);
    cudaFree(dev_i);
    cudaFree(dev_T_do);
    cudaFree(dev_TT);
    cudaFree(dev_Q1);
    cudaFree(dev_Q2);







    return cudaStatus;
}

