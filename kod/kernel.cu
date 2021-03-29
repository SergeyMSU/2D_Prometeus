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
#define TVD_ false
#define TVQ_ true
#define kor_Sol true

#define sss 500000000

using namespace std;

cudaError_t addWithCuda(void);


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
    K.print_konectiviti_short();
    K.print_point();
    K.print_cell();


    return cudaStatus;
}

