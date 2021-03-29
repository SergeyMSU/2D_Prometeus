#pragma once
#ifndef CELL_H
#define CELL_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>

#define M_inf  2.2 // 0.7 // 0.4 // 0.8
#define M_infty  M_inf
#define phi_0  1.627 // 17.0 // 4.878 //1.627
#define alpha  45.0
#define M_alf  12.927 // 12.0
#define epsilon_ (1.0/M_alf)
#define step  70000
#define omega 0.0 //6 //1600
#define M_0  1.1

#define kk_ 196.0
#define betta 0.0 //4.2426
#define chi 2.0
#define rr_0  0.01005127

#define kurant  0.2

//#define Max_E  6.0
//#define M_A  12.87
//#define M_infty  1.97
//#define kappa_  5.3743
//#define rr_0  0.00208527 // 0.00329761
//#define mu 0.2837

#define ga (5.0/3.0)          // Показатель адиабаты
#define ggg (5.0/3.0)
#define kv(x) ( (x)*(x) )
#define kvv(x,y,z)  (kv(x) + kv(y) + kv(z))

//#define U8(ro, p, u, v, w, bx, by, bz)  ( (p) / (ggg - 1.0) + 0.5 * (ro) * kvv(u,v,w) + kvv(bx,by,bz))
#define U8(ro, p, u, v, w, bx, by, bz)  ( (p) / (ggg - 1.0) + 0.5 * (ro) * kvv(u,v,w) + kvv(bx,by,bz) / 25.13274122871834590768)
#define skk(u,v,w,bx,by,bz) ( (u)*(bx) + (v)*(by) + (w)*(bz) )
#define g1 (ga - 1.0)
#define gg1 (ga - 1.0)
#define g2 (ga + 1.0)
#define gg2 (ga + 1.0)
#define gp ((g2/ga)/2.0)
#define gm ((g1/ga)/2.0)
#define gga ga
#define Omega 0.0
#define geo 0.00001


#define eps 10e-10
#define eps8 10e-8
#define pi 3.14159265358979323846
#define PI 3.14159265358979323846
#define cpi4 12.56637061435917295384
#define cpi8 25.13274122871834590768
#define spi4 ( 3.544907701811032 )
#define epsb 1e-6
#define eps_p 1e-6
#define eps_d 1e-3


struct Cell // Класс ячейка
{
	double ro[2];
	double p[2];
	double u[2];
	double v[2];
	double w[2];
	double Bx[2];
	double By[2];
	double Bz[2];
	double dx;               // Половина ширины ячейки
	double dy;               // Половина глубины ячейки
	double dz;               // Половина высоты ячейки
	double x;                // Координата центра ячейки
	double y;                // Координата центра ячейки
	double z;                // Координата центра ячейки
	int l;
	int r;
};



#endif // APPROXIMATION_H