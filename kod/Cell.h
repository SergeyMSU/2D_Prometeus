#pragma once
#ifndef CELL_H
#define CELL_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>

// Параметры под задачу. Только их можно менять
#define Kn_  242.944										// Число Кнудсена
#define Velosity_inf -2.54186							// Скорость плазмы на бесконечности
#define a_2 0.102584
#define sigma(x) (kv(1.0 - a_2 * log(x)))               // Диффернциальное сечение перезарядки
#define Distant 10										// a.u.
#define geo 0.00001										// Геометрическая точность
#define x_min_  -1000.0 // -2500.0 // -1300  //-2000                // -1500.0
#define x_max_ 650.0 // 450.0
#define y_max_ 1400.0 // 1600.0 //1840.0
#define chi_ 36.1059
#define n_p_LISM_ (3.0) 
#define n_H_LISM_ (1.0)

#define kurant  0.95


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