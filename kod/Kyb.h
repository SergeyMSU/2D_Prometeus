#pragma once
#include <vector>
#include "sensor.h"
#include <mutex>
#include "Cell.h"
using namespace std;
class Kyb;

class Kyb
{
public:
	double x;                // Координата центра ячейки
	double y;                // Координата центра ячейки
	double dx;               
	double dy;               
	double ro;
	double p;
	double u;
	double v;
	double Q;
	double ro_H1;
	double p_H1;
	double u_H1;
	double v_H1;
	double ro_H2;
	double p_H2;
	double u_H2;
	double v_H2;
	double ro_H3;
	double p_H3;
	double u_H3;
	double v_H3;
	double ro_H4;
	double p_H4;
	double u_H4;
	double v_H4;
	double F_n;
	double F_u;
	double F_v;
	double F_T;
	double I_u;
	double I_v;
	double I_T;
	vector <Kyb*> sosed;
	vector <Kyb*> boandary_1; // Соседи справа
	vector <Kyb*> boandary_2; // Соседи справа
	vector <Kyb*> boandary_3; // Соседи справа
	vector <Kyb*> boandary_4; // Соседи справа
	int number;
	int size;                // Показывает на каком этапе мельчения находится эта ячейка, это позволит определить точно её размер (без потери точности)
	mutex mut;

	static int move_number;
	bool drob;               // Удобная булевская переменная для использования в разных функциях 


	Kyb(double x, double y);      /// Конструктор класса

	/// Копирует значения переменных из B в A
	friend void Copy(Kyb* A, Kyb* B);

	// Принадлежит ли точка ячейке?
	bool Belong_slow(const double& xx, const double& yy, const double& DX, const double& DY);
	// Функция сама находит размер ячейки по переданным глобальным размерам, поэтому работает достаточно медленно

	bool Belong_fast(const double& xx, const double& yy, const double& dx, const double& dy);
	// Быстрый аналог предыдущей функции, но нужно подавать размеры ячейки в неё

	void Setup_boandary(const double& DX, const double& DY);
	// Расфвсовать соседей по границам - с какой они стараны

private:
	void initialization(double, double, int);
};

