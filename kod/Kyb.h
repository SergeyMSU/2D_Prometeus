#pragma once
#include <vector>
using namespace std;
class Kyb;

class Kyb
{
public:
	double x;                // Координата центра ячейки
	double y;                // Координата центра ячейки
	double ro;
	double p;
	double u;
	double v;
	double Q;
	vector <Kyb*> sosed;
	int number;
	int size;                // Показывает на каком этапе мельчения находится эта ячейка, это позволит определить точно её размер (без потери точности)

	static int move_number;
	bool drob;               // Удобная булевская переменная для использования в разных функциях 


	Kyb(double x, double y);      /// Конструктор класса

private:
	void initialization(double, double, int);
};

