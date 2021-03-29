#pragma once
#include <vector>
#include "Kyb.h"
#include "Cell.h"
using namespace std;
struct Cell;

class Konstruktor
{
public:
	vector <Kyb*> all_Kyb;
	int N;
	int M;
	double DX;            // размеры €чеек при создании сетки (€чеек с размером 1)
	double DY; 
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	Kyb* G1;
	Kyb* G2;
	Kyb* G3;
	Kyb* G4;
	vector <Kyb*> Sosed_1;  // Ёти вектора нужны на этапе дроблени€ сетки - они просто ускор€ют работу
	vector <Kyb*> Sosed_2;
	vector <Kyb*> Sosed_3;
	vector <Kyb*> Sosed_4;
	vector <Kyb*> Sosed_5;
	vector <Kyb*> Sosed_6;


	Konstruktor(int, int, double, double, double);      ///  онструктор класса

	virtual ~Konstruktor();


	void print_konectiviti_short(void);
	/// ѕоказывает между какими €чейками есть св€зи (дл€ “екѕлота)
	/// не включает в себ€ граничные €чейки G1 - G4



private:
	void initialization(int a, int b, double x1, double x2, double y2);
	void New_design(void);
	void konectiviti(void);
};

