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
	double y_max;
	Kyb* G1;
	Kyb* G2;
	Kyb* G3;
	Kyb* G4;
	vector <Kyb*> Sosed_1;  // Ёти вектора нужны на этапе дроблени€ сетки - они просто ускор€ют работу
	vector <Kyb*> Sosed_2;
	vector <Kyb*> Sosed_3;
	vector <Kyb*> Sosed_4;


	Konstruktor(int a, int b, double xmin, double xmax, double ymax);      ///  онструктор класса

	virtual ~Konstruktor();

	bool sosed_or_not(Kyb* A, Kyb* B);


	void print_konectiviti_short(void);
	/// ѕоказывает между какими €чейками есть св€зи (дл€ “екѕлота)
	/// не включает в себ€ граничные €чейки G1 - G4
	void print_point(void);
	void print_cell(void);


	void droblenie2_hand(Kyb* A);
	void Drobim(double x0, double y0, double r1, double r2);
	void Drobim(double x1, double x2, double r1);

	int get_size_conektiv(void);

	void Save_setka(string name);


	// ћеханика
	void initial_condition();
	void read_Cuda_massiv(double* ro, double* p, double* u, double* v, double* QQ);
	void print_Tecplot(void);



private:
	void initialization(int a, int b, double xmin, double xmax, double ymax);
	void New_design(void);
	void konectiviti(void);
	void number(void);
};

