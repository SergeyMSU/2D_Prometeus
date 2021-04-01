#pragma once
#include <vector>
#include "Kyb.h"
#include "Cell.h"
#include "sensor.h"
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
	double sqv_1;
	double sqv_2;
	double sqv_3;
	double sqv_4;
	double sqv_5;
	double sum_s;
	int Number1;
	int Number2;
	int Number3;
	int Number4;
	int AllNumber;
	Kyb* G1;
	Kyb* G2;
	Kyb* G3;
	Kyb* G4;
	vector <Kyb*> Sosed_1;  // Ёти вектора нужны на этапе дроблени€ сетки - они просто ускор€ют работу
	vector <Kyb*> Sosed_2;
	vector <Kyb*> Sosed_3;
	vector <Kyb*> Sosed_4;
	vector <Kyb*> Boandary_1;  // ячейки - граничащие с правой границей
	vector <Kyb*> Boandary_2;  // ячейки - граничащие с левой границей
	vector <Kyb*> Boandary_3;  // ячейки - граничащие с верхней границей
	vector<Sensor*> Sensors;


	Konstruktor(int a, int b, double xmin, double xmax, double ymax);      ///  онструктор класса

	virtual ~Konstruktor();

	bool sosed_or_not(Kyb* A, Kyb* B);


	// Ѕлок вывода геометрической информации по сетки

	void print_konectiviti_short(void);
	/// ѕоказывает между какими €чейками есть св€зи (дл€ “екѕлота)
	/// не включает в себ€ граничные €чейки G1 - G4
	void print_point(void);
	void print_cell(void);



	void droblenie2_hand(Kyb* A);
	void Drobim(double x0, double y0, double r1, double r2);
	void Drobim(double x1, double x2, double r1);
	int get_size_conektiv(void);



	// Ѕлок сохранени€ и загрузки сетки

	void Save_setka(string name);
	void Save_setka_multifluid(string name);
	void Download_setka(string name);
	void Download_setka_multifluid(string name);


	// √азова€ данамика
	void initial_condition();
	void read_Cuda_massiv(double* ro, double* p, double* u, double* v, double* QQ);
	void read_Cuda_massiv(double* ro, double* p, double* u, double* v, double* QQ,//
							double* ro_H1, double* p_H1, double* u_H1, double* v_H1,//
							double* ro_H2, double* p_H2, double* u_H2, double* v_H2, //
							double* ro_H3, double* p_H3, double* u_H3, double* v_H3, //
							double* ro_H4, double* p_H4, double* u_H4, double* v_H4);
	void print_Tecplot(void);
	void print_Tecplot_multifluid(void);


	// ћонте- арло
	void Velosity_initial(Sensor* s, double& Vx, double& Vy, double& Vz);
	void Velosity_initial2(Sensor* s, double& Vx, double& Vy, double& Vz);
	void Change_Velosity(Sensor* s, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, double& X, double& Y, double& Z, const double& cp);
	void M_K(void);  // основна€ функци€ запуска чатиц
	void Fly_exchenge(Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Kyb* ind, const double& mu);
	bool Flying_exchange(double& KSI, double& Vx, double& Vy, double& Vz, double& X, double& Y,//
		double& Z, Kyb*& next, Kyb* head, Kyb* prev, const double& mu, double& I_do);
	bool Peresechenie(const double& x0, const double& y0, const double& dx, const double& dy, const double& x, const double& y, const double& z, //
		const double& Vx, const double& Vy, const double& Vz, int& mode, double& t);



	Kyb* Konstruktor::Belong_point(int b, const double& x, const double& y);  // »щет какой граничной €чейке принадлежит точка
	// b - это номер границы   Boandary_b


	// Ѕлок проекции решени€ на ось x

	/// ‘ункци€, проектирующа€ решение на ось x
	/// <param name=""><ничего не принимает>
	/// <returns><возвращает вектор новых €чеек, центры которых наход€тс€ на оси>
	vector <Kyb*> Get_projection(void);
	double linear_funk(const double& x1, const double& y1, const double& x2, const double& y2, const double& t);



private:
	void initialization(int a, int b, double xmin, double xmax, double ymax);
	void New_design(void);
	void konectiviti(void);
	void number(void);

	void peresich(const double& y, const double& z, const double& Vy, const double& Vz, const double& R, double& t1, double& t2);
	// Ќаходит пересечение со сферой, вспомогательна€ функци€

	double minplus(const double& x, const double& y);
	double polar_angle(const double& x, const double& y);
};

