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
	double DX;            // размеры ячеек при создании сетки (ячеек с размером 1)
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
	vector <Kyb*> Sosed_1;  // Эти вектора нужны на этапе дробления сетки - они просто ускоряют работу
	vector <Kyb*> Sosed_2;
	vector <Kyb*> Sosed_3;
	vector <Kyb*> Sosed_4;
	vector <Kyb*> Boandary_1;  // Ячейки - граничащие с правой границей
	vector <Kyb*> Boandary_2;  // Ячейки - граничащие с левой границей
	vector <Kyb*> Boandary_3;  // Ячейки - граничащие с верхней границей
	vector<Sensor*> Sensors;


	Konstruktor(int a, int b, double xmin, double xmax, double ymax);      /// Конструктор класса

	virtual ~Konstruktor();

	bool sosed_or_not(Kyb* A, Kyb* B);


	// Блок вывода геометрической информации по сетки

	void print_konectiviti_short(void);
	/// Показывает между какими ячейками есть связи (для ТекПлота)
	/// не включает в себя граничные ячейки G1 - G4
	void print_point(void);
	void print_cell(void);



	void droblenie2_hand(Kyb* A);
	void Drobim(double x0, double y0, double r1, double r2);
	void Drobim(double x1, double x2, double r1);
	int get_size_conektiv(void);



	// Блок сохранения и загрузки сетки

	void Save_setka(string name);
	void Save_setka_multifluid(string name);
	void Download_setka(string name);
	void Download_setka_multifluid(string name);


	// Газовая данамика
	void initial_condition();
	void read_Cuda_massiv(double* ro, double* p, double* u, double* v, double* QQ);
	void read_Cuda_massiv(double* ro, double* p, double* u, double* v, double* QQ,//
							double* ro_H1, double* p_H1, double* u_H1, double* v_H1,//
							double* ro_H2, double* p_H2, double* u_H2, double* v_H2, //
							double* ro_H3, double* p_H3, double* u_H3, double* v_H3, //
							double* ro_H4, double* p_H4, double* u_H4, double* v_H4);
	void print_Tecplot(void);
	void print_Tecplot_multifluid(void);


	// Монте-Карло
	void M_K_training(void); // Подготовка к Монте-Карло
	// Функция должны быть обязательно запущена, если вы собираетесь использовать Монте-Карло в программе.
	// Запуск только один раз, но после дробления!!!


	void Velosity_initial(Sensor* s, double& Vx, double& Vy, double& Vz);
	void Velosity_initial2(Sensor* s, double& Vx, double& Vy, double& Vz);
	void Change_Velosity(Sensor* s, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, double& X, double& Y, double& Z, const double& cp);
	void M_K(void);  // основная функция запуска чатиц
	void Fly_exchenge(Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Kyb* ind, const double& mu);
	bool Flying_exchange(double& KSI, double& Vx, double& Vy, double& Vz, double& X, double& Y,//
		double& Z, Kyb*& next, Kyb* head, Kyb* prev, const double& mu, double& I_do);
	bool Peresechenie(const double& x0, const double& y0, const double& dx, const double& dy, const double& x, const double& y, const double& z, //
		const double& Vx, const double& Vy, const double& Vz, int& mode, double& t);



	Kyb* Belong_point(int b, const double& x, const double& y);  // Ищет какой граничной ячейке принадлежит точка
	// b - это номер границы   Boandary_b


	// Блок проекции решения на ось x

	/// Функция, проектирующая решение на ось x
	/// <param name=""><ничего не принимает>
	/// <returns><возвращает вектор новых ячеек, центры которых находятся на оси>
	vector <Kyb*> Get_projection(void);
	double linear_funk(const double& x1, const double& y1, const double& x2, const double& y2, const double& t);



private:
	void initialization(int a, int b, double xmin, double xmax, double ymax);
	void New_design(void);
	void konectiviti(void);
	void number(void);

	void peresich(const double& y, const double& z, const double& Vy, const double& Vz, const double& R, double& t1, double& t2);
	// Находит пересечение со сферой, вспомогательная функция

	double minplus(const double& x, const double& y);
	double polar_angle(const double& x, const double& y);

	double Velosity_1(const double& u, const double& cp);
	double Velosity_2(const double& u, const double& cp);
	double Velosity_3(const double& u, const double& cp);
};

