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
	double DX;            // ������� ����� ��� �������� ����� (����� � �������� 1)
	double DY; 
	double x_min;
	double x_max;
	double y_max;
	Kyb* G1;
	Kyb* G2;
	Kyb* G3;
	Kyb* G4;
	vector <Kyb*> Sosed_1;  // ��� ������� ����� �� ����� ��������� ����� - ��� ������ �������� ������
	vector <Kyb*> Sosed_2;
	vector <Kyb*> Sosed_3;
	vector <Kyb*> Sosed_4;


	Konstruktor(int a, int b, double xmin, double xmax, double ymax);      /// ����������� ������

	virtual ~Konstruktor();

	bool sosed_or_not(Kyb* A, Kyb* B);


	// ���� ������ �������������� ���������� �� �����

	void print_konectiviti_short(void);
	/// ���������� ����� ������ �������� ���� ����� (��� ��������)
	/// �� �������� � ���� ��������� ������ G1 - G4
	void print_point(void);
	void print_cell(void);



	void droblenie2_hand(Kyb* A);
	void Drobim(double x0, double y0, double r1, double r2);
	void Drobim(double x1, double x2, double r1);
	int get_size_conektiv(void);



	// ���� ���������� � �������� �����

	void Save_setka(string name);
	void Save_setka_multifluid(string name);
	void Download_setka(string name);
	void Download_setka_multifluid(string name);


	// ��������
	void initial_condition();
	void read_Cuda_massiv(double* ro, double* p, double* u, double* v, double* QQ);
	void read_Cuda_massiv(double* ro, double* p, double* u, double* v, double* QQ,//
							double* ro_H1, double* p_H1, double* u_H1, double* v_H1,//
							double* ro_H2, double* p_H2, double* u_H2, double* v_H2, //
							double* ro_H3, double* p_H3, double* u_H3, double* v_H3, //
							double* ro_H4, double* p_H4, double* u_H4, double* v_H4);
	void print_Tecplot(void);
	void print_Tecplot_multifluid(void);


	// ���� �������� ������� �� ��� x

	/// �������, ������������� ������� �� ��� x
	/// <param name=""><������ �� ���������>
	/// <returns><���������� ������ ����� �����, ������ ������� ��������� �� ���>
	vector <Kyb*> Get_projection(void);
	double linear_funk(const double& x1, const double& y1, const double& x2, const double& y2, const double& t);



private:
	void initialization(int a, int b, double xmin, double xmax, double ymax);
	void New_design(void);
	void konectiviti(void);
	void number(void);
};

