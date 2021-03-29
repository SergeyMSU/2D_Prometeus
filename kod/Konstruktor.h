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


	void print_konectiviti_short(void);
	/// ���������� ����� ������ �������� ���� ����� (��� ��������)
	/// �� �������� � ���� ��������� ������ G1 - G4
	void print_point(void);
	void print_cell(void);


	void droblenie2_hand(Kyb* A);
	void Drobim(double x0, double y0, double r1, double r2);



private:
	void initialization(int a, int b, double xmin, double xmax, double ymax);
	void New_design(void);
	void konectiviti(void);
	void number(void);
};

