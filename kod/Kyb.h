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
	double x;                // ���������� ������ ������
	double y;                // ���������� ������ ������
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
	vector <Kyb*> boandary_1; // ������ ������
	vector <Kyb*> boandary_2; // ������ ������
	vector <Kyb*> boandary_3; // ������ ������
	vector <Kyb*> boandary_4; // ������ ������
	int number;
	int size;                // ���������� �� ����� ����� ��������� ��������� ��� ������, ��� �������� ���������� ����� � ������ (��� ������ ��������)
	mutex mut;

	static int move_number;
	bool drob;               // ������� ��������� ���������� ��� ������������� � ������ �������� 


	Kyb(double x, double y);      /// ����������� ������

	/// �������� �������� ���������� �� B � A
	friend void Copy(Kyb* A, Kyb* B);

	// ����������� �� ����� ������?
	bool Belong_slow(const double& xx, const double& yy, const double& DX, const double& DY);
	// ������� ���� ������� ������ ������ �� ���������� ���������� ��������, ������� �������� ���������� ��������

	bool Belong_fast(const double& xx, const double& yy, const double& dx, const double& dy);
	// ������� ������ ���������� �������, �� ����� �������� ������� ������ � ��

	void Setup_boandary(const double& DX, const double& DY);
	// ����������� ������� �� �������� - � ����� ��� �������

private:
	void initialization(double, double, int);
};

