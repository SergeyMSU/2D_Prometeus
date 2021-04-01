#pragma once
#include <vector>
#include "sensor.h"
using namespace std;
class Kyb;

class Kyb
{
public:
	double x;                // ���������� ������ ������
	double y;                // ���������� ������ ������
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
	vector <Kyb*> sosed;
	int number;
	int size;                // ���������� �� ����� ����� ��������� ��������� ��� ������, ��� �������� ���������� ����� � ������ (��� ������ ��������)

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

private:
	void initialization(double, double, int);
};

