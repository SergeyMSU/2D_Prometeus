#pragma once
#include <vector>
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
	vector <Kyb*> sosed;
	int number;
	int size;                // ���������� �� ����� ����� ��������� ��������� ��� ������, ��� �������� ���������� ����� � ������ (��� ������ ��������)

	static int move_number;
	bool drob;               // ������� ��������� ���������� ��� ������������� � ������ �������� 


	Kyb(double x, double y);      /// ����������� ������

private:
	void initialization(double, double, int);
};

