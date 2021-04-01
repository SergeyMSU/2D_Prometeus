#include "Kyb.h"
using namespace std;
int Kyb::move_number = 0;

Kyb::Kyb(double x, double y)
{
	this->initialization(x, y, this->move_number);
	this->move_number++;
}

void Kyb::initialization(double x, double y, int nn)
{
	this->x = x;
	this->y = y;
	this->number = nn;
	this->drob = false;
	this->size = 1;             // Исходный размер по-умолчанию
	this->ro = 0.0;
	this->p = 0.0;
	this->u = 0.0;
	this->v = 0.0;
	this->Q = 0.0;
	this->F_n = 0.0;
	this->F_u = 0.0;
	this->F_v = 0.0;
	this->F_T = 0.0;
	this->I_u = 0.0;
	this->I_v = 0.0;
	this->I_T = 0.0;
}


void Copy(Kyb * A, Kyb * B)
{
	A->ro = B->ro;
	A->p = B->p;
	A->u = B->u;
	A->v = B->v;
	A->Q = B->Q;

	A->ro_H1 = B->ro_H1;
	A->p_H1 = B->p_H1;
	A->u_H1 = B->u_H1;
	A->v_H1 = B->v_H1;

	A->ro_H2 = B->ro_H2;
	A->p_H2 = B->p_H2;
	A->u_H2 = B->u_H2;
	A->v_H2 = B->v_H2;

	A->ro_H3 = B->ro_H3;
	A->p_H3 = B->p_H3;
	A->u_H3 = B->u_H3;
	A->v_H3 = B->v_H3;

	A->ro_H4 = B->ro_H4;
	A->p_H4 = B->p_H4;
	A->u_H4 = B->u_H4;
	A->v_H4 = B->v_H4;

	return;
}

bool Kyb::Belong_slow(const double& xx, const double& yy, const double& DX, const double& DY)
{
	double dx = (DX / pow(2, this->size - 1)) / 2.0;   // Половина длины ячейки
	double dy = (DY / pow(2, this->size - 1)) / 2.0;   // Половина ширины ячейки

	if (xx >= this->x - dx && xx <= this->x + dx)
	{
		if (yy >= this->y - dy && yy <= this->y + dy)
		{
			return true;
		}
	}

	return false;
}

bool Kyb::Belong_fast(const double& xx, const double& yy, const double& dx, const double& dy)
{

	if (xx >= this->x - dx && xx <= this->x + dx)
	{
		if (yy >= this->y - dy && yy <= this->y + dy)
		{
			return true;
		}
	}

	return false;
}

void Kyb::Setup_boandary(const double& DX, const double& DY)
{
	double dx1 = (DX / pow(2, this->size - 1)) / 2.0;   // Половина длины ячейки
	double dy1 = (DY / pow(2, this->size - 1)) / 2.0;   // Половина ширины ячейки
	double dx2, dy2;
	double n1, n2;

	for (auto& i : this->sosed)
	{
		dx2 = (DX / pow(2, i->size - 1)) / 2.0;   // Половина длины ячейки
		dy2 = (DY / pow(2, i->size - 1)) / 2.0;   // Половина ширины ячейки
		if (fabs(fabs(this->x - i->x) - dx1 - dx2) < geo)
		{
			n1 = (i->x - this->x);
			if (n1 > 0)
			{
				this->boandary_1.push_back(i);
			}
			else
			{
				this->boandary_2.push_back(i);
			}
		}
		else if (fabs(fabs(this->y - i->y) - dy1 - dy2) < geo)
		{
			n2 = (i->y - this->y);
			if (n2 > 0)
			{
				this->boandary_3.push_back(i);
			}
			else
			{
				this->boandary_4.push_back(i);
			}
		}
	}
}

