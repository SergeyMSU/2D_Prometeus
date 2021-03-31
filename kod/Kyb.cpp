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
	this->size = 1;             // Èñõîäíûé ğàçìåğ ïî-óìîë÷àíèş
	this->ro = 0.0;
	this->p = 0.0;
	this->u = 0.0;
	this->v = 0.0;
	this->Q = 0.0;
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

