#include "Kyb.h"
using namespace std;
int Kyb::move_number = 0;

Kyb::Kyb(double x, double y)
{
	this->initialization(x, y, this->move_number);
	this->move_number++;
}

void Kyb::initialization(double x , double y, int nn)
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

}

