#include "Konstruktor.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>


using namespace std;

Konstruktor::Konstruktor(int a, int b,  double x1, double x2, double y2)
{
	this->initialization(a, b, x1, x2, y2);
	this->New_design();                      // —оздаЄт начальную сетку (без св€зей)
	this->konectiviti();                     // —в€зывает начальные €чейки с их сосед€ми
}


Konstruktor::~Konstruktor()
{
	for (auto& i : this->all_Kyb)
	{
		delete i;
	}
	this->all_Kyb.clear();
}

void Konstruktor::initialization(int a, int b, double x1, double x2, double y2)
{
	this->N = a;
	this->M = b;
	this->DX = 0.0;
	this->DY = 0.0;
	this->x_min = x1;
	this->x_max = x2;
	this->y_max = y2;
	this->y_min = 0.0;
	this->Sosed_1.reserve(128);
	this->Sosed_2.reserve(128);
	this->Sosed_3.reserve(128);
	this->Sosed_4.reserve(128);
	this->Sosed_5.reserve(128);
	this->Sosed_6.reserve(128);
	this->all_Kyb.reserve(1000000);
}

void Konstruktor::New_design(void)
{
	double dx = (this->x_max - this->x_min) / (this->N - 1);
	double dy = (this->y_max) / this->M;
	this->DX = dx;
	this->DY = dy;

	//cout << "z = " << z_min + dz / 2.0 + dz * k << endl;
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < N; i++)
		{
			double y = y_min + j * dy;
			double x = x_min + i * dx;
			auto C = new Kyb(x, y);
			C->dx = dx;
			C->dy = dy;
			this->all_Kyb.push_back(C);
		}
	}
	
}

void Konstruktor::konectiviti(void)
{

	auto C1 = new Kyb(x_max + this->DX, 0.0);
	C1->number = -1;
	C1->size = 0;
	this->G1 = C1;

	auto C2 = new Kyb(x_min - this->DX, 0.0);
	C2->number = -2;
	C2->size = 0;
	this->G2 = C2;

	auto C3 = new Kyb(0.0, y_max + this->DY);
	C3->number = -3;
	C3->size = 0;
	this->G3 = C3;

	auto C4 = new Kyb(0.0, y_min - this->DY);
	C4->number = -4;
	C4->size = 0;
	this->G4 = C4;


	for (auto& i : this->all_Kyb)
	{
		int n = i->number % this->N;
		int m = (i->number - n) / this->N;
		// n,m,k - глобальный номер €чейки
		//cout << n << " " << m << " " << k << endl;
		if (i->number >= 0)
		{
			if (n > 0)
			{
				i->sosed.push_back(this->all_Kyb[m * this->N + (n - 1)]);
			}
			else
			{
				i->sosed.push_back(C2);
			}
			if (n < this->N - 1)
			{
				i->sosed.push_back(this->all_Kyb[m * this->N + (n + 1)]);
			}
			else
			{
				i->sosed.push_back(C1);
			}
			if (m > 0)
			{
				i->sosed.push_back(this->all_Kyb[(m - 1) * this->N + n]);
			}
			else
			{
				i->sosed.push_back(C4);
			}
			if (m < M - 1)
			{
				i->sosed.push_back(this->all_Kyb[(m + 1) * this->N + n]);
			}
			else
			{
				i->sosed.push_back(C3);
			}
		}
	}
}

void Konstruktor::print_konectiviti_short(void)
{
	int ll = 0;
	for (auto& i : this->all_Kyb)
	{
		for (auto& j : i->sosed)
		{
			if (j->number >= 0)
			{
				ll++;
			}
		}
	}
	ofstream fout;
	fout.open("Setka_konectiviti_short.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << 2 * ll;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->all_Kyb)
	{
		for (auto& j : i->sosed)
		{
			if (j->number >= 0)
			{
				fout << i->x << " " << i->y << endl;
				fout << (i->x + j->x) / 2.0 << " " << (i->y + j->y) / 2.0 << endl;
			}
		}
	}
	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}
}

