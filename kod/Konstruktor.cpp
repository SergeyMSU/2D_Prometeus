#include "Konstruktor.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>

#define geo  0.001

using namespace std;

Konstruktor::Konstruktor(int a, int b, int c, double x1, double x2, double y1, double y2, double z1, double z2)
{
	this->initialization(a, b, c, x1, x2, y1, y2, z1, z2);
	this->New_design();
	this->konectiviti();
}

Konstruktor::Konstruktor(string name)
{
	this->download_Setka(name);
	this->Sosed_1.reserve(128);
	this->Sosed_2.reserve(128);
	this->Sosed_3.reserve(128);
	this->Sosed_4.reserve(128);
	this->Sosed_5.reserve(128);
	this->Sosed_6.reserve(128);
}

Konstruktor::~Konstruktor()
{
	for (auto& i : this->all_Kyb)
	{
		delete i;
	}
	this->all_Kyb.clear();
}

void Konstruktor::initialization(int a, int b, int c, double x1, double x2, double y1, double y2, double z1, double z2)
{
	this->N = a;
	this->M = b;
	this->K = c;
	this->x_min = x1;
	this->x_max = x2;
	this->y_min = y1;
	this->y_max = y2;
	this->z_min = z1;
	this->z_max = z2;
	this->Sosed_1.reserve(128);
	this->Sosed_2.reserve(128);
	this->Sosed_3.reserve(128);
	this->Sosed_4.reserve(128);
	this->Sosed_5.reserve(128);
	this->Sosed_6.reserve(128);
	this->all_Kyb.reserve(50000000);
}

void Konstruktor::New_design(void)
{
	double dx = (this->x_max - this->x_min) / this->N;
	double dy = (this->y_max - this->y_min) / this->M;
	double dz = (this->z_max - this->z_min) / this->K;
	//cout << "dz = " << dz << endl;

	for (int k = 0; k < K; k++)
	{
		//cout << "z = " << z_min + dz / 2.0 + dz * k << endl;
		for (int j = 0; j < M; j++)
		{
			for (int i = 0; i < N; i++)
			{
				double x = x_min + dx / 2.0 + dx * i;
				double y = y_min + dy / 2.0 + dy * j;
				double z = z_min + dz / 2.0 + dz * k;
				auto C = new Kyb(x, y, z);
				C->dx = dx / 2.0;
				C->dy = dy / 2.0;
				C->dz = dz / 2.0;
				this->all_Kyb.push_back(C);
			}
		}
	}
}

void Konstruktor::konectiviti(void)
{
	double dx = (this->x_max - this->x_min) / this->N;
	double dy = (this->y_max - this->y_min) / this->M;
	double dz = (this->z_max - this->z_min) / this->K;

	auto C1 = new Kyb(x_max + dx/2.0, 0.0, 0.0);
	C1->dx = dx / 2.0;
	C1->dy = y_max - y_min;
	C1->dz = z_max - z_min;
	C1->number = -1;
	this->G1 = C1;

	auto C2 = new Kyb(x_min - dx / 2.0, 0.0, 0.0);
	C2->dx = dx / 2.0;
	C2->dy = y_max - y_min;
	C2->dz = z_max - z_min;
	C2->number = -2;
	this->G2 = C2;

	auto C3 = new Kyb(0.0, y_max + dy / 2.0, 0.0);
	C3->dx = x_max - x_min;
	C3->dy = dy / 2.0;
	C3->dz = z_max - z_min;
	C3->number = -3;
	this->G3 = C3;

	auto C4 = new Kyb(0.0, y_min - dy / 2.0, 0.0);
	C4->dx = x_max - x_min;
	C4->dy = dy / 2.0;
	C4->dz = z_max - z_min;
	C4->number = -4;
	this->G4 = C4;

	auto C5 = new Kyb(0.0, 0.0, z_max +  dz / 2.0);
	C5->dx = x_max - x_min;
	C5->dy = y_max - y_min;
	C5->dz = dz / 2.0;
	C5->number = -5;
	this->G5 = C5;

	auto C6 = new Kyb(0.0, 0.0, z_min - dz / 2.0);
	C6->dx = x_max - x_min;
	C6->dy = y_max - y_min;
	C6->dz = dz / 2.0;
	C6->number = -6;
	C6->move_number -= 6;
	this->G6 = C6;


	for (auto& i : this->all_Kyb)
	{
		int n = i->number % N;
		int L = this->N * this->M;
		int m = (i->number % L - n) / N;
		int k = (int)(i->number / L);
		// n,m,k - глобальный номер ячейки
		//cout << n << " " << m << " " << k << endl;
		if (i->number >= 0)
		{
			if (n > 0)
			{
				i->sosed.push_back(this->all_Kyb[k * L + m * N + (n - 1)]);
			}
			else
			{
				i->sosed.push_back(C2);
			}
			if (n < N - 1)
			{
				i->sosed.push_back(this->all_Kyb[k * L + m * N + (n + 1)]);
			}
			else
			{
				i->sosed.push_back(C1);
			}
			if (m > 0)
			{
				i->sosed.push_back(this->all_Kyb[k * L + (m - 1) * N + n]);
			}
			else
			{
				i->sosed.push_back(C4);
			}
			if (m < M - 1)
			{
				i->sosed.push_back(this->all_Kyb[k * L + (m + 1) * N + n]);
			}
			else
			{
				i->sosed.push_back(C3);
			}
			if (k > 0)
			{
				i->sosed.push_back(this->all_Kyb[(k - 1) * L + m * N + n]);
			}
			else
			{
				i->sosed.push_back(C6);
			}
			if (k < K - 1)
			{
				i->sosed.push_back(this->all_Kyb[(k + 1) * L + m * N + n]);
			}
			else
			{
				i->sosed.push_back(C5);
			}
		}
	}
}

void Konstruktor::print_my(int num)
{
	ofstream fout;
	fout.open("Yacheika.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\", \"Z\", \"B\"  ZONE T= \"HP\"" << endl;
	fout << this->all_Kyb[num]->x << " " << this->all_Kyb[num]->y << " " << this->all_Kyb[num]->z << " " << 0 << endl;
	for (auto& i : this->all_Kyb[num]->sosed)
	{
		fout << i->x << " " << i->y << " " << i->z << " " << 1 << endl;
		fout << i->x - i->dx << " " << i->y - i->dy << " " << i->z - i->dz << " " << 2 << endl;
		fout << i->x + i->dx << " " << i->y - i->dy << " " << i->z - i->dz << " " << 2 << endl;
		fout << i->x - i->dx << " " << i->y + i->dy << " " << i->z - i->dz << " " << 2 << endl;
		fout << i->x - i->dx << " " << i->y - i->dy << " " << i->z + i->dz << " " << 2 << endl;
		fout << i->x - i->dx << " " << i->y + i->dy << " " << i->z + i->dz << " " << 2 << endl;
		fout << i->x + i->dx << " " << i->y - i->dy << " " << i->z + i->dz << " " << 2 << endl;
		fout << i->x + i->dx << " " << i->y + i->dy << " " << i->z - i->dz << " " << 2 << endl;
		fout << i->x + i->dx << " " << i->y + i->dy << " " << i->z + i->dz << " " << 2 << endl;
	}
}

void Konstruktor::print_my(int num, int num2)
{
	ofstream fout;
	fout.open("Yacheika.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\", \"Z\", \"B\"  ZONE T= \"HP\"" << endl;
	fout << this->all_Kyb[num]->x << " " << this->all_Kyb[num]->y << " " << this->all_Kyb[num]->z << " " << 0 << endl;
	fout << this->all_Kyb[num2]->x << " " << this->all_Kyb[num2]->y << " " << this->all_Kyb[num2]->z << " " << 0 << endl;
	for (auto& i : this->all_Kyb[num]->sosed)
	{
		fout << i->x << " " << i->y << " " << i->z << " " << 1 << endl;
		fout << i->x - i->dx << " " << i->y - i->dy << " " << i->z - i->dz << " " << 2 << endl;
		fout << i->x + i->dx << " " << i->y - i->dy << " " << i->z - i->dz << " " << 2 << endl;
		fout << i->x - i->dx << " " << i->y + i->dy << " " << i->z - i->dz << " " << 2 << endl;
		fout << i->x - i->dx << " " << i->y - i->dy << " " << i->z + i->dz << " " << 2 << endl;
		fout << i->x - i->dx << " " << i->y + i->dy << " " << i->z + i->dz << " " << 2 << endl;
		fout << i->x + i->dx << " " << i->y - i->dy << " " << i->z + i->dz << " " << 2 << endl;
		fout << i->x + i->dx << " " << i->y + i->dy << " " << i->z - i->dz << " " << 2 << endl;
		fout << i->x + i->dx << " " << i->y + i->dy << " " << i->z + i->dz << " " << 2 << endl;
	}

	for (auto& i : this->all_Kyb[num2]->sosed)
	{
		fout << i->x << " " << i->y << " " << i->z << " " << 1 << endl;
		fout << i->x - i->dx << " " << i->y - i->dy << " " << i->z - i->dz << " " << 2 << endl;
		fout << i->x + i->dx << " " << i->y - i->dy << " " << i->z - i->dz << " " << 2 << endl;
		fout << i->x - i->dx << " " << i->y + i->dy << " " << i->z - i->dz << " " << 2 << endl;
		fout << i->x - i->dx << " " << i->y - i->dy << " " << i->z + i->dz << " " << 2 << endl;
		fout << i->x - i->dx << " " << i->y + i->dy << " " << i->z + i->dz << " " << 2 << endl;
		fout << i->x + i->dx << " " << i->y - i->dy << " " << i->z + i->dz << " " << 2 << endl;
		fout << i->x + i->dx << " " << i->y + i->dy << " " << i->z - i->dz << " " << 2 << endl;
		fout << i->x + i->dx << " " << i->y + i->dy << " " << i->z + i->dz << " " << 2 << endl;
	}
}

void Konstruktor::print_konectiviti(void)
{
	int ll = 0;
	for (auto& i : this->all_Kyb)
	{
		ll = ll + i->sosed.size();
	}
	ofstream fout;
	fout.open("Setka_konectiviti.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\", \"Z\"  ZONE T= \"HP\", N=  " << this->all_Kyb.size() + 6;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->all_Kyb)
	{
		fout << i->x << " " << i->y << " " << i->z << endl;
	}
	fout << this->G1->x << " " << this->G1->y << " " << this->G1->z << endl;
	fout << this->G2->x << " " << this->G2->y << " " << this->G2->z << endl;
	fout << this->G3->x << " " << this->G3->y << " " << this->G3->z << endl;
	fout << this->G4->x << " " << this->G4->y << " " << this->G4->z << endl;
	fout << this->G5->x << " " << this->G5->y << " " << this->G5->z << endl;
	fout << this->G6->x << " " << this->G6->y << " " << this->G6->z << endl;

	for (auto& i : this->all_Kyb)
	{
		if (i->number >= 0)
		{
			for (auto& j : i->sosed)
			{
				if (j->number >= 0)
				{
					fout << i->number + 1 << " " << j->number + 1 << endl;
				}
				else if (j->number == -6)
				{
					fout << i->number + 1 << " " << this->all_Kyb.size() + 6 << endl;
				}
				else if (j->number == -5)
				{
					fout << i->number + 1 << " " << this->all_Kyb.size() + 5 << endl;
				}
				else if (j->number == -4)
				{
					fout << i->number + 1 << " " << this->all_Kyb.size() + 4<< endl;
				}
				else if (j->number == -3)
				{
					fout << i->number + 1 << " " << this->all_Kyb.size() + 3 << endl;
				}
				else if (j->number == -2)
				{
					fout << i->number + 1 << " " << this->all_Kyb.size() + 2 << endl;
				}
				else if (j->number == -1)
				{
					fout << i->number + 1 << " " << this->all_Kyb.size() + 1 << endl;
				}
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
	fout << " VARIABLES = \"X\", \"Y\", \"Z\"  ZONE T= \"HP\", N=  " << 2 * ll;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->all_Kyb)
	{
		for (auto& j : i->sosed)
		{
			if (j->number >= 0)
			{
				fout << i->x << " " << i->y << " " << i->z << endl;
				fout << (i->x + j->x) / 2.0 << " " << (i->y + j->y) / 2.0 << " " << (i->z + j->z) / 2.0 << endl;
			}
		}
	}
	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}
}

void Konstruktor::print_konectiviti2(void)
{
	int ll = 0;
	for (auto& i : this->all_Kyb)
	{
		ll = ll + i->sosed.size();
	}
	ofstream fout;
	fout.open("Setka_konectiviti_2.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\", \"Z\"  ZONE T= \"HP\", N=  " << 2 * ll;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->all_Kyb)
	{
		for (auto& j : i->sosed)
		{
			fout << i->x << " " << i->y << " " << i->z << endl;
			if (j -> number >= 0)
			{
				fout << (i->x + j->x)/2.0 << " " << (i->y + j->y) / 2.0 << " " << (i->z + j->z) / 2.0 << endl;
			}
			else
			{
				fout << j->x << " " << j->y << " " << j->z << endl;
			}
		}
	}
	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}
}

int Konstruktor::get_size_conektiv(void)
{
	int ll = 0;
	for (auto& i : this->all_Kyb)
	{
		ll = ll + i->sosed.size();
	}
	return ll;
}

void Konstruktor::save_Setka(string name)
{
	ofstream fout, fout2;
	fout2.open("info_file.txt");
	fout2 << "M_inf = " << " " << M_inf << "  M_alf = " << M_alf << " step = " << step << "  M_0 = " << M_0 << "  kk = " << kk_ << "  chi = " << chi << "  rr_0 = " << rr_0 << endl;
	fout2.close();
	string name_ = "Golikov_Setka_file_" + name + ".txt";
	fout.open(name_);
	cout << " Save_  Yacheek vsego  " << this->all_Kyb.size() << endl;
	fout << this->all_Kyb.size() << endl;
	for (auto& i : this->all_Kyb)
	{
		fout << i->dx << " " << i->dy << " " << i->dz << " " << i->x << " " << i->y << " " << i->z << " " << //
			i->ro << " " << i->p << " " << i->u << " " << i->v << " " << i->w << " " << i->Bx << " " //
			<< i->By << " " << i->Bz << " " << i->Q << " " << i->number << endl;
	}
	fout << this->G1->dx << " " << this->G1->dy << " " << this->G1->dz << " " << this->G1->x << " " << this->G1->y << " " << this->G1->z << " " << //
		this->G1->ro << " " << this->G1->p << " " << this->G1->u << " " << this->G1->v << " " << this->G1->w << " " << this->G1->Bx << " " //
		<< this->G1->By << " " << this->G1->Bz << " " << this->G1->Q << " " << this->G1->number << endl;
	fout << this->G2->dx << " " << this->G2->dy << " " << this->G2->dz << " " << this->G2->x << " " << this->G2->y << " " << this->G2->z << " " << //
		this->G2->ro << " " << this->G2->p << " " << this->G2->u << " " << this->G2->v << " " << this->G2->w << " " << this->G2->Bx << " " //
		<< this->G2->By << " " << this->G2->Bz << " " << this->G2->Q << " " << this->G2->number << endl;
	
	fout << this->G3->dx << " " << this->G3->dy << " " << this->G3->dz << " " << this->G3->x << " " << this->G3->y << " " << this->G3->z << " " << //
		this->G3->ro << " " << this->G3->p << " " << this->G3->u << " " << this->G3->v << " " << this->G3->w << " " << this->G3->Bx << " " //
		<< this->G3->By << " " << this->G3->Bz << " " << this->G3->Q << " " << this->G3->number << endl;
	
	fout << this->G4->dx << " " << this->G4->dy << " " << this->G4->dz << " " << this->G4->x << " " << this->G4->y << " " << this->G4->z << " " << //
		this->G4->ro << " " << this->G4->p << " " << this->G4->u << " " << this->G4->v << " " << this->G4->w << " " << this->G4->Bx << " " //
		<< this->G4->By << " " << this->G4->Bz << " " << this->G4->Q << " " << this->G4->number << endl;
	
	fout << this->G5->dx << " " << this->G5->dy << " " << this->G5->dz << " " << this->G5->x << " " << this->G5->y << " " << this->G5->z << " " << //
		this->G5->ro << " " << this->G5->p << " " << this->G5->u << " " << this->G5->v << " " << this->G5->w << " " << this->G5->Bx << " " //
		<< this->G5->By << " " << this->G5->Bz << " " << this->G5->Q << " " << this->G5->number << endl;
	
	fout << this->G6->dx << " " << this->G6->dy << " " << this->G6->dz << " " << this->G6->x << " " << this->G6->y << " " << this->G6->z << " " << //
		this->G6->ro << " " << this->G6->p << " " << this->G6->u << " " << this->G6->v << " " << this->G6->w << " " << this->G6->Bx << " " //
		<< this->G6->By << " " << this->G6->Bz << " " << this->G6->Q << " " << this->G6->number;

	for (auto& i : this->all_Kyb)
	{
		fout << endl;
		fout << i->sosed.size();
		for (auto& j : i->sosed)
		{
			fout << " " <<  j->number;
		}
	}
}

void Konstruktor::download_Setka(string name)
{
	// Функция не работает, надо проверять
	ifstream fin;
	fin.open(name);
	if (!fin)
	{
		cout << "Could not open file! Check name file again!" << endl;
		exit(-1);
	}
	cout << "File is open!" << endl;
	int N = 0;
	fin >> N;
	double x, y, z, dx, dy, dz, ro, p, u, v, w, Bx, By, Bz, Q = 0.0;
	int num = 0;
	cout << "1 / 4" << endl;

	for (int i = 0; i < N; i++)
	{
		fin >> dx >> dy >> dz >> x >> y >> z >> ro >> p >> u >> v >> w >> Bx >> By >> Bz >> Q >>  num;
		auto C = new Kyb(x, y, z);
		C->dx = dx;
		C->dy = dy;
		C->dz = dz;
		C->number = num;
		C->x = x;
		C->y = y;
		C->z = z;
		C->ro = ro;
		C->p = p;
		C->u = u;
		C->v = v;
		C->w = w;
		C->Bx = Bx;
		C->By = By;
		C->Bz = Bz;
		C->Q = Q;
		this->all_Kyb.push_back(C);
	}
	cout << "2 / 4" << endl;

	fin >> dx >> dy >> dz >> x >> y >> z >> ro >> p >> u >> v >> w >> Bx >> By >> Bz  >> Q >> num;
	G1 = new Kyb(x, y, z);
	G1->dx = dx;G1->dy = dy;G1->dz = dz;
	G1->number = num;G1->x = x;G1->y = y;G1->z = z;G1->ro = ro;G1->p = p;G1->u = u;G1->v = v;
	G1->w = w;G1->Bx = Bx;G1->By = By;G1->Bz = Bz;
	fin >> dx >> dy >> dz >> x >> y >> z >> ro >> p >> u >> v >> w >> Bx >> By >> Bz >> Q >>  num;
	G2 = new Kyb(x, y, z);
	G2->dx = dx; G2->dy = dy; G2->dz = dz;
	G2->number = num; G2->x = x; G2->y = y; G2->z = z; G2->ro = ro; G2->p = p; G2->u = u; G2->v = v;
	G2->w = w; G2->Bx = Bx; G2->By = By; G2->Bz = Bz;
	fin >> dx >> dy >> dz >> x >> y >> z >> ro >> p >> u >> v >> w >> Bx >> By >> Bz >> Q >> num;
	G3 = new Kyb(x, y, z);
	G3->dx = dx; G3->dy = dy; G3->dz = dz;
	G3->number = num; G3->x = x; G3->y = y; G3->z = z; G3->ro = ro; G3->p = p; G3->u = u; G3->v = v;
	G3->w = w; G3->Bx = Bx; G3->By = By; G3->Bz = Bz;
	fin >> dx >> dy >> dz >> x >> y >> z >> ro >> p >> u >> v >> w >> Bx >> By >> Bz >> Q >> num;
	G4 = new Kyb(x, y, z);
	G4->dx = dx; G4->dy = dy; G4->dz = dz;
	G4->number = num; G4->x = x; G4->y = y; G4->z = z; G4->ro = ro; G4->p = p; G4->u = u; G4->v = v;
	G4->w = w; G4->Bx = Bx; G4->By = By; G4->Bz = Bz;
	fin >> dx >> dy >> dz >> x >> y >> z >> ro >> p >> u >> v >> w >> Bx >> By >> Bz >> Q >> num;
	G5 = new Kyb(x, y, z);
	G5->dx = dx; G5->dy = dy; G5->dz = dz;
	G5->number = num; G5->x = x; G5->y = y; G5->z = z; G5->ro = ro; G5->p = p; G5->u = u; G5->v = v;
	G5->w = w; G5->Bx = Bx; G5->By = By; G5->Bz = Bz;
	fin >> dx >> dy >> dz >> x >> y >> z >> ro >> p >> u >> v >> w >> Bx >> By >> Bz >> Q >> num;
	G6 = new Kyb(x, y, z);
	G6->dx = dx; G6->dy = dy; G6->dz = dz;
	G6->number = num; G6->x = x; G6->y = y; G6->z = z; G6->ro = ro; G6->p = p; G6->u = u; G6->v = v;
	G6->w = w; G6->Bx = Bx; G6->By = By; G6->Bz = Bz;

	int j = 0;
	int nn = 0;
	cout << "3 / 4" << endl;
	for (auto& i : this->all_Kyb)
	{
		fin >> j;
		/*if (j != 6)
		{
			cout << "error:  ne 6 yacheek  " << j << endl;
		}*/
		for (int k = 0; k < j; k++)
		{
			fin >> nn;
			if (nn == -1)
			{
				i->sosed.push_back(G1);
			}
			else if (nn == -2)
			{
				i->sosed.push_back(G2);
			}
			else if (nn == -3)
			{
				i->sosed.push_back(G3);
			}
			else if (nn == -4)
			{
				i->sosed.push_back(G4);
			}
			else if (nn == -5)
			{
				i->sosed.push_back(G5);
			}
			else if (nn == -6)
			{
				i->sosed.push_back(G6);
			}
			else
			{
				i->sosed.push_back(this->all_Kyb[nn]);
			}
		}
	}

	this->number();
	cout << "4 / 4" << endl;
	fin.close();
}

void Konstruktor::generate_Cuda_massiv(Cell* C, int* s)
{
	//C = new Cell[this->all_Kyb.size()];
	int nn = this->get_size_conektiv();
	//s = new int[nn];
	int c = 0;
	for (auto& i : this->all_Kyb)
	{
		for (auto& j : i->sosed)
		{
			s[c] = j->number;
			c++;
		}
	}

	int ll = 0;
	for (int i = 0; i < this->all_Kyb.size(); i++)
	{
		C[i].x = this->all_Kyb[i]->x;
		C[i].y = this->all_Kyb[i]->y;
		C[i].z = this->all_Kyb[i]->z;
		C[i].dx = this->all_Kyb[i]->dx;
		C[i].dy = this->all_Kyb[i]->dy;
		C[i].dz = this->all_Kyb[i]->dz;
		C[i].ro[0] = this->all_Kyb[i]->ro;
		C[i].p[0] = this->all_Kyb[i]->p;
		C[i].u[0] = this->all_Kyb[i]->u;
		C[i].v[0] = this->all_Kyb[i]->v;
		C[i].w[0] = this->all_Kyb[i]->w;
		C[i].Bx[0] = this->all_Kyb[i]->Bx;
		C[i].By[0] = this->all_Kyb[i]->By;
		C[i].Bz[0] = this->all_Kyb[i]->Bz;
		C[i].l = ll;
		C[i].r = ll + this->all_Kyb[i]->sosed.size() - 1;
		ll = ll + this->all_Kyb[i]->sosed.size();
	}
	return;
}

void Konstruktor::read_Cuda_massiv(double* ro, double* p, double* u, double* v, double* w, double* bx, double* by, double* bz, double* QQ)
{
	for (int i = 0; i < this->all_Kyb.size(); i++)
	{
		this->all_Kyb[i]->ro = ro[i];
		this->all_Kyb[i]->p =  p[i];
		this->all_Kyb[i]->u =  u[i];
		this->all_Kyb[i]->v =  v[i];
		this->all_Kyb[i]->w =  w[i];
		this->all_Kyb[i]->Bx = bx[i];
		this->all_Kyb[i]->By = by[i];
		this->all_Kyb[i]->Bz = bz[i];
		this->all_Kyb[i]->Q = QQ[i];
	}
}

void Konstruktor::Generate_sosed_for_TVD(int* s1, int* s2)
{
	int km = -1;
	double n1, n2, n3, v1, v2, v3, d, m, sk;
	int ni, nn;

	// Надо определить нужно ли вообще в этой ячейке делать ТВД процедуру
	for (auto& i : this->all_Kyb)
	{
		i->drob = true;
		for (auto& j : i->sosed)
		{
			if (j->number < 0 || fabs(i->dx - j->dx) > 0.001)
			{
				i->drob = false;
				break;
			}
		}
	}

	for (auto& i : this->all_Kyb)
	{
		if (i->drob == true)
		{
			for (auto& j : i->sosed)
			{
				for (auto& k : j->sosed)
				{
					if (k->number < 0 || fabs(j->dx - k->dx) > 0.001)
					{
						i->drob = false;
						break;
					}
				}
				if (i->drob == false)
				{
					break;
				}
			}
		}
	}

	// Определили, теперь занимаемся добавлением соседей
	int number = 0;
	for (auto& i : this->all_Kyb)
	{
		for (auto& j : i->sosed)
		{
			km++;
			if (i->drob == false)
			{
				s1[km] = -1;
				s2[km] = -1;
				continue;
			}
			number++;
			n1 = i->x - j->x;
			n2 = i->y - j->y;
			n3 = i->z - j->z;
			d = sqrt(kv(n1) + kv(n2) + kv(n3));
			n1 = n1 / d;
			n2 = n2 / d;
			n3 = n3 / d;
			ni = -1;
			nn = -1;
			m = 2.0;
			for (auto& k : i->sosed)
			{
				ni++;
				if (k->number != j->number)
				{
					v1 = i->x - k->x;
					v2 = i->y - k->y;
					v3 = i->z - k->z;
					d = sqrt(kv(v1) + kv(v2) + kv(v3));
					v1 = v1 / d;
					v2 = v2 / d;
					v3 = v3 / d;
					sk = skk(n1, n2, n3, v1, v2, v3);
					if (sk < m)
					{
						m = sk;
						nn = ni;
					}
				}
			}
			if (nn == -1)
			{
				cout << "error  bhgvshgvc2344343" << endl;
			}
			s1[km] = i->sosed[nn]->number;
		}
	}

	km = -1;
	for (auto& i : this->all_Kyb)
	{
		for (auto& j : i->sosed)
		{
			km++;
			if (i->drob == false)
			{
				s1[km] = -1;
				s2[km] = -1;
				continue;
			}
			
			n1 = j->x - i->x;
			n2 = j->y - i->y;
			n3 = j->z - i->z;
			d = sqrt(kv(n1) + kv(n2) + kv(n3));
			n1 = n1 / d;
			n2 = n2 / d;
			n3 = n3 / d;
			ni = -1;
			nn = -1;
			m = 2.0;
			for (auto& k : j->sosed)
			{
				ni++;
				if (k->number != i->number)
				{
					v1 = j->x - k->x;
					v2 = j->y - k->y;
					v3 = j->z - k->z;
					d = sqrt(kv(v1) + kv(v2) + kv(v3));
					v1 = v1 / d;
					v2 = v2 / d;
					v3 = v3 / d;
					sk = skk(n1, n2, n3, v1, v2, v3);
					if (sk < m)
					{
						m = sk;
						nn = ni;
					}
				}
			}
			if (nn == -1)
			{
				cout << "error  bhgvshgvc2344343" << endl;
			}
			s2[km] = j->sosed[nn]->number;
		}
	}

	cout << "TVD vvedeno dly " << number << " graney" << endl;
}

double Konstruktor::polar_angle(double x, double y)
{
	if (x * x + y * y < 0.0000001)
	{
		return 0.0;
	}
	else if (x < 0)
	{
		return atan(y / x) + 1.0 * PI;
	}
	else if (x > 0 && y >= 0)
	{
		return atan(y / x);
	}
	else if (x > 0 && y < 0)
	{
		return atan(y / x) + 2.0 * PI;
	}
	else if (y > 0 && x >= 0 && x <= 0)
	{
		return PI / 2.0;
	}
	else if (y < 0 && x >= 0 && x <= 0)
	{
		return  3.0 * PI / 2.0;
	}
	return 0.0;
}

void Konstruktor::spherical_skorost(double x, double y, double z, double Vx, double Vy, double Vz, double& Vr, double& Vphi, double& Vtheta)
{
	double r_1 = sqrt(x * x + y * y + z * z);
	double the_1 = acos(z / r_1);
	double phi_1 = this->polar_angle(x, y);

	Vr = Vx * sin(the_1) * cos(phi_1) + Vy * sin(the_1) * sin(phi_1) + Vz * cos(the_1);
	Vtheta = Vx * cos(the_1) * cos(phi_1) + Vy * cos(the_1) * sin(phi_1) - Vz * sin(the_1);
	Vphi = -Vx * sin(phi_1) + Vy * cos(phi_1);
}

void Konstruktor::dekard_skorost(double x, double y, double z, double Vr, double Vphi, double Vtheta, double& Vx, double& Vy, double& Vz)
{
	double r_2 = sqrt(x * x + y * y + z * z);
	double the_2 = acos(z / r_2);
	double phi_2 = this->polar_angle(x, y);

	if (sqrt(x * x + y * y) < 0.000001)
	{
		Vx = 0.0;
		Vy = 0.0;
		Vz = 0.0;
	}
	else
	{
		Vx = Vr * sin(the_2) * cos(phi_2) + Vtheta * cos(the_2) * cos(phi_2) - Vphi * sin(phi_2);
		Vy = Vr * sin(the_2) * sin(phi_2) + Vtheta * cos(the_2) * sin(phi_2) + Vphi * cos(phi_2);
		Vz = Vr * cos(the_2) - Vtheta * sin(the_2);
	}
}	

void Konstruktor::filling(void)
{

	double MM = kk_ /chi;
	double V_E = chi;
	double ro_E = MM / (4.0 * pi * V_E * rr_0 * rr_0);
	double P_E = ro_E * V_E * V_E / (ggg * M_0 * M_0);   // Мах другой в давлении
	double B_E = sqrt(kk_)/(M_alf * rr_0);

	for (auto& i : this->all_Kyb)
	{
		double dist = sqrt(i->x * i->x + i->y * i->y + i->z * i->z);
		//double dist2 = sqrt(kv(i->x + 0.8) + i->y * i->y + i->z * i->z);
		//double dist3 = kv(i->x + 1.8)/kv(2.9) + kv(i->y)/kv(2.9)  + kv(i->z)/kv(2.9);
		double dist3 = kv(i->x + 1.08) / kv(2.4) + kv(i->y) / kv(2.0) + kv(i->z) / kv(2.0);
		if (dist < 0.00005)
		{
			i->ro = 0.0;
			i->p = 0.0;
			i->u = 0.0;
			i->v = 0.0;
			i->w = 0.0;
			i->Bx = 0.0;
			i->By = 0.0;
			i->Bz = 0.0;
			i->Q = 0.0;
		}
		else if (dist3 <= 1.0)
		{
			/*i->ro = (1.0 / (phi_0 * phi_0 * r_0 * r_0)) * R * R / (dist * dist);
			i->p = 0.0204284 * pow(R / dist, 2.0 * ggg);
			i->u = phi_0 * i->x / dist;
			i->v = phi_0 * i->y / dist;
			i->w = phi_0 * i->z / dist;
			double BE = spi4 * epsilon_ * 2.0 * (1.0 / (dist/184.0));
			double the = acos(i->z / dist);
			double AA, BB, CC;
			this->dekard_skorost(i->x, i->y, i->z, 0.0, BE * sin(the), 0.0, AA, BB, CC);
			i->Bx = AA;
			i->By = BB;
			i->Bz = CC;
			i->Q = i->ro;*/
			/*i->ro = 0.25 / (dist*dist);
			i->p = ( 4.0 * 4.0/(ggg * ( 0.25 / (0.1 * 0.1) ) * M_0 * M_0 ) ) * pow(0.1 / dist, 2.0 * ggg);
			i->u = 4.0 * i->x / dist;
			i->v = 4.0 * i->y / dist;
			i->w = 4.0 * i->z / dist;
			double BE = spi4 * epsilon_ * 2.0 * (1.0 / dist);
			double the = acos(i->z / dist);
			double AA, BB, CC;
			this->dekard_skorost(i->x, i->y, i->z, 0.0, BE * sin(the), 0.0, AA, BB, CC);
			i->Bx = AA;
			i->By = BB;
			i->Bz = CC;
			i->Q = i->ro;*/
			i->ro = ro_E * rr_0 * rr_0 / (dist * dist);
			i->p = P_E * pow(rr_0 / dist, 2.0 * ggg);
			i->u = V_E * i->x / dist;
			i->v = V_E * i->y / dist;
			i->w = V_E * i->z / dist;
			double BE = B_E * (rr_0 / dist);
			double the = acos(i->z / dist);
			double AA, BB, CC;
			double BR = -B_E * kv((rr_0 / dist));    // Br
			this->dekard_skorost(i->x, i->y, i->z, BR, BE * sin(the), 0.0, AA, BB, CC);
			/*i->Bx = AA;
			i->By = BB;
			i->Bz = CC;*/
			i->Bx = 0.0;
			i->By = 0.0;
			i->Bz = 0.0;
			i->Q = i->ro;
		}
		else
		{
			i->ro = 1.0;
			i->p = 1.0/(ggg);
			i->u = -M_inf; //-1.0;
			i->v = 0.0;
			i->w = 0.0;
			i->Bx = -betta *cos(0.5235);
			i->By = -betta * sin(0.5235);
			i->Bz = 0.0;
			i->Q = 100.0;
		}
	}
}

void Konstruktor::filling_mini(void)
{
	double r_0 = 100.0 / 184.0;
	double R = 100.0;
	for (auto& i : this->all_Kyb)
	{
		double dist = sqrt(i->x * i->x + i->y * i->y + i->z * i->z);
		//double dist2 = sqrt(kv(i->x + 532.0) + kv(i->y - 100) + i->z * i->z);
		double dist3 = sqrt(kv(i->x + 155.0) + kv(i->y) + i->z * i->z);
		if (dist < 70)
		{
			i->ro = 0.0;
			i->p = 0.0;
			i->u = 0.0;
			i->v = 0.0;
			i->w = 0.0;
			i->Bx = 0.0;
			i->By = 0.0;
			i->Bz = 0.0;

		}
		else if (dist < 128.0 )
		{
			i->ro = (1.0 / (phi_0 * phi_0 * r_0 * r_0)) * R * R / (dist * dist);
			i->p = 0.0204284 * pow(R / dist, 2.0 * ggg);
			i->u = phi_0 * i->x / dist;
			i->v = phi_0 * i->y / dist;
			i->w = phi_0 * i->z / dist;
			i->Bx = 0.0;
			i->By = 0.0;
			i->Bz = 0.0;
			i->Q = i->ro;
		}

		/*if (dist2 < 110.0)
		{
			i->u = -0.2;
			i->Bx = 0.0;
			i->By = 0.0;
			i->Bz = 0.0;
		}*/
	}
}

void Konstruktor::filling_G_D(void)
{
	double r_0 = 0.05;
	for (auto& i : this->all_Kyb)
	{
		double dist = sqrt(i->x * i->x + i->y * i->y + i->z * i->z);
		if (dist < 0.05)
		{
			i->ro = 0.0;
			i->p = 0.0;
			i->u = 0.0;
			i->v = 0.0;
			i->w = 0.0;
			i->Bx = 0.0;
			i->By = 0.0;
			i->Bz = 0.0;

		}
		else if (dist <= 0.1)
		{
			i->ro = (((ggg + 1.0))/(ggg + 3.0) ) * pow(1.0/r_0 , 2.0) * pow(r_0 / dist, 2.0);
			i->p = (1.0/(ggg * M_0 * M_0)) * i->ro * pow(r_0 / dist, 2.0 * ggg);
			i->u = 1.0 * i->x / dist;
			i->v = 1.0 * i->y / dist;
			i->w = 1.0 * i->z / dist;
			i->Bx = 0.0;
			i->By = 0.0;
			i->Bz = 0.0;
		}
		else
		{
			i->ro = omega;
			i->p = 1.0;
			i->u = 0.0;
			i->v = 0.0;
			i->w = 0.0;
			i->Bx = 0.0;
			i->By = 0.0;
			i->Bz = 0.0;
		}
	}
}

void Konstruktor::print_Tecplot_z(double z, double T, string nam)
{
	double r_o = 0.25320769;
	ofstream fout;
	string name_f = "sd_z_" + to_string(z) + "__" + to_string(T) + "_" + nam + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Ro\", \"P\", \"P_all\", \"Vx\", \"Vy\", \"Vz\",  \"Bx\", \"By\", \"Bz\", \"BB\", \"Max\", \"Alf\",\"Q\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->all_Kyb)
	{
		if (fabs(i->z - z) <= i->dz)
		{
			double Max = 0.0;
			double Alf = 10000.0;
			double QQ = 0.0;
			if (i->ro > 0.000000000001)
			{
				QQ = i->Q / i->ro;
				Max = sqrt(kvv(i->u, i->v, i->w) / (ggg * i->p / i->ro));
				if (kvv(i->Bx, i->By, i->Bz) > 0)
				{
					Alf = sqrt(kvv(i->u, i->v, i->w) / (kvv(i->Bx, i->By, i->Bz) / (cpi4 * i->ro)));
				}
			}
			
			fout << i->x * r_o << " " << i->y * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->y * r_o * i->y * r_o) << " "  << i->ro << " " << i->p << " " //
				<< i->p + kvv(i->Bx, i->By, i->Bz) /cpi8 << " " << //
				i->u << " " << i->v << " " << i->w << //
				" " << i->Bx << " " << i->By << " " << i->Bz << " " << sqrt( kvv(i->Bx, i->By, i->Bz) ) <<  " " << Max << " " << Alf << " " << QQ << endl;
		}
	}
	fout.close();
}

void Konstruktor::print_Tecplot_z_20(double z, double T, string nam)
{
	double r_o = 0.25320769;
	ofstream fout;
	string name_f = "sd_z_" + to_string(z) + "__" + to_string(T) + "_" + nam + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Ro\", \"P\", \"P_all\", \"Vx\", \"Vy\", \"Vz\", \"VV\", \"Bx\", \"By\", \"Bz\", \"BB\", \"Max\", \"Alf\",\"Q\", \"jjj\", \"jx\", \"jy\",\"jz\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->all_Kyb)
	{
		if (fabs(i->z - z) <= i->dz)
		{
			double Max = 0.0;
			double Alf = 10000.0;
			double QQ = 0.0;
			if (i->ro > 0.000000000001)
			{
				QQ = i->Q / i->ro;
				Max = sqrt(kvv(i->u, i->v, i->w) / (ggg * i->p / i->ro));
				if (kvv(i->Bx, i->By, i->Bz) > 0)
				{
					Alf = sqrt(kvv(i->u, i->v, i->w) / (kvv(i->Bx, i->By, i->Bz) / (cpi4 * i->ro)));
				}
			}

			fout << i->x * r_o << " " << i->y * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->y * r_o * i->y * r_o) << " " << i->ro << " " << i->p << " " //
				<< i->p + kvv(i->Bx, i->By, i->Bz) / cpi8 << " " << //
				i->u << " " << i->v << " " << i->w << " " << sqrt(kvv(i->u, i->v, i->w)) << //
				" " << i->Bx << " " << i->By << " " << i->Bz << " " << sqrt(kvv(i->Bx, i->By, i->Bz)) << " " << Max << " " << Alf << " " << QQ << " " <<  sqrt(kv(i->jx) + kv(i->jy) + kv(i->jz)) << " " << i->jx << " " << i->jy << " " << i->jz << endl;

			//// Симметричные условия
			//fout << i->x * r_o << " " << -i->y * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->y * r_o * i->y * r_o) << " " << i->ro << " " << i->p << " " //
			//	<< i->p + kvv(i->Bx, i->By, i->Bz) / cpi8 << " " << //
			//	i->u << " " << -i->v << " " << i->w << " " << sqrt(kvv(i->u, i->v, i->w)) << //
			//	" " << -i->Bx << " " << i->By << " " << -i->Bz << " " << sqrt(kvv(i->Bx, i->By, i->Bz)) << " " << Max << " " << Alf << " " << QQ << " " << sqrt(kv(i->jx) + kv(i->jy) + kv(i->jz)) << " " << -i->jx << " " << i->jy << " " << -i->jz << endl;
		}
	}
	fout.close();
}

void Konstruktor::print_Tecplot_z_j(double z, double T, string nam)
{
	double r_o = 0.25320769;
	ofstream fout;
	string name_f = "sd_z_" + to_string(z) + "__" + to_string(T) + "_" + nam + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\",  \"Rho\", \"jjj\", \"jx\", \"jy\",\"jz\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->all_Kyb)
	{
		if ((fabs(i->z - z) <= i->dz)&&(i->j_ == true) )
		{
			fout << i->x * r_o << " " << i->y * r_o << " "  << sqrt(i->x * r_o * i->x * r_o + i->y * r_o * i->y * r_o) << " " <<  i->ro << " " << sqrt(kv(i->jx) + kv(i->jy) + kv(i->jz)) << " " <<   i->jx << " " << i->jy << " " << i->jz << endl;
		}
	}
	fout.close();
}

void Konstruktor::print_Tecplot_x(double x, double T, string nam)
{
	double r_o = 0.25320769;
	ofstream fout;
	string name_f = "sd_x_" + to_string(x) + "__" + to_string(T) + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"Y\", \"Z\", \"r\", \"Ro\", \"P\", \"P_all\", \"Vx\", \"Vy\", \"Vz\",  \"Bx\", \"By\", \"Bz\", \"BB\", \"Max\", \"Alf\", \"Q\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->all_Kyb)
	{
		if (fabs(i->x - x) <= i->dx)
		{
			double Max = 0.0;
			double Alf = 10000.0;
			double QQ = 0.0;
			if (i->ro > 0.000000000001)
			{
				QQ = i->Q / i->ro;
				Max = sqrt(kvv(i->u, i->v, i->w) / (ggg * i->p / i->ro));
				if (kvv(i->Bx, i->By, i->Bz) > 0.0001)
				{
					Alf = sqrt(kvv(i->u, i->v, i->w) / (kvv(i->Bx, i->By, i->Bz) / (cpi4 * i->ro)));
				}
			}

			fout << i->y * r_o << " " << i->z * r_o << " " << sqrt(i->z * r_o * i->z * r_o + i->y * r_o * i->y * r_o) << " " << i->ro << " " << i->p << " " //
				<< i->p + kvv(i->Bx, i->By, i->Bz) / cpi8 << " " << //
				i->u << " " << i->v << " " << i->w << //
				" " << i->Bx << " " << i->By << " " << i->Bz << " " << sqrt(kvv(i->Bx, i->By, i->Bz)) << " " << Max << " " << Alf << " " << QQ << endl;
		}
	}
	fout.close();
}

void Konstruktor::print_Tecplot_x_20(double x, double T, string nam)
{
	double r_o = 0.25320769;
	ofstream fout;
	string name_f = "sd_x_" + to_string(x) + "__" + to_string(T) + "_" + nam + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"Y\", \"Z\", \"r\", \"Ro\", \"P\", \"P_all\", \"Vx\", \"Vy\", \"Vz\", \"VV\", \"Bx\", \"By\", \"Bz\", \"BB\", \"Max\", \"Alf\",\"Q\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->all_Kyb)
	{
		if (fabs(i->x - x) <= i->dx)
		{
			double Max = 0.0;
			double Alf = 10000.0;
			double QQ = 0.0;
			if (i->ro > 0.000000000001)
			{
				QQ = i->Q / i->ro;
				Max = sqrt(kvv(i->u, i->v, i->w) / (ggg * i->p / i->ro));
				if (kvv(i->Bx, i->By, i->Bz) > 0.0001)
				{
					Alf = sqrt(kvv(i->u, i->v, i->w) / (kvv(i->Bx, i->By, i->Bz) / (cpi4 * i->ro)));
				}
			}

			fout << i->y * r_o << " " << i->z * r_o << " " << sqrt(i->z * r_o * i->z * r_o + i->y * r_o * i->y * r_o) << " " << i->ro << " " << i->p << " " //
				<< i->p + kvv(i->Bx, i->By, i->Bz) / cpi8 << " " << //
				i->u << " " << i->v << " " << i->w << " " <<  sqrt(kvv(i->u, i->v, i->w)) << //
				" " << i->Bx << " " << i->By << " " << i->Bz << " " << sqrt(kvv(i->Bx, i->By, i->Bz)) << " " << Max << " " << Alf << " " << QQ << endl;

			fout << -i->y * r_o << " " << i->z * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->y * r_o * i->y * r_o) << " " << i->ro << " " << i->p << " " //
				<< i->p + kvv(i->Bx, i->By, i->Bz) / cpi8 << " " << //
				i->u << " " << -i->v << " " << i->w << " " << sqrt(kvv(i->u, i->v, i->w)) << //
				" " << -i->Bx << " " << i->By << " " << -i->Bz << " " << sqrt(kvv(i->Bx, i->By, i->Bz)) << " " << Max << " " << Alf << " " << QQ << endl;

			fout << i->y * r_o << " " << -i->z * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->z * r_o * i->z * r_o) << " " << i->ro << " " << i->p << " " //
				<< i->p + kvv(i->Bx, i->By, i->Bz) / cpi8 << " " << //
				i->u << " " << i->v << " " << -i->w << " " << sqrt(kvv(i->u, i->v, i->w)) << //
				" " << i->Bx << " " << i->By << " " << -i->Bz << " " << sqrt(kvv(i->Bx, i->By, i->Bz)) << " " << Max << " " << Alf << " " << QQ << endl;

			fout << -i->y * r_o << " " << -i->z * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->z * r_o * i->z * r_o) << " " << i->ro << " " << i->p << " " //
				<< i->p + kvv(i->Bx, i->By, i->Bz) / cpi8 << " " << //
				i->u << " " << -i->v << " " << -i->w << " " << sqrt(kvv(i->u, i->v, i->w)) << //
				" " << -i->Bx << " " << i->By << " " << i->Bz << " " << sqrt(kvv(i->Bx, i->By, i->Bz)) << " " << Max << " " << Alf << " " << QQ << endl;
		}
	}
	fout.close();
}

void Konstruktor::print_Tecplot_y(double y, double T, string nam)
{
	double r_o = 0.25320769;
	ofstream fout;
	string name_f = "sd_y_" + to_string(y) + "__" + to_string(T) + "_" + nam + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"Z\", \"X\", \"r\", \"Ro\", \"P\", \"P_all\", \"Vx\", \"Vy\", \"Vz\",  \"Bx\", \"By\", \"Bz\", \"BB\", \"Max\", \"Alf\", \"Q\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->all_Kyb)
	{
		if (fabs(i->y - y) <= i->dy)
		{
			double Max = 0.0;
			double Alf = 10000.0;
			double QQ = 0.0;
			if (i->ro > 0.000000000001)
			{
				QQ = i->Q / i->ro;
				Max = sqrt(kvv(i->u, i->v, i->w) / (ggg * i->p / i->ro));
				if (kvv(i->Bx, i->By, i->Bz) > 0)
				{
					Alf = sqrt(kvv(i->u, i->v, i->w) / (kvv(i->Bx, i->By, i->Bz) / (cpi4 * i->ro)));
				}
			}

			fout << i->z * r_o << " " << i->x * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->z * r_o * i->z * r_o) << " " << i->ro << " " << i->p << " " //
				<< i->p + kvv(i->Bx, i->By, i->Bz) / cpi8 << " " << //
				i->u << " " << i->v << " " << i->w << //
				" " << i->Bx << " " << i->By << " " << i->Bz << " " << sqrt(kvv(i->Bx, i->By, i->Bz)) << " " << Max << " " << Alf  << " " << QQ << endl;
		}
	}
	fout.close();
}

void Konstruktor::print_Tecplot_y_20(double y, double T, string nam)
{
	double r_o = 0.25320769;
	ofstream fout;
	string name_f = "sd_y_" + to_string(y) + "__" + to_string(T) + "_" + nam + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"x\", \"Z\", \"r\", \"Ro\", \"P\", \"P_all\", \"Vx\", \"Vy\", \"Vz\", \"VV\",  \"Bx\", \"By\", \"Bz\", \"BB\", \"Max\", \"Alf\",\"Q\", \"jjj\", \"jx\", \"jy\",\"jz\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->all_Kyb)
	{
		if (fabs(i->y - y) <= i->dy)
		{
			double Max = 0.0;
			double Alf = 10000.0;
			double QQ = 0.0;
			if (i->ro > 0.000000000001)
			{
				QQ = i->Q / i->ro;
				Max = sqrt(kvv(i->u, i->v, i->w) / (ggg * i->p / i->ro));
				if (kvv(i->Bx, i->By, i->Bz) > 0)
				{
					Alf = sqrt(kvv(i->u, i->v, i->w) / (kvv(i->Bx, i->By, i->Bz) / (cpi4 * i->ro)));
				}
			}

			fout << i->x * r_o << " " << i->z * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->z * r_o * i->z * r_o) << " " << i->ro << " " << i->p << " " //
				<< i->p + kvv(i->Bx, i->By, i->Bz) / cpi8 << " " << //
				i->u << " " << i->v << " " << i->w << " " << sqrt(kvv(i->u, i->v, i->w)) << //
				" " << i->Bx << " " << i->By << " " << i->Bz << " " << sqrt(kvv(i->Bx, i->By, i->Bz)) << " " << Max << " " << Alf << " " << QQ << " " << sqrt(kv(i->jx) + kv(i->jy) + kv(i->jz)) << " " << i->jx << " " << i->jy << " " << i->jz << endl;

			//fout << i->x * r_o << " " << -i->z * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->z * r_o * i->z * r_o) << " " << i->ro << " " << i->p << " " //
			//	<< i->p + kvv(i->Bx, i->By, i->Bz) / cpi8 << " " << //
			//	i->u << " " << i->v << " " << -i->w << " " << sqrt(kvv(i->u, i->v, i->w)) << //
			//	" " << i->Bx << " " << i->By << " " << -i->Bz << " " << sqrt(kvv(i->Bx, i->By, i->Bz)) << " " << Max << " " << Alf << " " << QQ << " " << sqrt(kv(i->jx) + kv(i->jy) + kv(i->jz)) << " " << -i->jx << " " << -i->jy << " " << i->jz << endl;
		
		}
	}
	fout.close();
}

void Konstruktor::print_Tecplot_y_j(double y, double T, string nam)
{
	double r_o = 0.25320769;
	ofstream fout;
	string name_f = "sd_y_" + to_string(y) + "__" + to_string(T) + "_" + nam + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES =  \"X\", \"Z\", \"r\",  \"Rho\", \"jjj\", \"jx\", \"jy\",\"jz\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->all_Kyb)
	{
		if ((fabs(i->y - y) <= i->dy) && (i->j_ == true) )
		{
			fout << i->x * r_o << " " << i->z * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->y * r_o * i->y * r_o) << " " << i->ro << " " << sqrt(kv(i->jx) + kv(i->jy) + kv(i->jz)) << " " << i->jx << " " << i->jy << " " << i->jz << endl;
		}
	}
	fout.close();
}

void Konstruktor::print_some_point()
{
	ofstream fout;
	string name_f =  "Contakt.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"Z\", \"PP\",  ZONE T = \"HP\"" << endl;
	for (auto& i : this->all_Kyb)
	{
		if (i->Q > 50 && i->Q < 70)
		{
			fout << i->x << " " << i->y << " " << i->z << " " << i->p + kvv(i->Bx, i->By, i->Bz) / cpi8 << endl;
		}
	}
}

void Konstruktor::print_3D()
{
	ofstream fout;
	string name_f = "3D.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"Z\", \"Ro\", \"Vx\",\"Vy\",\"Vz\", ZONE T = \"HP\"" << endl;
	srand(123);
	for (auto& i : this->all_Kyb)
	{
		if (rand() % 60 == 0)
		{
			fout << i->x << " " << i->y << " " << i->z << " " << i->ro << " " <<  i->u << " " << i->v << " " << i->w  << endl;
		}
	}
}

void Konstruktor::print_3D_20(string nam)
{
	double r_o = 0.25320769;
	ofstream fout;
	string name_f = "3D_" + nam + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"Z\", \"Ro\", \"Vx\",\"Vy\",\"Vz\",\"Bx\",\"By\",\"Bz\" ZONE T = \"HP\"" << endl;
	srand(123);
	for (auto& i : this->all_Kyb)
	{
		if (rand() % 50 == 0)
		{
			fout << i->x * r_o << " " << i->y * r_o << " " << i->z * r_o << " " << i->ro << " " << i->u << " " << i->v << " " << i->w << " " << i->Bx << " " << i->By << " " << i->Bz << endl;

			fout << i->x * r_o << " " << -i->y * r_o << " " << i->z * r_o << " " << i->ro << " " << i->u << " " << -i->v << " " << i->w << " " << -i->Bx << " " << i->By << " " << -i->Bz << endl;

			fout << i->x * r_o << " " << i->y * r_o << " " << -i->z * r_o << " " << i->ro << " " << i->u << " " << i->v << " " << -i->w << " " << i->Bx << " " << i->By << " " << -i->Bz << endl;

			fout << i->x * r_o << " " << -i->y * r_o << " " << -i->z * r_o << " " << i->ro << " " << i->u << " " << -i->v << " " << -i->w << " " << -i->Bx << " " << i->By << " " << i->Bz << endl;
		}
	}
}

bool Konstruktor::sosed_or_not(Kyb* A, Kyb* B)
{
	bool a = (fabs(fabs(A->x - B->x) - A->dx - B->dx) < geo);
	bool b = (fabs(fabs(A->y - B->y) - A->dy - B->dy) < geo);
	bool c = (fabs(fabs(A->z - B->z) - A->dz - B->dz) < geo);

	int golos = 0;
	if (a == true)
	{
		golos++;
	}
	if (b == true)
	{
		golos++;
	}
	if (c == true)
	{
		golos++;
	}
	if (golos != 1)
	{
		return false;
	}

	bool d = fabs(A->x - B->x) < A->dx + B->dx;
	bool e = fabs(A->y - B->y) < A->dy + B->dy;
	bool g = fabs(A->z - B->z) < A->dz + B->dz;


	if ((a) && (e) && (g))
	{
		return true;
	}
	else if ((d) && (b) && (g))
	{
		return true;
	}
	else if ((d) && (e) && (c) )
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Konstruktor::find_sosed(Kyb* A)
{
	A->sosed.clear();
	for (auto& i : this->all_Kyb)
	{
		if (i->number != A->number)
		{
			if (sosed_or_not(A, i))
			{
				A->sosed.push_back(i);
			}
		}
	}
}

void Konstruktor::number(void)
{
	int n = 0;
	for (auto& i : this->all_Kyb)
	{
		if (i->number >= 0)
		{
			i->number = n;
			n++;
		}
	}
}
  
void Konstruktor::droblenie(Kyb* A, int NN)
{
	double dx_1 = 2.0 * A->dx / NN;
	double dy_1 = 2.0 * A->dy / NN;
	double dz_1 = 2.0 * A->dz / NN;
	double x_min = A->x - A->dx;
	double y_min = A->y - A->dy;
	double z_min = A->z - A->dz;

	vector <Kyb*> NEW;
	vector <Kyb*> Push;
	vector <Kyb*> Push2;
	Push.push_back(A);
	NEW.push_back(A);
	for (auto& k : A->sosed)
	{
		Push.push_back(k);
		Push2.push_back(k);
	}

	for (auto& k : A->sosed)  // Удаляю эту ячейку у соседей
	{
		vector <Kyb*> SS;
		for (auto& j : k->sosed)
		{
			if (j->number != A->number)
			{
				SS.push_back(j);
			}
		}
		k->sosed.clear();
		for (auto& j : SS)
		{
			k->sosed.push_back(j);
		}
	}


	A->sosed.clear();

	for (int k = 0; k < NN; k++)
	{
		//cout << "z = " << z_min + dz / 2.0 + dz * k << endl;
		for (int j = 0; j < NN; j++)
		{
			for (int i = 0; i < NN; i++)
			{
				double x = x_min + dx_1 / 2.0 + dx_1 * i;
				double y = y_min + dy_1 / 2.0 + dy_1 * j;
				double z = z_min + dz_1 / 2.0 + dz_1 * k;
				if ((k == 0) && (j == 0) && (i == 0))
				{
					A->x = x;
					A->y = y;
					A->z = z;
					A->dx = dx_1 / 2.0;
					A->dy = dy_1 / 2.0;
					A->dz = dz_1 / 2.0;
					for (auto& kk : Push2)
					{
						if (sosed_or_not(kk, A))
						{
							kk->sosed.push_back(A);
						}
					}
				}
				else
				{
					auto C = new Kyb(x, y, z);
					C->dx = dx_1 / 2.0;
					C->dy = dy_1 / 2.0;
					C->dz = dz_1 / 2.0;
					this->all_Kyb.push_back(C);
					Push.push_back(C);
					NEW.push_back(C);
					for (auto& kk : Push2)
					{
						if (sosed_or_not(kk, C))
						{
							kk->sosed.push_back(C);
						}
					}
				}
			}
		}
	}

	Push2.clear();

	this->number();

	for (auto& i : NEW)
	{
		i->sosed.clear();
		for (auto& j : Push)
		{
			if (i->number != j->number)
			{
				if (sosed_or_not(j, i))
				{
					i->sosed.push_back(j);
				}
			}
		}
	}
}

void Konstruktor::droblenie_fast(Kyb* A, int NN)
{
	double dx_1 = 2.0 * A->dx / NN;
	double dy_1 = 2.0 * A->dy / NN;
	double dz_1 = 2.0 * A->dz / NN;
	double x_min = A->x - A->dx;
	double y_min = A->y - A->dy;
	double z_min = A->z - A->dz;

	vector <Kyb*> NEW;
	vector <Kyb*> Push;
	vector <Kyb*> Push2;
	Push.push_back(A);
	NEW.push_back(A);
	for (auto& k : A->sosed)
	{
		Push.push_back(k);
		Push2.push_back(k);
	}

	for (auto& k : A->sosed)  // Удаляю эту ячейку у соседей
	// Т.е. чтобы в списке соседей у соседей А не было самой А
	{
		vector <Kyb*> SS;
		for (auto& j : k->sosed)
		{
			if (j->number != A->number)
			{
				SS.push_back(j);
			}
		}
		k->sosed.clear();
		for (auto& j : SS)
		{
			k->sosed.push_back(j);
		}
	}


	A->sosed.clear();

	for (int k = 0; k < NN; k++)
	{
		//cout << "z = " << z_min + dz / 2.0 + dz * k << endl;
		for (int j = 0; j < NN; j++)
		{
			for (int i = 0; i < NN; i++)
			{
				double x = x_min + dx_1 / 2.0 + dx_1 * i;
				double y = y_min + dy_1 / 2.0 + dy_1 * j;
				double z = z_min + dz_1 / 2.0 + dz_1 * k;
				if ((k == 0) && (j == 0) && (i == 0))
				{
					A->x = x;
					A->y = y;
					A->z = z;
					A->dx = dx_1 / 2.0;
					A->dy = dy_1 / 2.0;
					A->dz = dz_1 / 2.0;
					for (auto& kk : Push2)
					{
						if (sosed_or_not(kk, A))
						{
							kk->sosed.push_back(A);
						}
					}
				}
				else
				{
					auto C = new Kyb(x, y, z);
					C->dx = dx_1 / 2.0;
					C->dy = dy_1 / 2.0;
					C->dz = dz_1 / 2.0;
					this->all_Kyb.push_back(C);
					Push.push_back(C);
					NEW.push_back(C);
					if ((k == 0) || (j == 0) || (i == 0) || (k == NN - 1) || (j == NN - 1) || (i == NN - 1))
					{
						for (auto& kk : Push2)
						{
							if (sosed_or_not(kk, C))
							{
								kk->sosed.push_back(C);
							}
						}
					}
				}
			}
		}
	}

	Push2.clear();

	this->number();

	for (auto& i : NEW)
	{
		i->sosed.clear();
		for (auto& j : Push)
		{
			if (i->number != j->number)
			{
				if (sosed_or_not(j, i))
				{
					i->sosed.push_back(j);
				}
			}
		}
	}
}

void Konstruktor::droblenie2_hand(Kyb* A)
{
	double dx_1 = A->dx;
	double dy_1 = A->dy;
	double dz_1 = A->dz;
	double xx = A->x;
	double yy = A->y;
	double zz = A->z;

	for (auto& k : A->sosed)
	{
		if (k->number >= 0)
		{
			if (k->x > A->x && fabs(fabs(k->x - A->x) - dx_1 - k->dx) < geo)
			{
				Sosed_1.push_back(k);
			}
			else if (k->x < A->x && fabs(fabs(k->x - A->x) - dx_1 - k->dx) < geo)
			{
				Sosed_2.push_back(k);
			}
			else if (k->y > A->y && fabs(fabs(k->y - A->y) - dy_1 - k->dy) < geo)
			{
				Sosed_3.push_back(k);
			}
			else if (k->y < A->y && fabs(fabs(k->y - A->y) - dy_1 - k->dy) < geo)
			{
				Sosed_4.push_back(k);
			}
			else if (k->z > A->z && fabs(fabs(k->z - A->z) - dz_1 - k->dz) < geo)
			{
				Sosed_5.push_back(k);
			}
			else if (k->z < A->z && fabs(fabs(k->z - A->z) - dz_1 - k->dz) < geo)
			{
				Sosed_6.push_back(k);
			}
			else
			{
				cout << "Error 1086 hfgdhdgdd " << k->x << " " << k->y << " " << k->z << " " << k->dx << " " << k->dy << " " << k->dz << //
					" " << A->x << " " << A->y << " " << A->z << " " << A->dx << " " << A->dy << " " << A->dz << endl;
			}
		}
	}

	/*if (Sosed_6.size() == 0)
	{
		cout << "error  Sosed_6.size()" << endl;
		exit(-1);
	}*/

	for (auto& k : A->sosed)  // Удаляю эту ячейку у соседей
	// Т.е. чтобы в списке соседей у соседей А не было самой А
	{
		if (k->number >= 0)
		{
			int im = -1;
			bool ff = false;
			for (auto& j : k->sosed)
			{
				im++;
				if (j == A)
				{
					ff = true;
					break;
				}
			}
			if (ff = false)
			{
				cout << "erfer  1123" << endl;
			}
			auto iter = k->sosed.begin(); // указатель на первый элемент
			k->sosed.erase(iter + im);   // удаляем элемент
		}
	}

	A->sosed.clear();


	double ddx_1 = dx_1 / 2.0;
	double ddy_1 = dy_1 / 2.0;
	double ddz_1 = dz_1 / 2.0;

	double x = xx - ddx_1;
	double y = yy - ddy_1;
	double z = zz - ddz_1;
	A->x = x;
	A->y = y;
	A->z = z;
	A->dx = ddx_1;
	A->dy = ddy_1;
	A->dz = ddz_1;
	if (this->x_min > x)
	{
		this->x_min = x;
	}
	if (this->y_min > y)
	{
		this->y_min = y;
	}
	if (this->z_min > z)
	{
		this->z_min = z;
	}

	x = xx + ddx_1;
	y = yy - ddy_1;
	z = zz - ddz_1;
	auto C2 = new Kyb(x, y, z);
	C2->dx = A->dx;
	C2->dy = A->dy;
	C2->dz = A->dz;
	C2->ro = A->ro; C2->p = A->p; C2->u = A->u; C2->v = A->v; C2->w = A->w; C2->Bx = A->Bx; C2->By = A->By; C2->Bz = A->Bz; C2->Q = A->Q;
	this->all_Kyb.push_back(C2);

	x = xx - ddx_1;
	y = yy + ddy_1;
	z = zz - ddz_1;
	auto C3 = new Kyb(x, y, z);
	C3->dx = A->dx;
	C3->dy = A->dy;
	C3->dz = A->dz;
	C3->ro = A->ro; C3->p = A->p; C3->u = A->u; C3->v = A->v; C3->w = A->w; C3->Bx = A->Bx; C3->By = A->By; C3->Bz = A->Bz; C3->Q = A->Q;
	this->all_Kyb.push_back(C3);

	x = xx + ddx_1;
	y = yy + ddy_1;
	z = zz - ddz_1;
	auto C4 = new Kyb(x, y, z);
	C4->dx = A->dx;
	C4->dy = A->dy;
	C4->dz = A->dz;
	C4->ro = A->ro; C4->p = A->p; C4->u = A->u; C4->v = A->v; C4->w = A->w; C4->Bx = A->Bx; C4->By = A->By; C4->Bz = A->Bz; C4->Q = A->Q;
	this->all_Kyb.push_back(C4);

	x = xx - ddx_1;
	y = yy - ddy_1;
	z = zz + ddz_1;
	auto C5 = new Kyb(x, y, z);
	C5->dx = A->dx;
	C5->dy = A->dy;
	C5->dz = A->dz;
	C5->ro = A->ro; C5->p = A->p; C5->u = A->u; C5->v = A->v; C5->w = A->w; C5->Bx = A->Bx; C5->By = A->By; C5->Bz = A->Bz; C5->Q = A->Q;
	this->all_Kyb.push_back(C5);

	x = xx + ddx_1;
	y = yy - ddy_1;
	z = zz + ddz_1;
	auto C6 = new Kyb(x, y, z);
	C6->dx = A->dx;
	C6->dy = A->dy;
	C6->dz = A->dz;
	C6->ro = A->ro; C6->p = A->p; C6->u = A->u; C6->v = A->v; C6->w = A->w; C6->Bx = A->Bx; C6->By = A->By; C6->Bz = A->Bz; C6->Q = A->Q;
	this->all_Kyb.push_back(C6);

	x = xx - ddx_1;
	y = yy + ddy_1;
	z = zz + ddz_1;
	auto C7 = new Kyb(x, y, z);
	C7->dx = A->dx;
	C7->dy = A->dy;
	C7->dz = A->dz;
	C7->ro = A->ro; C7->p = A->p; C7->u = A->u; C7->v = A->v; C7->w = A->w; C7->Bx = A->Bx; C7->By = A->By; C7->Bz = A->Bz; C7->Q = A->Q;
	this->all_Kyb.push_back(C7);

	x = xx + ddx_1;
	y = yy + ddy_1;
	z = zz + ddz_1;
	auto C8 = new Kyb(x, y, z);
	C8->dx = A->dx;
	C8->dy = A->dy;
	C8->dz = A->dz;
	C8->ro = A->ro; C8->p = A->p; C8->u = A->u; C8->v = A->v; C8->w = A->w; C8->Bx = A->Bx; C8->By = A->By; C8->Bz = A->Bz; C8->Q = A->Q;
	this->all_Kyb.push_back(C8);
	if (this->x_max < x)
	{
		this->x_max = x;
	}
	if (this->y_max < y)
	{
		this->y_max = y;
	}
	if (this->z_max < z)
	{
		this->z_max = z;
	}

	A->sosed.push_back(C2); A->sosed.push_back(C3); A->sosed.push_back(C5);

	C2->sosed.push_back(A); C2->sosed.push_back(C4); C2->sosed.push_back(C6);

	C3->sosed.push_back(A); C3->sosed.push_back(C4); C3->sosed.push_back(C7);

	C4->sosed.push_back(C2); C4->sosed.push_back(C3); C4->sosed.push_back(C8);

	C5->sosed.push_back(A); C5->sosed.push_back(C6); C5->sosed.push_back(C7);

	C6->sosed.push_back(C2); C6->sosed.push_back(C5); C6->sosed.push_back(C8);

	C7->sosed.push_back(C3); C7->sosed.push_back(C5); C7->sosed.push_back(C8);

	C8->sosed.push_back(C4); C8->sosed.push_back(C6); C8->sosed.push_back(C7);

	for (auto& k : Sosed_1)
	{
		if (sosed_or_not(k, C2))
		{
			k->sosed.push_back(C2);
			C2->sosed.push_back(k);
		}
		if (sosed_or_not(k, C4))
		{
			k->sosed.push_back(C4);
			C4->sosed.push_back(k);
		}
		if (sosed_or_not(k, C6))
		{
			k->sosed.push_back(C6);
			C6->sosed.push_back(k);
		}
		if (sosed_or_not(k, C8))
		{
			k->sosed.push_back(C8);
			C8->sosed.push_back(k);
		}
	}

	for (auto& k : Sosed_2)
	{
		if (sosed_or_not(k, A))
		{
			k->sosed.push_back(A);
			A->sosed.push_back(k);
		}
		if (sosed_or_not(k, C3))
		{
			k->sosed.push_back(C3);
			C3->sosed.push_back(k);
		}
		if (sosed_or_not(k, C5))
		{
			k->sosed.push_back(C5);
			C5->sosed.push_back(k);
		}
		if (sosed_or_not(k, C7))
		{
			k->sosed.push_back(C7);
			C7->sosed.push_back(k);
		}
	}

	for (auto& k : Sosed_3)
	{
		if (sosed_or_not(k, C3))
		{
			k->sosed.push_back(C3);
			C3->sosed.push_back(k);
		}
		if (sosed_or_not(k, C4))
		{
			k->sosed.push_back(C4);
			C4->sosed.push_back(k);
		}
		if (sosed_or_not(k, C7))
		{
			k->sosed.push_back(C7);
			C7->sosed.push_back(k);
		}
		if (sosed_or_not(k, C8))
		{
			k->sosed.push_back(C8);
			C8->sosed.push_back(k);
		}
	}

	for (auto& k : Sosed_4)
	{
		if (sosed_or_not(k, A))
		{
			k->sosed.push_back(A);
			A->sosed.push_back(k);
		}
		if (sosed_or_not(k, C2))
		{
			k->sosed.push_back(C2);
			C2->sosed.push_back(k);
		}
		if (sosed_or_not(k, C5))
		{
			k->sosed.push_back(C5);
			C5->sosed.push_back(k);
		}
		if (sosed_or_not(k, C6))
		{
			k->sosed.push_back(C6);
			C6->sosed.push_back(k);
		}
	}

	for (auto& k : Sosed_5)
	{
		if (sosed_or_not(k, C5))
		{
			k->sosed.push_back(C5);
			C5->sosed.push_back(k);
		}
		if (sosed_or_not(k, C6))
		{
			k->sosed.push_back(C6);
			C6->sosed.push_back(k);
		}
		if (sosed_or_not(k, C7))
		{
			k->sosed.push_back(C7);
			C7->sosed.push_back(k);
		}
		if (sosed_or_not(k, C8))
		{
			k->sosed.push_back(C8);
			C8->sosed.push_back(k);
		}
	}

	for (auto& k : Sosed_6)
	{
		if (sosed_or_not(k, A))
		{
			k->sosed.push_back(A);
			A->sosed.push_back(k);
		}
		if (sosed_or_not(k, C2))
		{
			k->sosed.push_back(C2);
			C2->sosed.push_back(k);
		}
		if (sosed_or_not(k, C3))
		{
			k->sosed.push_back(C3);
			C3->sosed.push_back(k);
		}
		if (sosed_or_not(k, C4))
		{
			k->sosed.push_back(C4);
			C4->sosed.push_back(k);
		}
	}

	// Здесь нужно добавить границу, если надо)
	if (true)
	{
		if (sosed_or_not(G1, A))
		{
			A->sosed.push_back(G1);
		}
		if (sosed_or_not(G2, A))
		{
			A->sosed.push_back(G2);
		}
		if (sosed_or_not(G3, A))
		{
			A->sosed.push_back(G3);
		}
		if (sosed_or_not(G4, A))
		{
			A->sosed.push_back(G4);
		}
		if (sosed_or_not(G5, A))
		{
			A->sosed.push_back(G5);
		}
		if (sosed_or_not(G6, A))
		{
			A->sosed.push_back(G6);
		}


		if (sosed_or_not(G1, C2))
		{
			C2->sosed.push_back(G1);
		}
		if (sosed_or_not(G2, C2))
		{
			C2->sosed.push_back(G2);
		}
		if (sosed_or_not(G3, C2))
		{
			C2->sosed.push_back(G3);
		}
		if (sosed_or_not(G4, C2))
		{
			C2->sosed.push_back(G4);
		}
		if (sosed_or_not(G5, C2))
		{
			C2->sosed.push_back(G5);
		}
		if (sosed_or_not(G6, C2))
		{
			C2->sosed.push_back(G6);
		}


		if (sosed_or_not(G1, C3))
		{
			C3->sosed.push_back(G1);
		}
		if (sosed_or_not(G2, C3))
		{
			C3->sosed.push_back(G2);
		}
		if (sosed_or_not(G3, C3))
		{
			C3->sosed.push_back(G3);
		}
		if (sosed_or_not(G4, C3))
		{
			C3->sosed.push_back(G4);
		}
		if (sosed_or_not(G5, C3))
		{
			C3->sosed.push_back(G5);
		}
		if (sosed_or_not(G6, C3))
		{
			C3->sosed.push_back(G6);
		}


		if (sosed_or_not(G1, C4))
		{
			C4->sosed.push_back(G1);
		}
		if (sosed_or_not(G2, C4))
		{
			C4->sosed.push_back(G2);
		}
		if (sosed_or_not(G3, C4))
		{
			C4->sosed.push_back(G3);
		}
		if (sosed_or_not(G4, C4))
		{
			C4->sosed.push_back(G4);
		}
		if (sosed_or_not(G5, C4))
		{
			C4->sosed.push_back(G5);
		}
		if (sosed_or_not(G6, C4))
		{
			C4->sosed.push_back(G6);
		}


		if (sosed_or_not(G1, C5))
		{
			C5->sosed.push_back(G1);
		}
		if (sosed_or_not(G2, C5))
		{
			C5->sosed.push_back(G2);
		}
		if (sosed_or_not(G3, C5))
		{
			C5->sosed.push_back(G3);
		}
		if (sosed_or_not(G4, C5))
		{
			C5->sosed.push_back(G4);
		}
		if (sosed_or_not(G5, C5))
		{
			C5->sosed.push_back(G5);
		}
		if (sosed_or_not(G6, C5))
		{
			C5->sosed.push_back(G6);
		}


		if (sosed_or_not(G1, C6))
		{
			C6->sosed.push_back(G1);
		}
		if (sosed_or_not(G2, C6))
		{
			C6->sosed.push_back(G2);
		}
		if (sosed_or_not(G3, C6))
		{
			C6->sosed.push_back(G3);
		}
		if (sosed_or_not(G4, C6))
		{
			C6->sosed.push_back(G4);
		}
		if (sosed_or_not(G5, C6))
		{
			C6->sosed.push_back(G5);
		}
		if (sosed_or_not(G6, C6))
		{
			C6->sosed.push_back(G6);
		}


		if (sosed_or_not(G1, C7))
		{
			C7->sosed.push_back(G1);
		}
		if (sosed_or_not(G2, C7))
		{
			C7->sosed.push_back(G2);
		}
		if (sosed_or_not(G3, C7))
		{
			C7->sosed.push_back(G3);
		}
		if (sosed_or_not(G4, C7))
		{
			C7->sosed.push_back(G4);
		}
		if (sosed_or_not(G5, C7))
		{
			C7->sosed.push_back(G5);
		}
		if (sosed_or_not(G6, C7))
		{
			C7->sosed.push_back(G6);
		}


		if (sosed_or_not(G1, C8))
		{
			C8->sosed.push_back(G1);
		}
		if (sosed_or_not(G2, C8))
		{
			C8->sosed.push_back(G2);
		}
		if (sosed_or_not(G3, C8))
		{
			C8->sosed.push_back(G3);
		}
		if (sosed_or_not(G4, C8))
		{
			C8->sosed.push_back(G4);
		}
		if (sosed_or_not(G5, C8))
		{
			C8->sosed.push_back(G5);
		}
		if (sosed_or_not(G6, C8))
		{
			C8->sosed.push_back(G6);
		}
	}

	Sosed_1.clear();
	Sosed_2.clear();
	Sosed_3.clear();
	Sosed_4.clear();
	Sosed_5.clear();
	Sosed_6.clear();
	//this->number();
}

void Konstruktor::Drobim(double x1, double x2, double y1, double y2, double z1, double z2, int NN)
{
	for (auto& i : this->all_Kyb)
	{
		i->drob = false;
	}

	int ll = 0;
	for (auto& i : this->all_Kyb)
	{
		if ((i->x > x1) && (i->x < x2) && (i->y > y1) && (i->y < y2) && (i->z > z1) && (i->z < z2))
		{
			ll++;
			i->drob = true;
		}
	}
	
	int mm = this->all_Kyb.size();
	cout << ll << endl;
	for (int i = 0; i < mm; i++)
	{
		

		if (this->all_Kyb[i]->drob == true)
		{
			if (NN == 2)
			{
				this->droblenie2_hand(this->all_Kyb[i]);
			}
			else
			{
				this->droblenie_fast(this->all_Kyb[i], NN);
			}

			ll--;
			if (ll % 25000 == 0)
			{
				cout << ll << endl;
			}
		}
	}
	this->number();
}

void Konstruktor::Drobim(double x0, double y0, double z0, double r1, double r2, int NN, bool ff)
{
	for (auto& i : this->all_Kyb)
	{
		i->drob = false;
	}

	int ll = 0;
	double r;
	for (auto& i : this->all_Kyb)
	{
		r = sqrt(kv(i->x - x0) + kv(i->y - y0) + kv(i->z - z0));

		if ( (r > r1) && (r < r2) )
		{
			ll++;
			i->drob = true;
		}
	}

	int mm = this->all_Kyb.size();

	for (int i = 0; i < mm; i++)
	{
		if (this->all_Kyb[i]->drob == true)
		{
			if (NN == 2)
			{
				this->droblenie2_hand(this->all_Kyb[i]);
			}
			else 
			{
				this->droblenie_fast(this->all_Kyb[i], NN);
			}
			ll--;
			if (ll % 25000 == 0)
			{
				cout << ll << endl;
			}
		}
	}
	this->number();
}

void Konstruktor::Drobim(double x1, double x2, double r, int NN)
{
	for (auto& i : this->all_Kyb)
	{
		i->drob = false;
	}

	int ll = 0;
	double rr;
	for (auto& i : this->all_Kyb)
	{
		rr = sqrt(kv(i->y) + kv(i->z));
		if ( (i->x > x1) && (i->x < x2) && (rr < r) )
		{
				ll++;
				i->drob = true;
		}
	}

	int mm = this->all_Kyb.size();

	for (int i = 0; i < mm; i++)
	{
		if (this->all_Kyb[i]->drob == true)
		{
			if (NN == 2)
			{
				this->droblenie2_hand(this->all_Kyb[i]);
			}
			else
			{
				this->droblenie_fast(this->all_Kyb[i], NN);
			}
			ll--;
			if (ll % 25000 == 0)
			{
				cout << ll << endl;
			}
		}
	}
	this->number();
}

void Konstruktor::droblenie3(Kyb* A)
{
	double dx_1 = A->dx;
	double dy_1 = A->dy;
	double dz_1 = A->dz;
	double xx = A->x;
	double yy = A->y;
	double zz = A->z;


	double ddx_1 = dx_1 / 2.0;
	double ddy_1 = dy_1 / 2.0;
	double ddz_1 = dz_1 / 2.0;

	double x = xx - ddx_1;
	double y = yy - ddy_1;
	double z = zz - ddz_1;
	A->x = x;
	A->y = y;
	A->z = z;
	A->dx = ddx_1;
	A->dy = ddy_1;
	A->dz = ddz_1;
	if (this->x_min > x)
	{
		this->x_min = x;
	}
	if (this->y_min > y)
	{
		this->y_min = y;
	}
	if (this->z_min > z)
	{
		this->z_min = z;
	}

	x = xx + ddx_1;
	y = yy - ddy_1;
	z = zz - ddz_1;
	auto C2 = new Kyb(x, y, z);
	C2->dx = A->dx;
	C2->dy = A->dy;
	C2->dz = A->dz;
	C2->ro = A->ro; C2->p = A->p; C2->u = A->u; C2->v = A->v; C2->w = A->w; C2->Bx = A->Bx; C2->By = A->By; C2->Bz = A->Bz;
	this->all_Kyb.push_back(C2);

	x = xx - ddx_1;
	y = yy + ddy_1;
	z = zz - ddz_1;
	auto C3 = new Kyb(x, y, z);
	C3->dx = A->dx;
	C3->dy = A->dy;
	C3->dz = A->dz;
	C3->ro = A->ro; C3->p = A->p; C3->u = A->u; C3->v = A->v; C3->w = A->w; C3->Bx = A->Bx; C3->By = A->By; C3->Bz = A->Bz;
	this->all_Kyb.push_back(C3);

	x = xx + ddx_1;
	y = yy + ddy_1;
	z = zz - ddz_1;
	auto C4 = new Kyb(x, y, z);
	C4->dx = A->dx;
	C4->dy = A->dy;
	C4->dz = A->dz;
	C4->ro = A->ro; C4->p = A->p; C4->u = A->u; C4->v = A->v; C4->w = A->w; C4->Bx = A->Bx; C4->By = A->By; C4->Bz = A->Bz;
	this->all_Kyb.push_back(C4);

	x = xx - ddx_1;
	y = yy - ddy_1;
	z = zz + ddz_1;
	auto C5 = new Kyb(x, y, z);
	C5->dx = A->dx;
	C5->dy = A->dy;
	C5->dz = A->dz;
	C5->ro = A->ro; C5->p = A->p; C5->u = A->u; C5->v = A->v; C5->w = A->w; C5->Bx = A->Bx; C5->By = A->By; C5->Bz = A->Bz;
	this->all_Kyb.push_back(C5);

	x = xx + ddx_1;
	y = yy - ddy_1;
	z = zz + ddz_1;
	auto C6 = new Kyb(x, y, z);
	C6->dx = A->dx;
	C6->dy = A->dy;
	C6->dz = A->dz;
	C6->ro = A->ro; C6->p = A->p; C6->u = A->u; C6->v = A->v; C6->w = A->w; C6->Bx = A->Bx; C6->By = A->By; C6->Bz = A->Bz;
	this->all_Kyb.push_back(C6);

	x = xx - ddx_1;
	y = yy + ddy_1;
	z = zz + ddz_1;
	auto C7 = new Kyb(x, y, z);
	C7->dx = A->dx;
	C7->dy = A->dy;
	C7->dz = A->dz;
	C7->ro = A->ro; C7->p = A->p; C7->u = A->u; C7->v = A->v; C7->w = A->w; C7->Bx = A->Bx; C7->By = A->By; C7->Bz = A->Bz;
	this->all_Kyb.push_back(C7);

	x = xx + ddx_1;
	y = yy + ddy_1;
	z = zz + ddz_1;
	auto C8 = new Kyb(x, y, z);
	C8->dx = A->dx;
	C8->dy = A->dy;
	C8->dz = A->dz;
	C8->ro = A->ro; C8->p = A->p; C8->u = A->u; C8->v = A->v; C8->w = A->w; C8->Bx = A->Bx; C8->By = A->By; C8->Bz = A->Bz;
	this->all_Kyb.push_back(C8);
	if (this->x_max < x)
	{
		this->x_max = x;
	}
	if (this->y_max < y)
	{
		this->y_max = y;
	}
	if (this->z_max < z)
	{
		this->z_max = z;
	}
}

void Konstruktor::Delenie(double r1, double r2)
// Аккурантно, это функция не связывает ячейки, а просто дробит их
{
	for (auto& i : this->all_Kyb)
	{
		i->drob = false;
	}

	int ll = 0;
	double r;
	for (auto& i : this->all_Kyb)
	{
		r = sqrt(kv(i->x) + kv(i->y) + kv(i->z));
		if ((r > r1) && (r < r2))
		{
			ll++;
			i->drob = true;
		}
	}

	int mm = this->all_Kyb.size();

	for (int i = 0; i < mm; i++)
	{
		if (this->all_Kyb[i]->drob == true)
		{
			this->droblenie3(this->all_Kyb[i]);
			ll--;
			if (ll % 25000 == 0)
			{
				cout << ll << endl;
			}
		}
	}
	this->number();
}

void Konstruktor::Delenie(double x1, double x2, double y1, double y2, double z1, double z2)
{
	for (auto& i : this->all_Kyb)
	{
		i->drob = false;
	}

	int ll = 0;
	for (auto& i : this->all_Kyb)
	{
		if ((i->x > x1) && (i->x < x2) && (i->y > y1) && (i->y < y2) && (i->z > z1) && (i->z < z2))
		{
			ll++;
			i->drob = true;
		}
	}

	int mm = this->all_Kyb.size();
	cout << ll << endl;
	for (int i = 0; i < mm; i++)
	{
		if (this->all_Kyb[i]->drob == true)
		{
			this->droblenie3(this->all_Kyb[i]);

			ll--;
			if (ll % 25000 == 0)
			{
				cout << ll << endl;
			}
		}
	}
	this->number();
}

void Konstruktor::konect(double R)
{
	cout << "START" << endl;
	for (auto& i : this->all_Kyb)
	{
		i->sosed.clear();
	}

	int ll = this->all_Kyb.size();
	cout << ll << endl;


	int mm = this->all_Kyb.size();

		for (int i = 0; i < mm; i++)
		{
			for (auto& j : this->all_Kyb)
			{
				if (this->all_Kyb[i]->number != j->number)
				{
					if (fabs(this->all_Kyb[i]->x - j->x) <= R && fabs(this->all_Kyb[i]->y - j->y) <= R && fabs(this->all_Kyb[i]->z - j->z) <= R)
					{
						if (sosed_or_not(this->all_Kyb[i], j))
						{
							this->all_Kyb[i]->sosed.push_back(j);
							j->sosed.push_back(this->all_Kyb[i]);
						}
					}
				}
			}
				ll--;
				if (ll % 5000 == 0)
				{
					cout << ll << endl;
				}
		}

	for (auto& i : this->all_Kyb)
	{
		if (i->x > this->x_max - geo)
		{
			i->sosed.push_back(G1);
		}
		if (i->x < this->x_min + geo)
		{
			i->sosed.push_back(G2);
		}
		if (i->y > this->y_max - geo)
		{
			i->sosed.push_back(G3);
		}
		if (i->y < this->y_min + geo)
		{
			i->sosed.push_back(G4);
		}
		if (i->z > this->z_max - geo)
		{
			i->sosed.push_back(G5);
		}
		if (i->z < this->z_min + geo)
		{
			i->sosed.push_back(G6);
		}
	}
}

bool Konstruktor::get_square(Kyb* A, Kyb* B)
{
	// Аналог девайс функции для проверки правильности построения сетки
	if (fabs(fabs(A->x - B->x) - A->dx - B->dx) < 0.0004)
	{
		return true;
	}
	else if (fabs(fabs(A->y - B->y) - A->dy - B->dy) < 0.0004)
	{
		return true;
	}
	else if (fabs(fabs(A->z - B->z) - A->dz - B->dz) < 0.0004)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Konstruktor::count_j(void)
{
	for (auto& i : this->all_Kyb)
	{
		if (i->sosed.size() == 6)   // Будем считать только для ячеек, у которых 6 соседей (т.е. одинаковый размер ячейки у соседей)
		{
			bool jk = true;
			for (auto& j : i->sosed)
			{
				if (fabs(j->dx - i->dx) > 10 * geo && j->number >= 0)
				{
					jk = false;
				}
			}
			if (jk)
			{
				i->j_ = true;
				double Bx1 = 0.0, Bx2 = 0.0, Bx3 = 0.0, Bx4 = 0.0, Bx5 = 0.0, Bx6 = 0.0, By1 = 0.0, By2 = 0.0, By3 = 0.0, By4 = 0.0, By5 = 0.0, By6 = 0.0, Bz1 = 0.0, Bz2 = 0.0, Bz3 = 0.0, Bz4 = 0.0, Bz5 = 0.0, Bz6 = 0.0;
				int n1, n2, n3, n4, n5, n6;
				bool a1 = false, a2 = false, a3 = false, a4 = false, a5 = false, a6 = false;
				for (auto& j : i->sosed)
				{
					if (j->x > i->x && fabs(fabs(j->x - i->x) - i->dx - j->dx) < geo)
					{
						a1 = true;
						n1 = j->number;
						Bx1 = j->Bx;
						By1 = j->By;
						Bz1 = j->Bz;
					}
					else if (j->x < i->x && fabs(fabs(j->x - i->x) - i->dx - j->dx) < geo)
					{
						a2 = true;
						n2 = j->number;
						Bx2 = j->Bx;
						By2 = j->By;
						Bz2 = j->Bz;
					}
					else if (j->y > i->y && fabs(fabs(j->y - i->y) - i->dy - j->dy) < geo)
					{
						n3 = j->number;
						a3 = true;
						Bx3 = j->Bx;
						By3 = j->By;
						Bz3 = j->Bz;
					}
					else if (j->y < i->y && fabs(fabs(j->y - i->y) - i->dy - j->dy) < geo)
					{
						a4 = true;
						n4 = j->number;
						Bx4 = j->Bx;
						By4 = j->By;
						Bz4 = j->Bz;
					}
					else if (j->z > i->z && fabs(fabs(j->z - i->z) - i->dz - j->dz) < geo)
					{
						a5 = true;
						n5 = j->number;
						Bx5 = j->Bx;
						By5 = j->By;
						Bz5 = j->Bz;
					}
					else if (j->z < i->z && fabs(fabs(j->z - i->z) - i->dz - j->dz) < geo)
					{
						a6 = true;
						n6 = j->number;
						Bx6 = j->Bx;
						By6 = j->By;
						Bz6 = j->Bz;
					}
					else
					{
						cout << "Error 2585 in count_j" << j->x << " " << j->y << " " << j->z << " " << j->dx << " " << j->dy << " " << j->dz << //
							" " << i->x << " " << i->y << " " << i->z << " " << i->dx << " " << i->dy << " " << i->dz << endl;
						exit(-1);
					}
				}
				if (a1 == false || a2 == false || a3 == false || a4 == false || a5 == false || a6 == false)
				{
					cout << "Error 2610 in count_j" << endl;
					exit(-1);
				}
				if (n1 < 0)
				{
					Bx1 = i->Bx;
					By1 = i->By;
					Bz1 = i->Bz;
				}
				if (n2 < 0)
				{
					Bx2 = i->Bx;
					By2 = i->By;
					Bz2 = i->Bz;
				}
				if (n3 < 0)
				{
					Bx3 = i->Bx;
					By3 = i->By;
					Bz3 = i->Bz;
				}
				if (n4 < 0)
				{
					Bx4 = -i->Bx;
					By4 = i->By;
					Bz4 = -i->Bz;
				}
				if (n5 < 0)
				{
					Bx5 = i->Bx;
					By5 = i->By;
					Bz5 = i->Bz;
				}
				if (n6 < 0)
				{
					Bx6 = i->Bx;
					By6 = i->By;
					Bz6 = -i->Bz;
				}
				i->jx = (Bz3 - Bz4) / (4.0 * i->dy) - (By5 - By6) / (4.0 * i->dz);
				i->jy = (Bx5 - Bx6) / (4.0 * i->dz) - (Bz1 - Bz2) / (4.0 * i->dx);
				i->jz = (By1 - By2) / (4.0 * i->dx) - (Bx3 - Bx4) / (4.0 * i->dy);
			}
		}
	}
}