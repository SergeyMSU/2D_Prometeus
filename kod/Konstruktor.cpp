#include "Konstruktor.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <mutex>


using namespace std;

Konstruktor::Konstruktor(int a, int b,  double xmin, double xmax, double ymax)
{
	this->initialization(a, b, xmin, xmax, ymax);
	this->New_design();                      // Создаёт начальную сетку (без связей)
	this->konectiviti();                     // Связывает начальные ячейки с их соседями
}

Konstruktor::~Konstruktor()
{
	for (auto& i : this->all_Kyb)
	{
		delete i;
	}
	this->all_Kyb.clear();
}

void Konstruktor::initialization(int a, int b, double xmin, double xmax, double ymax)
{
	this->N = a;
	this->M = b;
	this->DX = 0.0;
	this->DY = 0.0;
	this->x_min = xmin;
	this->x_max = xmax;
	this->y_max = ymax;
	this->Sosed_1.reserve(128);
	this->Sosed_2.reserve(128);
	this->Sosed_3.reserve(128);
	this->Sosed_4.reserve(128);
	this->all_Kyb.reserve(1000000);
}

void Konstruktor::New_design(void)
{
	double dx = (this->x_max - this->x_min) / (this->N);
	double dy = (this->y_max) / this->M;
	this->DX = dx;
	this->DY = dy;

	double xx = this->x_min + dx / 2.0;
	double yy = dy / 2.0;

	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < N; i++)
		{
			double y = yy + j * dy;
			double x = xx + i * dx;
			auto C = new Kyb(x, y);
			this->all_Kyb.push_back(C);
		}
	}
	
}

void Konstruktor::konectiviti(void)
{

	auto C1 = new Kyb(this->x_max + this->DX, this->y_max / 2.0);
	C1->number = -1;
	C1->size = 0;
	this->G1 = C1;

	auto C2 = new Kyb(x_min - this->DX, this->y_max / 2.0);
	C2->number = -2;
	C2->size = 0;
	this->G2 = C2;

	auto C3 = new Kyb((this->x_max + this->x_min) / 2.0, y_max + this->DY);
	C3->number = -3;
	C3->size = 0;
	this->G3 = C3;

	auto C4 = new Kyb((this->x_max + this->x_min) / 2.0,  - this->DY);
	C4->number = -4;
	C4->size = 0;
	this->G4 = C4;


	for (auto& i : this->all_Kyb)
	{
		int n = i->number % this->N;
		int m = (i->number - n) / this->N;
		// n,m,k - глобальный номер ячейки
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
			ll++;
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
			fout << i->x << " " << i->y << endl;
			fout << (i->x + j->x) / 2.0 << " " << (i->y + j->y) / 2.0 << endl;
		}
	}
	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}

	fout.close();
}

void Konstruktor::print_cell(void)
{
	int ll = this->all_Kyb.size();
	ofstream fout;
	fout.open("Setka_all_cell.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << 3 * ll;
	fout << " , E= " << 2 * ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->all_Kyb)
	{
		double dx_1 = (this->DX / pow(2, i->size - 1)) / 2.0;   // Половина длины ячейки
		double dy_1 = (this->DY / pow(2, i->size - 1)) / 2.0;   // Половина ширины ячейки
		fout << i->x + dx_1 << " " << i->y - dy_1 << endl;
		fout << i->x - dx_1 << " " << i->y - dy_1 << endl;
		fout << i->x - dx_1 << " " << i->y + dy_1 << endl;
	}
	for (int i = 0; i < ll; i++)
	{
		fout << 3 * i + 1 << " " << 3 * i + 2 << endl;
		fout << 3 * i + 2 << " " << 3 * i + 3 << endl;
	}

	fout.close();
}

void Konstruktor::print_point(void)
{
	ofstream fout;
	fout.open("Setka_point.txt");
	for (auto& i : this->all_Kyb)
	{
		fout << i->x << " " << i->y << endl;
	}
	fout.close();
}

void Konstruktor::print_Tecplot(void)
{
	double r_o = 1.0;    // Размер расстояния
	double ro_o = 1.0;   // Размер плотности
	double p_o = 1.0;   // Размер давления
	double u_o = 1.0;   // Размер давления

	auto res = this->Get_projection();


	ofstream fout;
	string name_f = "2D_tecplot.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Ro\", \"P\", \"Vx\", \"Vy\", \"Max\",\"Q\",\"F_n\",\"F_u\",\"F_v\",\"F_T\", \"I_u\",\"I_v\",\"I_T\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->all_Kyb)
	{
			double Max = 0.0;
			double QQ = 0.0;
			if (i->ro > 0.000000000001)
			{
				QQ = i->Q / i->ro;
				Max = sqrt(kvv(i->u, i->v, 0.0) / (ggg * i->p / i->ro));
			}

			fout << i->x * r_o << " " << i->y * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->y * r_o * i->y * r_o) << //
				" " << i->ro * ro_o << " " << i->p * p_o << " " //
				<< i->u * u_o << " " << i->v * u_o << " " << Max << " "  << QQ << " " << i->F_n * ro_o << " " << //
				i->F_u * u_o << " " << i->F_v * u_o << " " << i->F_T << " " <<//
				i->I_u << " " << i->I_v << " " << i->I_T << " " << endl;
		
	}
	for (auto& i : res)
	{
		double Max = 0.0;
		double QQ = 0.0;
		if (i->ro > 0.000000000001)
		{
			QQ = i->Q / i->ro;
			Max = sqrt(kvv(i->u, i->v, 0.0) / (ggg * i->p / i->ro));
		}

		fout << i->x * r_o << " " << i->y * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->y * r_o * i->y * r_o) << //
			" " << i->ro * ro_o << " " << i->p * p_o << " " //
			<< i->u * u_o << " " << i->v * u_o << " " << Max << " " << QQ << " " << i->F_n * ro_o << " " << //
			i->F_u * u_o << " " << i->F_v * u_o << " " << i->F_T << " " <<//
			i->I_u << " " << i->I_v << " " << i->I_T << " " << endl;

	}
	fout.close();

	name_f = "1D_tecplot.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Ro\", \"P\", \"Vx\", \"Vy\", \"Max\",\"Q\",\"F_n\",\"F_u\",\"F_v\",\"F_T\", \"I_u\",\"I_v\",\"I_T\", ZONE T = \"HP\"" << endl;

	for (auto& i : res)
	{
		double Max = 0.0;
		double QQ = 0.0;
		if (i->ro > 0.00000001)
		{
			QQ = i->Q / i->ro;
			Max = sqrt(kvv(i->u, i->v, 0.0) / (ggg * i->p / i->ro));
		}

		fout << i->x * r_o << //
			" " << i->ro * ro_o << " " << i->p * p_o << " " //
			<< i->u * u_o << " " << i->v * u_o << " " << Max << " " << QQ << " " << i->F_n * ro_o << " " << //
			i->F_u * u_o << " " << i->F_v * u_o << " " << i->F_T << " " <<//
			i->I_u << " " << i->I_v << " " << i->I_T << " " << endl;
	}
	fout.close();
}

void Konstruktor::print_Tecplot_multifluid(void)
{
	double r_o = 1.0;    // Размер расстояния
	double ro_o = 0.06;   // Размер плотности
	double ro_o_H = 0.18;   // Размер плотности
	double p_o = 1.0;   // Размер давления
	double u_o = 10.38;   // Размер давления

	auto res = this->Get_projection();

	ofstream fout;
	string name_f = "2D_tecplot_multifluid.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Ro\", \"P\", \"Vx\", \"Vy\", \"Max\",\"Q\"," << //
		"\"Ro_H1\", \"P_H1\", \"Vx_H1\", \"Vy_H1\",\"Ro_H2\", \"P_H2\", \"Vx_H2\", \"Vy_H2\",\"Ro_H3\", \"P_H3\", \"Vx_H3\", \"Vy_H3\"," << //
		"\"Ro_H4\", \"P_H4\", \"Vx_H4\", \"Vy_H4\", \"Ro_H\", " << "ZONE T = \"HP\"" << endl;
	for (auto& i : this->all_Kyb)
	{
		double Max = 0.0;
		double QQ = 0.0;
		if (i->ro > 0.000000000001)
		{
			QQ = i->Q / i->ro;
			Max = sqrt(kvv(i->u, i->v, 0.0) / (ggg * i->p / i->ro));
		}

		fout << i->x * r_o << " " << i->y * r_o << " " << sqrt(i->x * r_o * i->x * r_o + i->y * r_o * i->y * r_o) << //
			" " << i->ro * ro_o << " " << i->p * p_o << " " //
			<< i->u * u_o << " " << i->v * u_o << " " << Max << " " << QQ << " " << //
			i->ro_H1 * ro_o_H << " " << i->p_H1 * p_o << " "<< i->u_H1 * u_o << " " << i->v_H1 * u_o << " " << //
			i->ro_H2 * ro_o_H << " " << i->p_H2 * p_o << " " << i->u_H2 * u_o << " " << i->v_H2 * u_o << " " << //
			i->ro_H3 * ro_o_H << " " << i->p_H3 * p_o << " " << i->u_H3 * u_o << " " << i->v_H3 * u_o << " " << //
			i->ro_H4 * ro_o_H << " " << i->p_H4 * p_o << " " << i->u_H4 * u_o << " " << i->v_H4 * u_o << " " << //
			(i->ro_H1 + i->ro_H2 + i->ro_H3 + i->ro_H4) * ro_o_H  << endl;

	}
	for (auto& i : res)
	{
		double Max = 0.0;
		double QQ = 0.0;
		if (i->ro > 0.000000000001)
		{
			QQ = i->Q / i->ro;
			Max = sqrt(kvv(i->u, i->v, 0.0) / (ggg * i->p / i->ro));
		}

		fout << i->x * r_o << " " << 0.0 << " " << sqrt(i->x * r_o * i->x * r_o + i->y * r_o * i->y * r_o) << //
			" " << i->ro * ro_o << " " << i->p * p_o << " " //
			<< i->u * u_o << " " << i->v * u_o << " " << Max << " " << QQ << " " << //
			i->ro_H1 * ro_o_H << " " << i->p_H1 * p_o << " " << i->u_H1 * u_o << " " << i->v_H1 * u_o << " " << //
			i->ro_H2 * ro_o_H << " " << i->p_H2 * p_o << " " << i->u_H2 * u_o << " " << i->v_H2 * u_o << " " << //
			i->ro_H3 * ro_o_H << " " << i->p_H3 * p_o << " " << i->u_H3 * u_o << " " << i->v_H3 * u_o << " " << //
			i->ro_H4 * ro_o_H << " " << i->p_H4 * p_o << " " << i->u_H4 * u_o << " " << i->v_H4 * u_o << " " << //
			(i->ro_H1 + i->ro_H2 + i->ro_H3 + i->ro_H4) * ro_o_H << endl;

	}
	fout.close();

	name_f = "1D_tecplot_multifluid.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Ro\", \"P\", \"Vx\", \"Vy\", \"Max\",\"Q\"," << //
		"\"Ro_H1\", \"P_H1\", \"Vx_H1\", \"Vy_H1\",\"Ro_H2\", \"P_H2\", \"Vx_H2\", \"Vy_H2\",\"Ro_H3\", \"P_H3\", \"Vx_H3\", \"Vy_H3\"," << //
		"\"Ro_H4\", \"P_H4\", \"Vx_H4\", \"Vy_H4\", \"Ro_H\", " << "ZONE T = \"HP\"" << endl;

	for (auto& i : res)
	{
		double Max = 0.0;
		double QQ = 0.0;
		if (i->ro > 0.000000000001)
		{
			QQ = i->Q / i->ro;
			Max = sqrt(kvv(i->u, i->v, 0.0) / (ggg * i->p / i->ro));
		}

		fout << i->x * r_o << " " << i->ro * ro_o << " " << i->p * p_o << " " //
			<< i->u * u_o << " " << i->v * u_o << " " << Max << " " << QQ << " " << //
			i->ro_H1 * ro_o_H << " " << i->p_H1 * p_o << " " << i->u_H1 * u_o << " " << i->v_H1 * u_o << " " << //
			i->ro_H2 * ro_o_H << " " << i->p_H2 * p_o << " " << i->u_H2 * u_o << " " << i->v_H2 * u_o << " " << //
			i->ro_H3 * ro_o_H << " " << i->p_H3 * p_o << " " << i->u_H3 * u_o << " " << i->v_H3 * u_o << " " << //
			i->ro_H4 * ro_o_H << " " << i->p_H4 * p_o << " " << i->u_H4 * u_o << " " << i->v_H4 * u_o << " " << //
			(i->ro_H1 + i->ro_H2 + i->ro_H3 + i->ro_H4) * ro_o_H << endl;
	}
	fout.close();
}

bool Konstruktor::sosed_or_not(Kyb* A, Kyb* B)
{
	if (B->number < 0)
	{
		cout << "wiufygwuyrg24y824392342" << endl;
		exit(-1);
	}

	if (A->number >= 0)
	{
		double Adx = (this->DX / pow(2, A->size - 1)) / 2.0;
		double Bdx = (this->DX / pow(2, B->size - 1)) / 2.0;
		double Ady = (this->DY / pow(2, A->size - 1)) / 2.0;
		double Bdy = (this->DY / pow(2, B->size - 1)) / 2.0;

		bool a = (fabs(fabs(A->x - B->x) - Adx - Bdx) < geo);
		bool b = (fabs(fabs(A->y - B->y) - Ady - Bdy) < geo);

		int golos = 0;
		if (a == true)
		{
			golos++;
		}
		if (b == true)
		{
			golos++;
		}
		if (golos != 1)
		{
			return false;
		}

		bool d = fabs(A->x - B->x) < Adx + Bdx;
		bool e = fabs(A->y - B->y) < Ady + Bdy;


		if ((a) && (e))
		{
			return true;
		}
		else if ((d) && (b))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		double Bdx = (this->DX / pow(2, B->size - 1)) / 2.0;
		double Bdy = (this->DY / pow(2, B->size - 1)) / 2.0;
		if (A == G1)
		{
			if (B->x + Bdx + geo > this->x_max)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		if (A == G2)
		{
			if (B->x - Bdx - geo < this->x_min)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		if (A == G3)
		{
			if (B->y + Bdy + geo > this->y_max)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		if (A == G4)
		{
			if (B->y - Bdy - geo < 0.0)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
	return false;
}

void Konstruktor::droblenie2_hand(Kyb* A)
// Дробление одной ячейки
{
	double dx_1 = (this->DX / pow(2, A->size - 1)) / 2.0;   // Половина длины ячейки
	double dy_1 = (this->DY / pow(2, A->size - 1)) / 2.0;   // Половина ширины ячейки
	double xx = A->x;
	double yy = A->y;

	for (auto& k : A->sosed)
	{
		double dx_2 = (this->DX / pow(2, k->size - 1)) / 2.0;   // Половина длины ячейки
		double dy_2 = (this->DY / pow(2, k->size - 1)) / 2.0;   // Половина ширины ячейки
		double geo_toch = 0.001 * dx_2;

		if (k->number >= 0)
		{
			if (k->x > xx && fabs(fabs(k->x - xx) - dx_1 - dx_2) < geo_toch)  // Соседи справа
			{
				Sosed_1.push_back(k);
			}
			else if (k->x < xx && fabs(fabs(k->x - xx) - dx_1 - dx_2) < geo_toch)  // Соседи слева
			{
				Sosed_2.push_back(k);
			}
			else if (k->y > yy && fabs(fabs(k->y - yy) - dy_1 - dy_2) < geo_toch) // Соседи сверху
			{
				Sosed_3.push_back(k);
			}
			else if (k->y < yy && fabs(fabs(k->y - yy) - dy_1 - dy_2) < geo_toch) // Соседи снизу
			{
				Sosed_4.push_back(k);
			}
			else
			{
				cout << "Error 1086 hfgdhdgdwefwefwcecsvrgewecegve45veb6ed   "  << k->x << " " << k->y << " " << xx << " " << yy << //
					" " << dx_2 << " " << dx_1 << endl;
				cout << k->size << " " << this->DX << endl;
			}
		}
	}


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
				cout << "erfer  112ugweufg2t378e253dwed234255" << endl;
			}
			auto iter = k->sosed.begin(); // указатель на первый элемент
			k->sosed.erase(iter + im);   // удаляем элемент
		}
	}

	A->sosed.clear();  // Чистим соседей самой А


	double ddx_1 = dx_1 / 2.0;
	double ddy_1 = dy_1 / 2.0;

	double x = xx - ddx_1;
	double y = yy - ddy_1;
	A->x = x;
	A->y = y;
	A->size += 1;                    // Повысели размер ячейки на следующий

	x = xx + ddx_1;
	y = yy - ddy_1;
	auto C2 = new Kyb(x, y);
	C2->ro = A->ro; C2->p = A->p; C2->u = A->u; C2->v = A->v;  C2->Q = A->Q;
	C2->size = A->size;
	this->all_Kyb.push_back(C2);

	x = xx - ddx_1;
	y = yy + ddy_1;
	auto C3 = new Kyb(x, y);
	C3->ro = A->ro; C3->p = A->p; C3->u = A->u; C3->v = A->v; C3->Q = A->Q;
	C3->size = A->size;
	this->all_Kyb.push_back(C3);

	x = xx + ddx_1;
	y = yy + ddy_1;
	auto C4 = new Kyb(x, y);
	C4->ro = A->ro; C4->p = A->p; C4->u = A->u; C4->v = A->v;  C4->Q = A->Q;
	C4->size = A->size;
	this->all_Kyb.push_back(C4);


	A->sosed.push_back(C2); A->sosed.push_back(C3);

	C2->sosed.push_back(A); C2->sosed.push_back(C4); 

	C3->sosed.push_back(A); C3->sosed.push_back(C4); 

	C4->sosed.push_back(C2); C4->sosed.push_back(C3); 

	for (auto& k : Sosed_1)
	{
		if (this->sosed_or_not(k, C2))
		{
			k->sosed.push_back(C2);
			C2->sosed.push_back(k);
		}
		if (this->sosed_or_not(k, C4))
		{
			k->sosed.push_back(C4);
			C4->sosed.push_back(k);
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
		


	}

	Sosed_1.clear();
	Sosed_2.clear();
	Sosed_3.clear();
	Sosed_4.clear();
}

void Konstruktor::Drobim(double x0, double y0, double r1, double r2)
{
	for (auto& i : this->all_Kyb)
	{
		i->drob = false;
	}

	int ll = 0;
	double r;
	for (auto& i : this->all_Kyb)
	{
		r = sqrt(kv(i->x - x0) + kv(i->y - y0));

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
			this->droblenie2_hand(this->all_Kyb[i]);
			ll--;
			/*if (ll % 25000 == 0)
			{
				cout << ll << endl;
			}*/
		}
	}
	this->number();
}

void Konstruktor::Drobim(double x1, double x2, double r1)
{
	for (auto& i : this->all_Kyb)
	{
		i->drob = false;
	}

	int ll = 0;
	double r;
	for (auto& i : this->all_Kyb)
	{
		r = sqrt(kv(i->x) + kv(i->y));
		if ((i->x > x1) && (i->x < x2) && (i->y < r1) && (r > 0.8 * Distant))
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
			this->droblenie2_hand(this->all_Kyb[i]);
			ll--;
			/*if (ll % 25000 == 0)
			{
				cout << ll << endl;
			}*/
		}
	}
	this->number();
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

int Konstruktor::get_size_conektiv(void)
{
	int ll = 0;
	for (auto& i : this->all_Kyb)
	{
		ll = ll + i->sosed.size();
	}
	return ll;
}

void Konstruktor::initial_condition()
{
	for (auto& i : this->all_Kyb)
	{
		/*i->ro_H1 = 0.000001;
		i->p_H1 = 0.00000001;
		i->u_H1 = 0.0;
		i->v_H1 = 0.0;
		i->ro_H2 = 0.000001;
		i->p_H2 = 0.00000001;
		i->u_H2 = 0.0;
		i->v_H2 = 0.0;
		i->ro_H3 = 0.000001;
		i->p_H3 = 0.00000001;
		i->u_H3 = 0.0;
		i->v_H3 = 0.0;
		i->ro_H4 = 0.000001;
		i->p_H4 = 0.00000001;
		i->u_H4 = 0.0;
		i->v_H4 = 0.0;*/
		double dist = sqrt(i->x * i->x + i->y * i->y);
		double r_0 = 1.0;
		double ro = (389.988 * 389.988) / (chi_ * chi_);
		double P_E = ro * chi_ * chi_ / (ggg * 5.0 * 5.0);
		if (dist <= 1.0)
		{
			i->ro = ro / (dist * dist);
			i->p = P_E * pow(r_0 / dist, 2.0 * ggg);
			i->u = chi_ * i->x / dist;
			i->v = chi_ * i->y / dist;
			i->Q = ro * r_0 * r_0 / (dist * dist);

			i->ro_H1 = 0.0000001;
			i->p_H1 = 0.000000001;
			i->u_H1 = chi_ * i->x / dist;
			i->v_H1 = chi_ * i->y / dist;

			/*i->ro_H2 = 0.000001;
			i->p_H2 = 0.00000001;
			i->u_H2 = 0.0;
			i->v_H2 = 0.0;
			i->ro_H3 = 0.000001;
			i->p_H3 = 0.00000001;
			i->u_H3 = 0.0;
			i->v_H3 = 0.0;
			i->ro_H4 = 0.000001;
			i->p_H4 = 0.00000001;
			i->u_H4 = 0.0;
			i->v_H4 = 0.0;*/
		}
		else if (dist <= Distant)
		{
			i->ro = ro / (dist * dist);
			i->p = P_E * pow(r_0 / dist, 2.0 * ggg);
			i->u = chi_ * i->x / dist;
			i->v = chi_ * i->y / dist;
			i->Q = ro * r_0 * r_0 / (dist * dist);

			i->ro_H1 = 0.0000001; // 0.0012 / 0.18; //(2.0 * (sigma(chi_) * 389.988) / (2.0 * Kn_ * chi_))* (1.0 / dist - 0.1 / kv(dist));
			i->p_H1 = 0.000000001;
			i->u_H1 = chi_ * i->x / dist;
			i->v_H1 = chi_ * i->y / dist;
		}
	}
}

void Konstruktor::read_Cuda_massiv(double* ro, double* p, double* u, double* v, double* QQ)
{
	for (int i = 0; i < this->all_Kyb.size(); i++)
	{
		this->all_Kyb[i]->ro = ro[i];
		this->all_Kyb[i]->p = p[i];
		this->all_Kyb[i]->u = u[i];
		this->all_Kyb[i]->v = v[i];
		this->all_Kyb[i]->Q = QQ[i];
	}
}

void Konstruktor::read_Cuda_massiv(double* ro, double* p, double* u, double* v, double* QQ,//
	double* ro_H1, double* p_H1, double* u_H1, double* v_H1,//
	double* ro_H2, double* p_H2, double* u_H2, double* v_H2, //
	double* ro_H3, double* p_H3, double* u_H3, double* v_H3, //
	double* ro_H4, double* p_H4, double* u_H4, double* v_H4)
{
	for (int i = 0; i < this->all_Kyb.size(); i++)
	{
		this->all_Kyb[i]->ro = ro[i];
		this->all_Kyb[i]->p = p[i];
		this->all_Kyb[i]->u = u[i];
		this->all_Kyb[i]->v = v[i];
		this->all_Kyb[i]->Q = QQ[i];
		this->all_Kyb[i]->ro_H1 = ro_H1[i];
		this->all_Kyb[i]->p_H1 = p_H1[i];
		this->all_Kyb[i]->u_H1 = u_H1[i];
		this->all_Kyb[i]->v_H1 = v_H1[i];
		this->all_Kyb[i]->ro_H2 = ro_H2[i];
		this->all_Kyb[i]->p_H2 = p_H2[i];
		this->all_Kyb[i]->u_H2 = u_H2[i];
		this->all_Kyb[i]->v_H2 = v_H2[i];
		this->all_Kyb[i]->ro_H3 = ro_H3[i];
		this->all_Kyb[i]->p_H3 = p_H3[i];
		this->all_Kyb[i]->u_H3 = u_H3[i];
		this->all_Kyb[i]->v_H3 = v_H3[i];
		this->all_Kyb[i]->ro_H4 = ro_H4[i];
		this->all_Kyb[i]->p_H4 = p_H4[i];
		this->all_Kyb[i]->u_H4 = u_H4[i];
		this->all_Kyb[i]->v_H4 = v_H4[i];
	}
}

void Konstruktor::Save_setka(string name)
{
	int ll = this->all_Kyb.size();
	ofstream fout;
	fout.open(name);

	for (auto& i : this->all_Kyb)
	{
		fout << i->ro << " " << i->p << " " << i->u << " " << i->v << " " << i->Q << endl;
	}

	fout.close();
}

void Konstruktor::Save_setka_MK(string name)
{
	int ll = this->all_Kyb.size();
	ofstream fout;
	fout.open("MK_" + name);

	for (auto& i : this->all_Kyb)
	{
		fout << i->ro << " " << i->p << " " << i->u << " " << i->v << " " << i->Q << //
			" " << i->F_n << " " << i->F_u << " " << i->F_v << " " << i->F_T << " " << i->I_u << " " << i->I_v << " " << i->I_T << endl;
	}

	fout.close();
}

void Konstruktor::Download_setka_MK(string name)
{
	int ll = this->all_Kyb.size();
	ifstream fout;
	fout.open(name);

	for (auto& i : this->all_Kyb)
	{
		fout >> i->ro >> i->p >> i->u >> i->v >> i->Q >> i->F_n >> i->F_u >> i->F_v >> i->F_T >> i->I_u >> i->I_v >> i->I_T;
	}

	fout.close();
}

void Konstruktor::Save_setka_multifluid(string name)
{
	int ll = this->all_Kyb.size();
	ofstream fout;
	fout.open(name);

	for (auto& i : this->all_Kyb)
	{
		fout << i->ro << " " << i->p << " " << i->u << " " << i->v << " " << i->Q << " " << //
			i->ro_H1 << " " << i->p_H1 << " " << i->u_H1 << " " << i->v_H1 << " " << //
			i->ro_H2 << " " << i->p_H2 << " " << i->u_H2 << " " << i->v_H2 << " " << //
			i->ro_H3 << " " << i->p_H3 << " " << i->u_H3 << " " << i->v_H3 << " " << //
			i->ro_H4 << " " << i->p_H4 << " " << i->u_H4 << " " << i->v_H4 << endl;

	}

	fout.close();
}

void Konstruktor::Download_setka(string name)
{
	int ll = this->all_Kyb.size();
	ifstream fout;
	fout.open(name);

	for (auto& i : this->all_Kyb)
	{
		fout >> i->ro >> i->p >>  i->u >>  i->v >> i->Q;
	}

	fout.close();
}

void Konstruktor::Download_setka_multifluid(string name)
{
	int ll = this->all_Kyb.size();
	ifstream fout;
	fout.open(name);

	for (auto& i : this->all_Kyb)
	{
		fout >> i->ro >> i->p >> i->u >> i->v >> i->Q >> i->ro_H1 >> i->p_H1 >> i->u_H1 >> i->v_H1 >> //
			i->ro_H2 >> i->p_H2 >> i->u_H2 >> i->v_H2 >> i->ro_H3 >> i->p_H3 >> i->u_H3 >> i->v_H3 >> i->ro_H4 >> i->p_H4 >> i->u_H4 >> i->v_H4;
	}

	fout.close();
}

double Konstruktor::linear_funk(const double& x1, const double& y1, const double& x2, const double& y2, const double& t)
{
	double k = (y1 - y2) / (x1 - x2);
	double b = (y1 * x2 - y2 * x1) / (x2 - x1);
	return k * t + b;
}

vector <Kyb*> Konstruktor::Get_projection(void)
{
	vector <Kyb*> Boandary;
	vector <Kyb*> res;

	for (auto& i : this->all_Kyb)
	{
		for (auto& j : i->sosed)
		{
			if (j->number == -4)
			{
				Boandary.push_back(i);
				break;
			}
		}
	}

	sort(Boandary.begin(), Boandary.end(), [](Kyb* i, Kyb* j)
		{
			return (i->x < j->x);
		});

	double dy_1, dy_2;
	for (auto& i : Boandary)
	{
		dy_1 = (this->DY / pow(2, i->size - 1)) / 2.0;   // Половина ширины ячейки
		for (auto& j : i->sosed)
		{
			dy_2 = (this->DY / pow(2, j->size - 1)) / 2.0;   // Половина ширины ячейки
			if (fabs(fabs(j->y - i->y) - dy_1 - dy_2) < dy_1 * 0.001 && j->y > i->y)
			{
				auto C1 = new Kyb(i->x, 0.0);
				C1->ro = linear_funk(i->y, i->ro, j->y, j->ro, 0.0);
				C1->p = linear_funk(i->y, i->p, j->y, j->p, 0.0);
				if (C1->p <= 0.0)
				{
					C1->p = i->p;
				}
				C1->u = linear_funk(i->y, i->u, j->y, j->u, 0.0);
				C1->v = linear_funk(i->y, i->v, j->y, j->v, 0.0);
				C1->Q = linear_funk(i->y, i->Q, j->y, j->Q, 0.0);

				C1->ro_H1 = linear_funk(i->y, i->ro_H1, j->y, j->ro_H1, 0.0);
				C1->p_H1 = linear_funk(i->y, i->p_H1, j->y, j->p_H1, 0.0);
				C1->u_H1 = linear_funk(i->y, i->u_H1, j->y, j->u_H1, 0.0);
				C1->v_H1 = linear_funk(i->y, i->v_H1, j->y, j->v_H1, 0.0);

				C1->ro_H2 = linear_funk(i->y, i->ro_H2, j->y, j->ro_H2, 0.0);
				C1->p_H2 = linear_funk(i->y, i->p_H2, j->y, j->p_H2, 0.0);
				C1->u_H2 = linear_funk(i->y, i->u_H2, j->y, j->u_H2, 0.0);
				C1->v_H2 = linear_funk(i->y, i->v_H2, j->y, j->v_H2, 0.0);

				C1->ro_H3 = linear_funk(i->y, i->ro_H3, j->y, j->ro_H3, 0.0);
				C1->p_H3 = linear_funk(i->y, i->p_H3, j->y, j->p_H3, 0.0);
				C1->u_H3 = linear_funk(i->y, i->u_H3, j->y, j->u_H3, 0.0);
				C1->v_H3 = linear_funk(i->y, i->v_H3, j->y, j->v_H3, 0.0);

				C1->ro_H4 = linear_funk(i->y, i->ro_H4, j->y, j->ro_H4, 0.0);
				C1->p_H4 = linear_funk(i->y, i->p_H4, j->y, j->p_H4, 0.0);
				C1->u_H4 = linear_funk(i->y, i->u_H4, j->y, j->u_H4, 0.0);
				C1->v_H4 = linear_funk(i->y, i->v_H4, j->y, j->v_H4, 0.0);

				C1->F_n = linear_funk(i->y, i->F_n, j->y, j->F_n, 0.0);
				C1->F_u = linear_funk(i->y, i->F_u, j->y, j->F_u, 0.0);
				C1->F_v = linear_funk(i->y, i->F_v, j->y, j->F_v, 0.0);
				C1->F_T = linear_funk(i->y, i->F_T, j->y, j->F_T, 0.0);

				C1->I_u = linear_funk(i->y, i->I_u, j->y, j->I_u, 0.0);
				C1->I_v = linear_funk(i->y, i->I_v, j->y, j->I_v, 0.0);
				C1->I_T = linear_funk(i->y, i->I_T, j->y, j->I_T, 0.0);

				res.push_back(C1);
				break;
			}
		}

	}

	return res;
}

void Konstruktor::Velosity_initial(Sensor* s, double& Vx, double& Vy, double& Vz)
{
	double ksi1 = s->MakeRandom();
	double ksi2 = s->MakeRandom();
	double a = sqrt(-log(1.0 - ksi2));
	Vy = a * cos(2.0 * pi * ksi1);
	Vz = a * sin(2.0 * pi * ksi1);
	//cout << Vy << endl;
	double ksi3, ksi4, ksi5, ksi6;
	double z = 0;
	double p1 = fabs(Velosity_inf) * sqrtpi / (1.0 + fabs(Velosity_inf) * sqrtpi);

	do
	{
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();

		if (p1 > ksi3)
		{
			z = cos(pi * ksi5) * sqrt(-log(ksi4));
		}
		else
		{
			if (ksi4 <= 0.5)
			{
				z = -sqrt(-log(2.0 * ksi4));
			}
			else
			{
				z = sqrt(-log(2.0 * (1.0 - ksi4)));
			}
		}
	} while (fabs(z + Velosity_inf) / (fabs(Velosity_inf) + fabs(z)) <= ksi6 || z > -Velosity_inf);

	Vx = z + Velosity_inf;
	return;
}

void Konstruktor::Velosity_initial2(Sensor* s, double& Vx, double& Vy, double& Vz)
{
	double ksi1 = s->MakeRandom();
	double ksi2 = s->MakeRandom();
	double a = sqrt(-log(1.0 - ksi2));
	Vy = a * cos(2.0 * pi * ksi1);
	Vz = a * sin(2.0 * pi * ksi1);
	//cout << Vy << endl;
	double ksi3, ksi4, ksi5, ksi6;
	double z = 0;
	double p1 = 0.5 * fabs(Velosity_inf) * sqrtpi / (0.5 + 0.5 * fabs(Velosity_inf) * sqrtpi);

	do
	{
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();

		if (p1 > ksi3)
		{
			z = cos(pi * ksi5) * sqrt(-log(ksi4));
		}
		else
		{
			z = sqrt(-log(1.0 - ksi4));

		}
	} while (fabs(z + Velosity_inf) / (fabs(Velosity_inf) + fabs(z)) <= ksi6 || z < -Velosity_inf);

	Vx = z + Velosity_inf;
	if (Vx <= 0)
	{
		cout << "dfEEERR 32424442" << endl;
	}
	return;
}

void Konstruktor::Change_Velosity(Sensor* s, const double& Ur, const double& Uthe, const double& Uphi, //
	const double& Vr, const double& Vthe, const double& Vphi, double& X, double& Y, double& Z, const double& cp)
{
	double x = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double p4 = 0.5 * sqrtpi * x / (1.0 + 0.5 * sqrtpi * x);
	double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6;
	double om1, om2, om3, lo;
	double y1, y2, y3;
	double v1, v2, v3, u1, u2, u3;
	double uu, yy, vv, D, ko;
	do
	{
		ksi1 = s->MakeRandom();
		ksi2 = s->MakeRandom();
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();
		//cout << "sd " << endl;
		if (p4 < ksi1)
		{
			om1 = 1.0 - 2.0 * ksi4;
			om2 = sqrt(1.0 - kv(om1)) * cos(2.0 * pi * ksi5);
			om3 = sqrt(1.0 - kv(om1)) * sin(2.0 * pi * ksi5);
			// Более экономичный алгоритм
			/*do
			{
				om2 = 1.0 - 2.0 * s->MakeRandom();
				om3 = 1.0 - 2.0 * s->MakeRandom();
				D = kv(om2) + kv(om3);
			} while (D > 1);
			ko = sqrt((1.0 - kv(om1)) / D);
			om2 = om2 * ko;
			om3 = om3 * ko;*/

			lo = sqrt(-log(ksi2 * ksi3));
			y1 = lo * om1;
			y2 = lo * om2;
			y3 = lo * om3;
		}
		else
		{
			y1 = sqrt(-log(ksi2)) * cos(pi * ksi3);
			y2 = sqrt(-log(ksi4)) * cos(2.0 * pi * ksi5);
			y3 = sqrt(-log(ksi4)) * sin(2.0 * pi * ksi5);
		}
		v1 = y1 + Ur;
		v2 = y2 + Uthe;
		v3 = y3 + Uphi;
		u1 = Vr - v1;
		u2 = Vthe - v2;
		u3 = Vphi - v3;
		uu = sqrt(kvv(u1, u2, u3));
		yy = sqrt(kvv(y1, y2, y3));
	} while (((uu * sigma2(uu, cp)) / (sigma2(x, cp) * (x + yy))) <= ksi6);
	//cout << v2 << endl;
	X = v1;
	Y = v2;
	Z = v3;
}

void Konstruktor::M_K_training(void)
{
	// Блок загрузки датчиков случайных чисел
	ifstream fin2;
	fin2.open("rnd_Dima.dat");
	double d, a1, b1, c;
	for (int i = 0; i < 270; i++)
	{
		fin2 >> d >> a1 >> b1 >> c;
		auto s = new Sensor(a1, b1, c);
		this->Sensors.push_back(s);
	}
	fin2.close();


	this->sqv_1 = (2.54189 * pi * (kv(this->y_max) - kv(350.0)));
	this->sqv_2 = (y_max * sqrtpi * (this->x_max - this->x_min));
	this->sqv_3 = (0.0000282543 * pi * kv(this->y_max));
	this->sqv_4 = (2.54189 * pi * kv(350.0));
	this->sum_s = this->sqv_1 + this->sqv_2 + this->sqv_3 + this->sqv_4;
	this->Number1 = 2025000 * 20;
	this->Number2 = (647838 * 3);
	this->Number3 = (16200 * 3);
	this->Number4 = (2025000 * 400);
	this->AllNumber = ((this->Number1) + (this->Number2) + (this->Number3) + (this->Number4));

	cout << "ALL NUMBER = " << this->AllNumber << endl;

	for (auto& i : this->all_Kyb)
	{
		for (auto& j : i->sosed)
		{
			if (j->number == -1)
			{
				this->Boandary_1.push_back(i);
			}
			if (j->number == -2)
			{
				this->Boandary_2.push_back(i);
			}
			if (j->number == -3)
			{
				this->Boandary_3.push_back(i);
			}
		}

		i->Setup_boandary(this->DX, this->DY);
	}
}

void Konstruktor::M_K(void)
{
	mutex mut_1;

	for (auto& i : this->all_Kyb)
	{
		i->F_n = 0.0;
		i->F_u = 0.0;
		i->F_v = 0.0;
		i->F_T = 0.0;

		i->I_u = 0.0;
		i->I_v = 0.0;
		i->I_T = 0.0;
	}


#pragma omp parallel for
	for (int index = 0; index < 270; index++)
	{
		bool info = false;
		if (index == 0)
		{
			info = true;
		}
		double ksi1, ksi2, x_0, y_0, z_0, r_0, ksi3, ksi4, ksi5, ksi6, ksi7, phi, Vr, Vphi, Vx;
		bool t = false;
		int n, m;
		double k;
		double mu1, mu2, mu3, mu4;
		mu1 = ((this->sqv_1) / this->sum_s) * (1.0 * this->AllNumber / this->Number1);
		mu2 = ((this->sqv_2) / this->sum_s) * (1.0 * this->AllNumber / this->Number2);
		mu3 = ((this->sqv_3) / this->sum_s) * (1.0 * this->AllNumber / this->Number3);
		mu4 = ((this->sqv_4) / this->sum_s) * (1.0 * this->AllNumber / this->Number4);
		Sensor* sens = Sensors[index];
		mut_1.lock();
		cout << index << " potok  is  270" << endl;
		mut_1.unlock();

		for (int ii = 0; ii < Number1 / 270; ii++)  //
		{
			//cout << ii << endl;
			t = false;
			double a, b, c;
			Velosity_initial(sens, a, b, c);
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();

			r_0 = sqrt(kv(350.0) + ksi1 * (kv(y_max) - kv(350.0)));   // Пересмотреть функцию
			phi = ksi2 * 2.0 * pi;
			y_0 = r_0 * cos(phi);
			z_0 = r_0 * sin(phi);

			//cout << a << endl;


			Kyb* Point = Belong_point(1, this->x_max - geo, r_0);  // Находит ячейку, которой принадлежит точка
			// РАБОТАЕТ ТОЛЬКО ЕСЛИ РАЗМЕРЫ ГРАНИЧНЫХ ЯЧЕЕК - ИЗНАЧАЛЬНЫЕ DX DY (иначе нужно поменять функцию)


			Fly_exchenge(sens, x_max - geo, y_0, z_0, a, b, c, Point, mu1, ii);
		}
		for (int ii = 0; ii < Number2 / 270; ii++)  // С боковой поверхности
		{
			t = false;
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();
			ksi3 = sens->MakeRandom();
			ksi4 = sens->MakeRandom();
			ksi5 = sens->MakeRandom();
			ksi6 = sens->MakeRandom();
			ksi7 = sens->MakeRandom();

			x_0 = (x_min + geo) + ksi1 * (x_max - 2.0 * geo - x_min);
			phi = ksi2 * 2.0 * pi;
			Vphi = cos(2.0 * pi * ksi3) * sqrt(-log(1.0 - ksi4));
			Vx = Velosity_inf + sin(2.0 * pi * ksi5) * sqrt(-log(1.0 - ksi6));
			Vr = -sqrt(-log(ksi7));
			y_0 = (y_max - geo) * cos(phi);
			z_0 = (y_max - geo) * sin(phi);


			Kyb* Point = Belong_point(3, x_0, y_max - geo);

			Fly_exchenge(sens, x_0, y_0, z_0, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi, Point, mu2, ii);
		}
		for (int ii = 0; ii < Number3 / 270; ii++)
		{
			//cout << ii << endl;
			t = false;
			double a, b, c;
			Velosity_initial2(sens, a, b, c);
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();

			r_0 = sqrt(ksi1 * y_max * y_max);
			phi = ksi2 * 2.0 * pi;
			y_0 = r_0 * cos(phi);
			z_0 = r_0 * sin(phi);




			Kyb* Point = Belong_point(2, x_min + geo, r_0);

			Fly_exchenge(sens, x_min + geo, y_0, z_0, a, b, c, Point, mu3, ii);
		}
		for (int ii = 0; ii < Number4 / 270; ii++)
		{
			//cout << ii << endl;
			t = false;
			double a, b, c;
			Velosity_initial(sens, a, b, c);
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();

			r_0 = sqrt(ksi1 * kv(350.0));
			phi = ksi2 * 2.0 * pi;
			y_0 = r_0 * cos(phi);
			z_0 = r_0 * sin(phi);

			//cout << a << endl;


			Kyb* Point = Belong_point(1, x_max - geo, r_0);

			Fly_exchenge(sens, x_max - geo, y_0, z_0, a, b, c, Point, mu4, ii);
		}
	}

	for (auto& k: this->all_Kyb)
	{
		double dx = (DX / pow(2, k->size - 1)) / 2.0;   // Половина длины ячейки
		double dy = (DY / pow(2, k->size - 1)) / 2.0;   // Половина ширины ячейки

		double no = (1.0 * AllNumber * (pi * kv(k->y + dy) * (2.0 * dx) - pi * kv(k->y - dy) * (2.0 * dx)));

		if (k->F_n > 0)
		{
			k->I_u = k->I_u / k->F_n;
			k->I_v = k->I_v / k->F_n;
			k->I_T = k->I_T / k->F_n;

			k->F_u = k->F_u / k->F_n;
			k->F_v = k->F_v / k->F_n;
			k->F_T = (2.0 / 3.0) * (k->F_T / k->F_n - kvv(k->F_u, k->F_v, 0.0));
		}
		else
		{
			k->I_u = 0.0;
			k->I_v = 0.0;
			k->I_T = 0.0;

			k->F_n = 0.0;
			k->F_u = 0.0;
			k->F_v = 0.0;
			k->F_T = 0.0;
		}

		k->F_n = sum_s * k->F_n / no;

		k->I_u = -k->I_u;
		k->I_v = -k->I_v;
		k->I_T = -(1.0 / 2.0) * k->I_T;
	}
}

void Konstruktor::Fly_exchenge(Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Kyb* ind, const double& mu, int ii)
{
	Kyb* next = nullptr;
	Kyb* prev = nullptr;
	Kyb* head = ind;
	double X = x_0, Y = y_0, Z = z_0;
	double KSI = -log(1.0 - sens->MakeRandom());
	double I_do = 0.0;
	int per = 0;

	double Ur, Uphi, Utheta;
	double Vr, Vphi, Vtheta;
	double uu, vv, ww;

	do
	{	
		/*if (ii == 0)
		{
			cout << X << " " << sqrt(kv(Y) + kv(Z)) << " " << 3 << endl;
			cout << head->x << " " << head->y << " " << 4 << endl;
		}*/
		if (Flying_exchange(KSI, Vx, Vy, Vz, X, Y, Z, next, head, prev, mu, I_do, ii) == false)
		{
			break;
		}
		if (KSI < 0.0)
		{
			per = 1;
			KSI = -log(1.0 - sens->MakeRandom());
			I_do = 0.0;
			prev = head;
			double sk = sqrt(head->p / head->ro);
			double alpha = this->polar_angle(Y, Z);
			double uuu, vvv;
			uuu = head->u;
			vvv = head->v;

			this->Change_Velosity(sens, uuu / sk, vvv * cos(alpha) / sk, vvv * sin(alpha) / sk,//
				Vx / sk, Vy / sk, Vz / sk, uu, vv, ww, sk);  // здесь подаются  r, theta, phi
			uu = uu * sk;
			vv = vv * sk;
			ww = ww * sk;

			//np1[head] += mu * (Vx - uu);
			//np1[head] += mu * (kvv(Vx, Vy, Vz) - kvv(uu, vv, ww));

			Vx = uu;
			Vy = vv;
			Vz = ww;
		}
		else
		{
			per = 0;
			prev = head;
			head = next;
		}
	} while (true);
	return;
}

bool Konstruktor::Flying_exchange(double& KSI, double& Vx, double& Vy, double& Vz, double& X, double& Y,//
	double& Z, Kyb*& next, Kyb* head, Kyb* prev, const double& mu, double& I_do, int ii)
	// Vx, Vy, Vz - скорость атома водорода
	// X,Y,Z - координата атома
	// head - текущая ячейки
{
	double dx = (this->DX / pow(2, head->size - 1)) / 2.0;   // Половина длины ячейки
	double dy = (this->DY / pow(2, head->size - 1)) / 2.0;   // Половина ширины ячейки
	double y0 = head->y;
	double x0 = head->x;
	double x, uz, uz_M, uz_E, t1, t2;// y, z, r;
	int mode = 0;
	double time = 1000000000;
	int step = 0;
	double cp = sqrt(head->p / head->ro);
	bool chenge_was = false;
	double u1, u2, u3;
	double al, be;
	double vx = head->u;
	double vy = head->v;
	double ro = head->ro;


	while(head->Belong_fast(X, sqrt(kv(Y) + kv(Z)), dx, dy) == false)
	{
		double alpha = polar_angle(Y, Z);
		double yy, zz;
		yy = y0 * cos(alpha);
		zz = y0 * sin(alpha);
		double nn = sqrt(kv(X - x0) + kv(Y - yy) + kv(Z - zz));
		double geot = dx / 9000.0;
		X = X - geot * (X - x0) / nn;
		Y = Y - geot * (Y - yy) / nn;
		Z = Z - geot * (Z - zz) / nn;
	}

	while (Peresechenie(x0, y0, dx, dy, X, Y, Z, Vx, Vy, Vz, mode, time) == false)
	{
		step++;
		if (step > 6)
		{
			cout << "Error  1605323534526853234245436" << endl;
			return false;
		}
		double alpha = polar_angle(Y, Z);
		double yy, zz;
		yy = y0 * cos(alpha);
		zz = y0 * sin(alpha);
		double nn = sqrt(kv(X - x0) + kv(Y - yy) + kv(Z - zz));
		double geot = dx / 100.0;
		X = X - geot * (X - x0) / nn;
		Y = Y - geot * (Y - yy) / nn;
		Z = Z - geot * (Z - zz) / nn;
	}

	//cout << mode << endl;

	double l = sqrt(kvv(time * Vx, time * Vy, time * Vz));
	//double alpha = polar_angle(Y + 0.5 * time * Vy, Z + 0.5 * time * Vz);
	double alpha = polar_angle(Y, Z);
	double u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
	uz = Velosity_1(u, cp);
	double sig = Kn_ * sqrt(kvv(Vx, Vy, Vz)) / (ro * uz * sigma(uz));
	double I = I_do + l / sig;

	if (I < KSI)  // Не произошла перезарядка
	{
		I_do = I;
		alpha = polar_angle(Y + 0.5 * time * Vy, Z + 0.5 * time * Vz);
	}
	else  // Перезарядка была
	{
		chenge_was = true;
		double ksi = (KSI - I_do) * sig;
		time = ksi / sqrt(kvv(Vx, Vy, Vz));
		next = head;
		I_do = 0.0;
		KSI = -1.0;
		alpha = polar_angle(Y + 0.5 * time * Vy, Z + 0.5 * time * Vz);  // Гарантирует расчёт угла в середине пути
	}


	double dalpha = fabs(polar_angle(Y + time * Vy, Z + time * Vz) - polar_angle(Y, Z));
	if (fabs(vy) > 0.01 * sqrt(kv(vx) + kv(vy)) && dalpha > pi / 90.0)
	{
		double tt;
		int n = ((int)(dalpha / (pi / 90.0)) + 1);
		double dt = time / (1.0 * n);
		for (int i = 0; i < n; i++)
		{
			tt = (i + 0.5) * dt;
			alpha = polar_angle(Y + tt * Vy, Z + tt * Vz);
			u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
			uz = Velosity_1(u, cp);
			uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi);
			uz_E = Velosity_3(u, cp);
			u1 = vx - Vx;
			u2 = vy * cos(alpha) - Vy;
			u3 = vy * sin(alpha) - Vz;
			double skalar = Vx * u1 + Vy * u2 + Vz * u3;


			head->mut.lock();
			head->F_n += dt * mu;
			head->F_u += dt * Vx * mu;
			head->F_v += dt * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			head->F_T += dt * kvv(Vx, Vy, Vz) * mu;


			head->I_u += mu * dt * uz_M * uz * sigma(uz_M) * u1 / u;
			head->I_v += mu * dt * uz_M * uz * sigma(uz_M) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;

			head->I_T += mu * dt * (0.5 * (3.0 * kv(cp) + 2.0 * kv(u)) * uz_E * sigma(uz_E) + 2.0 * uz_M * uz * sigma(uz_M) * skalar / u);
			head->mut.unlock();
		}
	}
	else
	{
		//Суммируем значения в источники
		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		uz = Velosity_1(u, cp);
		uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi);
		uz_E = Velosity_3(u, cp);
		u1 = vx - Vx;
		u2 = vy * cos(alpha) - Vy;
		u3 = vy * sin(alpha) - Vz;
		double skalar = Vx * u1 + Vy * u2 + Vz * u3;


		head->mut.lock();
		head->F_n += time * mu;
		head->F_u += time * Vx * mu;
		head->F_v += time * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		head->F_T += time * kvv(Vx, Vy, Vz) * mu;


		head->I_u += mu * time * uz_M * uz * sigma(uz_M) * u1 / u;
		head->I_v += mu * time * uz_M * uz * sigma(uz_M) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;

		head->I_T += mu * time * (0.5 * (3.0 * kv(cp) + 2.0 * kv(u)) * uz_E * sigma(uz_E) + 2.0 * uz_M * uz * sigma(uz_M) * skalar / u);
		head->mut.unlock();
	}



	// Передвигаем координату
	X = X + time * Vx;
	Y = Y + time * Vy;
	Z = Z + time * Vz;


	// Меняем ячейку правильно
	if (chenge_was == false)
	{
		if (mode == 1)
		{
			X = X + geo;
			if (head->boandary_1.size() == 1)
			{
				if (head->boandary_1[0]->number < 0)
				{
					return false;
				}
				else
				{
					next = head->boandary_1[0];
					return true;
				}
			}
			else
			{
				for (auto& i : head->boandary_1)
				{
					if (i->Belong_slow(X + geo, sqrt(kv(Y) + kv(Z)), this->DX, this->DY))
					{
						next = i;
						return true;
					}
				}
				cout << "Ne nashol soseda   123143254rf3dwedwwefw" << endl;
				/*for (auto& i : head->boandary_1)
				{
					cout << i->x << " " << i->y << " " << 0 << endl;
				}
				cout << X << " " << sqrt(kv(Y) + kv(Z)) << " " << 1 << endl;
				cout << head->x << " " << head->y << " " << 2 << endl;
				cout << head->boandary_1.size() << endl;
				cout << ii << endl;
				exit(-1);*/
				return false;
			}
		}
		if (mode == 2)
		{
			X = X - geo;
			//cout << "Size = " << head->boandary_2.size() << endl;
			if (head->boandary_2.size() == 1)
			{
				if (head->boandary_2[0]->number < 0)
				{
					return false;
				}
				else
				{
					next = head->boandary_2[0];
					return true;
				}
			}
			else
			{
				for (auto& i : head->boandary_2)
				{
					if (i->Belong_slow(X - geo, sqrt(kv(Y) + kv(Z)), this->DX, this->DY))
					{
						next = i;
						return true;
					}
				}

				cout << "Ne nashol soseda   wceferver34r3x4rc343r4" << endl;
				/*for (auto& i : head->boandary_2)
				{
					cout << i->x << " " << i->y << " " << 0 << endl;
				}
				cout << X << " " << sqrt(kv(Y) + kv(Z)) << " " << 1 << endl;
				cout << head->x << " " << head->y << " " << 2 << endl;
				cout << head->boandary_2.size() << endl;
				cout << ii << endl;
				exit(-1);*/
				return false;
			}
		}
		if (mode == 3)
		{
			double rr = sqrt(kv(Y) + kv(Z));
			Y = Y + geo * Y / rr;
			Z = Z + geo * Z / rr;
			if (head->boandary_3.size() == 1)
			{
				if (head->boandary_3[0]->number < 0)
				{
					return false;
				}
				else
				{
					next = head->boandary_3[0];
					return true;
				}
			}
			else
			{
				for (auto& i : head->boandary_3)
				{
					if (i->Belong_slow(X, sqrt(kv(Y) + kv(Z)) + geo, this->DX, this->DY))
					{
						next = i;
						return true;
					}
				}
				cout << "Ne nashol soseda   ukuijhgfer435353" << endl;
				/*for (auto& i : head->boandary_3)
				{
					cout << i->x << " " << i->y << " " << 0  << endl;
				}
				cout << X << " " << sqrt(kv(Y) + kv(Z)) << " " << 1 << endl;
				cout << head->x << " " << head->y << " " << 2 << endl;
				cout << head->boandary_3.size() << endl;
				cout << ii << endl;
				exit(-1);*/
				return false;
			}
		}
		if (mode == 4)
		{
			double rr = sqrt(kv(Y) + kv(Z));
			Y = Y - geo * Y / rr;
			Z = Z - geo * Z / rr;
			if (head->boandary_4.size() == 1)
			{
				if (head->boandary_4[0]->number < 0)
				{
					return false;
				}
				else
				{
					next = head->boandary_4[0];
					return true;
				}
			}
			else
			{
				for (auto& i : head->boandary_4)
				{
					if (i->Belong_slow(X, sqrt(kv(Y) + kv(Z)) - geo, this->DX, this->DY))
					{
						next = i;
						return true;
					}
				}
				cout << "Ne nashol soseda   wfcrfewv232345y6786v5vt5" << endl;
				/*for (auto& i : head->boandary_4)
				{
					cout << i->x << " " << i->y << " " << 0 << endl;
				}
				cout << X << " " << sqrt(kv(Y) + kv(Z)) << " " << 1 << endl;
				cout << head->x << " " << head->y << " " << 2 << endl;
				cout << head->boandary_4.size() << endl;
				cout << ii << endl;
				exit(-1);*/
				return false;
			}
		}
	}
	else
	{
		return true;
	}
	return false;
}

bool Konstruktor::Peresechenie(const double& x0, const double& y0, const double& dx, const double& dy, const double& x, const double& y, const double& z, //
	const double& Vx, const double& Vy, const double& Vz, int& mode, double& t)
{
	mode = 0;
	double t1 = -1.0, t2 = -1.0, t3 = -1.0, t4 = -1.0, t5 = -1.0, t6 = -1.0;
	if (fabs(Vx) > 0.0000001)
	{
		t1 = (x0 + dx - x) / Vx;
		t2 = (x0 - dx - x) / Vx;
	}
	if (kv(Vy) + kv(Vz) > 0.0000001)
	{
		this->peresich(y, z, Vy, Vz, y0 + dy, t3, t4);
		this->peresich(y, z, Vy, Vz, y0 - dy, t5, t6);
		t3 = minplus(t3, t4);
		t4 = minplus(t5, t6);
	}
	t = 100000000.0;
	if (t > t1 && t1 > 0.000000001)
	{
		t = t1;
		mode = 1;
	}
	if (t > t2 && t2 > 0.000000001)
	{
		t = t2;
		mode = 2;
	}
	if (t > t3 && t3 > 0.000000001)
	{
		t = t3;
		mode = 3;
	}
	if (t > t4 && t4 > 0.000000001)
	{
		t = t4;
		mode = 4;
	}

	if (fabs(Vx) > 0.0000001 && t1 < 0.000000001 && t2 < 0.000000001)
	{
		return false;
	}

	if (mode == 0)
	{
		return false;
	}
	return true;
}

void Konstruktor::peresich(const double& y, const double& z, const double& Vy, const double& Vz, const double& R, double& t1, double& t2)
{
	double b = (2.0 * y * Vy + 2.0 * z * Vz);
	double a = (kv(Vy) + kv(Vz));
	double D = b * b - 4.0 * (kv(y) + kv(z) - kv(R)) * a;
	if (D < 0)
	{
		t1 = -1.0;
		t2 = -1.0;
		return;
	}
	D = sqrt(D);
	t1 = (-b + D) / (2.0 * a);
	t2 = (-b - D) / (2.0 * a);
}


Kyb* Konstruktor::Belong_point(int b, const double& x, const double& y)
{
	if (b == 1)
	{
		for (auto& i : this->Boandary_1)
		{
			if (i->Belong_slow(x, y, this->DX, this->DY))
			{
				return i;
			}
		}
	}
	else if (b == 2)
	{
		for (auto& i : this->Boandary_2)
		{
			if (i->Belong_slow(x, y, this->DX, this->DY))
			{
				return i;
			}
		}
	}
	else if (b == 3)
	{
		for (auto& i : this->Boandary_3)
		{
			if (i->Belong_slow(x, y, this->DX, this->DY))
			{
				return i;
			}
		}
	}

	cout << "ERRORRRORORfireubfvwcefrvjywgkvygdcwkug324h324334" << endl;
	return nullptr;
}

double Konstruktor::minplus(const double& x, const double& y)
{
	if (x < 0.00000000001 && y < 0.00000000001)
	{
		return -1.0;
	}
	else if (x < 0.00000000001)
	{
		return y;
	}
	else if (y < 0.00000000001)
	{
		return x;
	}
	else if (x < y)
	{
		return x;
	}
	else
	{
		return y;
	}
}

double Konstruktor::polar_angle(const double& x, const double& y)
{
	if (x < 0)
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

double Konstruktor::Velosity_1(const double& u, const double& cp)
{
	if (u < 0.001)
	{
		return 2.0 * cp / sqrtpi + 2.0 * u * u / (3.0 * cp * sqrtpi) - u * u * u * u / (15.0 * cp * cp * cp * sqrtpi);
	}
	else
	{
		return  exp(-u * u / kv(cp)) * cp / sqrtpi + (u + kv(cp) / (2.0 * u)) * erf(u / cp);
	}
}

double Konstruktor::Velosity_2(const double& u, const double& cp)  // Считает на совсем скорость, а только её числитель (см. статью)
{
	if (u < 0.001)
	{
		return (8.0 / 3.0) * kv(cp) * kv(cp) * pi * u + (8.0 / 15.0) * kv(cp) * pi * u * u * u - (4.0 / 105.0) * pi * kv(u) * kv(u) * u;
	}
	else
	{
		return  cp * cp * cp * pi * (exp(-u * u / kv(cp)) * cp * u * 2.0 * (kv(cp) + 2.0 * kv(u)) +//
			sqrtpi * (4.0 * kv(u) * kv(u) + 4.0 * cp * cp * kv(u) - kv(cp) * kv(cp)) * erf(u / cp)) / (4.0 * u * u);
	}
}

double Konstruktor::Velosity_3(const double& u, const double& cp)
{
	if (u < 0.001)
	{
		return 8.0 * cp / (3.0 * sqrtpi) + 8.0 * u * u / (9.0 * cp * sqrtpi) - 44.0 * u * u * u * u / (135.0 * cp * cp * cp * sqrtpi);
	}
	else
	{
		return  exp(-u * u / kv(cp)) * cp * (5.0 * kv(cp) + 2.0 * kv(u)) / (sqrtpi * (3.0 * kv(cp) + 2.0 * kv(u))) +//
			(4.0 * kv(u) * kv(u) + 12.0 * cp * cp * kv(u) + 3.0 * kv(cp) * kv(cp)) * erf(u / cp) / (2.0 * u * (3.0 * kv(cp) + 2.0 * kv(u)));
	}
}