#include "Konstruktor.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>


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
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << 4 * ll;
	fout << " , E= " << 4 * ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->all_Kyb)
	{
		double dx_1 = (this->DX / pow(2, i->size - 1)) / 2.0;   // Половина длины ячейки
		double dy_1 = (this->DY / pow(2, i->size - 1)) / 2.0;   // Половина ширины ячейки
		fout << i->x - dx_1 << " " << i->y - dy_1 << endl;
		fout << i->x - dx_1 << " " << i->y + dy_1 << endl;
		fout << i->x + dx_1 << " " << i->y + dy_1 << endl;
		fout << i->x + dx_1 << " " << i->y - dy_1 << endl;
	}
	for (int i = 0; i < ll; i++)
	{
		fout << 4 * i + 1 << " " << 4 * i + 2 << endl;
		fout << 4 * i + 2 << " " << 4 * i + 3 << endl;
		fout << 4 * i + 3 << " " << 4 * i + 4 << endl;
		fout << 4 * i + 1 << " " << 4 * i + 4 << endl;
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


	ofstream fout;
	string name_f = "2D_tecplot.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Ro\", \"P\", \"Vx\", \"Vy\", \"Max\",\"Q\", ZONE T = \"HP\"" << endl;
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
				<< i->u * u_o << " " << i->v * u_o << " " << Max << " "  << QQ << endl;
		
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
			if (ll % 25000 == 0)
			{
				cout << ll << endl;
			}
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
		if ((i->x > x1) && (i->x < x2) && (i->y < r1) && (r > 0.95 * Distant))
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
			if (ll % 25000 == 0)
			{
				cout << ll << endl;
			}
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
		double dist = sqrt(i->x * i->x + i->y * i->y);
		double r_0 = 1.0;
		double ro = (389.988 * 389.988) / (chi_ * chi_);
		double P_E = ro * chi_ * chi_ / (ggg * 0.25 * 0.25);
		if (dist > Distant * 1.5)
		{
			i->ro = 1.0;
			i->p = 1.0;
			i->u = Velosity_inf;
			i->v = 0.0;
			i->Q = 100.0;
		}
		else
		{
			i->ro = ro / (dist * dist);
			i->p = P_E * pow(r_0 / dist, 2.0 * ggg);
			i->u = chi_ * i->x / dist;
			i->v = chi_ * i->y / dist;
			i->Q = ro * r_0 * r_0 / (dist * dist);
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