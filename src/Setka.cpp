#include "Setka.h"
#include "voro++.hh"
#include <vector>
#include <string>
#include <mutex>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <omp.h>
#include <iomanip>      // std::setprecision

using namespace std;
using namespace voro;

Setka::Setka(string name)
{
    cout << "Start  read  " << name << endl;
    this->initialization();
    this->Download_setka(name);
    this->renumber();

    this->Initialization_do_MHD();
    //this->Calc_normal();
    this->Reconstruct_medium2();
    cout << "End  read  " << name << endl;
}

Setka::Setka()
{
	this->initialization();

	// Управляющие переменные
	int sl_ = 15;                   // сколько слоёв точек, распределённые по радиусу    15
	int n_x = 33;                   // 33
	int n_y = 33;
	int n_z = 33;

	double dx = (x_max - x_min) / (n_x + 1);
	double dy = (y_max - y_min) / (n_y + 1);
	double dz = (z_max - z_min) / (n_z + 1);


	All_Cell.reserve(1000000);

	// Сначала загружаем ячейки, сферически распределённые 
	ifstream fout;
	fout.open("1000_sphere.txt");

	double x, y, z;
	double r, R;
	Cell* A = nullptr;

	for (int i = 0; i < 1000; i++)
	{
		fout >> x >> y >> z;
		for (int j = 1; j <= sl_; j++)                                 // Сколько слоёв?
		{
			r = sqrt(x * x + y * y + z * z);
			R = sl_L + j * (sl_R - sl_L)/(sl_ + 1);
			A = new Cell(x * R / r, y * R / r, z * R / r);
			A->type = C_base;
			this->All_Cell.push_back(A);
		}
		
	}

    double dh = (sl_R - sl_L) / (sl_ + 1) / 2.0;

	for (int i = 1; i <= n_x; i++)
	{
		for (int j = 1; j <= n_y; j++)
		{
			for (int k = 1; k <= n_z; k++)
			{
				x = this->x_min + i * (this->x_max - this->x_min) / (n_x + 1);
				y = this->y_min + j * (this->y_max - this->y_min) / (n_y + 1);
				z = this->z_min + k * (this->z_max - this->z_min) / (n_z + 1);

                int dr = 1;
                double xx, yy, zz;


                if (fabs(x) < 8.0 && fabs(y) < 6.0 && fabs(z) < 6.0)
                {
                    dr = 2;
                }
                if (fabs(x) < 5.0 && fabs(y) < 5.0 && fabs(z) < 5.0)
                {
                    dr = 4;
                }
                if (x > -4 && x < 3.0 && fabs(y) < 4.0 && fabs(z) < 4.0)
                {
                    dr = 6;
                }

                for (int ik = 1; ik <= dr; ik++)
                {
                    for (int jk = 1; jk <= dr; jk++)
                    {
                        for (int kk = 1; kk <= dr; kk++)
                        {
                            xx = x - dx / 2.0 + ik * dx / dr - dx / (2.0 * dr);
                            yy = y - dy / 2.0 + jk * dy / dr - dy / (2.0 * dr);
                            zz = z - dz / 2.0 + kk * dz / dr - dz / (2.0 * dr);

                            if (xx * xx + yy * yy + zz * zz > kv(sl_R + dh))
                            {
                                A = new Cell(xx, yy, zz);
                                A->type = C_base;
                                this->All_Cell.push_back(A);
                            }
                        }
                    }
                }
					
				
			}
		}
	}


	// Нумерация ячеек
	int ik = 0;
	for (auto& i : this->All_Cell)
	{
		i->number = ik;
		ik++;
	}

	cout << "All_points = " << ik - 1 << endl;


}

void Setka::renumber(void)
{
    cout << "renumber" << endl;

    int ik = 0;
    for (auto& i : this->All_Cell)
    {
        i->number = ik;
        ik++;
    }

    ik = 0;
    for (auto& i : this->All_Couple)
    {
        i->number = ik;
        ik++;
    }
}

void Setka::initialization(void)
{
	auto A = new Cell(0.0, 0.0, 0.0);
	A->type = C_wall_x_max;
	A->number = -1;
	this->Wall1 = A;

	A = new Cell(0.0, 0.0, 0.0);
	A->type = C_wall_x_min;
	A->number = -2;
	this->Wall2 = A;

	A = new Cell(0.0, 0.0, 0.0);
	A->type = C_wall_y_max;
	A->number = -3;
	this->Wall3 = A;

	A = new Cell(0.0, 0.0, 0.0);
	A->type = C_wall_y_min;
	A->number = -4;
	this->Wall4 = A;

	A = new Cell(0.0, 0.0, 0.0);
	A->type = C_wall_z_max;
	A->number = -5;
	this->Wall5 = A;

	A = new Cell(0.0, 0.0, 0.0);
	A->type = C_wall_z_min;
	A->number = -6;
	this->Wall6 = A;

	A = new Cell(0.0, 0.0, 0.0);
	A->type = C_sphere;
	A->number = -7;
	this->Wall7 = A;
}

void Setka::Print_points(void)
{
	cout << "Print_points" << endl;
	ofstream fout;
	string name_f = "all_points.txt";
	fout.open(name_f);

	for (auto& i : this->All_Cell)
	{
		fout << i->Center->x << " " << i->Center->y << " " << i->Center->z << endl;
	}
}

void Setka::Cut_Plane_long_z(void)
{
	double x, y, z, rsq, r;
	vector<int> neigh, f_vert, ney, edge, vershin;
	vector<double> X, Y, Ver, Versh;

	double nx, ny, nz, nn;

	cout << "Cut_Plane_long " << endl;

	ofstream fout;
	string name_f = "setka.txt";
	fout.open(name_f);

	int N = 0;
	int E = 0;

	container con(-8.1, 5.1, -6.1, 6.1, -6.1, 6.1, 26, 26, 26, false, false, false, 40);

	for (auto& i : this->All_Cell)
	{
		con.put(i->number, i->Center->x, i->Center->y, i->Center->z);
	}

	voronoicell_neighbor c;
	c_loop_all clo(con);
	if (clo.start()) do if (con.compute_cell(c, clo)) {
		clo.pos(x, y, z);
		//cout << clo.pid() << endl;
		// Весь код далее для печати сетки в файл Текплота. Сложности в том, что нужно по порядку расположить все узлы и грани сетки (Текплот так работает)
		c.nplane(0.0, 0.0, -2.0 * z + 0.00001, -77);

		c.neighbors(ney); // Получаем список номеров соседей по порядку
		int k = -1;  // номер нужного соседа
		// Ищем номер соседа "-7"
		int kk = -1;
		for (int j = 0; j < ney.size(); j++)
		{
			kk++;
			if (ney[j] == -77)
			{
				k = kk;
				break;
			}
		}

		if (k == -1)
		{
			continue;
		}
		//cout << "K = " << k << endl;
		//ищем нужные вершины:
		c.face_vertices(f_vert);
		int m = 0;
		int n = 0;
		int r = 0;
		while (r < k)
		{
			n = f_vert[m];
			m = m + n + 1;
			r++;
		}

		//cout << "m = " << m << endl;

		for (int j = m + 1; j < m + 1 + f_vert[m] - 1; j++)
		{
			vershin.push_back(N + f_vert[j] + 1);  // Вектор номеров нужных вершин
			vershin.push_back(N + f_vert[j + 1] + 1);  // Вектор номеров нужных вершин
		}

		vershin.push_back(N + f_vert[m + 1 + f_vert[m] - 1] + 1);
		vershin.push_back(N + f_vert[m + 1] + 1);

		c.vertices(Ver);

		for (int j = 0; j < Ver.size(); j = j + 3)
		{
			Versh.push_back(Ver[j] + x);  // Вектор всех вершин
			Versh.push_back(Ver[j + 1] + y);  // Вектор всех вершин
			Versh.push_back(Ver[j + 2] + z);  // Вектор всех вершин
		}

		N += Ver.size() / 3;
		E += f_vert[m];


	} while (clo.inc());

	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"Z\"  ZONE T= \"HP\", N=" << N << ", E=" << E << ", F=FEPOINT, ET=LINESEG " << endl;

	int jkl = 0;
	for (int j = 0; j < Versh.size(); j++)
	{
		jkl++;
		fout << Versh[j] << " ";
		if (jkl % 3 == 0)
		{
			fout << endl;
		}
	}

	jkl = 0;
	for (int j = 0; j < vershin.size(); j++)
	{
		jkl++;
		fout << vershin[j] << " ";
		if (jkl % 2 == 0)
		{
			fout << endl;
		}
	}
}

void Setka::Construct_start(void)
{
	cout << "Construct_start" << endl;
	double x, y, z, dd, r, nx, ny, nz;
	vector<int> ney;
	vector<double> Surf;
    vector<double> norm;
	int idd = 0;

	container con(this->x_min - 1.0 * RR_, this->x_max + 1.0 * RR_,   this->y_min - 1.0 * RR_, this->y_max + 1.0 * RR_, //
		this->z_min - 1.0 * RR_, this->z_max + 1.0 * RR_, 36, 36, 36, false, false, false, 40);

	for (auto& i : this->All_Cell)
	{
        if (i->include_ == true)
        {
            i->Grans.clear();
            con.put(i->number, i->Center->x, i->Center->y, i->Center->z);
        }
	}

	voronoicell_neighbor c;
	c_loop_all clo(con);
	if (clo.start()) do if (con.compute_cell(c, clo)) {
		clo.pos(x, y, z);
		idd = clo.pid();
		r = sqrt(kv(x) + kv(y) + kv(z));
		dd = r - (this->sl_L - 0.01 * RR_);
		nx = x / r;
		ny = y / r;
		nz = z / r;
		// Далее вручную разрезаем ячейку особыми гранями (стенками и границами)
		c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
		c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
		c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
		c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
		c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
		c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
		c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, - 2.0 * dd * nz, -7);

		c.neighbors(ney); // Получаем список номеров соседей по порядку
		c.face_areas(Surf);   // Площади граней
        c.normals(norm);

        double VV;
        VV = c.volume();  // Записали объём ячейки
        this->All_Cell[idd]->set_Volume(VV);

		for (int i = 0; i < ney.size(); i++)
		{
			auto A = new Gran();
			A->S = Surf[i];              // Площадь грани
            A->n1 = norm[i * 3];
            A->n2 = norm[i * 3 + 1];
            A->n3 = norm[i * 3 + 2];

			if (ney[i] >= 0)
			{
				A->Sosed = this->All_Cell[ney[i]];
			}
			else if (ney[i] == -1)
			{
				A->Sosed = this->Wall1;
			}
			else if (ney[i] == -2)
			{
				A->Sosed = this->Wall2;
			}
			else if (ney[i] == -3)
			{
				A->Sosed = this->Wall3;
			}
			else if (ney[i] == -4)
			{
				A->Sosed = this->Wall4;
			}
			else if (ney[i] == -5)
			{
				A->Sosed = this->Wall5;
			}
			else if (ney[i] == -1)
			{
				A->Sosed = this->Wall1;
			}
			else if (ney[i] == -6)
			{
				A->Sosed = this->Wall6;
			}
			else if (ney[i] == -7)
			{
				A->Sosed = this->Wall7;
			}
			else
			{
				cout << "ERROR  395  Setka  hfjkhgefyuevy34u54" << endl;
			}

			this->All_Cell[idd]->Grans.push_back(A);
		}

	} while (clo.inc());
}

void Setka::Add_couple_start(void)
{
    cout << "Add_couple_start" << endl;
    int Ni = 150;    // По z - оси
    int Nj = 200;
    double dd = 0.05;
    double r = 2.61;
    double phi = 0.0;
    
    double x, y, z, x2, y2;
    int kl = 0;
    for (int i = 0; i < Ni; i++)
    {
        for (int j = 0; j < Nj; j++)
        {
            kl++;
            phi = j * 2.0 * pi / Nj;
            x = r * cos(phi);
            y = r * sin(phi);
            x2 = (r + dd) * cos(phi);
            y2 = (r + dd) * sin(phi);
            z = -5.0 + i * 10.0 / (Ni - 1);

            auto A = new Cell(x, y, z);
            A->type = C_base;
            A->couple_ = true;
            A->include_ = true;
            A->mgd_ = true;
            this->All_Cell.push_back(A);

            auto B = new Cell(x2, y2, z);
            B->type = C_base;
            B->couple_ = true;
            B->include_ = true;
            B->mgd_ = true;
            this->All_Cell.push_back(B);

            auto Par = new Couple(A, B, dd);
            A->Par = Par;
            B->Par = Par;
            Par->Resolve();

            this->All_Couple.push_back(Par);
        }
    }

    cout << "All pars = " << kl << "    " << kl * 2.0 << endl;

    this->renumber();
}

void Setka::Disable_cells(void)
{
    cout << "Disable_cells" << endl;
    double x1, y1, z1, x2, y2, z2;
    for (auto& i : this->All_Couple)
    {
        i->A1->Can_neighbours.clear();
        i->A2->Can_neighbours.clear();
        i->A1->reconstruct_ = true;
        i->A2->reconstruct_ = true;

        i->A1->Center->get(x1, y1, z1);
        for (auto& j : i->A1->Grans)
        {
            if (j->Sosed->couple_ == false && j->Sosed->type == C_base)
            {
                j->Sosed->Center->get(x2, y2, z2);
                if (kv(x2 - x1) + kv(y2 - y1) + kv(z2 - z1) < kv(3.0 * i->dist))
                {
                    if (j->Sosed->include_ == true)
                    {
                        j->Sosed->include_ = false;
                        j->Sosed->reconstruct_ = true;
                        for (auto& k : j->Sosed->Grans)
                        {
                            k->Sosed->reconstruct_ = true;
                        }
                    }
                }
            }
        }

        i->A2->Center->get(x1, y1, z1);
        for (auto& j : i->A2->Grans)
        {
            if (j->Sosed->couple_ != true)
            {
                j->Sosed->Center->get(x2, y2, z2);
                if (kv(x2 - x1) + kv(y2 - y1) + kv(z2 - z1) < kv(3.0 * i->dist))
                {
                    if (j->Sosed->include_ == true && j->Sosed->type == C_base)
                    {
                        j->Sosed->include_ = false;
                        j->Sosed->reconstruct_ = true;
                        for (auto& k : j->Sosed->Grans)
                        {
                            k->Sosed->reconstruct_ = true;
                        }
                    }
                }
            }
        }
    }
}

void Setka::Save_setka(string name)
{
    cout << "Save setka  " << name << endl;
    ofstream fout;
    string namef = "Setka_" + name + ".txt";
    fout.open(namef);

    fout << this->All_Cell.size() << endl;
    for (auto& i : this->All_Cell)
    {
        fout << i->Center->x << " " << i->Center->y << " " << i->Center->z << endl;
        fout << i->type << " " << i->couple_ << " " << i->include_ << endl;
    }

    fout << this->All_Couple.size() << endl;
    for (auto& i : this->All_Couple)
    {
        fout << i->dist << " " << i->A1->number << " " << i->A2->number << endl;
    }

    for (auto& i : this->All_Cell)
    {
        if (i->include_ == false)
        {
            i->Grans.clear();
        }
        fout << i->Grans.size() << endl;
        for (auto& j : i->Grans)
        {
            fout << j->Sosed->number << endl;
        }

        fout << i->Grans_standart.size() << endl;
        for (auto& j : i->Grans_standart)
        {
            fout << j->number << endl;
        }
    }

    this->Save_MHD("MHD_" + name + ".txt");
}

void Setka::Download_setka(string name)
{
    ifstream fout;
    string namef = "Setka_" + name + ".txt";
    fout.open(namef);

    int N = 0;
    double x, y, z;
    int type;
    fout >> N;
    for (int i = 0; i < N; i++)
    {
        fout >> x >> y >> z;
        auto A = new Cell(x, y, z);
        fout >> type >>  A->couple_ >> A->include_;
        A->type = static_cast<Cell_type>(type);
        this->All_Cell.push_back(A);
    }

    int nk = 0;
    for (auto& i : this->All_Cell)
    {
        i->number = nk;
        nk++;
    }

    fout >> N;
    double d;
    int a, b;
    for (int i = 0; i < N; i++)
    {
        fout >> d >> a >> b;
        auto A = new Couple(this->All_Cell[a], this->All_Cell[b], d);
        this->All_Cell[a]->Par = A;
        this->All_Cell[b]->Par = A;
        A->Resolve();
        this->All_Couple.push_back(A);
    }

    nk = 0;
    for (auto& i : this->All_Couple)
    {
        i->number = nk;
        nk++;
    }

    for (auto& i : this->All_Cell)
    {
        int ni;
        fout >> ni;
        for (int j = 0; j < ni; j++)
        {
            fout >> a;
            if (a >= 0)
            {
                    i->Candidates.push_back(this->All_Cell[a]);
            }
            else if (a == -1)
                {
                    i->Candidates.push_back(this->Wall1);
                }
            else if (a == -2)
                {
                    i->Candidates.push_back(this->Wall2);
                }
            else if (a == -3)
                {
                    i->Candidates.push_back(this->Wall3);
                }
            else if (a == -4)
                {
                    i->Candidates.push_back(this->Wall4);
                }
            else if (a == -5)
                {
                    i->Candidates.push_back(this->Wall5);
                }
            else if (a == -6)
                {
                    i->Candidates.push_back(this->Wall6);
                }
            else if (a == -7)
                {
                     i->Candidates.push_back(this->Wall7);
                }
        }

        fout >> ni;
        for (int j = 0; j < ni; j++)
        {
            fout >> a;

            if (a >= 0)
            {
                i->Grans_standart.push_back(this->All_Cell[a]);
            }
            else if (a == -1)
            {
                i->Grans_standart.push_back(this->Wall1);
            }
            else if (a == -2)
            {
                i->Grans_standart.push_back(this->Wall2);
            }
            else if (a == -3)
            {
                i->Grans_standart.push_back(this->Wall3);
            }
            else if (a == -4)
            {
                i->Grans_standart.push_back(this->Wall4);
            }
            else if (a == -5)
            {
                i->Grans_standart.push_back(this->Wall5);
            }
            else if (a == -6)
            {
                i->Grans_standart.push_back(this->Wall6);
            }
            else if (a == -7)
            {
                i->Grans_standart.push_back(this->Wall7);
            }
        }
    }

    this->Download_MHD("MHD_" + name + ".txt");

    this->Construct_fast();
}

void Setka::include_cells(void)
{
    int kj = 0;
    for (auto& i : this->All_Cell)
    {
        if (i->include_ == false && i->couple_ == false)
        {
            for (auto& j : i->Grans)
            {
                delete j;
            }
            i->Grans.clear();

            bool b1 = true;
            bool b2 = false;
            double c1, c2, c3;
            double x1, x2, x3;
            double dist;
            i->Center->get(x1, x2, x3);
            for (auto& j : i->Grans_standart)
            {
                if (j->include_ == true && j->type == C_base)
                {
                    b2 = true;
                    if (j->couple_ == true)
                    {
                        j->Center->get(c1, c2, c3);
                        dist = sqrt(kv(c1 - x1) + kv(c2 - x2) + kv(c3 - x3));
                        if (dist <= j->Par->d_sosed * 1.5)
                        {
                            b1 = false;
                            break;
                        }
                    }

                    for (auto& k : j->Grans)
                    {
                        if (k->Sosed->couple_ == true)
                        {
                            k->Sosed->Center->get(c1, c2, c3);
                            dist = sqrt(kv(c1 - x1) + kv(c2 - x2) + kv(c3 - x3));
                            if (dist <= k->Sosed->Par->d_sosed * 1.5)
                            {
                                b1 = false;
                                break;
                            }
                        }
                    }
                }
            }

            if (b1 == true && b2 == true)
            {
                kj++;
                i->include_ = true;
                i->i_init_mgd = true;
                i->near_par = true;
                for (auto& j : i->Grans_standart)
                {
                    if (j->include_ == true)
                    {
                        j->Candidates.push_back(i);
                        i->Candidates.push_back(j);
                    }
                }
            }
        }
    }

    cout << "Vklucheno  " << kj << "  yacheek" << endl;
}

void Setka::exclude_cells(void)
{
    cout << "EXCLUDE" << endl;
    int kj = 0;
    double x, y, z;
    double x2, y2, z2;
    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true && i->couple_ != true)
        {
            for (auto& j : i->Grans)
            {
                if (j->Sosed->include_ == true && j->Sosed->couple_ == true && j->Sosed->type == C_base)
                {
                    i->Center->get(x, y, z);
                    j->Sosed->Center->get(x2, y2, z2);
                    if (sqrt(kv(x - x2) + kv(y - y2) + kv(z - z2)) <= j->Sosed->Par->d_sosed * 1.5)
                    {
                        kj++;
                        i->include_ = false;                    // Выключаем ячейку
                        i->i_include_candidate = false;
                        break;
                    }
                }
            }

        }
    }

    cout << "Viklucheno  " << kj << "  yacheek" << endl;
}

void Setka::save_soseds()
{
    cout << "Save soseds" << endl;
    for (auto& i : this->All_Cell)
    {
        for (auto& j : i->Grans)
        {
            i->Grans_standart.push_back(j->Sosed);
        }
    }
}

void Setka::Construct_fast(void)
{

#pragma omp parallel for
    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true)
        {
            double x, y, z, r, dd, nx, ny, nz;
            double A, B, C;
            voronoicell_neighbor c;
            vector<int> ney;
            ney.reserve(30);
            vector<double> Surf;
            Surf.reserve(30);
            vector<double> norm;
            norm.reserve(3);
            i->Center->get(x, y, z);

            for (auto& j : i->Grans)
            {
                delete j;
            }

            i->Grans.clear();



            c.init(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0);

            for (auto& j : i->Candidates)
            {
                if (j->type == C_base)
                {
                    A = j->Center->x - x;
                    B = j->Center->y - y;
                    C = j->Center->z - z;
                    c.nplane(A, B, C, (A * A + B * B + C * C), j->number);
                }
                else if (j == this->Wall1)
                {
                    c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
                }
                else if (j == this->Wall2)
                {
                    c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
                }
                else if (j == this->Wall3)
                {
                    c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
                }
                else if (j == this->Wall4)
                {
                    c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
                }
                else if (j == this->Wall5)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
                }
                else if (j == this->Wall6)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
                }
                else if (j == this->Wall7)
                {
                    r = sqrt(kv(x) + kv(y) + kv(z));
                    dd = r - (this->sl_L - 0.01 * RR_);
                    nx = x / r;
                    ny = y / r;
                    nz = z / r;
                    c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, -2.0 * dd * nz, -7);
                }
                else
                {
                    cout << "ERROR  395  Setka  wcewexdwadfeacefhfjkhgefyuevy34u54" << endl;
                }
            }



            c.neighbors(ney); // Получаем список номеров соседей по порядку
            c.face_areas(Surf);   // Площади граней
            c.normals(norm);

            double VV;
            VV = c.volume();  // Записали объём ячейки
            i->set_Volume(VV);

            double nn1, nn2, nn3;
            // Вычисляем движение для пар
            if (i->couple_ == true)
            {
                i->move1 = i->move2 = i->move3 = 0.0;
                int kkk;
                for (int j = 0; j < ney.size(); j++)
                {
                    if (ney[j] >= 0)
                    {
                        kkk = 1;
                        if (this->All_Cell[ney[j]]->couple_ == true)
                        {
                            kkk = 3;
                        }
                        nn1 = (this->All_Cell[ney[j]]->Center->x - x) / 2.0;
                        nn2 = (this->All_Cell[ney[j]]->Center->y - y) / 2.0;
                        nn3 = (this->All_Cell[ney[j]]->Center->z - z) / 2.0;
                        i->move1 += Surf[j] * kkk * (nn1) / 3.0;
                        i->move2 += Surf[j] * kkk * (nn2) / 3.0;
                        i->move3 += Surf[j] * kkk * (nn3) / 3.0;
                    }
                }

                i->move1 = i->move1 / VV;
                i->move2 = i->move2 / VV;
                i->move3 = i->move3 / VV;
            }


            for (int j = 0; j < ney.size(); j++)
            {
                auto A = new Gran();
                A->S = Surf[j];              // Площадь грани
                A->n1 = norm[j * 3];
                A->n2 = norm[j * 3 + 1];
                A->n3 = norm[j * 3 + 2];


                if (ney[j] >= 0)
                {
                    A->Sosed = this->All_Cell[ney[j]];
                    A->c1 = (x + this->All_Cell[ney[j]]->Center->x) / 2.0;   // Вектор центройда грани
                    A->c2 = (y + this->All_Cell[ney[j]]->Center->y) / 2.0;
                    A->c3 = (z + this->All_Cell[ney[j]]->Center->z) / 2.0;
                }
                else if (ney[j] == -1)
                {
                    A->Sosed = this->Wall1;
                }
                else if (ney[j] == -2)
                {
                    A->Sosed = this->Wall2;
                }
                else if (ney[j] == -3)
                {
                    A->Sosed = this->Wall3;
                }
                else if (ney[j] == -4)
                {
                    A->Sosed = this->Wall4;
                }
                else if (ney[j] == -5)
                {
                    A->Sosed = this->Wall5;
                }
                else if (ney[j] == -6)
                {
                    A->Sosed = this->Wall6;
                }
                else if (ney[j] == -7)
                {
                    A->Sosed = this->Wall7;
                }
                else
                {
                    cout << "ERROR  395  Setka  hfjkhgefyuevy34u54" << endl;
                }

                i->Grans.push_back(A);
            }

        }
    }
}

void Setka::Reconstruct_fast(bool wide)
{
    // wide = false - перестроиваем парные и их соседей.
    // wide = true - перестриваем парные и их соседей до 2 уровня включительно

#pragma omp parallel for
    for (auto& i : this->All_Cell)
    {
        bool bbb = false;
        if (wide == false)
        {
            bbb = i->near_par || i->couple_;
        }
        else
        {
            bbb = i->near_par || i->near_par_2 || i->couple_;
        }


        if (i->include_ == true && (bbb == true))
        {
            i->Candidates.clear();
            for (auto& j : i->Grans)
            {
                i->Candidates.push_back(j->Sosed);
            }
        }

    }

#pragma omp parallel for
    for (auto& i : this->All_Cell)
    {

        bool bbb = false;
        if (wide == false)
        {
            bbb = i->near_par || i->couple_;
        }
        else
        {
            bbb = i->near_par || i->near_par_2 || i->couple_;
        }

        if (i->include_ == true && (bbb == true))
        {
            double x, y, z, r, dd, nx, ny, nz;
            double A, B, C;
            voronoicell_neighbor c;
            vector<int> ney;
            ney.reserve(30);
            vector<double> Surf;
            Surf.reserve(30);
            vector<double> norm;
            norm.reserve(3);
            i->Center->get(x, y, z);

            for (auto& j : i->Grans)
            {
                delete j;
            }

            i->Grans.clear();



            c.init(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0);

            for (auto& j : i->Candidates)
            {
                if (j->type == C_base)
                {
                    A = j->Center->x - x;
                    B = j->Center->y - y;
                    C = j->Center->z - z;
                    c.nplane(A, B, C, (A * A + B * B + C * C), j->number);
                }
                else if (j == this->Wall1)
                {
                    c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
                }
                else if (j == this->Wall2)
                {
                    c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
                }
                else if (j == this->Wall3)
                {
                    c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
                }
                else if (j == this->Wall4)
                {
                    c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
                }
                else if (j == this->Wall5)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
                }
                else if (j == this->Wall6)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
                }
                else if (j == this->Wall7)
                {
                    r = sqrt(kv(x) + kv(y) + kv(z));
                    dd = r - (this->sl_L - 0.01 * RR_);
                    nx = x / r;
                    ny = y / r;
                    nz = z / r;
                    c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, -2.0 * dd * nz, -7);
                }
                else
                {
                    cout << "ERROR  395  Setka  wcewexdwadfeacefhfjkhgefyuevy34u54" << endl;
                }
            }



            c.neighbors(ney); // Получаем список номеров соседей по порядку
            c.face_areas(Surf);   // Площади граней
            c.normals(norm);

            double VV;
            VV = c.volume();  // Записали объём ячейки
            i->set_Volume(VV);

            double nn1, nn2, nn3;
            // Вычисляем движение для пар
            if (i->couple_ == true)
            {
                i->move1 = i->move2 = i->move3 = 0.0;
                int kkk;
                for (int j = 0; j < ney.size(); j++)
                {
                    if (ney[j] >= 0)
                    {
                        kkk = 1;
                        if (this->All_Cell[ney[j]]->couple_ == true)
                        {
                            kkk = 3;
                        }
                        nn1 = (this->All_Cell[ney[j]]->Center->x - x) / 2.0;
                        nn2 = (this->All_Cell[ney[j]]->Center->y - y) / 2.0;
                        nn3 = (this->All_Cell[ney[j]]->Center->z - z) / 2.0;
                        i->move1 += Surf[j] * kkk * (nn1) / 3.0;
                        i->move2 += Surf[j] * kkk * (nn2) / 3.0;
                        i->move3 += Surf[j] * kkk * (nn3) / 3.0;
                    }
                }

                i->move1 = i->move1 / VV;
                i->move2 = i->move2 / VV;
                i->move3 = i->move3 / VV;
            }


            for (int j = 0; j < ney.size(); j++)
            {
                auto A = new Gran();
                A->S = Surf[j];              // Площадь грани
                A->n1 = norm[j * 3];
                A->n2 = norm[j * 3 + 1];
                A->n3 = norm[j * 3 + 2];


                if (ney[j] >= 0)
                {
                    A->Sosed = this->All_Cell[ney[j]];
                    A->c1 = (x + this->All_Cell[ney[j]]->Center->x) / 2.0;   // Вектор центройда грани
                    A->c2 = (y + this->All_Cell[ney[j]]->Center->y) / 2.0;
                    A->c3 = (z + this->All_Cell[ney[j]]->Center->z) / 2.0;
                }
                else if (ney[j] == -1)
                {
                    A->Sosed = this->Wall1;
                }
                else if (ney[j] == -2)
                {
                    A->Sosed = this->Wall2;
                }
                else if (ney[j] == -3)
                {
                    A->Sosed = this->Wall3;
                }
                else if (ney[j] == -4)
                {
                    A->Sosed = this->Wall4;
                }
                else if (ney[j] == -5)
                {
                    A->Sosed = this->Wall5;
                }
                else if (ney[j] == -6)
                {
                    A->Sosed = this->Wall6;
                }
                else if (ney[j] == -7)
                {
                    A->Sosed = this->Wall7;
                }
                else
                {
                    cout << "ERROR  395  Setka  hfjkhgefyuevy34u54" << endl;
                }

                i->Grans.push_back(A);
            }

        }
    }
}

void Setka::Reconstruct_medium(bool wide)
{
    // wide = false - перестроиваем парные и их соседей.
    // wide = true - перестриваем парные и их соседей до 2 уровня включительно

    // Добавляем кандидатов для ячеек
    // При этом для парных добавляем соседей-соседей без повторения
#pragma omp parallel for
    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true)
        {
            bool bbb = false;
            if (wide == false)
            {
                bbb = i->near_par || i->couple_;
            }
            else
            {
                bbb = i->near_par || i->near_par_2 || i->couple_;
            }

            if (bbb == false)
            {
                continue;
            }

            i->Candidates.clear();

            if (i->couple_ == false)
            {
                for (auto& j : i->Grans)
                {
                    if (j->Sosed->number != i->number && j->Sosed->include_ == true)
                    {
                        i->Candidates.push_back(j->Sosed);
                    }
                }
            }
            else
            {
                for (auto& j : i->Grans)
                {
                    if (j->Sosed->number != i->number && j->Sosed->include_ == true)
                    {
                        i->Candidates.push_back(j->Sosed);
                    }

                    for (auto& k : j->Sosed->Grans)
                    {
                        if (k->Sosed->number != i->number && k->Sosed->include_ == true)
                        {
                            bool bnm = false;
                            for (auto& bi : i->Candidates)
                            {
                                if (bi->number == k->Sosed->number)
                                {
                                    bnm = true;
                                    break;
                                }
                            }
                            if (bnm == false)
                            {
                                i->Candidates.push_back(k->Sosed);
                            }
                        }
                    }

                }
            }
        }
    }

    // Вычисляем парные ячейки, при этом добавляем кандидатов в обычные от парных
#pragma omp parallel for
    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true && i->couple_ == true)
        {
            double x, y, z, r, dd, nx, ny, nz;
            double A, B, C;
            voronoicell_neighbor c;
            vector<int> ney;
            ney.reserve(30);
            vector<double> Surf;
            Surf.reserve(30);
            vector<double> norm;
            norm.reserve(3);
            i->Center->get(x, y, z);

            for (auto& j : i->Grans)
            {
                delete j;
            }

            i->Grans.clear();



            c.init(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0);

            for (auto& j : i->Candidates)
            {
                if (j->type == C_base)
                {
                    A = j->Center->x - x;
                    B = j->Center->y - y;
                    C = j->Center->z - z;
                    c.nplane(A, B, C, (A * A + B * B + C * C), j->number);
                }
                else if (j == this->Wall1)
                {
                    c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
                }
                else if (j == this->Wall2)
                {
                    c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
                }
                else if (j == this->Wall3)
                {
                    c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
                }
                else if (j == this->Wall4)
                {
                    c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
                }
                else if (j == this->Wall5)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
                }
                else if (j == this->Wall6)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
                }
                else if (j == this->Wall7)
                {
                    r = sqrt(kv(x) + kv(y) + kv(z));
                    dd = r - (this->sl_L - 0.01 * RR_);
                    nx = x / r;
                    ny = y / r;
                    nz = z / r;
                    c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, -2.0 * dd * nz, -7);
                }
                else
                {
                    cout << "ERROR  395  Setka  wcewexdwadfeacefhfjkhgefyuevy34u54" << endl;
                }
            }



            c.neighbors(ney); // Получаем список номеров соседей по порядку
            c.face_areas(Surf);   // Площади граней
            c.normals(norm);

            double VV;
            VV = c.volume();  // Записали объём ячейки
            i->set_Volume(VV);

            double nn1, nn2, nn3;
            // Вычисляем движение для пар
            if (true)
            {
                i->move1 = i->move2 = i->move3 = 0.0;
                int kkk;
                for (int j = 0; j < ney.size(); j++)
                {
                    if (ney[j] >= 0)
                    {
                        kkk = 1;
                        if (this->All_Cell[ney[j]]->couple_ == true)
                        {
                            kkk = 3;
                        }
                        nn1 = (this->All_Cell[ney[j]]->Center->x - x) / 2.0;
                        nn2 = (this->All_Cell[ney[j]]->Center->y - y) / 2.0;
                        nn3 = (this->All_Cell[ney[j]]->Center->z - z) / 2.0;
                        i->move1 += Surf[j] * kkk * (nn1) / 3.0;
                        i->move2 += Surf[j] * kkk * (nn2) / 3.0;
                        i->move3 += Surf[j] * kkk * (nn3) / 3.0;
                    }
                }

                i->move1 = i->move1 / VV;
                i->move2 = i->move2 / VV;
                i->move3 = i->move3 / VV;
            }


            for (int j = 0; j < ney.size(); j++)
            {
                auto A = new Gran();
                A->S = Surf[j];              // Площадь грани
                A->n1 = norm[j * 3];
                A->n2 = norm[j * 3 + 1];
                A->n3 = norm[j * 3 + 2];


                if (ney[j] >= 0)
                {
                    // Блок добавления текущей ячейки как соседа к обычной, которая есть сосед текущей парной
                    if (true)
                    {
                        if (this->All_Cell[ney[j]]->couple_ == false)
                        {
                            bool bnt = false;
                            for (auto& bi : this->All_Cell[ney[j]]->Candidates)
                            {
                                if (bi->number == i->number)
                                {
                                    bnt = true;
                                    break;
                                }
                            }

                            if (bnt == false)
                            {
                                this->All_Cell[ney[j]]->Candidates.push_back(i);
                            }
                        }
                    }

                    A->Sosed = this->All_Cell[ney[j]];
                    A->c1 = (x + this->All_Cell[ney[j]]->Center->x) / 2.0;   // Вектор центройда грани
                    A->c2 = (y + this->All_Cell[ney[j]]->Center->y) / 2.0;
                    A->c3 = (z + this->All_Cell[ney[j]]->Center->z) / 2.0;
                }
                else if (ney[j] == -1)
                {
                    A->Sosed = this->Wall1;
                }
                else if (ney[j] == -2)
                {
                    A->Sosed = this->Wall2;
                }
                else if (ney[j] == -3)
                {
                    A->Sosed = this->Wall3;
                }
                else if (ney[j] == -4)
                {
                    A->Sosed = this->Wall4;
                }
                else if (ney[j] == -5)
                {
                    A->Sosed = this->Wall5;
                }
                else if (ney[j] == -6)
                {
                    A->Sosed = this->Wall6;
                }
                else if (ney[j] == -7)
                {
                    A->Sosed = this->Wall7;
                }
                else
                {
                    cout << "ERROR  395  Setka  hfjkhgefyuevy34u54" << endl;
                }

                i->Grans.push_back(A);
            }

        }
    }

    // Теперь перестраиваем обычные ячейки
#pragma omp parallel for
    for (auto& i : this->All_Cell)
    {
        if (i->include_ == false || i->couple_ == true)
        {
            continue;
        }

        bool bbb = false;
        if (wide == false)
        {
            bbb = i->near_par;
        }
        else
        {
            bbb = i->near_par || i->near_par_2;
        }

        if (bbb == true)
        {
            double x, y, z, r, dd, nx, ny, nz;
            double A, B, C;
            voronoicell_neighbor c;
            vector<int> ney;
            ney.reserve(30);
            vector<double> Surf;
            Surf.reserve(30);
            vector<double> norm;
            norm.reserve(3);
            i->Center->get(x, y, z);

            for (auto& j : i->Grans)
            {
                delete j;
            }

            i->Grans.clear();



            c.init(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0);

            for (auto& j : i->Candidates)
            {
                if (j->type == C_base)
                {
                    A = j->Center->x - x;
                    B = j->Center->y - y;
                    C = j->Center->z - z;
                    c.nplane(A, B, C, (A * A + B * B + C * C), j->number);
                }
                else if (j == this->Wall1)
                {
                    c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
                }
                else if (j == this->Wall2)
                {
                    c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
                }
                else if (j == this->Wall3)
                {
                    c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
                }
                else if (j == this->Wall4)
                {
                    c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
                }
                else if (j == this->Wall5)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
                }
                else if (j == this->Wall6)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
                }
                else if (j == this->Wall7)
                {
                    r = sqrt(kv(x) + kv(y) + kv(z));
                    dd = r - (this->sl_L - 0.01 * RR_);
                    nx = x / r;
                    ny = y / r;
                    nz = z / r;
                    c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, -2.0 * dd * nz, -7);
                }
                else
                {
                    cout << "ERROR  395  Setka  wcewexdwadfeacefhfjkhgefyuevy34u54" << endl;
                }
            }



            c.neighbors(ney); // Получаем список номеров соседей по порядку
            c.face_areas(Surf);   // Площади граней
            c.normals(norm);

            double VV;
            VV = c.volume();  // Записали объём ячейки
            i->set_Volume(VV);

            double nn1, nn2, nn3;
            // Вычисляем движение для пар
            if (i->couple_ == true)
            {
                i->move1 = i->move2 = i->move3 = 0.0;
                int kkk;
                for (int j = 0; j < ney.size(); j++)
                {
                    if (ney[j] >= 0)
                    {
                        kkk = 1;
                        if (this->All_Cell[ney[j]]->couple_ == true)
                        {
                            kkk = 3;
                        }
                        nn1 = (this->All_Cell[ney[j]]->Center->x - x) / 2.0;
                        nn2 = (this->All_Cell[ney[j]]->Center->y - y) / 2.0;
                        nn3 = (this->All_Cell[ney[j]]->Center->z - z) / 2.0;
                        i->move1 += Surf[j] * kkk * (nn1) / 3.0;
                        i->move2 += Surf[j] * kkk * (nn2) / 3.0;
                        i->move3 += Surf[j] * kkk * (nn3) / 3.0;
                    }
                }

                i->move1 = i->move1 / VV;
                i->move2 = i->move2 / VV;
                i->move3 = i->move3 / VV;
            }


            for (int j = 0; j < ney.size(); j++)
            {
                auto A = new Gran();
                A->S = Surf[j];              // Площадь грани
                A->n1 = norm[j * 3];
                A->n2 = norm[j * 3 + 1];
                A->n3 = norm[j * 3 + 2];


                if (ney[j] >= 0)
                {
                    A->Sosed = this->All_Cell[ney[j]];
                    A->c1 = (x + this->All_Cell[ney[j]]->Center->x) / 2.0;   // Вектор центройда грани
                    A->c2 = (y + this->All_Cell[ney[j]]->Center->y) / 2.0;
                    A->c3 = (z + this->All_Cell[ney[j]]->Center->z) / 2.0;
                }
                else if (ney[j] == -1)
                {
                    A->Sosed = this->Wall1;
                }
                else if (ney[j] == -2)
                {
                    A->Sosed = this->Wall2;
                }
                else if (ney[j] == -3)
                {
                    A->Sosed = this->Wall3;
                }
                else if (ney[j] == -4)
                {
                    A->Sosed = this->Wall4;
                }
                else if (ney[j] == -5)
                {
                    A->Sosed = this->Wall5;
                }
                else if (ney[j] == -6)
                {
                    A->Sosed = this->Wall6;
                }
                else if (ney[j] == -7)
                {
                    A->Sosed = this->Wall7;
                }
                else
                {
                    cout << "ERROR  395  Setka  hfjkhgefyuevy34u54" << endl;
                }

                i->Grans.push_back(A);
            }

        }
    }

}

bool Setka::Reconstruct_medium2(bool wide)
{
    // Перестроение ВСЕЙ (не только соседних с парами) сетки!!!
    // Добавляем соседей-соседей и строим

    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true)
        {
            i->Candidates.clear();

            for (auto& j : i->Grans)
            {
                if (j->Sosed->number != i->number && j->Sosed->include_ == true)
                {
                    i->Candidates.push_back(j->Sosed);
                    j->Sosed->Candidates.push_back(i);
                }
            }

            for (auto& jj : i->Grans)
            {
                for (auto& j : jj->Sosed->Grans)
                {
                    if (j->Sosed->number != i->number && j->Sosed->include_ == true)
                    {
                        i->Candidates.push_back(j->Sosed);
                    }
                }
            }

            
        }
    }

    // Вычисляем парные ячейки, при этом добавляем кандидатов в обычные от парных
    bool bkl = true;
#pragma omp parallel for
    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true)
        {
            double x, y, z, r, dd, nx, ny, nz;
            double A, B, C;
            voronoicell_neighbor c;
            vector<int> ney;
            ney.reserve(30);
            vector<double> Surf;
            Surf.reserve(30);
            vector<double> norm;
            norm.reserve(3);
            i->Center->get(x, y, z);
            int nkl = i->Grans.size();
            for (auto& j : i->Grans)
            {
                delete j;
            }

            i->Grans.clear();



            c.init(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0);


            for (auto& j : i->Candidates)
            {
                if (j->type == C_base)
                {
                    A = j->Center->x - x;
                    B = j->Center->y - y;
                    C = j->Center->z - z;
                    c.nplane(A, B, C, (A * A + B * B + C * C), j->number);
                }
                else if (j == this->Wall1)
                {
                    c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
                }
                else if (j == this->Wall2)
                {
                    c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
                }
                else if (j == this->Wall3)
                {
                    c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
                }
                else if (j == this->Wall4)
                {
                    c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
                }
                else if (j == this->Wall5)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
                }
                else if (j == this->Wall6)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
                }
                else if (j == this->Wall7)
                {
                    r = sqrt(kv(x) + kv(y) + kv(z));
                    dd = r - (this->sl_L - 0.01 * RR_);
                    nx = x / r;
                    ny = y / r;
                    nz = z / r;
                    c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, -2.0 * dd * nz, -7);
                }
                else
                {
                    cout << "ERROR  3ewd  Setka  wcewdwedewdewazxczxdwadfeacefhfjkhgefyuevy34u54" << endl;
                }
            }



            c.neighbors(ney); // Получаем список номеров соседей по порядку
            c.face_areas(Surf);   // Площади граней
            c.normals(norm);

            double VV;
            VV = c.volume();  // Записали объём ячейки
            i->set_Volume(VV);

            double nn1, nn2, nn3;
            // Вычисляем движение для пар
            if (i->couple_ == true)
            {
                i->move1 = i->move2 = i->move3 = 0.0;
                int kkk;
                for (int j = 0; j < ney.size(); j++)
                {
                    if (ney[j] >= 0)
                    {
                        kkk = 1;
                        if (this->All_Cell[ney[j]]->couple_ == true)
                        {
                            kkk = 3;
                        }
                        nn1 = (this->All_Cell[ney[j]]->Center->x - x) / 2.0;
                        nn2 = (this->All_Cell[ney[j]]->Center->y - y) / 2.0;
                        nn3 = (this->All_Cell[ney[j]]->Center->z - z) / 2.0;
                        i->move1 += Surf[j] * kkk * (nn1) / 3.0;
                        i->move2 += Surf[j] * kkk * (nn2) / 3.0;
                        i->move3 += Surf[j] * kkk * (nn3) / 3.0;
                    }
                }

                i->move1 = i->move1 / VV;
                i->move2 = i->move2 / VV;
                i->move3 = i->move3 / VV;
            }


            for (int j = 0; j < ney.size(); j++)
            {
                auto A = new Gran();
                A->S = Surf[j];              // Площадь грани
                A->n1 = norm[j * 3];
                A->n2 = norm[j * 3 + 1];
                A->n3 = norm[j * 3 + 2];


                if (ney[j] >= 0)
                {
                    A->Sosed = this->All_Cell[ney[j]];
                    A->c1 = (x + this->All_Cell[ney[j]]->Center->x) / 2.0;   // Вектор центройда грани
                    A->c2 = (y + this->All_Cell[ney[j]]->Center->y) / 2.0;
                    A->c3 = (z + this->All_Cell[ney[j]]->Center->z) / 2.0;
                }
                else if (ney[j] == -1)
                {
                    A->Sosed = this->Wall1;
                }
                else if (ney[j] == -2)
                {
                    A->Sosed = this->Wall2;
                }
                else if (ney[j] == -3)
                {
                    A->Sosed = this->Wall3;
                }
                else if (ney[j] == -4)
                {
                    A->Sosed = this->Wall4;
                }
                else if (ney[j] == -5)
                {
                    A->Sosed = this->Wall5;
                }
                else if (ney[j] == -6)
                {
                    A->Sosed = this->Wall6;
                }
                else if (ney[j] == -7)
                {
                    A->Sosed = this->Wall7;
                }
                else
                {
                    cout << "ERROR  395  Setka  hfjkhgefyuevy34u54" << endl;
                }

                i->Grans.push_back(A);
            }

            if (nkl != i->Grans.size())
            {
                bkl = false;
            }
        }
    }
    return bkl;
}

bool Setka::Reconstruct_medium3(bool wide)
{
    // wide = false - перестроиваем парные и их соседей.
    // wide = true - перестриваем парные и их соседей до 2 уровня включительно

    // Добавляем кандидатов для ячеек
    // Добавляем соседей-соседей по алгоритму без повторений, но ВСЁ ТОЛЬКО для ячеек близких к разрыву

    // Обнулим маячок
    bool bbb;

    //cout << "A1" << endl;
    for (auto& i : this->All_Cell)
    {
        i->i_bbb = false;
        if (i->include_ == true)
        {
            bbb = false;
            if (wide == false)
            {
                bbb = i->near_par || i->couple_;
            }
            else
            {
                bbb = i->near_par || i->near_par_2 || i->couple_;
            }
            i->i_bbb = bbb;
            if (bbb == false)
            {
                continue;
            }

            for (auto& j : i->Grans)
            {
                if (j->Sosed->type == C_base && j->Sosed->include_ == true)
                {
                    i->Candidates.push_back(j->Sosed);
                }
            }
        }
    }

    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true)
        {
            for (auto& j : i->Grans)
            {
                j->Sosed->i_sosed = false;
                for (auto& k : j->Sosed->Grans)
                {
                    k->Sosed->i_sosed = false;
                }

                if (j->Sosed->i_boundary_1 == true)
                {
                    i->i_boundary_2 = true;
                }
            }

            for (auto& j : i->Grans)
            {
                j->Sosed->i_sosed = true;
            }

            for (auto& j : i->Grans)
            {
                if (j->Sosed->type == C_base)
                {
                    for (auto& k : j->Sosed->Grans)
                    {
                        if (k->Sosed->type == C_base && k->Sosed->i_sosed == false && k->Sosed->include_ == true)
                        {
                            k->Sosed->i_sosed = true;
                            if (k->Sosed->i_bbb == true)
                            {
                                k->Sosed->Candidates.push_back(i);
                            }
                        }
                    }
                }
            }
        }
    }
    //cout << "A2" << endl;

    // Вычисляем парные ячейки, при этом добавляем кандидатов в обычные от парных
    //cout << "B1" << endl;
    bool bkl = true;
#pragma omp parallel for
    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true)
        {
            if (i->i_bbb == false)
            {
                continue;
            }

            double x, y, z, r, dd, nx, ny, nz;
            double A, B, C;
            voronoicell_neighbor c;
            vector<int> ney;
            ney.reserve(30);
            vector<double> Surf;
            Surf.reserve(30);
            vector<double> norm;
            norm.reserve(3);
            i->Center->get(x, y, z);
            int bnl = i->Grans.size();
            int bnl2 = 0;
            for (auto& j : i->Grans)
            {
                bnl2 = bnl2 + j->Sosed->number;
                delete j;
            }

            i->Grans.clear();



            c.init(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0);

            for (auto& j : i->Candidates)
            {

                if (j->type == C_base)
                {
                    if (j->include_ == false)
                    {
                        cout << "no include sfinefrlerf" << endl;
                    }

                    A = j->Center->x - x;
                    B = j->Center->y - y;
                    C = j->Center->z - z;
                    c.nplane(A, B, C, (A * A + B * B + C * C), j->number);
                }
            }


            if (i->i_boundary_1 || i->i_boundary_2)
            {
                c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
                c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
                c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
                c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
                c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
                c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
            }




            c.neighbors(ney); // Получаем список номеров соседей по порядку
            c.face_areas(Surf);   // Площади граней
            c.normals(norm);

            double VV;
            VV = c.volume();  // Записали объём ячейки
            i->set_Volume(VV);

            double nn1, nn2, nn3;
            // Вычисляем движение для пар
            if (i->couple_ == true)
            {
                i->move1 = i->move2 = i->move3 = 0.0;
                int kkk;
                int normk = 0;
                int nk = 0;
                for (int j = 0; j < ney.size(); j++)
                {
                    if (ney[j] >= 0)
                    {
                        kkk = 1;
                        if (this->All_Cell[ney[j]]->couple_ == true)
                        {
                            kkk = 1;   // 3
                        }
                        nn1 = (this->All_Cell[ney[j]]->Center->x - x) / 2.0;
                        nn2 = (this->All_Cell[ney[j]]->Center->y - y) / 2.0;
                        nn3 = (this->All_Cell[ney[j]]->Center->z - z) / 2.0;
                        i->move1 += Surf[j] * kkk * (nn1) / 3.0;
                        i->move2 += Surf[j] * kkk * (nn2) / 3.0;
                        i->move3 += Surf[j] * kkk * (nn3) / 3.0;
                        normk = normk + kkk;
                        nk++;
                    }
                }

                i->move1 = i->move1 / VV * nk / normk;
                i->move2 = i->move2 / VV * nk / normk;
                i->move3 = i->move3 / VV * nk / normk;
            }

            int bnk2 = 0;
            for (int j = 0; j < ney.size(); j++)
            {
                auto A = new Gran();
                A->S = Surf[j];              // Площадь грани
                A->n1 = norm[j * 3];
                A->n2 = norm[j * 3 + 1];
                A->n3 = norm[j * 3 + 2];


                if (ney[j] >= 0)
                {
                    A->Sosed = this->All_Cell[ney[j]];
                    A->c1 = (x + this->All_Cell[ney[j]]->Center->x) / 2.0;   // Вектор центройда грани
                    A->c2 = (y + this->All_Cell[ney[j]]->Center->y) / 2.0;
                    A->c3 = (z + this->All_Cell[ney[j]]->Center->z) / 2.0;
                }
                else if (ney[j] == -1)
                {
                    A->Sosed = this->Wall1;
                }
                else if (ney[j] == -2)
                {
                    A->Sosed = this->Wall2;
                }
                else if (ney[j] == -3)
                {
                    A->Sosed = this->Wall3;
                }
                else if (ney[j] == -4)
                {
                    A->Sosed = this->Wall4;
                }
                else if (ney[j] == -5)
                {
                    A->Sosed = this->Wall5;
                }
                else if (ney[j] == -6)
                {
                    A->Sosed = this->Wall6;
                }
                else if (ney[j] == -7)
                {
                    A->Sosed = this->Wall7;
                }
                else
                {
                    cout << "ERROR  395  Setka  hfjkhgefyuevy34u54" << endl;
                }
                bnk2 = bnk2 + A->Sosed->number;
                i->Grans.push_back(A);
            }

            i->Candidates.clear();
            if (bnl != i->Grans.size() || bnl2 != bnk2)
            {
                bkl = false;
            }

        }
    }
    //cout << "B2" << endl;
    return bkl;
}

bool Setka::Reconstruct_medium4(bool wide)
{
    // wide = false - перестроиваем парные и их соседей.
    // wide = true - перестриваем парные и их соседей до 2 уровня включительно

    // Добавляем кандидатов для ячеек
    // Добавляем соседей-соседей по алгоритму без повторений, но ВСЁ ТОЛЬКО для ячеек близких к разрыву

    // Обнулим маячок
    bool bbb;

    //cout << "A1" << endl;

    // Не уверен, что здесь нужно это очищать
    /*for (auto& i : this->All_Cell)
    {
        i->Candidates.clear();
    }*/

    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true)
        {
            for (auto& j : i->Grans)
            {
                if (j->Sosed->type == C_base && j->Sosed->include_ == true)
                {
                    i->Candidates.push_back(j->Sosed);
                }
            }
        }
    }

    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true)
        {
            for (auto& j : i->Grans)
            {
                j->Sosed->i_sosed = false;
                for (auto& k : j->Sosed->Grans)
                {
                    k->Sosed->i_sosed = false;
                }

                if (j->Sosed->i_boundary_1 == true)
                {
                    i->i_boundary_2 = true;
                }
            }

            for (auto& j : i->Grans)
            {
                j->Sosed->i_sosed = true;
            }

            for (auto& j : i->Grans)
            {
                if (j->Sosed->type == C_base)
                {
                    for (auto& k : j->Sosed->Grans)
                    {
                        if (k->Sosed->type == C_base && k->Sosed->i_sosed == false && k->Sosed->include_ == true)
                        {
                            k->Sosed->i_sosed = true;
                            k->Sosed->Candidates.push_back(i);
                        }
                    }
                }
            }
        }
    }
    //cout << "A2" << endl;

    // Вычисляем парные ячейки, при этом добавляем кандидатов в обычные от парных
    //cout << "B1" << endl;
    bool bkl = true;
#pragma omp parallel for
    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true)
        {
            double x, y, z, r, dd, nx, ny, nz;
            double A, B, C;
            voronoicell_neighbor c;
            vector<int> ney;
            ney.reserve(30);
            vector<double> Surf;
            Surf.reserve(30);
            vector<double> norm;
            norm.reserve(3);
            i->Center->get(x, y, z);
            int bnl = i->Grans.size();
            int bnl2 = 0;
            for (auto& j : i->Grans)
            {
                bnl2 = bnl2 + j->Sosed->number;
                delete j;
            }

            i->Grans.clear();



            c.init(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0);

            for (auto& j : i->Candidates)
            {

                if (j->type == C_base)
                {
                    if (j->include_ == false)
                    {
                        cout << "no include sfinefrlerf" << endl;
                    }

                    A = j->Center->x - x;
                    B = j->Center->y - y;
                    C = j->Center->z - z;
                    c.nplane(A, B, C, (A * A + B * B + C * C), j->number);
                }
            }


            if (true)//(i->i_boundary_1 || i->i_boundary_2)
            {
                c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
                c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
                c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
                c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
                c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
                c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);

                r = sqrt(kv(x) + kv(y) + kv(z));
                dd = r - (this->sl_L - 0.01 * RR_);
                nx = x / r;
                ny = y / r;
                nz = z / r;
                c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, -2.0 * dd * nz, -7);
            }




            c.neighbors(ney); // Получаем список номеров соседей по порядку
            c.face_areas(Surf);   // Площади граней
            c.normals(norm);

            double VV;
            VV = c.volume();  // Записали объём ячейки
            i->set_Volume(VV);

            double nn1, nn2, nn3;
            // Вычисляем движение для пар
            if (i->couple_ == true)
            {
                i->move1 = i->move2 = i->move3 = 0.0;
                int kkk;
                for (int j = 0; j < ney.size(); j++)
                {
                    if (ney[j] >= 0)
                    {
                        kkk = 1;
                        if (this->All_Cell[ney[j]]->couple_ == true)
                        {
                            kkk = 3;
                        }
                        nn1 = (this->All_Cell[ney[j]]->Center->x - x) / 2.0;
                        nn2 = (this->All_Cell[ney[j]]->Center->y - y) / 2.0;
                        nn3 = (this->All_Cell[ney[j]]->Center->z - z) / 2.0;
                        i->move1 += Surf[j] * kkk * (nn1) / 3.0;
                        i->move2 += Surf[j] * kkk * (nn2) / 3.0;
                        i->move3 += Surf[j] * kkk * (nn3) / 3.0;
                    }
                }

                i->move1 = i->move1 / VV;
                i->move2 = i->move2 / VV;
                i->move3 = i->move3 / VV;
            }

            int bnk2 = 0;
            for (int j = 0; j < ney.size(); j++)
            {
                auto A = new Gran();
                A->S = Surf[j];              // Площадь грани
                A->n1 = norm[j * 3];
                A->n2 = norm[j * 3 + 1];
                A->n3 = norm[j * 3 + 2];


                if (ney[j] >= 0)
                {
                    A->Sosed = this->All_Cell[ney[j]];
                    A->c1 = (x + this->All_Cell[ney[j]]->Center->x) / 2.0;   // Вектор центройда грани
                    A->c2 = (y + this->All_Cell[ney[j]]->Center->y) / 2.0;
                    A->c3 = (z + this->All_Cell[ney[j]]->Center->z) / 2.0;
                }
                else if (ney[j] == -1)
                {
                    A->Sosed = this->Wall1;
                }
                else if (ney[j] == -2)
                {
                    A->Sosed = this->Wall2;
                }
                else if (ney[j] == -3)
                {
                    A->Sosed = this->Wall3;
                }
                else if (ney[j] == -4)
                {
                    A->Sosed = this->Wall4;
                }
                else if (ney[j] == -5)
                {
                    A->Sosed = this->Wall5;
                }
                else if (ney[j] == -6)
                {
                    A->Sosed = this->Wall6;
                }
                else if (ney[j] == -7)
                {
                    A->Sosed = this->Wall7;
                }
                else
                {
                    cout << "ERROR  395  Setka  hfjkhgefyuevy34u54" << endl;
                }
                bnk2 = bnk2 + A->Sosed->number;
                i->Grans.push_back(A);
            }

            i->Candidates.clear();
            if (bnl != i->Grans.size() || bnl2 != bnk2)
            {
                bkl = false;
            }

        }
    }
    //cout << "B2" << endl;
    return bkl;
}

void Setka::Reconstruct_couple(bool t1)
{
    // Будем перестраивать только включённые ячейки с установленным флагом по реконструкции (те, которые вблизи пар).
    //cout << "Reconstruct_couple  " << t1 << endl;

    // Сначала работает с парными ячейками, так как их в любом случае надо перестраивать
    if (t1 == false)
    {
#pragma omp parallel for
        for (auto& i : this->All_Cell)
        {
            if (i->couple_ == true)
            {
                i->Can_neighbours.clear();
                for (auto& j : i->Grans)
                {
                    if (j->Sosed->number != i->number && j->Sosed->include_ == true)
                    {
                        i->Can_neighbours.insert(make_pair(j->Sosed->number, j->Sosed));
                    }

                    for (auto& k : j->Sosed->Grans)
                    {
                        if (k->Sosed->number != i->number && k->Sosed->include_ == true)
                        {
                            i->Can_neighbours.insert(make_pair(k->Sosed->number, k->Sosed));
                        }
                    }
                    
                }
            }
        }
    }

#pragma omp parallel for
    for (auto& i : this->All_Cell)
    {
        if (i->couple_ == true)
        {
            double x, y, z, r, dd, nx, ny, nz;
            double A, B, C;
            voronoicell_neighbor c;
            vector<int> ney;
            ney.reserve(30);
            vector<double> Surf;
            Surf.reserve(30);
            vector<double> norm;
            norm.reserve(3);
            i->Center->get(x, y, z);

            for (auto& j : i->Grans)
            {
                delete j;
            }

            i->Grans.clear();



            c.init(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0);

            for (auto& j : i->Can_neighbours)
            {
                if (j.second->type == C_base)
                {
                    A = j.second->Center->x - x;
                    B = j.second->Center->y - y;
                    C = j.second->Center->z - z;
                    c.nplane(A, B, C, (A * A + B * B + C * C), j.second->number);
                }
                else if (j.second == this->Wall1)
                {
                    c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
                }
                else if (j.second == this->Wall2)
                {
                    c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
                }
                else if (j.second == this->Wall3)
                {
                    c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
                }
                else if (j.second == this->Wall4)
                {
                    c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
                }
                else if (j.second == this->Wall5)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
                }
                else if (j.second == this->Wall6)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
                }
                else if (j.second == this->Wall7)
                {
                    r = sqrt(kv(x) + kv(y) + kv(z));
                    dd = r - (this->sl_L - 0.01 * RR_);
                    nx = x / r;
                    ny = y / r;
                    nz = z / r;
                    c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, -2.0 * dd * nz, -7);
                }
                else
                {
                    cout << "ERROR  395  Setka  wcewexdwadfeacefhfjkhgefyuevy34u54" << endl;
                }
            }



            c.neighbors(ney); // Получаем список номеров соседей по порядку
            c.face_areas(Surf);   // Площади граней
            c.normals(norm);

            double VV;
            VV = c.volume();  // Записали объём ячейки
            i->set_Volume(VV);

            double nn1, nn2, nn3;
            // Вычисляем движение для пар
            if (true)
            {
                i->move1 = i->move2 = i->move3 = 0.0;
                int kkk;
                for (int j = 0; j < ney.size(); j++)
                {
                    if (ney[j] >= 0)
                    {
                        kkk = 1;
                        if (this->All_Cell[ney[j]]->couple_ == true)
                        {
                            kkk = 3;
                        }
                        nn1 = (this->All_Cell[ney[j]]->Center->x - x) / 2.0;
                        nn2 = (this->All_Cell[ney[j]]->Center->y - y) / 2.0;
                        nn3 = (this->All_Cell[ney[j]]->Center->z - z) / 2.0;
                        i->move1 += Surf[j] * kkk * (nn1) / 3.0;
                        i->move2 += Surf[j] * kkk * (nn2) / 3.0;
                        i->move3 += Surf[j] * kkk * (nn3) / 3.0;
                    }
                }

                i->move1 = i->move1 / VV;
                i->move2 = i->move2 / VV;
                i->move3 = i->move3 / VV;
            }


            for (int j = 0; j < ney.size(); j++)
            {
                auto A = new Gran();
                A->S = Surf[j];              // Площадь грани
                A->n1 = norm[j * 3];
                A->n2 = norm[j * 3 + 1];
                A->n3 = norm[j * 3 + 2];


                if (ney[j] >= 0)
                {
                    A->Sosed = this->All_Cell[ney[j]];
                    A->Sosed->reconstruct_ = true;
                    A->c1 = (x + this->All_Cell[ney[j]]->Center->x) / 2.0;   // Вектор центройда грани
                    A->c2 = (y + this->All_Cell[ney[j]]->Center->y) / 2.0;
                    A->c3 = (z + this->All_Cell[ney[j]]->Center->z) / 2.0;
                }
                else if (ney[j] == -1)
                {
                    A->Sosed = this->Wall1;
                }
                else if (ney[j] == -2)
                {
                    A->Sosed = this->Wall2;
                }
                else if (ney[j] == -3)
                {
                    A->Sosed = this->Wall3;
                }
                else if (ney[j] == -4)
                {
                    A->Sosed = this->Wall4;
                }
                else if (ney[j] == -5)
                {
                    A->Sosed = this->Wall5;
                }
                else if (ney[j] == -6)
                {
                    A->Sosed = this->Wall6;
                }
                else if (ney[j] == -7)
                {
                    A->Sosed = this->Wall7;
                }
                else
                {
                    cout << "ERROR  395  Setka  hfjkhgefyuevy34u54" << endl;
                }

                i->Grans.push_back(A);
            }

        }
    }

    // Сначала в списке ПАР расширяем список кандидатов на соседство всевозможными соседями
 
    if (t1 == false)
    {
#pragma omp parallel for
        for (auto& i : this->All_Cell)
        {
            if (i->include_ == true && i->reconstruct_ == true && i->couple_ == false)
            {
                i->Can_neighbours.clear();
                for (auto& j : i->Grans)
                {
                    if (j->Sosed->number != i->number && j->Sosed->include_ == true)
                    {
                        i->Can_neighbours.insert(make_pair(j->Sosed->number, j->Sosed));
                    }

                    for (auto& k : j->Sosed->Grans)
                    {
                        if (k->Sosed->number != i->number && k->Sosed->include_ == true)
                        {
                            i->Can_neighbours.insert(make_pair(k->Sosed->number, k->Sosed));
                        }
                    }
                }
            }
        }
    }


    //for(int ii = 0; ii < this->All_Cell.size(); ii++)
#pragma omp parallel for
    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true && i->reconstruct_ == true && i->couple_ == false)
        {
            double x, y, z, r, dd, nx, ny, nz;
            double A, B, C;
            voronoicell_neighbor c;
            vector<int> ney;
            ney.reserve(30);
            vector<double> Surf;
            Surf.reserve(30);
            vector<double> norm;
            norm.reserve(3);
            i->Center->get(x, y, z);


            for (auto& j : i->Grans)
            {
               delete j;
            }

            i->Grans.clear();
            
            

            c.init(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0);

            for (auto& j : i->Can_neighbours)
            {
                if (j.second->type == C_base)
                {
                    A = j.second->Center->x - x;
                    B = j.second->Center->y - y;
                    C = j.second->Center->z - z;
                    c.nplane(A, B, C, (A * A + B * B + C * C), j.second->number);
                }
                else if (j.second == this->Wall1)
                {
                    c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
                }
                else if (j.second == this->Wall2)
                {
                    c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
                }
                else if (j.second == this->Wall3)
                {
                    c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
                }
                else if (j.second == this->Wall4)
                {
                    c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
                }
                else if (j.second == this->Wall5)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
                }
                else if (j.second == this->Wall6)
                {
                    c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
                }
                else if (j.second == this->Wall7)
                {
                    r = sqrt(kv(x) + kv(y) + kv(z));
                    dd = r - (this->sl_L - 0.01 * RR_);
                    nx = x / r;
                    ny = y / r;
                    nz = z / r;
                    c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, -2.0 * dd * nz, -7);
                }
                else
                {
                    cout << "ERROR  395  Setka  wcewexdwadfeacefhfjkhgefyuevy34u54" << endl;
                }
            }

            

            c.neighbors(ney); // Получаем список номеров соседей по порядку
            c.face_areas(Surf);   // Площади граней
            c.normals(norm);

            double VV;
            VV = c.volume();  // Записали объём ячейки
            i->set_Volume(VV);

            

            // Блок где находился геометрический центр каждой грани (очень долго работает, поэтому удалили)
            if (false)//(i->couple_ == true)
            {
                // Нужно вычислять сдвиг ячейки только для парных ячеек!
                i->move1 = i->move2 = i->move3 = 0.0;
                vector<double> versh;
                versh.reserve(c.p * 3);
                vector<Point*> ver;
                ver.reserve(c.p);
                vector<int> face_vert;
                face_vert.reserve(c.p * 3);

                c.vertices(versh);
                for (int ik = 0; ik < c.p * 3; ik = ik + 3)
                {
                    auto P = new Point(versh[ik], versh[ik + 1], versh[ik + 2]);
                    //cout << P->x << " " << P->y << " " << P->z << endl;
                    ver.push_back(P);
                }


                c.face_vertices(face_vert);
                double m1, m2, m3;
                int ik = 0;
                int g = -1;
            

                while (ik < face_vert.size())
                {
                    g++;
                    m1 = 0.0;
                    m2 = 0.0;
                    m3 = 0.0;
                    for (int kk = 0; kk < face_vert[ik]; kk++)
                    {
                        m1 = m1 + ver[face_vert[ik + kk + 1]]->x;
                        m2 = m2 + ver[face_vert[ik + kk + 1]]->y;
                        m3 = m3 + ver[face_vert[ik + kk + 1]]->z;
                    }

                    i->move1 += Surf[g] * (m1 / face_vert[ik]) / 3.0;
                    i->move2 += Surf[g] * (m2 / face_vert[ik]) / 3.0;
                    i->move3 += Surf[g] * (m3 / face_vert[ik]) / 3.0;
                    ik = ik + face_vert[ik] + 1;
                }

                i->move1 /= 0.1;
                i->move2 /= 0.1;
                i->move3 /= 0.1;

                //cout << i->move1 << " " << i->move2 << " " << i->move3 << endl;

                for (auto& vb : ver)
                {
                    delete vb;
                }

            }

            
            for (int j = 0; j < ney.size(); j++)
            {
                auto A = new Gran();
                A->S = Surf[j];              // Площадь грани
                A->n1 = norm[j * 3];
                A->n2 = norm[j * 3 + 1];
                A->n3 = norm[j * 3 + 2];
                

                if (ney[j] >= 0)
                {
                    A->Sosed = this->All_Cell[ney[j]];
                    A->c1 = (x + this->All_Cell[ney[j]]->Center->x) / 2.0;   // Вектор центройда грани
                    A->c2 = (y + this->All_Cell[ney[j]]->Center->y) / 2.0;
                    A->c3 = (z + this->All_Cell[ney[j]]->Center->z) / 2.0;
                }
                else if (ney[j] == -1)
                {
                    A->Sosed = this->Wall1;
                }
                else if (ney[j] == -2)
                {
                    A->Sosed = this->Wall2;
                }
                else if (ney[j] == -3)
                {
                    A->Sosed = this->Wall3;
                }
                else if (ney[j] == -4)
                {
                    A->Sosed = this->Wall4;
                }
                else if (ney[j] == -5)
                {
                    A->Sosed = this->Wall5;
                }
                else if (ney[j] == -6)
                {
                    A->Sosed = this->Wall6;
                }
                else if (ney[j] == -7)
                {
                    A->Sosed = this->Wall7;
                }
                else
                {
                    cout << "ERROR  395  Setka  hfjkhgefyuevy34u54" << endl;
                }

                i->Grans.push_back(A);
            }

        }
    }

}

void Setka::Cut_Surface(void)
{
    double x, y, z, r, dd, nx, ny, nz;
    double A, B, C;
    voronoicell_neighbor c;

    vector<int> neigh, f_vert, ney, edge, vershin;
    vector<double> X, Y, Ver, Versh;

    cout << "Cut_Surface " << endl;

    ofstream fout;
    string name_f = "Cut_Surface.txt";
    fout.open(name_f);

    int N = 0;
    int E = 0;

    for (auto& i : this->All_Cell)
    {
        if (i->include_ == false)
        {
            continue;
        }

        if (i->couple_ == false)
        {
            continue;
        }

        if (i->Par->A1 != i)
        {
            continue;
        }

        x = i->Center->x;
        y = i->Center->y;
        z = i->Center->z;

        c.init(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0);

        for (auto& j : i->Grans)
        {
            if (j->Sosed->type == C_base)
            {
                A = j->Sosed->Center->x - x;
                B = j->Sosed->Center->y - y;
                C = j->Sosed->Center->z - z;
                c.nplane(A, B, C, (A * A + B * B + C * C), j->Sosed->number);
            }
        }

        // В этой функции можно резать лишнее, т.к. эта функция печати, в других лучше сделать по-другому
        r = sqrt(kv(x) + kv(y) + kv(z));
        dd = r - (this->sl_L - 0.01 * RR_);
        nx = x / r;
        ny = y / r;
        nz = z / r;
        // Далее вручную разрезаем ячейку особыми гранями (стенками и границами)
        c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
        c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
        c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
        c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
        c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
        c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);

        // Весь код далее для печати сетки в файл Текплота. Сложности в том, что нужно по порядку расположить все узлы и грани сетки (Текплот так работает)
        //c.nplane(0.0, 0.0, -2.0 * z + 0.00001 * RR_, -77);

        int sosed_couple = i->Par->A2->number;


        c.neighbors(ney); // Получаем список номеров соседей по порядку
        int k = -1;  // номер нужного соседа
        // Ищем номер соседа "-7"
        int kk = -1;
        for (int j = 0; j < ney.size(); j++)
        {
            kk++;
            if (ney[j] == sosed_couple)
            {
                k = kk;
                break;
            }
        }

        if (k == -1)
        {
            continue;
        }
        //cout << "K = " << k << endl;
        //ищем нужные вершины:
        c.face_vertices(f_vert);
        int m = 0;
        int n = 0;
        int r = 0;
        while (r < k)
        {
            n = f_vert[m];
            m = m + n + 1;
            r++;
        }

        //cout << "m = " << m << endl;

        for (int j = m + 1; j < m + 1 + f_vert[m] - 1; j++)
        {
            vershin.push_back(N + f_vert[j] + 1);  // Вектор номеров нужных вершин
            vershin.push_back(N + f_vert[j + 1] + 1);  // Вектор номеров нужных вершин
        }

        vershin.push_back(N + f_vert[m + 1 + f_vert[m] - 1] + 1);
        vershin.push_back(N + f_vert[m + 1] + 1);

        c.vertices(Ver);

        for (int j = 0; j < Ver.size(); j = j + 3)
        {
            Versh.push_back(Ver[j] + x);  // Вектор всех вершин
            Versh.push_back(Ver[j + 1] + y);  // Вектор всех вершин
            Versh.push_back(Ver[j + 2] + z);  // Вектор всех вершин
        }

        N += Ver.size() / 3;
        E += f_vert[m];
    }

    fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"Z\"  ZONE T= \"HP\", N=" << N << ", E=" << E << ", F=FEPOINT, ET=LINESEG " << endl;

    int jkl = 0;
    for (int j = 0; j < Versh.size(); j++)
    {
        jkl++;
        fout << Versh[j] << " ";
        if (jkl % 3 == 0)
        {
            fout << endl;
        }
    }

    jkl = 0;
    for (int j = 0; j < vershin.size(); j++)
    {
        jkl++;
        fout << vershin[j] << " ";
        if (jkl % 2 == 0)
        {
            fout << endl;
        }
    }
}

void Setka::set_normal(void)
{
#pragma omp parallel for
    for (auto& i : this->All_Couple)
    {
        double n1, n2, n3;
        i->get_normal(n1, n2, n3);
        i->n1 = n1;
        i->n2 = n2;
        i->n3 = n3;
    }
}

void Setka::move_par(void)
{

#pragma omp parallel for
    for (auto& i : this->All_Couple)
    {
        double m1, m2, m3;
        double n1, n2, n3;
        n1 = i->n1;
        n2 = i->n2;
        n3 = i->n3;
        m1 = i->A1->move1;
        m2 = i->A1->move2;
        m3 = i->A1->move3;
        m1 += i->A2->move1;
        m2 += i->A2->move2;
        m3 += i->A2->move3;
        m1 /= 2.0;
        m2 /= 2.0;
        m3 /= 2.0;
        double sk = m1 * n1 + m2 * n2 + m3 * n3;
        m1 = m1 - sk * n1;
        m2 = m2 - sk * n2;
        m3 = m3 - sk * n3;
        i->move(m1 * 0.0007, m2 * 0.0007, m3 * 0.0007);   // 0.001 с Коэффициент нужен для замедления движения, иначе оно слишком большое
    }

}

void Setka::Calc_normal(void)
{

#pragma omp parallel for
    for (auto& ii : this->All_Couple)
    {
        double v1, v2, v3;
        v1 = ii->n1;
        v2 = ii->n2;
        v3 = ii->n3;
        //cout << "do = " << v1 << " " << v2 << " " << v3 << endl;
        ii->n1 = 0.0;
        ii->n2 = 0.0;
        ii->n3 = 0.0;
        map <int, Couple*> Nei;
        vector <Point> NN;
        NN.reserve(9);
        double x, y, z;
        double a, b, c, d, e, f;
        double n1, n2, n3, nn;

        // Сначала ищем всех соседей
        for (auto& j : ii->A1->Grans)
        {
            if (j->Sosed->couple_ == true)
            {
                if (j->Sosed->Par != ii)
                {
                    Nei.insert(make_pair(j->Sosed->Par->number, j->Sosed->Par));
                }
            }
        }

        for (auto& j : ii->A2->Grans)
        {
            if (j->Sosed->couple_ == true)
            {
                if (j->Sosed->Par != ii)
                {
                    Nei.insert(make_pair(j->Sosed->Par->number, j->Sosed->Par));
                }
            }
        }


        int S = Nei.size();

        if (S < 1)
        {
            cout << "Net par!!!  703  sihfbkgvsdueqx13e3r2eda" << endl;
            exit(-3);
        }

        if (S == 1)
        {
            for (auto& kk : Nei)
            {
                kk.second->get_normal(x, y, z);
                ii->n1 = x;
                ii->n2 = y;
                ii->n3 = z;
            }
            continue;
        }

        for (auto& kk : Nei)
        {
            kk.second->get_centr(x, y, z);
            NN.push_back(Point(x, y, z));
        }

        if (S == 2)
        {
            S = 3;
            ii->get_centr(x, y, z);
            NN.push_back(Point(x, y, z));
        }

        // Пробегаемся по всем парам - соседям
        for (unsigned short int i = 0; i < S - 2; i++)
        {
            for (unsigned short int j = i + 1; j < S - 1; j++)
            {
                for (unsigned short int k = j + 1; k < S; k++)
                {
                    a = NN[j].x - NN[i].x;
                    b = NN[k].x - NN[i].x;
                    c = NN[j].y - NN[i].y;
                    d = NN[k].y - NN[i].y;
                    e = NN[j].z - NN[i].z;
                    f = NN[k].z - NN[i].z;
                    n1 = c * f - e * d;
                    n2 = -a * f + e * b;
                    n3 = a * d - b * c;
                    nn = sqrt(kv(n1) + kv(n2) + kv(n3));
                    n1 = n1 / nn;
                    n2 = n2 / nn;
                    n3 = n3 / nn;

                    /*if (i == 0 && j == 1 && k == 2)
                    {
                        if (v1 * n1 + v2 * n2 + v3 * n3 >= 0.0)
                        {
                            v1 = n1;
                            v2 = n2;
                            v3 = n3;
                        }
                        else
                        {
                            v1 = -n1;
                            v2 = -n2;
                            v3 = -n3;
                        }
                    }*/

                    if (v1 * n1 + v2 * n2 + v3 * n3 >= 0.0)
                    {
                        ii->n1 += n1;
                        ii->n2 += n2;
                        ii->n3 += n3;
                    }
                    else
                    {
                        ii->n1 -= n1;
                        ii->n2 -= n2;
                        ii->n3 -= n3;
                    }
                    
                }
            }
        }

        nn = sqrt(kv(ii->n1) + kv(ii->n2) + kv(ii->n3));
        ii->n1 = ii->n1 / nn;
        ii->n2 = ii->n2 / nn;
        ii->n3 = ii->n3 / nn;

        /*if (v1 * n1 + v2 * n2 + v3 * n3 < 0.5)
        {
            cout << v1 << " " << v2 << " " << v3 << endl;
            cout << n1 << " " << n2 << " " << n3 << endl;
            cout << S << endl;
            for (unsigned short int i = 0; i < S; i++)
            {
                cout << NN[i].x << " " << NN[i].y << " " << NN[i].z << endl;
            }
            exit(-1);
        }*/

        //cout << "posle = " << ii->n1 << " " << ii->n2 << " " << ii->n3 << endl;

        //cout << ii->n1 << " - " << ii->n2 << " - " << ii->n3 << endl;
        //cout << v1 << " = " << v2 << " = " << v3 << " " << endl;

    }

    for (auto& ii : this->All_Couple)
    {
        ii->orient();
    }
}

void Setka::Calc_normal2(void)
{

//#pragma omp parallel for
    for (auto& ii : this->All_Couple)
    {
        ii->n1 = 0.0;
        ii->n2 = 0.0;
        ii->n3 = 0.0;
        map <int, Couple*> Nei;
        vector <Point> NN;
        NN.reserve(9);
        double x, y, z;
        double a, b, c, d, e, f;
        double n1, n2, n3, nn;

        // Сначала ищем всех соседей
        for (auto& j : ii->A1->Grans)
        {
            if (j->Sosed->couple_ == true)
            {
                if (j->Sosed->Par != ii)
                {
                    Nei.insert(make_pair(j->Sosed->Par->number, j->Sosed->Par));
                }
            }
        }

        for (auto& j : ii->A2->Grans)
        {
            if (j->Sosed->couple_ == true)
            {
                if (j->Sosed->Par != ii)
                {
                    Nei.insert(make_pair(j->Sosed->Par->number, j->Sosed->Par));
                }
            }
        }


        int S = Nei.size();

        if (S < 1)
        {
            cout << "Net par!!!  703  sihfbkgvsdueqx13e3r2eda" << endl;
            exit(-3);
        }


        for (auto& kk : Nei)
        {
            kk.second->get_normal(x, y, z);
            NN.push_back(Point(x, y, z));
        }


        // Пробегаемся по всем парам - соседям
        for (auto& kk : NN)
        {
            ii->n1 = ii->n1 + kk.x;
            ii->n2 = ii->n2 + kk.y;
            ii->n3 = ii->n3 + kk.z;
        }

        nn = sqrt(kv(ii->n1) + kv(ii->n2) + kv(ii->n3));
        ii->n1 = ii->n1 / nn;
        ii->n2 = ii->n2 / nn;
        ii->n3 = ii->n3 / nn;
    }

    for (auto& ii : this->All_Couple)
    {
        ii->orient();
    }
}

void Setka::Cut_Plane_z(double R)
{
	double x, y, z, r, dd, nx, ny, nz;
	double A, B, C;
	voronoicell_neighbor c;

	vector<int> neigh, f_vert, ney, edge, vershin;
	vector<double> X, Y, Ver, Versh;

	cout << "Cut_Plane_z " << endl;

	ofstream fout;
	string name_f = "setka.txt";
	fout.open(name_f);

	int N = 0;
	int E = 0;

	for (auto& i : this->All_Cell)
	{
		if (fabs(i->Center->z) > R)
		{
			continue;
		}

        if (i->include_ == false)
        {
            continue;
        }

        /*if (i->couple_ == false)
        {
            continue;
        }*/

		x = i->Center->x;
		y = i->Center->y;
		z = i->Center->z;

		c.init(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0);

		for (auto& j : i->Grans)
		{
			if (j->Sosed->type == C_base)
			{
				A = j->Sosed->Center->x - x;
				B = j->Sosed->Center->y - y;
				C = j->Sosed->Center->z - z;
				c.nplane(A, B, C, (A * A + B * B + C * C), j->Sosed->number);
			}
		}

		// В этой функции можно резать лишнее, т.к. эта функция печати, в других лучше сделать по-другому
		r = sqrt(kv(x) + kv(y) + kv(z));
		dd = r - (this->sl_L - 0.01 * RR_);
		nx = x / r;
		ny = y / r;
		nz = z / r;
		// Далее вручную разрезаем ячейку особыми гранями (стенками и границами)
		c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
		c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
		c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
		c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
		c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
		c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
		c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, -2.0 * dd * nz, -7);

		// Весь код далее для печати сетки в файл Текплота. Сложности в том, что нужно по порядку расположить все узлы и грани сетки (Текплот так работает)
		c.nplane(0.0, 0.0, -2.0 * z + 0.00001 * RR_, -77);

		c.neighbors(ney); // Получаем список номеров соседей по порядку
		int k = -1;  // номер нужного соседа
		// Ищем номер соседа "-7"
		int kk = -1;
		for (int j = 0; j < ney.size(); j++)
		{
			kk++;
			if (ney[j] == -77)
			{
				k = kk;
				break;
			}
		}

		if (k == -1)
		{
			continue;
		}
		//cout << "K = " << k << endl;
		//ищем нужные вершины:
		c.face_vertices(f_vert);
		int m = 0;
		int n = 0;
		int r = 0;
		while (r < k)
		{
			n = f_vert[m];
			m = m + n + 1;
			r++;
		}

		//cout << "m = " << m << endl;

		for (int j = m + 1; j < m + 1 + f_vert[m] - 1; j++)
		{
			vershin.push_back(N + f_vert[j] + 1);  // Вектор номеров нужных вершин
			vershin.push_back(N + f_vert[j + 1] + 1);  // Вектор номеров нужных вершин
		}

		vershin.push_back(N + f_vert[m + 1 + f_vert[m] - 1] + 1);
		vershin.push_back(N + f_vert[m + 1] + 1);

		c.vertices(Ver);

		for (int j = 0; j < Ver.size(); j = j + 3)
		{
			Versh.push_back(Ver[j] + x);  // Вектор всех вершин
			Versh.push_back(Ver[j + 1] + y);  // Вектор всех вершин
			Versh.push_back(Ver[j + 2] + z);  // Вектор всех вершин
		}

		N += Ver.size() / 3;
		E += f_vert[m];
	}

	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"Z\"  ZONE T= \"HP\", N=" << N << ", E=" << E << ", F=FEPOINT, ET=LINESEG " << endl;

	int jkl = 0;
	for (int j = 0; j < Versh.size(); j++)
	{
		jkl++;
		fout << Versh[j] << " ";
		if (jkl % 3 == 0)
		{
			fout << endl;
		}
	}

	jkl = 0;
	for (int j = 0; j < vershin.size(); j++)
	{
		jkl++;
		fout << vershin[j] << " ";
		if (jkl % 2 == 0)
		{
			fout << endl;
		}
	}

}

void Setka::Cut_Plane_y(double R)
{
    double x, y, z, r, dd, nx, ny, nz;
    double A, B, C;
    voronoicell_neighbor c;

    vector<int> neigh, f_vert, ney, edge, vershin;
    vector<double> X, Y, Ver, Versh;

    cout << "Cut_Plane_y " << endl;

    ofstream fout;
    string name_f = "setka_y.txt";
    fout.open(name_f);

    int N = 0;
    int E = 0;

    for (auto& i : this->All_Cell)
    {
        if (fabs(i->Center->y) > R)
        {
            continue;
        }

        if (i->include_ == false)
        {
            continue;
        }

        /*if (i->couple_ == false)
        {
            continue;
        }*/

        x = i->Center->x;
        y = i->Center->y;
        z = i->Center->z;

        c.init(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0);

        for (auto& j : i->Grans)
        {
            if (j->Sosed->type == C_base)
            {
                A = j->Sosed->Center->x - x;
                B = j->Sosed->Center->y - y;
                C = j->Sosed->Center->z - z;
                c.nplane(A, B, C, (A * A + B * B + C * C), j->Sosed->number);
            }
        }

        // В этой функции можно резать лишнее, т.к. эта функция печати, в других лучше сделать по-другому
        r = sqrt(kv(x) + kv(y) + kv(z));
        dd = r - (this->sl_L - 0.01 * RR_);
        nx = x / r;
        ny = y / r;
        nz = z / r;
        // Далее вручную разрезаем ячейку особыми гранями (стенками и границами)
        c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
        c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
        c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
        c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
        c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
        c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
        c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, -2.0 * dd * nz, -7);

        // Весь код далее для печати сетки в файл Текплота. Сложности в том, что нужно по порядку расположить все узлы и грани сетки (Текплот так работает)
        c.nplane(0.0, -2.0 * y + 0.00001 * RR_, 0.0, -77);

        c.neighbors(ney); // Получаем список номеров соседей по порядку
        int k = -1;  // номер нужного соседа
        // Ищем номер соседа "-7"
        int kk = -1;
        for (int j = 0; j < ney.size(); j++)
        {
            kk++;
            if (ney[j] == -77)
            {
                k = kk;
                break;
            }
        }

        if (k == -1)
        {
            continue;
        }
        //cout << "K = " << k << endl;
        //ищем нужные вершины:
        c.face_vertices(f_vert);
        int m = 0;
        int n = 0;
        int r = 0;
        while (r < k)
        {
            n = f_vert[m];
            m = m + n + 1;
            r++;
        }

        //cout << "m = " << m << endl;

        for (int j = m + 1; j < m + 1 + f_vert[m] - 1; j++)
        {
            vershin.push_back(N + f_vert[j] + 1);  // Вектор номеров нужных вершин
            vershin.push_back(N + f_vert[j + 1] + 1);  // Вектор номеров нужных вершин
        }

        vershin.push_back(N + f_vert[m + 1 + f_vert[m] - 1] + 1);
        vershin.push_back(N + f_vert[m + 1] + 1);

        c.vertices(Ver);

        for (int j = 0; j < Ver.size(); j = j + 3)
        {
            Versh.push_back(Ver[j] + x);  // Вектор всех вершин
            Versh.push_back(Ver[j + 1] + y);  // Вектор всех вершин
            Versh.push_back(Ver[j + 2] + z);  // Вектор всех вершин
        }

        N += Ver.size() / 3;
        E += f_vert[m];
    }

    fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"Z\"  ZONE T= \"HP\", N=" << N << ", E=" << E << ", F=FEPOINT, ET=LINESEG " << endl;

    int jkl = 0;
    for (int j = 0; j < Versh.size(); j++)
    {
        jkl++;
        fout << Versh[j] << " ";
        if (jkl % 3 == 0)
        {
            fout << endl;
        }
    }

    jkl = 0;
    for (int j = 0; j < vershin.size(); j++)
    {
        jkl++;
        fout << vershin[j] << " ";
        if (jkl % 2 == 0)
        {
            fout << endl;
        }
    }

}

void Setka::Cut_Plane_z_Tecplot(double R)
{
	double x, y, z, r, dd, nx, ny, nz;
	double A, B, C;
	voronoicell_neighbor c;
    double X, Y, Z;
	vector<int> f_vert, ney, vershin;
	vector<double> Ver, Versh;

	cout << "Cut_Plane_z_Tecplot" << endl;

	ofstream fout;
	string name_f = "2D_tecplot_z_plane.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Rho\", \"P\", \"Vx\", \"Vy\", \"Vz\"," << //
		"\"Bx\", \"By\", \"Bz\", ZONE T = \"HP\"" << endl;


	for (auto& i : this->All_Cell)
	{
		if (fabs(i->Center->z) > R)
		{
			continue;
		}

        if (i->include_ == false)
        {
            continue;
        }

        X = Y = Z = 0.0;

		x = i->Center->x;
		y = i->Center->y;
		z = i->Center->z;

		c.init(-20.0, 20.0, -20.0, 20.0, -20.0, 20.0);

		for (auto& j : i->Grans)
		{
			if (j->Sosed->type == C_base)
			{
				A = j->Sosed->Center->x - x;
				B = j->Sosed->Center->y - y;
				C = j->Sosed->Center->z - z;
				c.nplane(A, B, C, (A * A + B * B + C * C), j->Sosed->number);
			}
		}

		// В этой функции можно резать лишнее, т.к. эта функция печати, в других лучше сделать по-другому
		r = sqrt(kv(x) + kv(y) + kv(z));
		dd = r - (this->sl_L - 0.01 * RR_);
		nx = x / r;
		ny = y / r;
		nz = z / r;
		// Далее вручную разрезаем ячейку особыми гранями (стенками и границами)
		c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
		c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
		c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
		c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
		c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
		c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
		c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, -2.0 * dd * nz, -7);

		// Весь код далее для печати сетки в файл Текплота. Сложности в том, что нужно по порядку расположить все узлы и грани сетки (Текплот так работает)
		c.nplane(0.0, 0.0, -2.0 * z + 0.00001 * RR_, -77);

		c.neighbors(ney); // Получаем список номеров соседей по порядку
		int k = -1;  // номер нужного соседа
		// Ищем номер соседа "-7"
		int kk = -1;
		for (int j = 0; j < ney.size(); j++)
		{
			kk++;
			if (ney[j] == -77)
			{
				k = kk;
				break;
			}
		}

		if (k == -1)
		{
			continue;
		}

		//ищем нужные вершины:
		c.face_vertices(f_vert);
		int m = 0;
		int n = 0;
		int r = 0;
		while (r < k)
		{
			n = f_vert[m];
			m = m + n + 1;
			r++;
		}

		//cout << "m = " << m << endl;
        c.vertices(Ver);
        int kkk = 0;
		for (int j = m + 1; j < m + 1 + f_vert[m]; j++)
		{
            kkk++;
            X = X + Ver[f_vert[j] * 3] + x;
            Y = Y + Ver[f_vert[j] * 3 + 1] + y;
            Z = Z + Ver[f_vert[j] * 3 + 2] + z;
		}

        X = X / kkk;
        Y = Y / kkk;
        Z = Z / kkk;

        fout << X << " " << Y << " " << sqrt(kv(X) + kv(Y)) << " " << i->par[0].ro << " " <<//
            i->par[0].p << " " << i->par[0].u << " " << i->par[0].v << " " << i->par[0].w << " " << //
            i->par[0].bx << " " << i->par[0].by << " " << i->par[0].bz << " " << endl;

	}

}

void Setka::Cut_Plane_y_Tecplot(double R)
{
    double x, y, z, r, dd, nx, ny, nz;
    double A, B, C;
    voronoicell_neighbor c;
    double X, Y, Z;
    vector<int> f_vert, ney, vershin;
    vector<double> Ver, Versh;

    cout << "Cut_Plane_y_Tecplot" << endl;

    ofstream fout;
    string name_f = "2D_tecplot_y_plane.txt";
    fout.open(name_f);
    fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Z\", \"r\", \"Rho\", \"P\", \"Vx\", \"Vy\", \"Vz\"," << //
        "\"Bx\", \"By\", \"Bz\", ZONE T = \"HP\"" << endl;


    for (auto& i : this->All_Cell)
    {
        if (fabs(i->Center->y) > R)
        {
            continue;
        }

        if (i->include_ == false)
        {
            continue;
        }

        X = Y = Z = 0.0;

        x = i->Center->x;
        y = i->Center->y;
        z = i->Center->z;

        c.init(-20.0, 20.0, -20.0, 20.0, -20.0, 20.0);

        for (auto& j : i->Grans)
        {
            if (j->Sosed->type == C_base)
            {
                A = j->Sosed->Center->x - x;
                B = j->Sosed->Center->y - y;
                C = j->Sosed->Center->z - z;
                c.nplane(A, B, C, (A * A + B * B + C * C), j->Sosed->number);
            }
        }

        // В этой функции можно резать лишнее, т.к. эта функция печати, в других лучше сделать по-другому
        r = sqrt(kv(x) + kv(y) + kv(z));
        dd = r - (this->sl_L - 0.01 * RR_);
        nx = x / r;
        ny = y / r;
        nz = z / r;
        // Далее вручную разрезаем ячейку особыми гранями (стенками и границами)
        c.nplane(2.0 * (x_max - x), 0.0, 0.0, -1);
        c.nplane(2.0 * (x_min - x), 0.0, 0.0, -2);
        c.nplane(0.0, 2.0 * (y_max - y), 0.0, -3);
        c.nplane(0.0, 2.0 * (y_min - y), 0.0, -4);
        c.nplane(0.0, 0.0, 2.0 * (z_max - z), -5);
        c.nplane(0.0, 0.0, 2.0 * (z_min - z), -6);
        c.nplane(-2.0 * dd * nx, -2.0 * dd * ny, -2.0 * dd * nz, -7);

        // Весь код далее для печати сетки в файл Текплота. Сложности в том, что нужно по порядку расположить все узлы и грани сетки (Текплот так работает)
        c.nplane(0.0, -2.0 * y + 0.00001 * RR_, 0.0, -77);

        c.neighbors(ney); // Получаем список номеров соседей по порядку
        int k = -1;  // номер нужного соседа
        // Ищем номер соседа "-7"
        int kk = -1;
        for (int j = 0; j < ney.size(); j++)
        {
            kk++;
            if (ney[j] == -77)
            {
                k = kk;
                break;
            }
        }

        if (k == -1)
        {
            continue;
        }

        //ищем нужные вершины:
        c.face_vertices(f_vert);
        int m = 0;
        int n = 0;
        int r = 0;
        while (r < k)
        {
            n = f_vert[m];
            m = m + n + 1;
            r++;
        }

        //cout << "m = " << m << endl;
        c.vertices(Ver);
        int kkk = 0;
        for (int j = m + 1; j < m + 1 + f_vert[m]; j++)
        {
            kkk++;
            X = X + Ver[f_vert[j] * 3] + x;
            Y = Y + Ver[f_vert[j] * 3 + 1] + y;
            Z = Z + Ver[f_vert[j] * 3 + 2] + z;
        }

        X = X / kkk;
        Y = Y / kkk;
        Z = Z / kkk;

        fout << X << " " << Z << " " << sqrt(kv(X) + kv(Y)) << " " << i->par[0].ro << " " <<//
            i->par[0].p << " " << i->par[0].u << " " << i->par[0].v << " " << i->par[0].w << " " << //
            i->par[0].bx << " " << i->par[0].by << " " << i->par[0].bz << " " << endl;

    }

}

void Setka::Initialization_do_MHD(void)
{
    cout << "Initialization_do_MHD" << endl;

    for (auto& i : this->All_Cell)
    {
        i->Candidates.clear();
        double V = 0.0;
        i->get_Volume(V);
        i->set_Volume(V);
        i->Center->renew();
        i->near_par_2 = false;
        i->near_par = false;
        if (i->include_ == true)
        {
            if (i->couple_ == false)
            {
                for (auto& kk : i->Grans)
                {
                    if (kk->Sosed->couple_ == true)
                    {
                        i->near_par = true;
                        break;
                    }
                }
            }
        }
    }

    for (auto& ii : this->All_Cell)
    {
        if (ii->include_ == true)
        {
            if (ii->couple_ == false && ii->near_par == false)
            {
                for (auto& kk : ii->Grans)
                {
                    if (kk->Sosed->near_par == true)
                    {
                        ii->near_par_2 = true;
                        break;
                    }
                }
            }
        }
    }

    for (auto& i : this->All_Cell)
    {
        for (auto& j : i->Grans)
        {
            if (j->Sosed->type != C_base)
            {
                i->i_boundary_1 = true;
            }
        }
    }

    for (auto& i : this->All_Couple)
    {
        i->Resolve();
        double n1, n2, n3;
        i->get_normal(n1, n2, n3);
        i->n1 = n1;
        i->n2 = n2;
        i->n3 = n3;
    }

    for (auto& i : this->All_Cell)
    {
        if (i->include_ == false)
        {
            i->mgd_ = false;
            continue;
        }
        else
        {
            i->mgd_ = true;
        }

        bool inner = false;
        bool intry = false;
        for (auto& j : i->Grans)
        {
            if (j->Sosed->type == C_sphere)
            {
                inner = true;
            }

            if (j->Sosed->type == C_wall_x_max)
            {
                intry = true;
            }
        }

        if (intry == true || inner == true)
        {
            i->mgd_ = false;
        }
    }

    for (auto& i : this->All_Cell)
    {
        if (i->include_ == true && i->mgd_ == true)
        {
            if (i->par[0].ro == 0.0)
            {
                cout << i->couple_ << " " << i->type << " " << i->Center->x << " " << i->Center->y << " " << i->Center->z << endl;
                cout << i->number << " " << i->Grans.size() << endl;
                exit(-1);
            }

            if (i->par[1].ro == 0.0)
            {
                cout << "2 --- " << i->couple_ << " " << i->type << " " << i->Center->x << " " << i->Center->y << " " << i->Center->z << endl;
                exit(-1);
            }
        }
    }

    for (auto& i : this->All_Cell)
    {
        if (i->i_init_mgd == true)
        {
            double ro = 0.0;
            double p = 0.0;
            double u = 0.0;
            double v = 0.0;
            double w = 0.0;
            double bx = 0.0;
            double by = 0.0;
            double bz = 0.0;
            int kk = 0;
            for (auto& j : i->Grans)
            {
                if (j->Sosed->type == C_base && j->Sosed->include_ == true && j->Sosed->i_init_mgd == false)
                {
                    kk++;
                    ro = ro + j->Sosed->par[0].ro;
                    p = p + j->Sosed->par[0].p;
                    u = u + j->Sosed->par[0].u;
                    v = v + j->Sosed->par[0].v;
                    w = w + j->Sosed->par[0].w;
                    bx = bx + j->Sosed->par[0].bx;
                    by = by + j->Sosed->par[0].by;
                    bz = bz + j->Sosed->par[0].bz;
                }
            }

            if (kk == 0)
            {
                cout << "ERROR  3844  init mgd" << endl;
                i->par[0].ro = 1.0;
                i->par[0].p = 1.0 / (ggg);
                i->par[0].u = -M_inf; //-1.0;
                i->par[0].v = 0.0;
                i->par[0].w = 0.0;
                i->par[0].bx = 0.0; // -betta * cos(0.5235);
                i->par[0].by = 0.0; // -betta * cos(0.5235);
                i->par[0].bz = 0.0;
            }
            else
            {
                i->par[0].ro = ro / kk;
                i->par[0].p = p / kk;
                i->par[0].u = u / kk;
                i->par[0].v = v / kk;
                i->par[0].w = w / kk;
                i->par[0].bx = bx / kk;
                i->par[0].by = by / kk;
                i->par[0].bz = bz / kk;
            }
        }

    }

    for (auto& i : this->All_Cell)
    {
        i->i_init_mgd = false;
    }

    this->Calc_normal();
    this->Calc_normal2();
    this->Calc_normal();
}

void Setka::Go_MHD(int times)
{
    cout << "Go_MHD   " << times << endl;
    int now1 = 1;
    int now2 = 0;
    double T[2];
    T[0] = T[1] = 0.00000001;
    mutex mut;

    for (int st = 0; st <= times; st++)
    {
        if (st % 100 == 0 && st > 0)
        {
            cout << "Go_MHD    " << st << " " << T[now2] << endl;
        }
        now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
        now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
        T[now2] = 100000000;

#pragma omp parallel for
        for(int ii = 0; ii < this->All_Cell.size(); ii++)
        {
            auto K = this->All_Cell[ii];
            if (K->mgd_ == false)
            {
                K->par[now2] = K->par[now1];                                                   
                continue;
            }

            double x, y, z, r;
            K->Center->get(x, y, z);
            r = sqrt(kv(x) + kv(y) + kv(z));
            double P[8] = { 0.0 };
            P[0] = P[1] = P[2] = P[3] = P[4] = P[5] = P[6] = P[7] = 0.0;
            double Potok[10] = { 0.0 };
            Potok[0] = Potok[1] = Potok[2] = Potok[3] = Potok[4] = Potok[5] = Potok[6] = Potok[7] = Potok[8] = Potok[9] = 0.0;
            double tmin = 1000;
			double S = 0.0; // Площадь грани
            double ro = K->par[now1].ro;
            double p = K->par[now1].p;
            double vx = K->par[now1].u;
            double vy = K->par[now1].v;
            double vz = K->par[now1].w;
            double bx = K->par[now1].bx;
            double by = K->par[now1].by;
            double bz = K->par[now1].bz;
            double Volume = 0.0;
            K->get_Volume(Volume);

            double ro2, p2, vx2, vy2, vz2, bx2, by2, bz2, sks, dist;

            double n1, n2, n3, nn;
            double x2, y2, z2;

            for (auto& i : K->Grans)
            {
                int mm = 3;                                                          // Метод  ---------------------------------------
                double ro2, p2, vx2, vy2, vz2, bx2, by2, bz2;
                
                if (r < 1.5)
                {
                    mm = 2;
                }
                
                i->Sosed->Center->get(x2, y2, z2);
                if (i->Sosed->type == C_base)
                {
                    ro2 = i->Sosed->par[now1].ro;
                    p2 = i->Sosed->par[now1].p;
                    vx2 = i->Sosed->par[now1].u;
                    vy2 = i->Sosed->par[now1].v;
                    vz2 = i->Sosed->par[now1].w;
                    bx2 = i->Sosed->par[now1].bx;
                    by2 = i->Sosed->par[now1].by;
                    bz2 = i->Sosed->par[now1].bz;
                }
                else
                {
                    ro2 = ro;
                    p2 = p;
                    vx2 = vx;
                    vy2 = vy;
                    vz2 = vz;
                    bx2 = bx;
                    by2 = by;
                    bz2 = bz;
                    if (i->Sosed->type == C_wall_x_min)
                    {
                        vx2 = -M_inf;
                    }

                    if (i->Sosed->type != C_sphere)
                    {
                        mm = 1;
                    }

                }

                S = i->S;
                n1 = i->n1;
                n2 = i->n2;
                n3 = i->n3;
                dist = sqrt(kv(x2 - x) + kv(y2 - y) + kv(z2 - z)) / 2.0;
                sks = n1 * (bx + bx2) / 2.0 + n2 * (by + by2) / 2.0 + n3 * (bz + bz2) / 2.0;
                Potok[8] = Potok[8] + sks * S;
                /*if (r < 1.4)
                {
                    mm = 1;
                }*/
                double PQ;
                double Vc;
                if (mm == 3 || mm == 2)
                {
                    tmin = min(tmin, HLLDQ_Korolkov(ro, 1.0, p, vx, vy, vz, bx, by, bz, ro2, 1.0, p2, vx2, vy2, vz2, bx2, by2, bz2, P, PQ, n1, n2, n3, dist, mm, Vc));
                }
                else
                {
                    tmin = min(tmin, HLLD_Alexashov(ro, p, vx, vy, vz, bx, by, bz, ro2, p2, vx2, vy2, vz2, bx2, by2, bz2, P, n1, n2, n3, dist, mm, Vc));
                }
                for (int k = 0; k < 8; k++)  // Суммируем все потоки в ячейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }

            }

            double ro3, p3, u3, v3, w3, bx3, by3, bz3;

            ro3 = ro - T[now1] * Potok[0] / Volume;
            if (ro3 <= 0.0)
            {
                printf("ERROR -  dssdbfhfshjskfutytqqazz\n");
                printf("%lf, %lf, %lf, %lf\n", x, y, z, ro3);
                ro3 = 0.0001;
            }
            u3 = (ro * vx - T[now1] * (Potok[1] + (bx / cpi4) * Potok[8]) / Volume) / ro3;
            v3 = (ro * vy - T[now1] * (Potok[2] + (by / cpi4) * Potok[8]) / Volume) / ro3;
            w3 = (ro * vz - T[now1] * (Potok[3] + (bz / cpi4) * Potok[8]) / Volume) / ro3;
            bx3 = (bx - T[now1] * (Potok[4] + vx * Potok[8]) / Volume);
            by3 = (by - T[now1] * (Potok[5] + vy * Potok[8]) / Volume);
            bz3 = (bz - T[now1] * (Potok[6] + vz * Potok[8]) / Volume);
            p3 = ((U8(ro, p, vx, vy, vz, bx, by, bz) - T[now1] * (Potok[7] + (skk(vx, vy, vz, bx, by, bz) / cpi4) * Potok[8])//
                / Volume) - 0.5 * ro3 * kvv(u3, v3, w3) - kvv(bx3, by3, bz3) / cpi8) * (ggg - 1.0);

            if (p3 <= 0)
            {
                p3 = 0.000001;
            }

            K->par[now2].ro = ro3;
            K->par[now2].p = p3;
            K->par[now2].u = u3;
            K->par[now2].v = v3;
            K->par[now2].w = w3;
            K->par[now2].bx = bx3;
            K->par[now2].by = by3;
            K->par[now2].bz = bz3;


            if (tmin < T[now2])
            {
                mut.lock();
                T[now2] = tmin;
                mut.unlock();
            }

        }

    }
}

void Setka::Culc_couple(int now1, const double& time)
{
#pragma omp parallel for
    for (auto& i : this->All_Couple)
    {
        auto K = i->A1;
        double ro = K->par[now1].ro;
        double p = K->par[now1].p;
        double vx = K->par[now1].u;
        double vy = K->par[now1].v;
        double vz = K->par[now1].w;
        double bx = K->par[now1].bx;
        double by = K->par[now1].by;
        double bz = K->par[now1].bz;
        double ro2, p2, vx2, vy2, vz2, bx2, by2, bz2;
        K = i->A2;
        ro2 = K->par[now1].ro;
        p2 = K->par[now1].p;
        vx2 = K->par[now1].u;
        vy2 = K->par[now1].v;
        vz2 = K->par[now1].w;
        bx2 = K->par[now1].bx;
        by2 = K->par[now1].by;
        bz2 = K->par[now1].bz;
        double n1, n2, n3;
        n1 = i->n1;
        n2 = i->n2;
        n3 = i->n3;

        double P[8] = { 0.0 };
        P[0] = P[1] = P[2] = P[3] = P[4] = P[5] = P[6] = P[7] = 0.0;
        double PQ = 0.0;
        double Vc;

        HLLDQ_Korolkov(ro, 1.0, p, vx, vy, vz, bx, by, bz, ro2, 1.0, p2, vx2, vy2, vz2, bx2, by2, bz2, P, PQ, n1, n2, n3, 1.0, 3, Vc, 0.0);
        
        Vc = Vc * time;
        //cout << Vc << endl;
        i->A1->Center->move(n1 * Vc, n2 * Vc, n3 * Vc);  // Сразу передвинули пару
        i->A2->Center->move(n1 * Vc, n2 * Vc, n3 * Vc);
    }
}

void Setka::Start_MHD(int times)
{
    // Основная функция
    // Здесь предумострены все движения пар, какие функции для этого нужно вызывать.
    // Также она максимально оптимизирована должна быть.

    // Важный момент 1:
    // Для того, чтобы найти скорость грани нужно знать где грань была в предыдущий момент времени
    // Поэтому в каждой ячейке есть список старых граней
    // В конце этой функции в каждой ячейке удаляется список старых-старых граней и он обновляется новыми для последующего использования

    this->Initialization_do_MHD();

    cout << "Start_MHD   " << times << endl;
    int now1 = 1;
    int now2 = 0;
    double T[2];
    T[0] = T[1] = 0.00000001;
    mutex mut;
    bool bkl = true;

    for (int st = 0; st < times; st++)
    {
        if (st % 5 == 0)
        {
            cout << "Start_MHD    " << st << " " << T[now2] << endl;
        }
        now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
        now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
        T[now2] = 100000000;

       
        
        this->Culc_couple(now1, T[now1]);    // Считаем движение пар из задачи о распаде разрыва и двигаем их!!!
        this->Reconstruct_medium3(true);
        this->move_par();
        this->Calc_normal();
        bkl = true;
        do
        {
            bkl = this->Reconstruct_medium3(true);
            //cout << "bkl = " << bkl << endl;
            if (st % 5 == 0)
            {
                cout << "bkl = " << bkl << endl;
            }
        } while (bkl == false);

        //this->Reconstruct_medium3(true);
        //this->Reconstruct_medium3(true);

        // Важно! Это перестроение должно обязатьльно быть, так как во время него создаются новые грани (физически новые переменные)
        // а старые грани в конце функции очистятся (указатели сотрутся)
        // а если мы не обновим грани, то они просто удаляется
        //cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;

#pragma omp parallel for
        for (int ii = 0; ii < this->All_Cell.size(); ii++)
        {
            auto K = this->All_Cell[ii];
            K->Candidates.clear();
            K->near_par = false;
            K->near_par_2 = false;
            K->i_boundary_1 = false;
            K->i_boundary_2 = false;
            K->i_include_candidate = false;
            if (K->include_ == false)
            {
                K->par[now2] = K->par[now1];
                continue;
            }

            /*if (K->couple_ == false)
            {
                K->reconstruct_ = false;
            }*/

            //cout << "ii = " << ii << "   type = " << K->type << endl;
            if (K->mgd_ == false)
            {
                K->par[now2] = K->par[now1];
                continue;
            }

            double x, y, z, r;
            double x_do, y_do, z_do;
            K->Center->get(x, y, z);
            K->Center->get_do(x_do, y_do, z_do);
            r = sqrt(kv(x) + kv(y) + kv(z));
            double P[8] = { 0.0 };
            P[0] = P[1] = P[2] = P[3] = P[4] = P[5] = P[6] = P[7] = 0.0;
            double Potok[10] = { 0.0 };
            Potok[0] = Potok[1] = Potok[2] = Potok[3] = Potok[4] = Potok[5] = Potok[6] = Potok[7] = Potok[8] = Potok[9] = 0.0;
            double tmin = 1000;
            double S = 0.0; // Площадь грани
            double ro = K->par[now1].ro;
            double p = K->par[now1].p;
            double vx = K->par[now1].u;
            double vy = K->par[now1].v;
            double vz = K->par[now1].w;
            double bx = K->par[now1].bx;
            double by = K->par[now1].by;
            double bz = K->par[now1].bz;
            double Volume = 0.0;
            double Volume_do = 0.0;
            K->get_Volume(Volume);
            K->get_Volume_do(Volume_do);

            double ro2, p2, vx2, vy2, vz2, bx2, by2, bz2, sks, dist;

            double n1, n2, n3, nn;
            double x2, y2, z2;
            double x2_do, y2_do, z2_do;

            for (auto& i : K->Grans)
            {
                int mm = 3;                                                          // Метод  ---------------------------------------
                double ro2, p2, vx2, vy2, vz2, bx2, by2, bz2;

                if (r < 1.5)
                {
                    mm = 2;
                }

                i->Sosed->Center->get(x2, y2, z2);
                i->Sosed->Center->get_do(x2_do, y2_do, z2_do);
                if (i->Sosed->type == C_base)
                {
                    if (i->Sosed->couple_ == true)
                    {
                        K->near_par = true;
                    }
                    ro2 = i->Sosed->par[now1].ro;
                    p2 = i->Sosed->par[now1].p;
                    vx2 = i->Sosed->par[now1].u;
                    vy2 = i->Sosed->par[now1].v;
                    vz2 = i->Sosed->par[now1].w;
                    bx2 = i->Sosed->par[now1].bx;
                    by2 = i->Sosed->par[now1].by;
                    bz2 = i->Sosed->par[now1].bz;
                }
                else
                {
                    K->i_boundary_1 = true;
                    ro2 = ro;
                    p2 = p;
                    vx2 = vx;
                    vy2 = vy;
                    vz2 = vz;
                    bx2 = bx;
                    by2 = by;
                    bz2 = bz;
                    if (i->Sosed->type == C_wall_x_min)
                    {
                        vx2 = -M_inf;
                    }

                    if (i->Sosed->type != C_sphere)
                    {
                        mm = 1;
                    }

                }

                S = i->S;
                n1 = i->n1;
                n2 = i->n2;
                n3 = i->n3;
                dist = sqrt(kv(x2 - x) + kv(y2 - y) + kv(z2 - z)) / 2.0;
                sks = n1 * (bx + bx2) / 2.0 + n2 * (by + by2) / 2.0 + n3 * (bz + bz2) / 2.0;
                Potok[8] = Potok[8] + sks * S;

                if (K->couple_ == false && i->Sosed->couple_ == true && dist < 0.05)
                {
                    K->i_include_candidate = true;                    // Кандидат на отлючение
                }
                // Вычисляем скорость движения грани
                double wv = 0.0;
                if (i->Sosed->type == C_base)
                {
                    // Старый вариант с дополнительным списком граней
                    /*if (K->Grans_do.size() > 0)
                    {
                        if (K->Grans_do.count(i->Sosed->number) == 1)
                        {
                            double dx, dy, dz;
                            dx = (i->c1 - K->Grans_do[i->Sosed->number]->c1) / T[now1];
                            dy = (i->c2 - K->Grans_do[i->Sosed->number]->c2) / T[now1];
                            dz = (i->c3 - K->Grans_do[i->Sosed->number]->c3) / T[now1];
                            wv = n1 * dx + n2 * dy + n3 * dz;
                        }
                    }*/
                    double dx, dy, dz;
                    dx = 0.5 * (x + x2 - x_do - x2_do);
                    dy = 0.5 * (y + y2 - y_do - y2_do);
                    dz = 0.5 * (z + z2 - z_do - z2_do);
                    wv = n1 * dx + n2 * dy + n3 * dz;

                }

                //cout << wv << endl;

                /*if (r < 1.4)
                {
                    mm = 1;
                }*/
                double PQ;
                if (mm == 3 || mm == 2)
                {
                    tmin = min(tmin, HLLDQ_Korolkov(ro, 1.0, p, vx, vy, vz, bx, by, bz, ro2, 1.0, p2, vx2, vy2, vz2, bx2, by2, bz2, P, PQ, n1, n2, n3, dist, mm, wv));
                }
                else
                {
                    tmin = min(tmin, HLLD_Alexashov(ro, p, vx, vy, vz, bx, by, bz, ro2, p2, vx2, vy2, vz2, bx2, by2, bz2, P, n1, n2, n3, dist, mm, wv));
                }
                for (int k = 0; k < 8; k++)  // Суммируем все потоки в ячейке
                {
                    Potok[k] = Potok[k] + P[k] * S;
                }

            }

            double ro3, p3, u3, v3, w3, bx3, by3, bz3;

            ro3 = ro * Volume_do / Volume - T[now1] * Potok[0] / Volume;
            if (ro3 <= 0.0)
            {
                printf("ERROR -  dssdbfhfshjskfutytqqazz\n");
                printf("x = %lf, y = %lf, z = %lf, ro3 = %lf, V_do = %lf, V = %lf, Potok = %lf, T = %lf\n", x, //
                    y, z, ro3, Volume_do, Volume, Potok[0], T[now1]);
                exit(-1);
                ro3 = 0.0001;
            }
            u3 = (ro * vx * Volume_do / Volume - T[now1] * (Potok[1] + (bx / cpi4) * Potok[8]) / Volume) / ro3;
            v3 = (ro * vy * Volume_do / Volume - T[now1] * (Potok[2] + (by / cpi4) * Potok[8]) / Volume) / ro3;
            w3 = (ro * vz * Volume_do / Volume - T[now1] * (Potok[3] + (bz / cpi4) * Potok[8]) / Volume) / ro3;
            bx3 = (bx * Volume_do / Volume - T[now1] * (Potok[4] + vx * Potok[8]) / Volume);
            by3 = (by * Volume_do / Volume - T[now1] * (Potok[5] + vy * Potok[8]) / Volume);
            bz3 = (bz * Volume_do / Volume - T[now1] * (Potok[6] + vz * Potok[8]) / Volume);
            p3 = ((U8(ro, p, vx, vy, vz, bx, by, bz) * Volume_do / Volume - T[now1] * (Potok[7] + (skk(vx, vy, vz, bx, by, bz) / cpi4) * Potok[8])//
                / Volume) - 0.5 * ro3 * kvv(u3, v3, w3) - kvv(bx3, by3, bz3) / cpi8) * (ggg - 1.0);

            if (p3 <= 0)
            {
                p3 = 0.000001;
            }

            K->par[now2].ro = ro3;
            K->par[now2].p = p3;
            K->par[now2].u = u3;
            K->par[now2].v = v3;
            K->par[now2].w = w3;
            K->par[now2].bx = bx3;
            K->par[now2].by = by3;
            K->par[now2].bz = bz3;


            if (tmin < T[now2])
            {
                mut.lock();
                T[now2] = tmin;
                mut.unlock();
            }

            K->Center->renew();
        }


        // Находим соседей второго уровня к граничным ячейкам
        for (auto& ii: this->All_Cell)
        {
            if (ii->include_ == true)
            {
                if (ii->couple_ == false && ii->near_par == false)
                {
                    for (auto& kk : ii->Grans)
                    {
                        if (kk->Sosed->near_par == true)
                        {
                            ii->near_par_2 = true;
                            break;
                        }
                    }
                }
            }
        }
        // Нужен также алгоритм пересмотра кандидатов на реконструкцию
        // + также нужен алгоритм для включения отключенных ячеек и выключения включённых
        
    }
}

void Setka::Rebuild(void)
{
    this->Initialization_do_MHD();
    for (auto& ii : this->All_Couple)
    {
        if (ii->A1->include_ == false)
        {
            continue;
        }
        ii->d_sosed = 0.0;
        double x, y, z;
        double a, b, c, d, e, f;
        double n1, n2, n3, nn;

        ii->get_centr(x, y, z);
        double x2, y2, z2;

        // Сначала ищем всех соседей
        for (auto& j : ii->A1->Grans)
        {
            if (j->Sosed->couple_ == true)
            {
                if (j->Sosed->Par != ii)
                {
                    j->Sosed->Par->get_centr(x2, y2, z2);
                    ii->d_sosed = max(ii->d_sosed, sqrt(kv(x - x2) + kv(y - y2) + kv(z - z2)));
                }
            }
        }

        for (auto& j : ii->A2->Grans)
        {
            if (j->Sosed->couple_ == true)
            {
                if (j->Sosed->Par != ii)
                {
                    j->Sosed->Par->get_centr(x2, y2, z2);
                    ii->d_sosed = max(ii->d_sosed, sqrt(kv(x - x2) + kv(y - y2) + kv(z - z2)));
                }
            }
        }
    }

    bool bkl = true;
    this->include_cells();
    this->Construct_start();
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);
    this->Reconstruct_medium2(true);
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);
    this->Reconstruct_medium2(true);
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);
    //this->Initialization_do_MHD();

    this->exclude_cells();
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);
    this->Reconstruct_medium2(true);
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);
    this->Reconstruct_medium2(true);
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);
    this->Initialization_do_MHD();
}

void Setka::Zapusk(void)
{
    this->Initialization_do_MHD();
    cout << " bbb  =  " << this->Reconstruct_medium3(true) << endl;
    bool bkl = true;
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);

    this->Rebuild();

    this->Initialization_do_MHD();
    this->Init();
    this->Go_MHD(5000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);

    this->Rebuild();

    this->Initialization_do_MHD();
    this->Init();
    cout << "SECTION 2" << endl;
    this->Go_MHD(5000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);

    this->Rebuild();

    this->Initialization_do_MHD();
    this->Init();
    cout << "SECTION 3" << endl;
    this->Go_MHD(5000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);

    this->Save_setka("vers_8");

    this->Rebuild();

    this->Initialization_do_MHD();
    this->Init();
    cout << "SECTION 4" << endl;
    this->Go_MHD(5000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);

    this->Rebuild();

    this->Initialization_do_MHD();
    this->Init();
    cout << "SECTION 5" << endl;
    this->Go_MHD(5000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);

    this->Save_setka("vers_9");

    this->Rebuild();

    this->Initialization_do_MHD();
    this->Init();
    cout << "SECTION 6" << endl;
    this->Go_MHD(5000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);

    this->Rebuild();

    this->Initialization_do_MHD();
    this->Init();
    cout << "SECTION 7" << endl;
    this->Go_MHD(5000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    this->Go_MHD(1000);
    this->Start_MHD(100);
    do
    {
        bkl = this->Reconstruct_medium4(true);
        cout << "reconstruct_4 = " << bkl << endl;
    } while (bkl == false);

    this->Save_setka("vers_10");

    this->Rebuild();

    this->Initialization_do_MHD();
    this->Init();
    cout << "SECTION 8" << endl;
    this->Go_MHD(5000);
    this->Start_MHD(300);

    //this->Rebuild();

    this->Initialization_do_MHD();
}

void Setka::Save_MHD(string name_f)
{
    ofstream fout;
    fout.open(name_f);
    double X, Y, Z;
    for (auto& i : this->All_Cell)
    {
        i->Center->get(X, Y, Z);
        fout << X << " " << Y << " " << Z << " " << i->par[0].ro << " " <<//
            i->par[0].p << " " << i->par[0].u << " " << i->par[0].v << " " << i->par[0].w << " " << //
            i->par[0].bx << " " << i->par[0].by << " " << i->par[0].bz << " " << endl;
    }
}

void Setka::Download_MHD(string name_f)
{
    cout << "Download_MHD  from  " << name_f << endl;
    ifstream fout;
    fout.open(name_f);
    double X, Y, Z;
    for (auto& i : this->All_Cell)
    {
        fout >> X >> Y >> Z >> i->par[0].ro >> i->par[0].p >> i->par[0].u >> i->par[0].v >> //
            i->par[0].w >> i->par[0].bx >> i->par[0].by >> i->par[0].bz;

        i->par[1] = i->par[0];
    }
}

void Setka::Init(void)
{

    double x, y, z, r;
    for (auto& i : this->All_Cell)
    {
        /*bool inner = false;
        bool intry = false;
        for (auto& j : i->Grans)
        {
            if (j->Sosed->type == C_sphere)
            {
                inner = true;
            }

            if (j->Sosed->type == C_wall_x_max)
            {
                intry = true;
            }
        }

        if (intry == true || inner == true)
        {
            i->mgd_ = false;
        }*/

        if (i->mgd_ == true)
        {
            continue;
        }

        i->Center->get(x, y, z);
        r = sqrt(kv(x) + kv(y) + kv(z));
        if (r < 1.3)
        {
            i->par[0].ro = 1.0 / (kv(chi) * kv(r));
            i->par[0].p = i->par[0].ro * kv(chi) / (ggg * kv(5.0));
            i->par[0].u = chi * x / r;
            i->par[0].v = chi * y / r;
            i->par[0].w = chi * z / r;
            double B = sqrt(4.0 * pi) / (MA * r);
            double the = acos(z / r);
            double AA, BB, CC;
            double BR = 0.0;    // Br
            this->dekard_skorost(x, y, z, BR, B * sin(the), 0.0, AA, BB, CC);
            i->par[0].bx = AA;
            i->par[0].by = BB;
            i->par[0].bz = CC;

            i->par[1] = i->par[0];
        }
        else
        {
            i->par[0].ro = 1.0;
            i->par[0].p = 1.0 / (ggg);
            i->par[0].u = -M_inf; //-1.0;
            i->par[0].v = 0.0;
            i->par[0].w = 0.0;
            i->par[0].bx = 0.0; // -betta * cos(0.5235);
            i->par[0].by = 0.0; // -betta * cos(0.5235);
            i->par[0].bz = 0.0;

            i->par[1] = i->par[0];
        }

    }
}

void Setka::dekard_skorost(double x, double y, double z, double Vr, double Vphi, double Vtheta, double& Vx, double& Vy, double& Vz)
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

double Setka::polar_angle(double x, double y)
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

double Setka::HLLD_Alexashov(const double& ro_L, const double& p_L, const double& v1_L, const double& v2_L, const double& v3_L,//
    const double& Bx_L, const double& By_L, const double& Bz_L, const double& ro_R, const double& p_R, const double& v1_R, const double& v2_R, const double& v3_R,//
    const double& Bx_R, const double& By_R, const double& Bz_R, double* P, const double& n1, const double& n2, const double& n3, const double& rad, int metod, //
    double& Vcon, const double& wv)
{
    //cout << "HLLD_Alexashov" << endl;
    int x0 = 0, x1 = 1, x2 = 2;
    double aco[3][3];
    int n_state = metod;
    //c-------  n_state=0   - one speed LAX
    //c-------  n_state=1   - two speed LAX (HLL,(Harten-Lax-van-Leer))
    //c-------  n_state=2   - two-state (3 speed) HLLC (Contact Discontinuity)
    //c-------  n_state=3   - multi-state (5 speed) HLLD (All Discontinuity)

    double FR[8], FL[8], dq[8];
    double FW[8], UL[8], UZ[8], UR[8];
    double UZL[8], UZR[8];
    double UZZL[8], UZZR[8];
    double vL[3], vR[3], bL[3], bR[3];
    double vzL[3], vzR[3], bzL[3], bzR[3];
    double vzzL[3], vzzR[3], bzzL[3], bzzR[3];
    double qv[3], qb[3];



    //double wv = 0.0;
    int n_disco = 0; // Для определения скоростей характеристик  1


    double r1 = ro_L;
    double u1 = v1_L;
    double v1 = v2_L;
    double w1 = v3_L;
    double p1 = p_L;
    double bx1 = Bx_L / spi4;
    double by1 = By_L / spi4;
    double bz1 = Bz_L / spi4;


    double r2 = ro_R;
    double u2 = v1_R;
    double v2 = v2_R;
    double w2 = v3_R;
    double p2 = p_R;
    double bx2 = Bx_R / spi4;
    double by2 = By_R / spi4;
    double bz2 = Bz_R / spi4;

    double ro = (r2 + r1) / x2;
    double au = (u2 + u1) / x2;
    double av = (v2 + v1) / x2;
    double aw = (w2 + w1) / x2;
    double ap = (p2 + p1) / x2;
    double abx = (bx2 + bx1) / x2;
    double aby = (by2 + by1) / x2;
    double abz = (bz2 + bz1) / x2;

    double al = n1;
    double be = n2;
    double ge = n3;

    double bk = abx * al + aby * be + abz * ge;
    double b2 = kv(abx) + kv(aby) + kv(abz);

    double  d = b2 - kv(bk);
    aco[0][0] = al;
    aco[1][0] = be;
    aco[2][0] = ge;

    if (d > 0.00001)
    {
        d = sqrt(d);
        aco[0][1] = (abx - bk * al) / d;
        aco[1][1] = (aby - bk * be) / d;
        aco[2][1] = (abz - bk * ge) / d;
        aco[0][2] = (aby * ge - abz * be) / d;
        aco[1][2] = (abz * al - abx * ge) / d;
        aco[2][2] = (abx * be - aby * al) / d;
    }
    else
    {
        double aix, aiy, aiz;
        if ((fabs(al) < fabs(be)) && (fabs(al) < fabs(ge)))
        {
            aix = x1;
            aiy = x0;
            aiz = x0;
        }
        else if (fabs(be) < fabs(ge))
        {
            aix = x0;
            aiy = x1;
            aiz = x0;
        }
        else
        {
            aix = x0;
            aiy = x0;
            aiz = x1;
        }
        double aik = aix * al + aiy * be + aiz * ge;
        d = sqrt(x1 - kv(aik));
        aco[0][1] = (aix - aik * al) / d;
        aco[1][1] = (aiy - aik * be) / d;
        aco[2][1] = (aiz - aik * ge) / d;
        aco[0][2] = (aiy * ge - aiz * be) / d;
        aco[1][2] = (aiz * al - aix * ge) / d;
        aco[2][2] = (aix * be - aiy * al) / d;
    }

    aco[0][0] = al;
    aco[1][0] = be;
    aco[2][0] = ge;

    //if (fabs(skk(aco[0][0], aco[1][0], aco[2][0], aco[0][1], aco[1][1], aco[2][1])) > 0.000001 || //
    //    fabs(skk(aco[0][0], aco[1][0], aco[2][0], aco[0][2], aco[1][2], aco[2][2])) > 0.000001 || //
    //    fabs(skk(aco[0][2], aco[1][2], aco[2][2], aco[0][1], aco[1][1], aco[2][1])) > 0.000001 || //
    //    fabs(kvv(aco[0][0], aco[1][0], aco[2][0]) - 1.0) > 0.000001 || fabs(kvv(aco[0][1], aco[1][1], aco[2][1]) - 1.0) > 0.000001 ||//
    //    fabs(kvv(aco[0][2], aco[1][2], aco[2][2]) - 1.0) > 0.000001)
    //{
    //    printf("Ne normal  174fdcdsaxes\n");
    //}


    for (int i = 0; i < 3; i++)
    {
        vL[i] = aco[0][i] * u1 + aco[1][i] * v1 + aco[2][i] * w1;
        vR[i] = aco[0][i] * u2 + aco[1][i] * v2 + aco[2][i] * w2;
        bL[i] = aco[0][i] * bx1 + aco[1][i] * by1 + aco[2][i] * bz1;
        bR[i] = aco[0][i] * bx2 + aco[1][i] * by2 + aco[2][i] * bz2;
    }



    double aaL = bL[0] / sqrt(r1);
    double b2L = kv(bL[0]) + kv(bL[1]) + kv(bL[2]);
    double b21 = b2L / r1;
    double cL = sqrt(ga * p1 / r1);
    double qp = sqrt(b21 + cL * (cL + 2.0 * aaL));
    double qm = sqrt(b21 + cL * (cL - 2.0 * aaL));
    double cfL = (qp + qm) / x2;
    double ptL = p1 + b2L / x2;

    double aaR = bR[0] / sqrt(r2);
    double b2R = kv(bR[0]) + kv(bR[1]) + kv(bR[2]);
    double b22 = b2R / r2;
    double cR = sqrt(ga * p2 / r2);
    qp = sqrt(b22 + cR * (cR + 2.0 * aaR));
    qm = sqrt(b22 + cR * (cR - 2.0 * aaR));
    double cfR = (qp + qm) / x2;
    double ptR = p2 + b2R / x2;

    double aC = (aaL + aaR) / x2;
    double b2o = (b22 + b21) / x2;
    double cC = sqrt(ga * ap / ro);
    qp = sqrt(b2o + cC * (cC + x2 * aC));
    qm = sqrt(b2o + cC * (cC - x2 * aC));
    double cfC = (qp + qm) / x2;
    double vC1 = (vL[0] + vR[0]) / x2;

    double SL, SR;

    if (true)
    {
        SL = min(vL[0], vR[0]) - max(cfL, cfR);
        SR = max(vL[0], vR[0]) + max(cfL, cfR);
        //cout << "SL,R  =  " << SL << " " << SR << endl;
    }
    else if (n_disco == 1)
    {
        SL = min((vL[0] - cfL), (vC1 - cfC));
        SR = max((vR[0] + cfR), (vC1 + cfC));
    }
    else if (n_disco == 0)
    {
        SL = min((vL[0] - cfL), (vR[0] - cfR));
        SR = max((vL[0] + cfL), (vR[0] + cfR));
    }
    else if (n_disco == 2)
    {
        double SL_1 = min((vL[0] - cfL), (vC1 - cfC));
        double SR_1 = max((vR[0] + cfR), (vC1 + cfC));
        double SL_2 = min((vL[0] - cfL), (vR[0] - cfR));
        double SR_2 = max((vL[0] + cfL), (vR[0] + cfR));
        double oo = 0.75;
        double oo1 = 1.0 - oo;
        SL = oo * SL_1 + oo1 * SL_2;
        SR = oo * SR_1 + oo1 * SR_2;
    }

    //cout << "SL,  SR   =  " << std::setprecision(16) << SL << " " << SR << endl;


    double suR = SR - vR[0];
    double suL = SL - vL[0];
    double SM = (suR * r2 * vR[0] - ptR + ptL - suL * r1 * vL[0]) / (suR * r2 - suL * r1);
    Vcon = SM;
    //cout << "SM   =  " << std::setprecision(16) << SM << endl;
    //cout << "SM = " << SM << endl;

    // double dsl = SL;
    // double dsc = SM;
    // double dsp = SR;

    //if ((SR < SL) || (SL > SM) || (SR < SM))
    //{
    //    printf("ERROR -  254 fghrvtrgr\n");
    //    printf("%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",//
    //        vL[0], vR[0], cfL, cfR, ro_L, ro_R, p_L, p_R, suR, suL);
    //}

    double UU = max(fabs(SL), fabs(SR));
    double time = krit * rad / UU;

    double TR0, TL0;
    if (n_state == 0)
    {
        TR0 = fabs(vL[0] + vR[0]) / x2 + cfC;
        TL0 = -TR0;
        SR = TR0;
        SL = TL0;
    }


    double upt1 = (kv(u1) + kv(v1) + kv(w1)) / 2.0;
    double sbv1 = u1 * bx1 + v1 * by1 + w1 * bz1;

    double upt2 = (kv(u2) + kv(v2) + kv(w2)) / 2.0;
    double sbv2 = u2 * bx2 + v2 * by2 + w2 * bz2;

    double e1 = p1 / g1 + r1 * upt1 + b2L / x2;
    double e2 = p2 / g1 + r2 * upt2 + b2R / x2;

    FL[0] = r1 * vL[0];
    FL[1] = r1 * vL[0] * vL[0] + ptL - kv(bL[0]);
    FL[2] = r1 * vL[0] * vL[1] - bL[0] * bL[1];
    FL[3] = r1 * vL[0] * vL[2] - bL[0] * bL[2];
    FL[4] = (e1 + ptL) * vL[0] - bL[0] * sbv1;
    FL[5] = 0.0;
    FL[6] = vL[0] * bL[1] - vL[1] * bL[0];
    FL[7] = vL[0] * bL[2] - vL[2] * bL[0];

    FR[0] = r2 * vR[0];
    FR[1] = r2 * vR[0] * vR[0] + ptR - kv(bR[0]);
    FR[2] = r2 * vR[0] * vR[1] - bR[0] * bR[1];
    FR[3] = r2 * vR[0] * vR[2] - bR[0] * bR[2];
    FR[4] = (e2 + ptR) * vR[0] - bR[0] * sbv2;
    FR[5] = 0.0;
    FR[6] = vR[0] * bR[1] - vR[1] * bR[0];
    FR[7] = vR[0] * bR[2] - vR[2] * bR[0];

    UL[0] = r1;
    UL[4] = e1;
    UR[0] = r2;
    UR[4] = e2;

    for (int i = 0; i < 3; i++)
    {

        UL[i + 1] = r1 * vL[i];
        UL[i + 5] = bL[i];
        UR[i + 1] = r2 * vR[i];
        UR[i + 5] = bR[i];
    }

    for (int ik = 0; ik < 8; ik++)
    {
        UZ[ik] = (SR * UR[ik] - SL * UL[ik] + FL[ik] - FR[ik]) / (SR - SL);
    }


    if (n_state <= 1)
    {

        for (int ik = 0; ik < 8; ik++)
        {
            dq[ik] = UR[ik] - UL[ik];
        }



        double TL = SL;
        double TR = SR;
        if (SL > wv)
        {
            TL = 0.0;
            for (int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv * UL[ik];
            }
        }
        else if ((SL <= wv) && (wv <= SR))
        {
            for (int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv * UZ[ik];
            }
        }
        else if (SR < wv)
        {
            TR = 0.0;
            for (int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv * UR[ik];
            }
        }
        else
        {
            printf("ERROR  329 87732, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", r1, r2, p1, p2, al, be, ge);
        }


        double a = TR * TL;
        double b = TR - TL;

        P[0] = (TR * FL[0] - TL * FR[0] + a * dq[0]) / b - FW[0];
        P[4] = (TR * FL[4] - TL * FR[4] + a * dq[4]) / b - FW[4];

        for (int ik = 1; ik < 4; ik++)
        {
            qv[ik - 1] = (TR * FL[ik] - TL * FR[ik] + a * dq[ik]) / b - FW[ik];
        }

        for (int ik = 5; ik < 8; ik++)
        {
            qb[ik - 5] = (TR * FL[ik] - TL * FR[ik] + a * dq[ik]) / b - FW[ik];
        }


        // Улучшения для BN
        if (false)
        {
            double SN = max(fabs(SL), fabs(SR)) / 2.0;       // Дмитрий Борисович здесь не делит на 2

            double wbn = 0.0;
            if (wv >= SR)
            {
                wbn = wv * bR[0];
            }
            else if (wv <= SL)
            {
                wbn = wv * bL[0];
            }
            else
            {
                wbn = wv * (bL[0] + bR[0]) / x2;
            }

            qb[0] = -SN * (bR[0] - bL[0]) - wbn;
        }
        // у Дмитрия Борисовича это улучшение только для HLLD схемы

        for (int ik = 0; ik < 3; ik++)
        {
            P[ik + 1] = aco[ik][0] * qv[0] + aco[ik][1] * qv[1] + aco[ik][2] * qv[2];
            P[ik + 5] = aco[ik][0] * qb[0] + aco[ik][1] * qb[1] + aco[ik][2] * qb[2];
            P[ik + 5] = spi4 * P[ik + 5];
        }

        double SWAP = P[4];
        P[4] = P[5];
        P[5] = P[6];
        P[6] = P[7];
        P[7] = SWAP;
        return time;
    }
    if (n_state == 2)
    {
        double suRm = suR / (SR - SM);
        double suLm = suL / (SL - SM);
        double rzR = r2 * suRm;
        double rzL = r1 * suLm;
        vzR[0] = SM;
        vzL[0] = SM;
        double ptzR = ptR + r2 * suR * (SM - vR[0]);
        double ptzL = ptL + r1 * suL * (SM - vL[0]);
        double ptz = (ptzR + ptzL) / x2;
        bzR[0] = UZ[5];
        bzL[0] = UZ[5];

        vzR[1] = UZ[2] / UZ[0];
        vzR[2] = UZ[3] / UZ[0];
        vzL[1] = vzR[1];
        vzL[2] = vzR[2];

        vzR[1] = vR[1] + UZ[5] * (bR[1] - UZ[6]) / suR / r2;
        vzR[2] = vR[2] + UZ[5] * (bR[2] - UZ[7]) / suR / r2;
        vzL[1] = vL[1] + UZ[5] * (bL[1] - UZ[6]) / suL / r1;
        vzL[2] = vL[2] + UZ[5] * (bL[2] - UZ[7]) / suL / r1;

        bzR[1] = UZ[6];
        bzR[2] = UZ[7];
        bzL[1] = bzR[1];
        bzL[2] = bzR[2];

        double sbvz = (UZ[5] * UZ[1] + UZ[6] * UZ[2] + UZ[7] * UZ[3]) / UZ[0];

        double ezR = e2 * suRm + (ptz * SM - ptR * vR[0] + UZ[5] * (sbv2 - sbvz)) / (SR - SM);
        double ezL = e1 * suLm + (ptz * SM - ptL * vL[0] + UZ[5] * (sbv1 - sbvz)) / (SL - SM);

        if (false)//(fabs(UZ[5]) < epsb)
        {
            vzR[1] = vR[1];
            vzR[2] = vR[2];
            vzL[1] = vL[1];
            vzL[2] = vL[2];
            bzR[1] = bR[1] * suRm;
            bzR[2] = bR[2] * suRm;
            bzL[1] = bL[1] * suLm;
            bzL[2] = bL[2] * suLm;
        }

        UZL[0] = rzL;
        UZL[4] = ezL;
        UZR[0] = rzR;
        UZR[4] = ezR;

        for (int ik = 0; ik < 3; ik++)
        {
            UZL[ik + 1] = vzL[ik] * rzL;
            UZL[ik + 5] = bzL[ik];
            UZR[ik + 1] = vzR[ik] * rzR;
            UZR[ik + 5] = bzR[ik];
        }

        if (SL > wv)
        {
            P[0] = FL[0] - wv * UL[0];
            P[4] = FL[4] - wv * UL[4];
            for (int ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FL[ik] - wv * UL[ik];
            }

            for (int ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FL[ik] - wv * UL[ik];
            }
        }
        if (SL <= wv && SM >= wv)
        {
            P[0] = FL[0] + SL * (rzL - r1) - wv * UZL[0];
            P[4] = FL[4] + SL * (ezL - e1) - wv * UZL[4];
            for (int ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            }

            for (int ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            }
        }
        if (SM < wv && SR >= wv)
        {
            P[0] = FR[0] + SR * (rzR - r2) - wv * UZR[0];
            P[4] = FR[4] + SR * (ezR - e2) - wv * UZR[4];
            for (int ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
            }
            for (int ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
            }
        }

        if (SR < wv)
        {
            P[0] = FR[0] - wv * UR[0];
            P[4] = FR[4] - wv * UR[4];
            for (int ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FR[ik] - wv * UR[ik];
            }
            for (int ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FR[ik] - wv * UR[ik];
            }
        }

        // Улучшения для BN
        // у Дмитрия Борисовича это улучшение только для HLLD схемы
        if (true)
        {
            double SN = max(fabs(SL), fabs(SR));  // У Дмитрия Борисовича пополам не делится

            double wbn = 0.0;
            if (wv >= SR)
            {
                wbn = wv * bR[0];
            }
            else if (wv <= SL)
            {
                wbn = wv * bL[0];
            }
            else
            {
                wbn = wv * (bL[0] + bR[0]) / x2;
            }

            qb[0] = -SN * (bR[0] - bL[0]) - wbn;
        }

        for (int ik = 0; ik < 3; ik++)
        {
            P[ik + 1] = aco[ik][0] * qv[0] + aco[ik][1] * qv[1] + aco[ik][2] * qv[2];
            P[ik + 5] = aco[ik][0] * qb[0] + aco[ik][1] * qb[1] + aco[ik][2] * qb[2];
            P[ik + 5] = spi4 * P[ik + 5];
        }

        double SWAP = P[4];
        P[4] = P[5];
        P[5] = P[6];
        P[6] = P[7];
        P[7] = SWAP;
        return time;
    }
    if (n_state == 3)
    {

        double ptz = (suR * r2 * ptL - suL * r1 * ptR + r1 * r2 * suR * suL * (vR[0] - vL[0])) / (suR * r2 - suL * r1);

        vzL[0] = SM;
        vzR[0] = SM;
        vzzL[0] = SM;
        vzzR[0] = SM;
        double ptzL = ptz;
        double ptzR = ptz;
        double ptzzL = ptz;
        double ptzzR = ptz;

        double suRm = suR / (SR - SM);
        double suLm = suL / (SL - SM);
        double rzR = r2 * suRm;
        double rzL = r1 * suLm;

        double bn = UZ[5];
        double bn2 = bn * bn;
        bzL[0] = bn;
        bzR[0] = bn;
        bzzL[0] = bn;
        bzzR[0] = bn;

        double ttR = r2 * suR * (SR - SM) - bn2;
        double tvR, tbR, tvL, tbL;
        if (fabs(ttR) <= 0.00000001)
        {
            tvR = x0;
            tbR = x0;
        }
        else
        {
            tvR = (SM - vR[0]) / ttR;
            tbR = (r2 * suR * suR - bn2) / ttR;
        }

        double ttL = r1 * suL * (SL - SM) - bn2;
        if (fabs(ttL) <= 0.00000001)
        {
            tvL = x0;
            tbL = x0;
        }
        else
        {
            tvL = (SM - vL[0]) / ttL;
            tbL = (r1 * suL * suL - bn2) / ttL;
        }

        vzL[1] = vL[1] - bn * bL[1] * tvL;
        vzL[2] = vL[2] - bn * bL[2] * tvL;
        vzR[1] = vR[1] - bn * bR[1] * tvR;
        vzR[2] = vR[2] - bn * bR[2] * tvR;

        bzL[1] = bL[1] * tbL;
        bzL[2] = bL[2] * tbL;
        bzR[1] = bR[1] * tbR;
        bzR[2] = bR[2] * tbR;

        double sbvL = bzL[0] * vzL[0] + bzL[1] * vzL[1] + bzL[2] * vzL[2];
        double sbvR = bzR[0] * vzR[0] + bzR[1] * vzR[1] + bzR[2] * vzR[2];

        double ezR = e2 * suRm + (ptz * SM - ptR * vR[0] + bn * (sbv2 - sbvR)) / (SR - SM);
        double ezL = e1 * suLm + (ptz * SM - ptL * vL[0] + bn * (sbv1 - sbvL)) / (SL - SM);

        double rzzR = rzR;
        double rzzL = rzL;
        double rzRs = sqrt(rzR);
        double rzLs = sqrt(rzL);
        double rzss = rzRs + rzLs;
        double rzps = rzRs * rzLs;

        double SZL = SM - fabs(bn) / rzLs;
        double SZR = SM + fabs(bn) / rzRs;

        int ibn = 0;
        double sbn;
        if (fabs(bn) > 0.000001)
        {
            sbn = 1.0 * sgn(bn);
            ibn = 1;
        }
        else
        {
            sbn = 0.0;
            ibn = 0;
            SZL = SM;
            SZR = SM;
        }

        //cout << "SZ  =  " <<  SZL << " " << SZR << endl;

        vzzL[1] = (rzLs * vzL[1] + rzRs * vzR[1] + sbn * (bzR[1] - bzL[1])) / rzss;
        vzzL[2] = (rzLs * vzL[2] + rzRs * vzR[2] + sbn * (bzR[2] - bzL[2])) / rzss;
        vzzR[1] = vzzL[1];
        vzzR[2] = vzzL[2];

        bzzL[1] = (rzLs * bzR[1] + rzRs * bzL[1] + sbn * rzps * (vzR[1] - vzL[1])) / rzss;
        bzzL[2] = (rzLs * bzR[2] + rzRs * bzL[2] + sbn * rzps * (vzR[2] - vzL[2])) / rzss;
        bzzR[1] = bzzL[1];
        bzzR[2] = bzzL[2];

        double sbzz = bzzL[0] * vzzL[0] + bzzL[1] * vzzL[1] + bzzL[2] * vzzL[2];

        double ezzR = ezR + rzRs * sbn * (sbvR - sbzz);
        double ezzL = ezL - rzLs * sbn * (sbvL - sbzz);

        UZL[0] = rzL;
        UZL[4] = ezL;
        UZR[0] = rzR;
        UZR[4] = ezR;

        for (int ik = 0; ik < 3; ik++)
        {
            UZL[ik + 1] = vzL[ik] * rzL;
            UZL[ik + 5] = bzL[ik];
            UZR[ik + 1] = vzR[ik] * rzR;
            UZR[ik + 5] = bzR[ik];
        }



        UZZL[0] = rzzL;
        UZZL[4] = ezzL;
        UZZR[0] = rzzR;
        UZZR[4] = ezzR;
        for (int ik = 0; ik < 3; ik++)
        {
            UZZL[ik + 1] = vzzL[ik] * rzzL;
            UZZL[ik + 5] = bzzL[ik];
            UZZR[ik + 1] = vzzR[ik] * rzzR;
            UZZR[ik + 5] = bzzR[ik];
        }



        if (SL > wv)
        {
            //cout << "1 zone" << endl;
            P[0] = FL[0] - wv * UL[0];
            P[4] = FL[4] - wv * UL[4];
            for (int ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FL[ik] - wv * UL[ik];
            }

            for (int ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FL[ik] - wv * UL[ik];
            }
        }

        if ((SL <= wv) && (SZL >= wv))
        {
            //cout << "2 zone" << endl;
            int ik = 0;
            P[ik] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            ik = 4;
            P[ik] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            for (int ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            }

            for (int ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            }
        }
        //c------ FZZ
        if (ibn == 1)
        {

            if ((SZL <= wv) && (SM >= wv))
            {
                //cout << "3 zone" << endl;
                int ik = 0;
                P[ik] = FL[ik] + SZL * (UZZL[ik] - UZL[ik]) + SL * (UZL[ik] - UL[ik]) - wv * UZZL[ik];
                ik = 4;
                P[ik] = FL[ik] + SZL * (UZZL[ik] - UZL[ik]) + SL * (UZL[ik] - UL[ik]) - wv * UZZL[ik];
                for (int ik = 1; ik < 4; ik++)
                {
                    qv[ik - 1] = FL[ik] + SZL * (UZZL[ik] - UZL[ik]) + SL * (UZL[ik] - UL[ik]) - wv * UZZL[ik];
                }

                for (int ik = 5; ik < 8; ik++)
                {
                    qb[ik - 5] = FL[ik] + SZL * (UZZL[ik] - UZL[ik]) + SL * (UZL[ik] - UL[ik]) - wv * UZZL[ik];
                }
            }

            if ((SM <= wv) && (SZR >= wv))
            {
                //cout << "4 zone" << endl;
                int ik = 0;
                P[ik] = FR[ik] + SZR * (UZZR[ik] - UZR[ik]) + SR * (UZR[ik] - UR[ik]) - wv * UZZR[ik];
                ik = 4;
                P[ik] = FR[ik] + SZR * (UZZR[ik] - UZR[ik]) + SR * (UZR[ik] - UR[ik]) - wv * UZZR[ik];
                for (int ik = 1; ik < 4; ik++)
                {
                    qv[ik - 1] = FR[ik] + SZR * (UZZR[ik] - UZR[ik]) + SR * (UZR[ik] - UR[ik]) - wv * UZZR[ik];
                }

                for (int ik = 5; ik < 8; ik++)
                {
                    qb[ik - 5] = FR[ik] + SZR * (UZZR[ik] - UZR[ik]) + SR * (UZR[ik] - UR[ik]) - wv * UZZR[ik];
                }
            }

        }
        //c------ 
        if ((SZR <= wv) && (SR >= wv))
        {
            //cout << "5 zone" << endl;
            int ik = 0;
            P[ik] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
            ik = 4;
            P[ik] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];

            for (int ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
            }

            for (int ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
            }
        }

        if (SR < wv)
        {
            //cout << "6 zone" << endl;
            P[0] = FR[0] - wv * UR[0];
            P[4] = FR[4] - wv * UR[4];
            for (int ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FR[ik] - wv * UR[ik];
            }


            for (int ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FR[ik] - wv * UR[ik];
            }
        }


        //c----- Bn
        //double SN = max(fabs(SL), fabs(SR));

        double SN = max(fabs(SL), fabs(SR)) / 2.0;    // У Дмитрия Борисовича пополам не делится

        double wbn = 0.0;
        if (wv >= SR)
        {
            wbn = wv * bR[0];
        }
        else if (wv <= SL)
        {
            wbn = wv * bL[0];
        }
        else
        {
            wbn = wv * (bL[0] + bR[0]) / x2;
        }

        qb[0] = -SN * (bR[0] - bL[0]) - wbn;

        //c-----


        for (int ik = 0; ik < 3; ik++)
        {
            P[ik + 1] = aco[ik][0] * qv[0] + aco[ik][1] * qv[1] + aco[ik][2] * qv[2];
            P[ik + 5] = aco[ik][0] * qb[0] + aco[ik][1] * qb[1] + aco[ik][2] * qb[2];
            P[ik + 5] = spi4 * P[ik + 5];
        }

        double SWAP = P[4];
        P[4] = P[5];
        P[5] = P[6];
        P[6] = P[7];
        P[7] = SWAP;
        return time;
    }

    return time;
}

double Setka::HLLDQ_Korolkov(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L, const double& v3_L,//
    const double& Bx_L, const double& By_L, const double& Bz_L, const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& v3_R,//
    const double& Bx_R, const double& By_R, const double& Bz_R, double* P, double& PQ, const double& n1, const double& n2,//
    const double& n3, const double& rad, int metod, double& Vcon, const double& wv)
{// Не работает, если скорость грани не нулевая
    //cout << "HLLDQ_Korolkov " << endl;
    double bx_L = Bx_L / spi4;
    double by_L = By_L / spi4;
    double bz_L = Bz_L / spi4;

    double bx_R = Bx_R / spi4;
    double by_R = By_R / spi4;
    double bz_R = Bz_R / spi4;

    double t1 = 0.0;
    double t2 = 0.0;
    double t3 = 0.0;

    double m1 = 0.0;
    double m2 = 0.0;
    double m3 = 0.0;

    double abx = (bx_R + bx_L) / 2.0;
    double aby = (by_R + by_L) / 2.0;
    double abz = (bz_R + bz_L) / 2.0;

    double al = n1;
    double be = n2;
    double ge = n3;

    double bk = abx * al + aby * be + abz * ge;
    double b2 = kv(abx) + kv(aby) + kv(abz);

    double  d = b2 - kv(bk);

    if (d > 0.00001)
    {
        d = sqrt(d);
        t1 = (abx - bk * al) / d;
        t2 = (aby - bk * be) / d;
        t3 = (abz - bk * ge) / d;
        m1 = (aby * ge - abz * be) / d;
        m2 = (abz * al - abx * ge) / d;
        m3 = (abx * be - aby * al) / d;
    }
    else
    {
        double aix, aiy, aiz;
        if ((fabs(al) < fabs(be)) && (fabs(al) < fabs(ge)))
        {
            aix = 1.0;
            aiy = 0.0;
            aiz = 0.0;
        }
        else if (fabs(be) < fabs(ge))
        {
            aix = 0.0;
            aiy = 1.0;
            aiz = 0.0;
        }
        else
        {
            aix = 0.0;
            aiy = 0.0;
            aiz = 1.0;
        }
        double aik = aix * al + aiy * be + aiz * ge;
        d = sqrt(1.0 - kv(aik));
        t1 = (aix - aik * al) / d;
        t2 = (aiy - aik * be) / d;
        t3 = (aiz - aik * ge) / d;
        m1 = (aiy * ge - aiz * be) / d;
        m2 = (aiz * al - aix * ge) / d;
        m3 = (aix * be - aiy * al) / d;
    }

    double u1, v1, w1, u2, v2, w2;
    u1 = v1_L * n1 + v2_L * n2 + v3_L * n3;
    v1 = v1_L * t1 + v2_L * t2 + v3_L * t3;
    w1 = v1_L * m1 + v2_L * m2 + v3_L * m3;
    u2 = v1_R * n1 + v2_R * n2 + v3_R * n3;
    v2 = v1_R * t1 + v2_R * t2 + v3_R * t3;
    w2 = v1_R * m1 + v2_R * m2 + v3_R * m3;

    double bn1, bt1, bm1, bn2, bt2, bm2;
    bn1 = bx_L * n1 + by_L * n2 + bz_L * n3;
    bt1 = bx_L * t1 + by_L * t2 + bz_L * t3;
    bm1 = bx_L * m1 + by_L * m2 + bz_L * m3;
    bn2 = bx_R * n1 + by_R * n2 + bz_R * n3;
    bt2 = bx_R * t1 + by_R * t2 + bz_R * t3;
    bm2 = bx_R * m1 + by_R * m2 + bz_R * m3;

    //cout << " = " << bt2 * bt2 + bm2 * bm2 << endl;

    double sqrtroL = sqrt(ro_L);
    double sqrtroR = sqrt(ro_R);
    double ca_L = bn1 / sqrtroL;
    double ca_R = bn2 / sqrtroR;
    double cL = sqrt(ggg * p_L / ro_L);
    double cR = sqrt(ggg * p_R / ro_R);

    double bb_L = kv(bx_L) + kv(by_L) + kv(bz_L);
    double bb_R = kv(bx_R) + kv(by_R) + kv(bz_R);

    double aL = (kv(bx_L) + kv(by_L) + kv(bz_L)) / ro_L;
    double aR = (kv(bx_L) + kv(by_L) + kv(bz_L)) / ro_L;

    double uu_L = (kv(v1_L) + kv(v2_L) + kv(v3_L)) / 2.0;
    double uu_R = (kv(v1_R) + kv(v2_R) + kv(v3_R)) / 2.0;

    double cfL = sqrt((ggg * p_L + bb_L + //
        sqrt(kv(ggg * p_L + bb_L) - 4.0 * ggg * p_L * kv(bn1))) / (2.0 * ro_L));
    double cfR = sqrt((ggg * p_R + bb_R + //
        sqrt(kv(ggg * p_R + bb_R) - 4.0 * ggg * p_R * kv(bn2))) / (2.0 * ro_R));


    double SL = min(u1, u2) - max(cfL, cfR);
    double SR = max(u1, u2) + max(cfL, cfR);
    //cout << "SL,  SR   =  " << std::setprecision(16) << SL << " " << SR << endl;
    double pTL = p_L + bb_L / 2.0;
    double pTR = p_R + bb_R / 2.0;

    double suR = (SR - u2);
    double suL = (SL - u1);

    double SM = (suR * ro_R * u2 - suL * ro_L * u1 - pTR + pTL) //
        / (suR * ro_R - suL * ro_L);

    Vcon = SM;

    //cout << "SM   =  " << std::setprecision(16) << SM << endl;

    double PTT = (suR * ro_R * pTL - suL * ro_L * pTR + ro_L * ro_R * suR * suL * (u2 - u1))//
        / (suR * ro_R - suL * ro_L);

    double UU = max(fabs(SL), fabs(SR));
    double time = krit * rad / UU;

    double FL[9], FR[9], UL[9], UR[9];

    double e1 = p_L / g1 + ro_L * uu_L + bb_L / 2.0;
    double e2 = p_R / g1 + ro_R * uu_R + bb_R / 2.0;


    FL[0] = ro_L * u1;
    FL[1] = ro_L * u1 * u1 + pTL - kv(bn1);
    FL[2] = ro_L * u1 * v1 - bn1 * bt1;
    FL[3] = ro_L * u1 * w1 - bn1 * bm1;
    FL[4] = (e1 + pTL) * u1 - bn1 * (u1 * bn1 + v1 * bt1 + w1 * bm1);
    //cout << uu_L << endl;
    FL[5] = 0.0;
    FL[6] = u1 * bt1 - v1 * bn1;
    FL[7] = u1 * bm1 - w1 * bn1;
    FL[8] = Q_L * u1;

    FR[0] = ro_R * u2;
    FR[1] = ro_R * u2 * u2 + pTR - kv(bn2);
    FR[2] = ro_R * u2 * v2 - bn2 * bt2;
    FR[3] = ro_R * u2 * w2 - bn2 * bm2;
    FR[4] = (e2 + pTR) * u2 - bn2 * (u2 * bn2 + v2 * bt2 + w2 * bm2);
    FR[5] = 0.0;
    FR[6] = u2 * bt2 - v2 * bn2;
    FR[7] = u2 * bm2 - w2 * bn2;
    FR[8] = Q_R * u2;

    UL[0] = ro_L;
    UL[1] = ro_L * u1;
    UL[2] = ro_L * v1;
    UL[3] = ro_L * w1;
    UL[4] = e1;
    UL[5] = bn1;
    UL[6] = bt1;
    UL[7] = bm1;
    UL[8] = Q_L;

    UR[0] = ro_R;
    UR[1] = ro_R * u2;
    UR[2] = ro_R * v2;
    UR[3] = ro_R * w2;
    UR[4] = e2;
    UR[5] = bn2;
    UR[6] = bt2;
    UR[7] = bm2;
    UR[8] = Q_R;

    double bn = (SR * UR[5] - SL * UL[5] + FL[5] - FR[5]) / (SR - SL);
    double bt = (SR * UR[6] - SL * UL[6] + FL[6] - FR[6]) / (SR - SL);
    double bm = (SR * UR[7] - SL * UL[7] + FL[7] - FR[7]) / (SR - SL);
    double bbn = bn * bn;

    double ro_LL = ro_L * (SL - u1) / (SL - SM);
    double ro_RR = ro_R * (SR - u2) / (SR - SM);
    double Q_LL = Q_L * (SL - u1) / (SL - SM);
    double Q_RR = Q_R * (SR - u2) / (SR - SM);

    if (metod == 2)   // HLLC  + mgd
    {
        double sbv1 = u1 * bn1 + v1 * bt1 + w1 * bm1;
        double sbv2 = u2 * bn2 + v2 * bt2 + w2 * bm2;

        double UZ0 = (SR * UR[0] - SL * UL[0] + FL[0] - FR[0]) / (SR - SL);
        double UZ1 = (SR * UR[1] - SL * UL[1] + FL[1] - FR[1]) / (SR - SL);
        double UZ2 = (SR * UR[2] - SL * UL[2] + FL[2] - FR[2]) / (SR - SL);
        double UZ3 = (SR * UR[3] - SL * UL[3] + FL[3] - FR[3]) / (SR - SL);
        double UZ4 = (SR * UR[4] - SL * UL[4] + FL[4] - FR[4]) / (SR - SL);
        double vzL, vzR, vLL, wLL, vRR, wRR, ppLR, btt1, bmm1, btt2, bmm2, ee1, ee2;


        double suRm = suR / (SR - SM);
        double suLm = suL / (SL - SM);
        double rzR = ro_R * suRm;
        double rzL = ro_L * suLm;

        double ptzR = pTR + ro_R * suR * (SM - u2);
        double ptzL = pTL + ro_L * suL * (SM - u1);
        double ptz = (ptzR + ptzL) / 2.0;


        vRR = UZ2 / UZ0;
        wRR = UZ3 / UZ0;
        vLL = vRR;
        wLL = wRR;

        /*vRR = v2 + bn * (bt2 - bt) / suR / ro_R;
        wRR = w2 + bn * (bm2 - bm) / suR / ro_R;
        vLL = v1 + bn * (bt1 - bt) / suL / ro_L;
        wLL = w1 + bn * (bm1 - bm) / suL / ro_L;*/

        btt2 = bt;
        bmm2 = bm;
        btt1 = btt2;
        bmm1 = bmm2;

        double sbvz = (bn * UZ1 + bt * UZ2 + bm * UZ3) / UZ0;

        ee2 = e2 * suRm + (ptz * SM - pTR * u2 + bn * (sbv2 - sbvz)) / (SR - SM);
        ee1 = e1 * suLm + (ptz * SM - pTL * u1 + bn * (sbv1 - sbvz)) / (SL - SM);

        /*if (fabs(bn) < 0.000001 )
        {
            vRR = v2;
            wRR = w2;
            vLL = v1;
            wLL = w1;
            btt2 = bt2 * suRm;
            bmm2 = bm2 * suRm;
            btt1 = bt1 * suLm;
            bmm1 = bm1 * suLm;
        }*/

        /*ppLR = (pTL + ro_L * (SL - u1) * (SM - u1) + pTR + ro_R * (SR - u2) * (SM - u2)) / 2.0;

        if (fabs(bn) < 0.000001)
        {
            vLL = v1;
            wLL = w1;
            vRR = v2;
            wRR = w2;

            btt1 = bt1 * (SL - u1) / (SL - SM);
            btt2 = bt2 * (SR - u2) / (SR - SM);

            bmm1 = bm1 * (SL - u1) / (SL - SM);
            bmm2 = bm2 * (SR - u2) / (SR - SM);

            ee1 = ((SL - u1) * e1 - pTL * u1 + ppLR * SM) / (SL - SM);
            ee2 = ((SR - u2) * e2 - pTL * u2 + ppLR * SM) / (SR - SM);
        }
        else
        {
            btt2 = btt1 = (SR * UR[6] - SL * UL[6] + FL[6] - FR[6]) / (SR - SL);
            bmm2 = bmm1 = (SR * UR[7] - SL * UL[7] + FL[7] - FR[7]) / (SR - SL);
            vLL = v1 + bn * (bt1 - btt1) / (ro_L * (SL - u1));
            vRR = v2 + bn * (bt2 - btt2) / (ro_R * (SR - u2));

            wLL = w1 + bn * (bm1 - bmm1) / (ro_L * (SL - u1));
            wRR = w2 + bn * (bm2 - bmm2) / (ro_R * (SR - u2));

            double sks1 = u1 * bn1 + v1 * bt1 + w1 * bm1 - SM * bn - vLL * btt1 - wLL * bmm1;
            double sks2 = u2 * bn2 + v2 * bt2 + w2 * bm2 - SM * bn - vRR * btt2 - wRR * bmm2;

            ee1 = ((SL - u1) * e1 - pTL * u1 + ppLR * SM + bn * sks1) / (SL - SM);
            ee2 = ((SR - u2) * e2 - pTR * u2 + ppLR * SM + bn * sks2) / (SR - SM);
        }*/


        double  ULL[9], URR[9], PO[9];
        ULL[0] = ro_LL;
        ULL[1] = ro_LL * SM;
        ULL[2] = ro_LL * vLL;
        ULL[3] = ro_LL * wLL;
        ULL[4] = ee1;
        ULL[5] = bn;
        ULL[6] = btt1;
        ULL[7] = bmm1;
        ULL[8] = Q_LL;

        URR[0] = ro_RR;
        URR[1] = ro_RR * SM;
        URR[2] = ro_RR * vRR;
        URR[3] = ro_RR * wRR;
        URR[4] = ee2;
        URR[5] = bn;
        URR[6] = btt2;
        URR[7] = bmm2;
        URR[8] = Q_RR;

        if (SL >= 0.0)
        {
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FL[i] - wv * UL[i];
            }
        }
        else if (SL < 0.0 && SM >= 0.0)
        {
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FL[i] + SL * ULL[i] - SL * UL[i] - wv * ULL[i];
            }
        }
        else if (SR > 0.0 && SM < 0.0)
        {
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FR[i] + SR * URR[i] - SR * UR[i] - wv * URR[i];
            }
        }
        else if (SR <= 0.0)
        {
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FR[i] - wv * UR[i];
            }
        }



        double SN = max(fabs(SL), fabs(SR));

        double wbn = 0.0;
        if (wv >= SR)
        {
            wbn = wv * bn2;
        }
        else if (wv <= SL)
        {
            wbn = wv * bn1;
        }
        else
        {
            wbn = wv * (bn1 + bn2) / 2.0;
        }

        PO[5] = -SN * (bn2 - bn1) - wbn;

        P[1] = n1 * PO[1] + t1 * PO[2] + m1 * PO[3];
        P[2] = n2 * PO[1] + t2 * PO[2] + m2 * PO[3];
        P[3] = n3 * PO[1] + t3 * PO[2] + m3 * PO[3];
        P[5] = spi4 * (n1 * PO[5] + t1 * PO[6] + m1 * PO[7]);
        P[6] = spi4 * (n2 * PO[5] + t2 * PO[6] + m2 * PO[7]);
        P[7] = spi4 * (n3 * PO[5] + t3 * PO[6] + m3 * PO[7]);
        P[0] = PO[0];
        P[4] = PO[4];
        PQ = PO[8];

        double SWAP = P[4];
        P[4] = P[5];
        P[5] = P[6];
        P[6] = P[7];
        P[7] = SWAP;
        return time;

    }
    else if (metod == 3)  // HLLD
    {

        double ttL = ro_L * suL * (SL - SM) - bbn;   // Может быть здесь не такой bn нужен !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        double ttR = ro_R * suR * (SR - SM) - bbn;

        double vLL, wLL, vRR, wRR, btt1, bmm1, btt2, bmm2;

        double PTTL = PTT;
        double PTTR = PTT;

        if (fabs(ttL) >= 0.00001)
        {
            vLL = v1 - bn * bt1 * (SM - u1) / ttL;
            wLL = w1 - bn * bm1 * (SM - u1) / ttL;
            btt1 = bt1 * (ro_L * suL * suL - bbn) / ttL;
            bmm1 = bm1 * (ro_L * suL * suL - bbn) / ttL;
        }
        else
        {
            vLL = v1;
            wLL = w1;
            btt1 = 0.0;
            bmm1 = 0.0;
            ro_LL = ro_L;
            PTTL = pTL;
            cout << "2252   jfhfghjieuye  LLL" << endl;
            cout << ro_L << " " << p_L << " " << ro_R << " " << p_R << " " << v1_L << " " << v1_R << " " << Bx_L << " " << Bx_R << endl;
            exit(-1);
        }

        if (fabs(ttR) >= 0.00001)
        {
            vRR = v2 - bn * bt2 * (SM - u2) / ttR;
            wRR = w2 - bn * bm2 * (SM - u2) / ttR;
            btt2 = bt2 * (ro_R * suR * suR - bbn) / ttR;
            bmm2 = bm2 * (ro_R * suR * suR - bbn) / ttR;
            //cout << "tbr = " << (ro_R * suR * suR - bbn) / ttR << endl;
            //cout << "bt2 = " << bt2 << endl;
        }
        else
        {
            vRR = v2;
            wRR = w2;
            btt2 = 0.0;
            bmm2 = 0.0;
            ro_RR = ro_R;
            PTTR = pTR;
            cout << "2272   jfhfghjieuye  RRR" << endl;
        }

        double eLL = (e1 * suL + PTTL * SM - pTL * u1 + bn * //              опять  bn  не тот !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ((u1 * bn1 + v1 * bt1 + w1 * bm1) - (SM * bn + vLL * btt1 + wLL * bmm1))) //
            / (SL - SM);
        double eRR = (e2 * suR + PTTR * SM - pTR * u2 + bn * //
            ((u2 * bn2 + v2 * bt2 + w2 * bm2) - (SM * bn + vRR * btt2 + wRR * bmm2))) //
            / (SR - SM);

        double sqrtroLL = sqrt(ro_LL);
        double sqrtroRR = sqrt(ro_RR);
        double SLL = SM - fabs(bn) / sqrtroLL;  //              опять  bn  не тот !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        double SRR = SM + fabs(bn) / sqrtroRR;  //              опять  bn  не тот !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        double idbn = 1.0;
        if (fabs(bn) > 0.000001)
        {
            idbn = 1.0 * sgn(bn);
        }
        else
        {
            idbn = 0.0;
            SLL = SM;
            SRR = SM;
            //cout << "idbn == 0" << endl;
        }

        double vLLL = (sqrtroLL * vLL + sqrtroRR * vRR + //
            idbn * (btt2 - btt1)) / (sqrtroLL + sqrtroRR);

        double wLLL = (sqrtroLL * wLL + sqrtroRR * wRR + //
            idbn * (bmm2 - bmm1)) / (sqrtroLL + sqrtroRR);

        double bttt = (sqrtroLL * btt2 + sqrtroRR * btt1 + //
            idbn * sqrtroLL * sqrtroRR * (vRR - vLL)) / (sqrtroLL + sqrtroRR);

        double bmmm = (sqrtroLL * bmm2 + sqrtroRR * bmm1 + //
            idbn * sqrtroLL * sqrtroRR * (wRR - wLL)) / (sqrtroLL + sqrtroRR);

        double eLLL = eLL - idbn * sqrtroLL * ((SM * bn + vLL * btt1 + wLL * bmm1) //
            - (SM * bn + vLLL * bttt + wLLL * bmmm));
        double eRRR = eRR + idbn * sqrtroRR * ((SM * bn + vRR * btt2 + wRR * bmm2) //
            - (SM * bn + vLLL * bttt + wLLL * bmmm));
        //cout << " = " << bn << " " << btt2 << " " << bmm2 << endl;
        //cout << "sbvr = " << (SM * bn + vRR * btt2 + wRR * bmm2) << endl;
        double  ULL[9], URR[9], ULLL[9], URRR[9];

        ULL[0] = ro_LL;
        ULL[1] = ro_LL * SM;
        ULL[2] = ro_LL * vLL;
        ULL[3] = ro_LL * wLL;
        ULL[4] = eLL;
        ULL[5] = bn;
        ULL[6] = btt1;
        ULL[7] = bmm1;
        ULL[8] = Q_LL;

        URR[0] = ro_RR;
        //cout << ro_RR << endl;
        URR[1] = ro_RR * SM;
        URR[2] = ro_RR * vRR;
        URR[3] = ro_RR * wRR;
        URR[4] = eRR;
        URR[5] = bn;
        URR[6] = btt2;
        URR[7] = bmm2;
        URR[8] = Q_RR;

        ULLL[0] = ro_LL;
        ULLL[1] = ro_LL * SM;
        ULLL[2] = ro_LL * vLLL;
        ULLL[3] = ro_LL * wLLL;
        ULLL[4] = eLLL;
        ULLL[5] = bn;
        ULLL[6] = bttt;
        ULLL[7] = bmmm;
        ULLL[8] = Q_LL;

        URRR[0] = ro_RR;
        URRR[1] = ro_RR * SM;
        URRR[2] = ro_RR * vLLL;
        URRR[3] = ro_RR * wLLL;
        URRR[4] = eRRR;
        URRR[5] = bn;
        URRR[6] = bttt;
        URRR[7] = bmmm;
        URRR[8] = Q_RR;

        double PO[9];

        if (SL >= 0.0)
        {
            //cout << "SL >= 0.0" << endl;
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FL[i] - wv * UL[i];
            }
        }
        else if (SL <= 0.0 && SLL >= 0.0)
        {
            //cout << "SL < 0.0 && SLL >= 0.0" << endl;
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FL[i] + SL * ULL[i] - SL * UL[i] - wv * ULL[i];
            }
            //cout << ULL[0] << endl;
        }
        else if (SLL < 0.0 && SM >= 0.0)
        {
            //cout << "SLL <= 0.0 && SM >= 0.0" << endl;
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FL[i] + SLL * ULLL[i] - (SLL - SL) * ULL[i] - SL * UL[i] - wv * ULLL[i];
            }
        }
        else if (SM <= 0.0 && SRR > 0.0)
        {
            //cout << "SM < 0.0 && SRR > 0.0" << endl;
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FR[i] + SRR * URRR[i] - (SRR - SR) * URR[i] - SR * UR[i] - wv * URRR[i];
            }
            //cout << "P4 = " << URRR[4] << endl;
        }
        else if (SR >= 0.0 && SRR <= 0.0)
        {
            //cout << "SR > 0.0 && SRR <= 0.0" << endl;
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FR[i] + SR * URR[i] - SR * UR[i] - wv * URR[i];
            }
            //cout << URR[0] << endl;
        }
        else if (SR <= 0.0)
        {
            //cout << "SR <= 0.0" << endl;
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FR[i] - wv * UR[i];;
            }
        }



        double SN = max(fabs(SL), fabs(SR))/2.0;
        double wbn = 0.0;
        if (wv >= SR)
        {
            wbn = wv * bn2;
        }
        else if (wv <= SL)
        {
            wbn = wv * bn1;
        }
        else
        {
            wbn = wv * (bn1 + bn2) / 2.0;
        }

        PO[5] = -SN * (bn2 - bn1) - wbn;

        P[1] = n1 * PO[1] + t1 * PO[2] + m1 * PO[3];
        P[2] = n2 * PO[1] + t2 * PO[2] + m2 * PO[3];
        P[3] = n3 * PO[1] + t3 * PO[2] + m3 * PO[3];
        P[5] = spi4 * (n1 * PO[5] + t1 * PO[6] + m1 * PO[7]);
        P[6] = spi4 * (n2 * PO[5] + t2 * PO[6] + m2 * PO[7]);
        P[7] = spi4 * (n3 * PO[5] + t3 * PO[6] + m3 * PO[7]);
        P[0] = PO[0];
        P[4] = PO[4];
        PQ = PO[8];

        double SWAP = P[4];
        P[4] = P[5];
        P[5] = P[6];
        P[6] = P[7];
        P[7] = SWAP;
        return time;
    }
    return time;
}

void Setka::HLLD_Test(void)
{
    double P[8] = { 0.0 };
    P[0] = P[1] = P[2] = P[3] = P[4] = P[5] = P[6] = P[7] = 0.0;

    double ro, p, vx, vy, vz, bx, by, bz, ro2, p2, vx2, vy2, vz2, bx2, by2, bz2;

    ro = 1.12;
    p = 1.35;
    vx = 0.11;
    vy = 0.34;
    vz = 0.13;
    bx = -0.1;
    by = 0.4;
    bz = 0.11;

    ro2 = 2.12;
    p2 = 0.65;
    vx2 = 0.16;
    vy2 = 0.14;
    vz2 = 0.23;
    bx2 = -1.7;
    by2 = 0.9;
    bz2 = 0.51;


    double Vc;
    HLLD_Alexashov(ro, p, vx, vy, vz, bx, by, bz, ro2, p2, vx2, vy2, vz2, bx2, by2, bz2, P, 1.0, 0.0, 0.0, 1.0, 2, Vc);

    for (int i = 0; i < 8; i++)
    {
        cout << "P[" << i << "] = " << std::setprecision(16) << P[i] << endl;
    }

    cout << endl;

    double PQ;
    HLLDQ_Korolkov(ro, 1.0, p, vx, vy, vz, bx, by, bz, ro2, 1.0, p2, vx2, vy2, vz2, bx2, by2, bz2, P, PQ, 1.0, 0.0, 0.0, 1.0, 2, Vc);

    for (int i = 0; i < 8; i++)
    {
        cout << "P[" << i << "] = " << std::setprecision(16) << P[i] << endl;
    }
}