#include "voro++.hh"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iterator>
#include "Help.h"
#include "Setka.h"
#include "Cell.h"
#include <unistd.h>

class Couple;


// Гид по примерам

// Пример 1: создаётся в ручную набор точек на плоскости, далее также в ручную создаётся ячейка и режатся плоскостями для каждой точки
// далее это всё печатается в текплот файл
// ОШИБКА! Это всё не совсем правильно сделано
// Аторы при посроении помещают ячейку в начало координат, а потом разрезают её плоскостями, поэтому нужно потом все узлы сдвинуть!
// Именно по этому не существует плоскости, с нулевым расстоянием, потому что там сама точка находится! Это прикольно, но слегка не удобно
// нужно переделать правильно! см. Пример 3

// Пример 2: как в предыдущём примере, но используется контейнер, ничего не режется в ручную. Точки берутся из файла
// также печатается результат в файл
// Также установили границу - цилиндр

// Пример 3: то же самое, что пример 1, но более правильно организован разрез

// Планы: 
// Хорошо бы научиться менять координаты точек, которые уже находятся в контейнере  (скорее всего это не очень правильно и удобно)

using namespace voro;
using namespace std;

double rnd() { return double(rand()) / RAND_MAX; }

int sign(double x)
{
	if (x > 0)
	{
		return 1;
	}
	else if (x < 0)
	{
		return -1;
	}
	return 0;
}


int main()
{
	if (true)
	{
		cout << "SLEEP" << endl;
		//sleep(7200);
		cout << "NOT SLEEP" << endl;

		Setka B = Setka("vers_28");  // 19  до bn

		B.Initialization_do_MHD();

		//B.Calc_normal();
		//B.Construct_start();
		bool bkl = true;
		if (false)
		{
			B.Add_par();
			B.renumber();
			B.Construct_start();
			B.Rebuild();
			B.Calc_normal();
			bool bkl = true;
			do
			{
				bkl = B.Reconstruct_medium4(true);
				cout << "bkl = " << bkl << endl;
			} while (bkl == false);
			B.Calc_normal2();
			bkl = true;
			do
			{
				bkl = B.Reconstruct_medium4(true);
				cout << "bkl = " << bkl << endl;
			} while (bkl == false);
			B.Calc_normal();
			bkl = true;
			do
			{
				bkl = B.Reconstruct_medium4(true);
				cout << "bkl = " << bkl << endl;
			} while (bkl == false);
		}

		for (int iii = 0; iii < 0; iii++)
		{
			B.Smooth();
			do
			{
				bkl = B.Reconstruct_medium4(true);
				cout << "bkl = " << bkl << endl;
			} while (bkl == false);
			B.Calc_normal();
			do
			{
				bkl = B.Reconstruct_medium4(true);
				cout << "bkl = " << bkl << endl;
			} while (bkl == false);
			B.Calc_normal2();
			do
			{
				bkl = B.Reconstruct_medium4(true);
				cout << "bkl = " << bkl << endl;
			} while (bkl == false);
			B.Calc_normal();
			do
			{
				bkl = B.Reconstruct_medium4(true);
				cout << "bkl = " << bkl << endl;
			} while (bkl == false);
		}

		if (false)
		{
			int yuyu = 0;
			for (auto& ii : B.All_Couple)
			{
				if (ii->A1->include_ == false)
				{
					continue;
				}

				map <int, Couple*> Nei;
				double x, y, z;
				double x2, y2, z2;

				ii->get_centr(x, y, z);

				// Сначала ищем всех соседей
				for (auto& j : ii->A1->Grans)
				{
					if (j->Sosed->couple_ == true && j->Sosed->include_ == true)
					{
						if (j->Sosed->Par != ii)
						{
							Nei.insert(make_pair(j->Sosed->Par->number, j->Sosed->Par));
						}
					}
				}

				for (auto& j : ii->A2->Grans)
				{
					if (j->Sosed->couple_ == true && j->Sosed->include_ == true)
					{
						if (j->Sosed->Par != ii)
						{
							Nei.insert(make_pair(j->Sosed->Par->number, j->Sosed->Par));
						}
					}
				}



				for (auto& i : Nei)
				{
					i.second->get_centr(x2, y2, z2);
					double a1, a2, a3;
					a1 = (x2 - x);
					a2 = (y2 - y);
					a3 = (z2 - z);
					double aa = sqrt(kv(a1) + kv(a2) + kv(a3));
					a1 = a1 / aa;
					a2 = a2 / aa;
					a3 = a3 / aa;
					if (fabs(a1 * ii->n1 + a2 * ii->n2 + a3 * ii->n3) > 0.52)
					{
						yuyu++;
						yuyu++;
						i.second->A1->include_ = false;
						i.second->A1->Candidates.clear();
						i.second->A1->Grans.clear();
						i.second->A2->include_ = false;
						i.second->A2->Candidates.clear();
						i.second->A2->Grans.clear();

						ii->A1->include_ = false;
						ii->A1->Candidates.clear();
						ii->A1->Grans.clear();
						ii->A2->include_ = false;
						ii->A2->Candidates.clear();
						ii->A2->Grans.clear();
						break;
					}
				}
				Nei.clear();
			}
				cout << "Par izolated  " << yuyu << endl;
			
		}


		//B.Construct_start();
		//B.Calc_normal();

		do
		{
			bkl = B.Reconstruct_medium4(true);
			cout << "bkl = " << bkl << endl;
		} while (bkl == false);

		/*for (auto& i : B.All_Cell)
		{
			if (i->Center->x < -4.0)
			{
				i->par[0].bx = 0.0;
				i->par[0].by = 0.0;
				i->par[0].bz = 0.0;
				i->par[1].bx = 0.0;
				i->par[1].by = 0.0;
				i->par[1].bz = 0.0;
			}
		}*/


		B.Zapusk();
		B.Save_setka("vers_29");
		//B.Construct_start();

		B.Cut_Plane_z(1.0);
		B.Cut_Plane_y(1.0);
		B.Cut_Surface();
		B.Cut_Plane_z_Tecplot(1.0);
		B.Cut_Plane_y_Tecplot(1.0);

		exit(-1);
	}

	Setka A = Setka();
	A.Construct_start();
	A.save_soseds();
	//A.HLLD_Test();
	//return 0;
	A.Add_couple_start();
	//A.Print_points();
	A.Construct_start();
	A.Disable_cells();
	A.Construct_start();
	/*A.Reconstruct_couple(); 
	A.Reconstruct_couple();
	A.Reconstruct_couple(true);*/
	A.set_normal();
	//A.move_par();
	/*A.Reconstruct_couple();
	A.Reconstruct_couple();
	A.Reconstruct_couple();*/
	A.Calc_normal();
	/*A.Reconstruct_couple();
	A.Reconstruct_medium2();*/

	A.Init();
	A.Initialization_do_MHD();

	A.Cut_Plane_z(1.0);
	A.Cut_Plane_y(1.0);
	A.Cut_Surface();
	//A.Cut_Plane_y_Tecplot(1.0);
	//A.Download_MHD("Vers_1.txt");
	//A.Go_MHD(30);
	//A.Calc_normal();                           // Проверить в самом начала, мне кажется там что-то не правильно
	//A.Start_MHD(100);
	//A.Zapusk();
	//A.Reconstruct_couple();
	A.Cut_Plane_z_Tecplot(1.0);
	A.Cut_Plane_y_Tecplot(1.0);
	//A.Save_MHD("Vers_2.txt");
	A.Cut_Plane_z(1.0);
	A.Cut_Plane_y(1.0);
	A.Cut_Surface();

	A.Save_setka("vers_1");

	/*for (auto& i : A.All_Cell[80000]->Grans)
	{
		cout << i->S << " " << i->Sosed->number << " " << i->n1 << " " << i->n2 << " " << i->n3 << " " <<  endl;
	}*/

	return 0;
	double x, y, z, rsq, r;
	vector<int> neigh, f_vert, ney, edge, vershin;
	vector<double> X, Y, Ver, Versh;

	X.push_back(-1.0);
	Y.push_back(0.0);

	X.push_back(0.0);
	Y.push_back(0.0);

	X.push_back(1.0);
	Y.push_back(0.0);

	X.push_back(-1.0);
	Y.push_back(1.0);

	X.push_back(0.0);
	Y.push_back(1.0);

	X.push_back(1.0);
	Y.push_back(1.0);

	X.push_back(0.5);
	Y.push_back(0.5);

	X.push_back(-0.5);
	Y.push_back(0.5);

	X.push_back(-1.0);
	Y.push_back(2.0);

	X.push_back(0.0);
	Y.push_back(0.5);

	X.push_back(0.0);
	Y.push_back(2.0);

	double nx, ny, nz, nn;

	cout << "START " << endl;

	ofstream fout;
	string name_f = "setka.txt";
	fout.open(name_f);

	ofstream fout2;
	name_f = "points.txt";
	fout2.open(name_f);

	int N = 0;
	int E = 0;

	for (int i = 0; i < X.size(); i++) // X.size()
	{
		fout2 << X[i] << " " << Y[i] << endl;
	}

	for (int i = 0; i < X.size(); i++) // X.size()
	{
		cout << "D = " << i << endl;
		voronoicell_neighbor c;
		c.init(-20.0, 20.0, -20.0, 20.0, -20.0, 20.0);

		for (int j = 0; j < X.size(); j++)
		{
			if (i == j)
			{
				continue;
			}

			// Странно организован разрез, но как есть. Здесь плоскость задаётся нормалью и квадратом расстояния до нуля (с коэффициентом, см код ниже)
			double A = (X[j] - X[i]);
			double B = (Y[j] - Y[i]);
			double C = 0.0;
			c.nplane(A, B, C, (A * A + B * B + C * C), j);
			// Вообще обязательно чтобы четвёртое число было квадратом первых трёх со знаком (какой плоскостью резать)
		}

		// Весь код далее для печати сетки в файл Текплота. Сложности в том, что нужно по порядку расположить все узлы и грани сетки (Текплот так работает)
		c.nplane(0.0, 0.0, 1.0, 0.01, -7);

		c.neighbors(ney); // Получаем список номеров соседей по порядку
		int k = -1;  // номер нужного соседа
		// Ищем номер соседа "-7" - скорее всего это последний сосед (надо проверить)  ответ - нет!!!
		for (int j = 0; j < ney.size(); j++)
		{
			k++;
			if (ney[j] == -7)
			{
				break;
			}
		}

		if (k == -1)
		{
			continue;
		}
		cout << "K = " << k << endl;
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

		cout << "m = " << m << endl;

		for (int j = m + 1; j < m + 1 + f_vert[m] - 1; j++)
		{
			vershin.push_back(N + f_vert[j] + 1);  // Вектор номеров нужных вершин
			vershin.push_back(N + f_vert[j + 1] + 1);  // Вектор номеров нужных вершин
		}

		vershin.push_back(N + f_vert[m + 1 + f_vert[m] - 1] + 1);
		vershin.push_back(N + f_vert[m + 1] + 1);

		c.vertices(Ver);

		for (int j = 0; j < Ver.size(); j++)
		{
			if (j % 3 == 0)
			{
				Versh.push_back(Ver[j] + X[i]);  // Вектор всех вершин
			}
			else if (j % 3 == 1)
			{
				Versh.push_back(Ver[j] + Y[i]);  // Вектор всех вершин
			}
			else if(j % 3 == 2)
			{
				Versh.push_back(Ver[j]);  // Вектор всех вершин
			}
		}

		N += Ver.size() / 3;
		E += f_vert[m];

		if (true)
		{
			printf("Total vertices      : %d\n", c.p);
			printf("Vertex positions    : "); c.output_vertices(); puts("");
			printf("Vertex orders       : "); c.output_vertex_orders(); puts("");
			printf("Max rad. sq. vertex : %g\n\n", 0.25 * c.max_radius_squared());
			cout << endl;
			printf("Total faces         : %d\n", c.number_of_faces());
			printf("Surface area        : %g\n", c.surface_area());
			printf("Face freq. table    : "); c.output_face_freq_table(); puts("");
			printf("Face orders         : "); c.output_face_orders(); puts("");
			printf("Face areas          : "); c.output_face_areas(); puts("");
			printf("Face normals        : "); c.output_normals(); puts("");
			printf("Face vertices       : "); c.output_face_vertices(); puts("\n");
			c.neighbors(ney);
			c.vertices(Ver);
			c.face_vertices(f_vert);
			for (int j = 0; j < Ver.size(); j++)
			{
				cout << Ver[j] << " ";
			}
			cout << endl;
			for (int j = 0; j < ney.size(); j++)
			{
				cout << ney[j] << " ";
			}
			cout << endl;
			for (int j = 0; j < f_vert.size(); j++)
			{
				cout << f_vert[j] << " ";
			}
			cout << endl;
		}

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