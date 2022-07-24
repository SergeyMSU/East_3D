#pragma once
#include <vector>
#include <map>
#include "Help.h"
using namespace std;
class Point;
class Gran;


class Couple
{
public:
	Cell* A1;
	Cell* A2;
	double dist;       // Расстояние между парными ячейками
	double n1;
	double n2;
	double n3;
	double dvig1 = 0.0;                   // Вектор сдвига для парных ячеек
	double dvig2 = 0.0;                   // Вектор сдвига для парных ячеек
	double dvig3 = 0.0;                   // Вектор сдвига для парных ячеек
	double mo1 = 0.0;                   // Вектор сдвига для парных ячеек
	double mo2 = 0.0;                   // Вектор сдвига для парных ячеек
	double mo3 = 0.0;                   // Вектор сдвига для парных ячеек
	double d_sosed;
	int number = -1;

	Couple(Cell* A1, Cell* A2, const double& dist);

	void Resolve(void); // пересчитывает расстояния между ячейками
	void get_centr(double& x, double& y, double& z);
	void get_normal(double& x, double& y, double& z);  // Вычисляет нормаль пары по центрам A1 и A2, а не просто берёт её.
	void orient(void); // Ориентировать пару по её нормали (меняет координаты пары)
	void move(const double& m1, const double& m2, const double& m3);


};

