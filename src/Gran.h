#pragma once

#include "Help.h"
#include "Cell.h"

using namespace std;

class Cell;

class Gran
{
public:
	Cell* Sosed;
	double S;                // площадь грани     нормаль к грани не нужна, т.к. её можно вычислить, используя координаты соседа
	double n1 = 0.0;         // Нормали - быстрее будет, если их вычислять заранее, однако если памяти будет не хватать, то можно убрать
	double n2 = 0.0;
	double n3 = 0.0;
	double c1 = 0.0;         // Перемещение грани относительно предыдущего положения 
	double c2 = 0.0;
	double c3 = 0.0;
	double dl = 0.0;
	bool check_ = false;

	Gran();

	void set_normal(const double& a, const double& b, const double& c);
	void get_normal(double& a, double& b, double& c);
};

