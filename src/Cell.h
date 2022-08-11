#pragma once
#include <vector>
#include <map>
#include <set>
#include "Help.h"
using namespace std;
class Point;
class Gran;
class Couple;

enum Cell_type  // Тип грани нужен для граничных условий
{
	C_base,                     // Базовая ячейка \ обычная   включая парные, это не важно
	C_wall_x_max,               // Мнимая ячейка - стена
	C_wall_x_min,               // Мнимая ячейка - стена
	C_wall_y_max,               // Мнимая ячейка - стена
	C_wall_y_min,               // Мнимая ячейка - стена
	C_wall_z_max,               // Мнимая ячейка - стена
	C_wall_z_min,               // Мнимая ячейка - стена
	C_sphere,                   // Мнимая ячейка - внутренняя сфера
};

struct Parametr
{
	double ro = 0.0;
	double p = 0.0;
	double u = 0.0;
	double v = 0.0;
	double w = 0.0;
	double bx = 0.0;
	double by = 0.0;
	double bz = 0.0;
	double phi = 0.0;
	double divB = 0.0;
	double sivi = 0.0;
};

class Cell
{
public:
	Point* Center = nullptr;                // Центр ячейки
	Couple* Par = nullptr;               // Пара ячейки (для парных ячеек)
	map <int, Cell*> Can_neighbours;     // Ячейки - соседи
	//map <int, Gran*> Grans_do;     // Грани до добавления новых (нужно для расчёта перемещения грани)
	vector <Cell*> Candidates;         //  Реальные грани с соседями и площадью грани и нормалью
	vector <Gran*> Grans;         //  Реальные грани с соседями и площадью грани и нормалью
	vector <Gran*> Grans_TVD;         //  Ячейки - противоположенные к каждой грани
	vector <Cell*> Grans_standart;         //  Ячейки - противоположенные к каждой грани
	// Эти грани задаются при первоначальном декартовом построении сетки, их не нужно никогда менять
	// Более того, если их "потерять" то восстановить уже не получится.
	// Введены для возможности возврата ячейки, которая была удалена из сетки.
	double Potok[9] = { 0.0 };

	//map <int, Gran*> mp;
	Parametr par[2];              // газодинамические параметры в ячейке
	int number;                   // Номер ячейки
	Cell_type type;
	double move1;                   // Вектор сдвига для парных ячеек
	double move2;                   // Вектор сдвига для парных ячеек
	double move3;                   // Вектор сдвига для парных ячеек

	mutex acces_TVD;
	mutex acces_POTOK;
	bool check_ = false;

	

	// Управляющие параметры
	int n_cop = 0;                // Какая она по счёту в паре
	bool mgd_ = true;             // Нужно ли считать законы сохранения и менять параметры в этой ячейке
	// Пока не знаю зачем это, это вроде нужно для всех включённых ячеек?!?!?!?
	bool couple_ = false;         // Является ли ячейка парной?  // очень часто нужно это знать везде
	// Главное что этот параметр задаётся и больше не может меняться!
	bool include_ = true;         // Включена ли ячейка сейчас в сетку?
	// Очень важны параметр, нужен для возможность включать\выключать ячейки
	bool reconstruct_ = false;    // Нужно ли перестраивать эту ячейку (в случае подозрения на изменение соседей)
	// Этот параметр должен перерасти в следующие два, пока не уверен, что он нужен
	bool near_par = false;        // явлюяется ли ячейка граничащей и парной ячейкой
	// Важно что этот параметр постоянно меняется и его нужно периодически обновлять!!!
	bool near_par_2 = false;      // сосед второго уровня с парной ячейкой, также динамический параметр

	bool i_sosed = false;         // Параметр для добавления соседей в эту ячейку (он просто используется в некоторых функциях)
	bool i_bbb = false;           // вспомогательный маяк
	bool i_boundary_1 = false;    // Рядом с границей
	bool i_boundary_2 = false;    // Сосед второго порядка к границе
	bool i_include_candidate = false;     // Кандидат на отключение
	bool i_init_mgd = false;
	bool i_delete = false;                 // Нужно ли удалить ячейку?

	bool extern_boundary = false;            // Граничит ли ячейка с границей (внешней ил внутренней, не важно.

	bool TVD_reconstruct = false;                 // Нужно ли пересчитывать ТВД?
	

	Cell(const double& x, const double& y, const double& z);  // Конструктор сам создаёт центр ячейки в памети и сохраняет указатель на него
	void set_Volume(const double& V);
	void get_Volume(double& V);
	void get_Volume_do(double& V);

	// Работа с ТВД
	bool Get_TVD(Gran* G, Gran*& A);
	int Get_TVD_Param(Gran* G, Parametr& p1, Parametr& p2, int now, int n_gran, bool bnkl = false);
	// Ищет снесённые значения из данной ячейки на грань G, при этом A - сосед "сзади", т.е. со стороны самой ячейки


private:
	double Volume;                // Объём ячейки  (вынесен сюда, для того чтобы автоматически сохранять предыдущее значение
	double Volume_do;                // Объём ячейки

}; 

