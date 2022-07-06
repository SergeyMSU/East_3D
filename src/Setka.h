#pragma once
#include "Help.h"
#include "Cell.h"
#include "voro++.hh"
#include <vector>
#include <string>
#include <mutex>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <list>

using namespace std;

class Cell;
class Couple;


class Setka
{
public:
	vector <Cell*> All_Cell;              // Вектор всех ячеек сетки
	vector <Couple*> All_Couple;              // Вектор всех ячеек сетки
	list <Cell*> All_Cell_off;             // Список выключенных ячеек

	// Все особые ячейки создаются в специальной функции inicializtion
	Cell* Wall1;                          // Особые ячейки - границы       -  C_wall_x_max
	Cell* Wall2;                          // C_wall_x_min
	Cell* Wall3;
	Cell* Wall4;
	Cell* Wall5;
	Cell* Wall6;
	Cell* Wall7;

	double x_min = -15.0;
	double x_max = 10.0;
	double y_min = -10.1;
	double y_max = 10.0;
	double z_min = -10.0;
	double z_max = 10.1;
	double sl_L = 0.6;                    // От какого радиуса слои включительно   0.3
	double sl_R = 1.5;

	Setka();                              // Начальное заполнение сетки точками


	void Initialization_do_MHD(void);

	void Construct_start(void);                 // Начальное построение сетки (с использованием контейнера), работает долго
	// Эту функцию нужно использовать только при первоначальном построении сетки, далешь нужно использовать более быстрые алгоритмы.
	// Она находит всех соседей для каждой ячейки (включая граничных соседей), далее записывает всё в "грани" (площадь и номер соседа)
	void Reconstruct_couple(bool t1 = false);
	// построение сетки с учётом имеющихся соседей у ячейки, также учитываются отключенные ячейки (они не строятся)
	// перестроение начинается от пар (эта функция следует после отключения ненужных ячеек
	// Также добавляются соседи-соседей как всегда
	// t1 - нужно ли для парных ячеек считать вектор сдвига (последующего перемещения по алгоритму Алексашова)?
	void Reconstruct_fast(bool wide = true);  // Быстрое перестроение сетки
	// Перестраивает сетку по предыдущим соседям, новых не добавляет
	// т.е. функция работает при небольшом движении, часто пользовать ей нельзя так как может потерять соседей.
	void Reconstruct_medium(bool wide = true);
	// перестраиваем через добавление соседей-соседей
	// для парных ячеек, при этом не парным добавляем только те, которые у них были до этого
	void Reconstruct_medium2(bool wide = true);
	void Reconstruct_medium3(bool wide = true);
	void renumber(void);

	void include_cells();  // включить ячейки
	void exclude_cells();  // выключить ячейки

	/// ------------------------------------------------------------------------------------------------------------------------
	///  Управление парами и всё что с этим связано
	/// ------------------------------------------------------------------------------------------------------------------------
	void Add_couple_start(void);                 // Начальное добавление пар в сетку (просто цилиндр из пар)
	void Disable_cells(void);                    // Отключение ячеек, слишком близких к парам
	// пробегаемся по всем соседям у пар и отключаем их в случае, если они слишком близко
	// "отключаем" - это просто заносим их в список отключённых и устанавливаем флажок, реально они ещё не отключены!
	void Calc_normal(void);                    // Вычисляет новые нормали у парных ячеек (по соседям)  ???????????????????????????????????????????????????????????????????????????????????????????
	// Также переориентирует пару
	
	void move_par(void);                       // Передвинуть пары (НЕ по нормали) в соответствии с функцией регуляризации  ??????????????????????????????????????????????????????????????????????
	void set_normal(void);                     // Просто вычисляет нормаль каждой пары и устанавливает её (для начальной установки нормалей)
	void Cut_Surface(void);            // Печать поверхности между парами


	// Вывод геометрии сетки
	void Print_points(void);              // Выводит все точки (просто три координаты каждой точки) - далее можно их рисовать в текплоте

	void Cut_Plane_long_z(void);            // Разрез всех ячееки плоскостью z = 0
	// Разрез проводится через помещение всех точек в контейнер Voro++

	void Cut_Plane_z(double R = 0.0);            // Разрез всех ячееки плоскостью z = 0
	// Быстрая функция, использует уже найденных соседей каждой ячейки
	// R - расстояние, за которым резать уже не надо.
	void Cut_Plane_y(double R = 0.0);            // Разрез всех ячееки плоскостью z = 0

	void Cut_Plane_z_Tecplot(double R = 0.0);    // Печать газодинамические переменных в файл в плоскости z = 0
	// сначала режет ячейки плоскостью, потом находит точку - центройд грани, и в этих координатах печатает газодинамические переменные

	void Cut_Plane_y_Tecplot(double R = 0.0);    // Печать газодинамических переменных в файл в плоскости y = 0





	// Преобразования координат

	void dekard_skorost(double x, double y, double z, double Vr, double Vphi, double Vtheta, double& Vx, double& Vy, double& Vz);
	double polar_angle(double x, double y);



	// МГД
	void Init(void);         // Инициализация газодинамических параметров
	void Go_MHD(int times);  // запускаем счёт МГД на стационарной сетке  times - шагов
	void Start_MHD(int times);   // запускаем программу на нестационарной сетке!!!  ОСНОВНАЯ ФУНКЦИЯ
	void Culc_couple(int now1, const double& time);    // Считаем движение пар из задачи о распаде разрыва и двигаем их!!!

	void Save_MHD(string name_f = "Save_MHD.txt");
	void Download_MHD(string name_f = "Save_MHD.txt");
	// Простые функции просто сохраняющие все МГД-переменные в файл, при этом структура сетки не сохраняется.

	double HLLD_Alexashov(const double& ro_L, const double& p_L, const double& v1_L, const double& v2_L, const double& v3_L,//
		const double& Bx_L, const double& By_L, const double& Bz_L, const double& ro_R, const double& p_R, const double& v1_R, const double& v2_R, const double& v3_R,//
		const double& Bx_R, const double& By_R, const double& Bz_R, double* P, const double& n1, const double& n2, const double& n3,//
		const double& rad, int metod, double& Vcon, const double& wv = 0.0);
	
	double HLLDQ_Korolkov(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L, const double& v3_L,//
		const double& Bx_L, const double& By_L, const double& Bz_L, const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& v3_R,//
		const double& Bx_R, const double& By_R, const double& Bz_R, double* P, double& PQ, const double& n1, const double& n2,//
		const double& n3, const double& rad, int metod, double& Vcon, const double& wv = 0.0);
	
	
	void HLLD_Test(void);



private:
	void initialization(void);            // Инициализируем особые ячейки-границы

};

