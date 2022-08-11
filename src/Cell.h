#pragma once
#include <vector>
#include <map>
#include <set>
#include "Help.h"
using namespace std;
class Point;
class Gran;
class Couple;

enum Cell_type  // ��� ����� ����� ��� ��������� �������
{
	C_base,                     // ������� ������ \ �������   ������� ������, ��� �� �����
	C_wall_x_max,               // ������ ������ - �����
	C_wall_x_min,               // ������ ������ - �����
	C_wall_y_max,               // ������ ������ - �����
	C_wall_y_min,               // ������ ������ - �����
	C_wall_z_max,               // ������ ������ - �����
	C_wall_z_min,               // ������ ������ - �����
	C_sphere,                   // ������ ������ - ���������� �����
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
	Point* Center = nullptr;                // ����� ������
	Couple* Par = nullptr;               // ���� ������ (��� ������ �����)
	map <int, Cell*> Can_neighbours;     // ������ - ������
	//map <int, Gran*> Grans_do;     // ����� �� ���������� ����� (����� ��� ������� ����������� �����)
	vector <Cell*> Candidates;         //  �������� ����� � �������� � �������� ����� � ��������
	vector <Gran*> Grans;         //  �������� ����� � �������� � �������� ����� � ��������
	vector <Gran*> Grans_TVD;         //  ������ - ����������������� � ������ �����
	vector <Cell*> Grans_standart;         //  ������ - ����������������� � ������ �����
	// ��� ����� �������� ��� �������������� ���������� ���������� �����, �� �� ����� ������� ������
	// ����� ����, ���� �� "��������" �� ������������ ��� �� ���������.
	// ������� ��� ����������� �������� ������, ������� ���� ������� �� �����.
	double Potok[9] = { 0.0 };

	//map <int, Gran*> mp;
	Parametr par[2];              // ���������������� ��������� � ������
	int number;                   // ����� ������
	Cell_type type;
	double move1;                   // ������ ������ ��� ������ �����
	double move2;                   // ������ ������ ��� ������ �����
	double move3;                   // ������ ������ ��� ������ �����

	mutex acces_TVD;
	mutex acces_POTOK;
	bool check_ = false;

	

	// ����������� ���������
	int n_cop = 0;                // ����� ��� �� ����� � ����
	bool mgd_ = true;             // ����� �� ������� ������ ���������� � ������ ��������� � ���� ������
	// ���� �� ���� ����� ���, ��� ����� ����� ��� ���� ���������� �����?!?!?!?
	bool couple_ = false;         // �������� �� ������ ������?  // ����� ����� ����� ��� ����� �����
	// ������� ��� ���� �������� ������� � ������ �� ����� ��������!
	bool include_ = true;         // �������� �� ������ ������ � �����?
	// ����� ����� ��������, ����� ��� ����������� ��������\��������� ������
	bool reconstruct_ = false;    // ����� �� ������������� ��� ������ (� ������ ���������� �� ��������� �������)
	// ���� �������� ������ ��������� � ��������� ���, ���� �� ������, ��� �� �����
	bool near_par = false;        // ��������� �� ������ ���������� � ������ �������
	// ����� ��� ���� �������� ��������� �������� � ��� ����� ������������ ���������!!!
	bool near_par_2 = false;      // ����� ������� ������ � ������ �������, ����� ������������ ��������

	bool i_sosed = false;         // �������� ��� ���������� ������� � ��� ������ (�� ������ ������������ � ��������� ��������)
	bool i_bbb = false;           // ��������������� ����
	bool i_boundary_1 = false;    // ����� � ��������
	bool i_boundary_2 = false;    // ����� ������� ������� � �������
	bool i_include_candidate = false;     // �������� �� ����������
	bool i_init_mgd = false;
	bool i_delete = false;                 // ����� �� ������� ������?

	bool extern_boundary = false;            // �������� �� ������ � �������� (������� �� ����������, �� �����.

	bool TVD_reconstruct = false;                 // ����� �� ������������� ���?
	

	Cell(const double& x, const double& y, const double& z);  // ����������� ��� ������ ����� ������ � ������ � ��������� ��������� �� ����
	void set_Volume(const double& V);
	void get_Volume(double& V);
	void get_Volume_do(double& V);

	// ������ � ���
	bool Get_TVD(Gran* G, Gran*& A);
	int Get_TVD_Param(Gran* G, Parametr& p1, Parametr& p2, int now, int n_gran, bool bnkl = false);
	// ���� �������� �������� �� ������ ������ �� ����� G, ��� ���� A - ����� "�����", �.�. �� ������� ����� ������


private:
	double Volume;                // ����� ������  (������� ����, ��� ���� ����� ������������� ��������� ���������� ��������
	double Volume_do;                // ����� ������

}; 

