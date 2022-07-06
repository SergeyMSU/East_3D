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
	vector <Cell*> All_Cell;              // ������ ���� ����� �����
	vector <Couple*> All_Couple;              // ������ ���� ����� �����
	list <Cell*> All_Cell_off;             // ������ ����������� �����

	// ��� ������ ������ ��������� � ����������� ������� inicializtion
	Cell* Wall1;                          // ������ ������ - �������       -  C_wall_x_max
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
	double sl_L = 0.6;                    // �� ������ ������� ���� ������������   0.3
	double sl_R = 1.5;

	Setka();                              // ��������� ���������� ����� �������


	void Initialization_do_MHD(void);

	void Construct_start(void);                 // ��������� ���������� ����� (� �������������� ����������), �������� �����
	// ��� ������� ����� ������������ ������ ��� �������������� ���������� �����, ������ ����� ������������ ����� ������� ���������.
	// ��� ������� ���� ������� ��� ������ ������ (������� ��������� �������), ����� ���������� �� � "�����" (������� � ����� ������)
	void Reconstruct_couple(bool t1 = false);
	// ���������� ����� � ������ ��������� ������� � ������, ����� ����������� ����������� ������ (��� �� ��������)
	// ������������ ���������� �� ��� (��� ������� ������� ����� ���������� �������� �����
	// ����� ����������� ������-������� ��� ������
	// t1 - ����� �� ��� ������ ����� ������� ������ ������ (������������ ����������� �� ��������� ����������)?
	void Reconstruct_fast(bool wide = true);  // ������� ������������ �����
	// ������������� ����� �� ���������� �������, ����� �� ���������
	// �.�. ������� �������� ��� ��������� ��������, ����� ���������� �� ������ ��� ��� ����� �������� �������.
	void Reconstruct_medium(bool wide = true);
	// ������������� ����� ���������� �������-�������
	// ��� ������ �����, ��� ���� �� ������ ��������� ������ ��, ������� � ��� ���� �� �����
	void Reconstruct_medium2(bool wide = true);
	void Reconstruct_medium3(bool wide = true);
	void renumber(void);

	void include_cells();  // �������� ������
	void exclude_cells();  // ��������� ������

	/// ------------------------------------------------------------------------------------------------------------------------
	///  ���������� ������ � �� ��� � ���� �������
	/// ------------------------------------------------------------------------------------------------------------------------
	void Add_couple_start(void);                 // ��������� ���������� ��� � ����� (������ ������� �� ���)
	void Disable_cells(void);                    // ���������� �����, ������� ������� � �����
	// ����������� �� ���� ������� � ��� � ��������� �� � ������, ���� ��� ������� ������
	// "���������" - ��� ������ ������� �� � ������ ����������� � ������������� ������, ������� ��� ��� �� ���������!
	void Calc_normal(void);                    // ��������� ����� ������� � ������ ����� (�� �������)  ???????????????????????????????????????????????????????????????????????????????????????????
	// ����� ��������������� ����
	
	void move_par(void);                       // ����������� ���� (�� �� �������) � ������������ � �������� �������������  ??????????????????????????????????????????????????????????????????????
	void set_normal(void);                     // ������ ��������� ������� ������ ���� � ������������� � (��� ��������� ��������� ��������)
	void Cut_Surface(void);            // ������ ����������� ����� ������


	// ����� ��������� �����
	void Print_points(void);              // ������� ��� ����� (������ ��� ���������� ������ �����) - ����� ����� �� �������� � ��������

	void Cut_Plane_long_z(void);            // ������ ���� ������ ���������� z = 0
	// ������ ���������� ����� ��������� ���� ����� � ��������� Voro++

	void Cut_Plane_z(double R = 0.0);            // ������ ���� ������ ���������� z = 0
	// ������� �������, ���������� ��� ��������� ������� ������ ������
	// R - ����������, �� ������� ������ ��� �� ����.
	void Cut_Plane_y(double R = 0.0);            // ������ ���� ������ ���������� z = 0

	void Cut_Plane_z_Tecplot(double R = 0.0);    // ������ ���������������� ���������� � ���� � ��������� z = 0
	// ������� ����� ������ ����������, ����� ������� ����� - �������� �����, � � ���� ����������� �������� ���������������� ����������

	void Cut_Plane_y_Tecplot(double R = 0.0);    // ������ ���������������� ���������� � ���� � ��������� y = 0





	// �������������� ���������

	void dekard_skorost(double x, double y, double z, double Vr, double Vphi, double Vtheta, double& Vx, double& Vy, double& Vz);
	double polar_angle(double x, double y);



	// ���
	void Init(void);         // ������������� ���������������� ����������
	void Go_MHD(int times);  // ��������� ���� ��� �� ������������ �����  times - �����
	void Start_MHD(int times);   // ��������� ��������� �� �������������� �����!!!  �������� �������
	void Culc_couple(int now1, const double& time);    // ������� �������� ��� �� ������ � ������� ������� � ������� ��!!!

	void Save_MHD(string name_f = "Save_MHD.txt");
	void Download_MHD(string name_f = "Save_MHD.txt");
	// ������� ������� ������ ����������� ��� ���-���������� � ����, ��� ���� ��������� ����� �� �����������.

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
	void initialization(void);            // �������������� ������ ������-�������

};

