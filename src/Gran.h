#pragma once

#include "Help.h"
#include "Cell.h"

using namespace std;

class Cell;

class Gran
{
public:
	Cell* Sosed;
	double S;                // ������� �����     ������� � ����� �� �����, �.�. � ����� ���������, ��������� ���������� ������
	double n1 = 0.0;         // ������� - ������� �����, ���� �� ��������� �������, ������ ���� ������ ����� �� �������, �� ����� ������
	double n2 = 0.0;
	double n3 = 0.0;
	double c1 = 0.0;         // ����������� ����� ������������ ����������� ��������� 
	double c2 = 0.0;
	double c3 = 0.0;
	double dl = 0.0;
	bool check_ = false;

	Gran();

	void set_normal(const double& a, const double& b, const double& c);
	void get_normal(double& a, double& b, double& c);
};

