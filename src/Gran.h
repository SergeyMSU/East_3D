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

	Gran();
};

