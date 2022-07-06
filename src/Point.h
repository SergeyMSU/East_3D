#pragma once

using namespace std;



class Point
{
public:
	
	double x;
	double y;
	double z;

	double x_do;
	double y_do;
	double z_do;

	Point();
	Point(const double& x, const double& y, const double& z);

	void set(const double& x, const double& y, const double& z);           // ������������� ��������
	void get(double& x, double& y, double& z);                             // �������� ��������
	void get_do(double& x, double& y, double& z);
	void move(const double& Vx, const double& Vy, const double& Vz);       // ����������� ����� �� ������
	void renew(void);                                                      // ������������ ������ ���������� ������

private:
};

