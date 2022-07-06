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

	void set(const double& x, const double& y, const double& z);           // Устанавливает значения
	void get(double& x, double& y, double& z);                             // Получает значения
	void get_do(double& x, double& y, double& z);
	void move(const double& Vx, const double& Vy, const double& Vz);       // Передвигает точку на вектор
	void renew(void);                                                      // Перезанисать старые координаты новыми

private:
};

