#include "Point.h"

Point::Point() : x(0.0), y(0.0), z(0.0)
{
	return;
}


Point::Point(const double& x, const double& y, const double& z) : x(x), y(y), z(z)
{
	return;
}

void Point::set(const double& x, const double& y, const double& z)
{
	this->x = x;
	this->y = y;
	this->z = z;
	return;
}

void Point::renew(void)
{
	this->x_do = this->x;
	this->y_do = this->y;
	this->z_do = this->z;
}

void Point::get(double& x, double& y, double& z)
{
	x = this->x;
	y = this->y;
	z = this->z;
	return;
}

void Point::get_do(double& x, double& y, double& z)
{
	x = this->x_do;
	y = this->y_do;
	z = this->z_do;
	return;
}

void Point::move(const double& Vx, const double& Vy, const double& Vz)
{
	this->x += Vx;
	this->y += Vy;
	this->z += Vz;
	return;
}
