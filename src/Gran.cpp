#include "Gran.h"

Gran::Gran()
{

}


void Gran::set_normal(const double& a, const double& b, const double& c)
{
	this->n1 = a;
	this->n2 = b;
	this->n3 = c;
}

void Gran::get_normal(double& a, double& b, double& c)
{
	a = this->n1;
	b = this->n2;
	c = this->n3;
}