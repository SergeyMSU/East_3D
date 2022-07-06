#include "Couple.h"

Couple::Couple(Cell* A1, Cell* A2, const double& dist) : n1(0.0), n2(0.0), n3(0.0)
{
	this->A1 = A1;
	this->A2 = A2;
	this->dist = dist;
}

void Couple::Resolve(void)
{
	double x1, y1, z1, x2, y2, z2;

	this->A1->Center->get(x1, y1, z1);
	this->A2->Center->get(x2, y2, z2);

	this->dist = sqrt(kv(x2 - x1) + kv(y2 - y1) + kv(z2 - z1));
}

void Couple::get_centr(double& x, double& y, double& z)
{
	x = (this->A1->Center->x + this->A2->Center->x) / 2.0;
	y = (this->A1->Center->y + this->A2->Center->y) / 2.0;
	z = (this->A1->Center->z + this->A2->Center->z) / 2.0;
}

void Couple::get_normal(double& x, double& y, double& z)
{
	double x1, x2, y1, y2, z1, z2;
	this->A1->Center->get(x1, y1, z1);
	this->A2->Center->get(x2, y2, z2);
	x = (x2 - x1) / this->dist;
	y = (y2 - y1) / this->dist;
	z = (z2 - z1) / this->dist;
}

void Couple::move(const double& m1, const double& m2, const double& m3)
{
	this->A1->Center->move(m1, m2, m3);
	this->A2->Center->move(m1, m2, m3);
}

void Couple::orient(void)
{
	double cx, cy, cz;
	this->get_centr(cx, cy, cz);
	this->A1->Center->set(cx - this->n1 * this->dist / 2.0, cy - this->n2 * this->dist / 2.0, cz - this->n3 * this->dist / 2.0);
	this->A2->Center->set(cx + this->n1 * this->dist / 2.0, cy + this->n2 * this->dist / 2.0, cz + this->n3 * this->dist / 2.0);
}
