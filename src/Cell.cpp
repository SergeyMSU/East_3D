#include "Cell.h"


Cell::Cell(const double& x, const double& y, const double& z)
{
	this->Center = new Point(x, y, z);
	this->Grans.reserve(15);
	this->Candidates.reserve(50);
}

void Cell::set_Volume(const double& V)
{
	this->Volume_do = this->Volume;
	this->Volume = V;
}

void Cell::get_Volume(double& V)
{
	V = this->Volume;
}

void Cell::get_Volume_do(double& V)
{
	V = this->Volume_do;
}