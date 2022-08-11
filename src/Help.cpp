
#include "Help.h"

double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& x3, const double& t3, const double& y)
// Главное значение с параметрами 2
// Строим линии между 1 и 2,  2 и 3, потом находим минмодом значение в y
{
	if (true)
	{
		double d = minmod((t1 - t2) / (x1 - x2), (t2 - t3) / (x2 - x3));
		return  (d * (y - x2) + t2);
	}
	else
	{
		// Новая процедура со сжатием
		double dUl = (t2 - t1) / (x2 - x1);
		double dUr = (t3 - t2) / (x3 - x2);
		return t2 + 0.5 * ((1.0 - eta_) * minmod(dUl, betta_ * dUr) + (1.0 + eta_) * minmod(betta_ * dUl, dUr)) * (y - x2);
	}

}

double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& y)
// Главное значение с параметрами 2
{
	double d = (t1 - t2) / (x1 - x2);
	return  (d * (y - x2) + t2);
}

double minmod(const double& x, const double& y)
{
	if (sign(x) + sign(y) == 0)
	{
		return 0.0;
	}
	else
	{
		return   ((sign(x) + sign(y)) / 2.0) * min(fabs(x), fabs(y));  ///minmod
		//return (2*x*y)/(x + y);   /// vanleer
	}
}

double sign(const double& x)
{
	if (x > 0)
	{
		return 1.0;
	}
	else if (x < 0)
	{
		return -1.0;
	}
	else
	{
		return 0.0;
	}
}