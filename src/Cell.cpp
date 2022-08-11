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

bool Cell::Get_TVD(Gran* G, Gran*& A)
{
	if (G == nullptr)
	{
		return false;
	}


	if (G->Sosed->type == C_base)
	{
		double v1, v2, v3;
		v1 = -G->n1;
		v2 = -G->n2;
		v3 = -G->n3;
		
		double u1, u2, u3, dd;
		double M = 0.0;
		for (auto& i: this->Grans)
		{
			if (i->Sosed->number == G->Sosed->number || i->Sosed->type != C_base)
			{
				continue;
			}

			i->get_normal(u1, u2, u3);

			dd = u1 * v1 + u2 * v2 + u3 * v3;
			if (dd > M)
			{
				M = dd;
				A = i;
			}
		}

		if (M < 0.5)
		{
			A = nullptr;
			return false;
		}
		else
		{
			return true;
		}
	}
	else
	{
		A = nullptr;
		return false;
	}

	return true;
}

int Cell::Get_TVD_Param(Gran* G, Parametr& p1, Parametr& p2, int now, int n_gran, bool bnkl)
{
	Gran* A1 = nullptr;
	double x, y, z;
	this->Center->get(x, y, z);
	double x1, y1, z1;
	G->Sosed->Center->get(x1, y1, z1);
	double x2, y2, z2;
	double a, b;
	a = sqrt(kvv(x1 - x, y1 - y, z1 - z));
	bool bbb = false;
	//cout << " A " << endl;
	if (this->Grans_TVD[n_gran] != nullptr)
	{
		if (this->TVD_reconstruct == false)
		{
			bbb = true;
			A1 = this->Grans_TVD[n_gran];
		}
		else
		{
			this->Grans_TVD[n_gran]->Sosed->Center->get(x2, y2, z2);
			double u1, u2, u3;
			double v1, v2, v3;
			u1 = G->n1;
			u2 = G->n2;
			u3 = G->n3;
			v1 = -this->Grans_TVD[n_gran]->n1;
			v2 = -this->Grans_TVD[n_gran]->n2;
			v3 = -this->Grans_TVD[n_gran]->n3;

			if (u1 * v1 + u2 * v2 + u3 * v3 > 0.9)
			{
				bbb = true;
				A1 = this->Grans_TVD[n_gran];
				//cout << "Uspeshno " << endl;
			}
		}
	}
	else
	{
		if (this->TVD_reconstruct == false)
		{
			goto afg;
		}
	}
	//cout << " B " << endl;

	if (bbb == false)
	{
		if (bnkl == true)
		{
			bbb = (this->Grans_TVD[n_gran] != nullptr);
			cout << bbb << "    " << x1 << " " << y1 << " " << z1 << "    " << x << " " << y << " " << z << endl;
			cout << x2 << " " << y2 << " " << z2 << endl;
			cout << this->number << endl << endl;
		}

		bbb = this->Get_TVD(G, A1);
		this->acces_TVD.lock();
		this->Grans_TVD[n_gran] = nullptr;
		this->acces_TVD.unlock();
		//cout << "NEusp " << endl;
	}
	//cout << " C " << endl;
	// —носим значени€ со стороны €чейки!
	if (bbb == true)
	{
		this->acces_TVD.lock();
		this->Grans_TVD[n_gran] = A1;
		this->acces_TVD.unlock();
		A1->Sosed->Center->get(x2, y2, z2);
		b = sqrt(kvv(x2 - x, y2 - y, z2 - z));

		p1.ro = linear(-(a/2.0 + b), A1->Sosed->par[now].ro, -a/2.0, this->par[now].ro, a/2.0, G->Sosed->par[now].ro, 0.0);
		p1.p = linear(-(a / 2.0 + b), A1->Sosed->par[now].p, -a / 2.0, this->par[now].p, a / 2.0, G->Sosed->par[now].p, 0.0);
		p1.u = linear(-(a / 2.0 + b), A1->Sosed->par[now].u, -a / 2.0, this->par[now].u, a / 2.0, G->Sosed->par[now].u, 0.0);
		p1.v = linear(-(a / 2.0 + b), A1->Sosed->par[now].v, -a / 2.0, this->par[now].v, a / 2.0, G->Sosed->par[now].v, 0.0);
		p1.w = linear(-(a / 2.0 + b), A1->Sosed->par[now].w, -a / 2.0, this->par[now].w, a / 2.0, G->Sosed->par[now].w, 0.0);
		p1.bx = linear(-(a / 2.0 + b), A1->Sosed->par[now].bx, -a / 2.0, this->par[now].bx, a / 2.0, G->Sosed->par[now].bx, 0.0);
		p1.by = linear(-(a / 2.0 + b), A1->Sosed->par[now].by, -a / 2.0, this->par[now].by, a / 2.0, G->Sosed->par[now].by, 0.0);
		p1.bz = linear(-(a / 2.0 + b), A1->Sosed->par[now].bz, -a / 2.0, this->par[now].bz, a / 2.0, G->Sosed->par[now].bz, 0.0);

		if (p1.ro <= 0.0)
		{
			p1.ro = this->par[now].ro;
		}
		if (p1.p <= 0.0)
		{
			p1.p = this->par[now].p;
		}
	}
	else
	{
		afg:
		p1.ro = this->par[now].ro;
		p1.p = this->par[now].p;
		p1.u = this->par[now].u;
		p1.v = this->par[now].v;
		p1.w = this->par[now].w;
		p1.bx = this->par[now].bx;
		p1.by = this->par[now].by;
		p1.bz = this->par[now].bz;
	}
	//cout << " D " << endl;
	// —носим с другой стороны
	G->Sosed->acces_TVD.lock();
	G->Sosed->Grans_TVD.resize(G->Sosed->Grans.size(), nullptr);
	G->Sosed->acces_TVD.unlock();
	n_gran = -1;
	Gran* G2 = nullptr;
	bool bn = false;
	for (auto& i : G->Sosed->Grans)
	{
		n_gran++;
		if (i->Sosed->number == this->number)
		{
			G2 = i;
			bn = true;
			break;
		}
	}
	//cout << " F " << endl;
	if (bn == true)
	{
		bbb = false;
		if (G->Sosed->Grans_TVD[n_gran] != nullptr)
		{
			if (G->Sosed->TVD_reconstruct == false)
			{
				bbb = true;
				A1 = G->Sosed->Grans_TVD[n_gran];
			}
			else
			{
				G->Sosed->Grans_TVD[n_gran]->Sosed->Center->get(x2, y2, z2);
				double u1, u2, u3;
				double v1, v2, v3;
				u1 = G2->n1;
				u2 = G2->n2;
				u3 = G2->n3;
				v1 = -G->Sosed->Grans_TVD[n_gran]->n1;
				v2 = -G->Sosed->Grans_TVD[n_gran]->n2;
				v3 = -G->Sosed->Grans_TVD[n_gran]->n3;
				
				if (u1 * v1 + u2 * v2 + u3 * v3 > 0.9)
				{
					bbb = true;
					A1 = G->Sosed->Grans_TVD[n_gran];
				}
			}
		}
		else
		{
			if (G->Sosed->TVD_reconstruct == false)
			{
				goto ak;
			}
		}

		if (bbb == false)
		{
			bbb = G->Sosed->Get_TVD(G2, A1);
			G->Sosed->acces_TVD.lock();
			G->Sosed->Grans_TVD[n_gran] = nullptr;
			G->Sosed->acces_TVD.unlock();
		}

		if (bbb == true)
		{
			G->Sosed->acces_TVD.lock();
			G->Sosed->Grans_TVD[n_gran] = A1;
			G->Sosed->acces_TVD.unlock();
			A1->Sosed->Center->get(x2, y2, z2);
			b = sqrt(kvv(x2 - x1, y2 - y1, z2 - z1));

			p2.ro = linear(-(a / 2.0 + b), A1->Sosed->par[now].ro, -a / 2.0, G->Sosed->par[now].ro, a / 2.0, this->par[now].ro, 0.0);
			p2.p = linear(-(a / 2.0 + b), A1->Sosed->par[now].p, -a / 2.0, G->Sosed->par[now].p, a / 2.0, this->par[now].p, 0.0);
			p2.u = linear(-(a / 2.0 + b), A1->Sosed->par[now].u, -a / 2.0, G->Sosed->par[now].u, a / 2.0, this->par[now].u, 0.0);
			p2.v = linear(-(a / 2.0 + b), A1->Sosed->par[now].v, -a / 2.0, G->Sosed->par[now].v, a / 2.0, this->par[now].v, 0.0);
			p2.w = linear(-(a / 2.0 + b), A1->Sosed->par[now].w, -a / 2.0, G->Sosed->par[now].w, a / 2.0, this->par[now].w, 0.0);
			p2.bx = linear(-(a / 2.0 + b), A1->Sosed->par[now].bx, -a / 2.0, G->Sosed->par[now].bx, a / 2.0, this->par[now].bx, 0.0);
			p2.by = linear(-(a / 2.0 + b), A1->Sosed->par[now].by, -a / 2.0, G->Sosed->par[now].by, a / 2.0, this->par[now].by, 0.0);
			p2.bz = linear(-(a / 2.0 + b), A1->Sosed->par[now].bz, -a / 2.0, G->Sosed->par[now].bz, a / 2.0, this->par[now].bz, 0.0);

			if (p2.ro <= 0.0)
			{
				p2.ro = G->Sosed->par[now].ro;
			}
			if (p2.p <= 0.0)
			{
				p2.p = G->Sosed->par[now].p;
			}
		}
		else
		{
			goto ak;
		}
	}
	else
	{
		ak:
		p2.ro = G->Sosed->par[now].ro;
		p2.p = G->Sosed->par[now].p;
		p2.u = G->Sosed->par[now].u;
		p2.v = G->Sosed->par[now].v;
		p2.w = G->Sosed->par[now].w;
		p2.bx = G->Sosed->par[now].bx;
		p2.by = G->Sosed->par[now].by;
		p2.bz = G->Sosed->par[now].bz;
		return -1;
	}

	return n_gran;
}