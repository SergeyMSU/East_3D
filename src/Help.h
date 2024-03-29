#pragma once

#include <vector>
#include <string>
#include <mutex>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include "voro++.hh"

#include "Cell.h"
#include "Point.h"
#include "Setka.h"
#include "Gran.h"
#include "Couple.h"

double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& x3, const double& t3, const double& y);
double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& y);
double minmod(const double& x, const double& y);
double sign(const double& x);


#define kv(x) ( (x)*(x) )
#define kurant  0.2
#define krit  kurant

#define ga (5.0/3.0)          // ���������� ��������
#define ggg (5.0/3.0)
#define kvv(x,y,z)  (kv(x) + kv(y) + kv(z))

#define U8(ro, p, u, v, w, bx, by, bz)  ( (p) / (ggg - 1.0) + 0.5 * (ro) * kvv(u,v,w) + kvv(bx,by,bz) / 25.13274122871834590768)
#define skk(u,v,w,bx,by,bz) ( (u)*(bx) + (v)*(by) + (w)*(bz) )

#define g1 (ga - 1.0)
#define gg1 (ga - 1.0)
#define g2 (ga + 1.0)
#define gg2 (ga + 1.0)
#define gp ((g2/ga)/2.0)
#define gm ((g1/ga)/2.0)
#define gga ga


#define eps 10e-10
#define eps8 10e-8
#define pi 3.14159265358979323846
#define PI 3.14159265358979323846
#define cpi4 12.56637061435917295384
#define cpi8 25.13274122871834590768
#define spi4 ( 3.544907701811032 )
#define epsb 1e-6
#define eps_p 1e-6
#define eps_d 1e-3


#define RR_ 1.0                   // ����������� ������ ������  (�.�. ��� ����� �������� ��������� ����� ��������)

#define chi 2.0  // 32.2305
#define MA 12.9345
#define eta 150.0
#define M_inf 0.05
#define betta 0.0

#define normB true           // ����� �� �������� ���������� ���������� ���������� ���� �� ��������?
#define TVD_ true
#define  tens_ 0.003            // ����������� ��������� �����������

#define eta_ 0.5    // �������� ��������
#define betta_ 1.0  // ������


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}



