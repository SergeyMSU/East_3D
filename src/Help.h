#pragma once

#include <vector>
#include <string>
#include <mutex>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "voro++.hh"

#include "Cell.h"
#include "Point.h"
#include "Setka.h"
#include "Gran.h"
#include "Couple.h"


#define kv(x) ( (x)*(x) )
#define kurant  0.2
#define krit  kurant

#define ga (5.0/3.0)          // ѕоказатель адиабаты
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


#define RR_ 1.0                   // ’арактерный размер задачи  (т.е. все сетка занимает несколько таких размеров)

#define chi 2.0
#define MA 12.0
#define M_inf 0.05
#define betta 4.2426


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}
