#ifndef __kernel__
#define __kernel__

#include "numerics.h"

double logplus(double x, double y, double rl);

double logminus(double x, double y, double rl);

double ptrate(double x, double y, double rl, double rl2, double mk, double fk);

double kernV(double x, double y, double rl, double rl2, double mk, double fk);

double kernA(double x, double y, double rl, double rl2, double mk, double fk);

double kern1(double x, double y, double rl, double rl2, double mk, double fk);

double kern2(double x, double y, double rl, double rl2, double mk, double fk);

double kernVV(double x, double y, double rl, double rl2, double mk);

double kernAA(double x, double y, double rl, double rl2, double mk);

double kern11(double x, double y, double rl, double rl2, double mk);
  
double kern22(double x, double y, double rl, double rl2, double mk);

double kernA1(double x, double y, double rl, double rl2, double mk);

double kern12(double x, double y, double rl, double rl2, double mk);



















#endif
