#ifndef chebyshev_hpp
#define chebyshev_hpp

#include <stdio.h>
#include <cmath>

void chebyshevpoints(double *cx, double *cy, int nx, int ny, double *F, double $
double chebyshevtrans(double a, double b, double x);
void fill_chebyshev(double *xx, int nx, double a, double b, double x);
void fill_chebysheva(double *xx, int nx);
void fill_chebyshevb(double *xx, int nx);
void fill_Fx(double *Fx, double *cx, int nx, double a, double b);
void fill_Fy(double *Fy, double *cy, int ny, double c, double d);
void multiplication(double *res, double *A, int na, int m, int nb, double *B);
void interpolation_tensor(double *T, double *TT, double *Fx, double *F, double $
void horscalar(double *res, double *x, double *T, int nx, int ny);
void verscalar(double *res, double *T, double *y, int nx, int ny);
double scalar(double *x, double *y, int n);
double chebyshevvalue(double x, double y, double *T, int nx, int ny, double a, $
void copy(double *a1, double *a2, int n);
void copy(int i, double *a, double *T, int n);
void copy(double *a, double *T, int n, int i);
double scalar(double *x, double *y, double *T, int nx, int ny);

#endif /* chebyshev_hpp */

