#ifndef chebyshev_hpp
#define chebyshev_hpp

#include <stdio.h>
#include <cmath>

void chebyshevpoints(double *cx, double *cy, int nx, int ny, double *F, double a, double b, double c, double d, double (*f)(double, double));
double chebyshevtrans(double a, double b, double x);
void fill_Fx(double *Fx, double *cx, int nx, double a, double b);
void fill_Fy(double *Fy, double *cy, int ny, double c, double d);
void multiplication(double *res, double *A, int na, int m, int nb, double *B);
double fscalar_ij(int i, int nx, double *Fx, double *F, double *Fy, int ny, int j);
double scalar_ij(int i, int nx, double *Fx, double *T, double *Fy, int ny, int j);
double scalar_aj(int nx, double *T, double *Fy, int ny, int j);
double scalar_bj(int nx, double *T, double *Fy, int ny, int j);
double scalar_ic(int i, int nx, double *Fx, double *T, int ny);
double scalar_id(int i, int nx, double *Fx, double *T, int ny);
double scalar_ac(int nx, double *T, int ny);
double scalar_ad(int nx, double *T, int ny);
double scalar_bc(int nx, double *T, int ny);
double scalar_bd(int nx, double *T, int ny);
void interpolation_tensor(double *T, double *Fx, double *F, double *Fy, int nx, int ny);
void Fill_TT(double *Fx, int nx, double *T, double *Fy, int ny, double *TT);

#endif /* chebyshev_hpp */
