#include "chebyshev.hpp"

#define PI 3.1415926535897932384626433832795

void chebyshevpoints(double *cx, double *cy, int nx, int ny, double *F, double a, double b, double c, double d, double (*f)(double, double)){
    double stepx=PI/nx;
    double anglex=stepx/2.0;
    double semisumx=(a+b)/2;
    double semidifx=(b-a)/2;
    double stepy=PI/ny;
    double angley=stepy/2.0;
    double semisumy=(c+d)/2;
    double semidify=(d-c)/2;
    for(int i=nx-1;i>=0;--i){
        cx[i]=cos(anglex);
        cx[i]*=semidifx;
        cx[i]+=semisumx;
        for(int j=ny-1;j>=0;--j){
            cy[j]=cos(angley);
            cy[j]*=semidify;
            cy[j]+=semisumy;
            F[i*ny+j]=f(cx[i], cy[j]);
            printf("f(%d,%d)=%lf\n", i, j, F[i*ny+j]);
            angley+=stepy;
        }
        angley=stepy/2.0;
        anglex+=stepx;
    }
}

double chebyshevtrans(double a, double b, double x){
    return (2*x-(b+a))/(b-a);
}

void fill_chebyshev(double *xx, int nx, double a, double b, double x){
    double trans=chebyshevtrans(a, b, x);
    xx[0]=1;
    xx[1]=trans;
    trans*=2;
    for(int i=2; i<nx; ++i)
        xx[i]=trans*xx[i-1]-xx[i-2];
}

void fill_chebysheva(double *xx, int nx){
    xx[0]=1;
    xx[1]=-1;
    for(int i=2; i<nx; ++i)
        xx[i]=-2*xx[i-1]-xx[i-2];
}

void fill_chebyshevb(double *xx, int nx){
    xx[0]=1;
    xx[1]=1;
    for(int i=2; i<nx; ++i)
        xx[i]=2*xx[i-1]-xx[i-2];
}


void fill_Fx(double *Fx, double *cx, int nx, double a, double b){
    for(int j=0; j<nx; ++j){
        Fx[j]=1;
        Fx[nx+j]=chebyshevtrans(a, b, cx[j]);
    }
    for(int j=0; j<nx; ++j){
        double trans=2*chebyshevtrans(a, b, cx[j]);
        for(int i=2; i<nx; ++i){
            Fx[i*nx+j]=trans*Fx[(i-1)*nx+j]-Fx[(i-2)*nx+j];
        }
    }
}

void fill_Fy(double *Fy, double *cy, int ny, double c, double d){
    for(int i=0; i<ny; ++i){
        Fy[i*ny]=1;
        Fy[i*ny+1]=chebyshevtrans(c, d, cy[i]);
    }
    for(int i=0; i<ny; ++i){
        double trans=2*chebyshevtrans(c, d, cy[i]);
        for(int j=2; j<ny; ++j){
            Fy[i*ny+j]=trans*Fy[i*ny+j-1]-Fy[i*ny+j-2];
        }
    }
}

void multiplication(double *res, double *A, int na, int m, int nb, double *B){
    for(int i=0; i<na; ++i){
        for(int j=0; j<nb; ++j){
            res[i*nb+j]=0;
            for(int k=0; k<m; ++k)
                res[i*nb+j]+=(A[i*m+k]*B[k*nb+j]);
        }
    }
}

void interpolation_tensor(double *T, double *TT, double *Fx, double *F, double *Fy, int nx, int ny){
    multiplication(TT, Fx, nx, nx, ny, F);
    multiplication(T, TT, nx, ny, ny, Fy);
    for(int i=0; i<nx; ++i){
        for(int j=0; j<ny; ++j){
            T[i*ny+j]/=nx;
            T[i*ny+j]/=ny;
            if(i)
                T[i*ny+j]*=2;
            if(j)
                T[i*ny+j]*=2;
        }
    }
    /*for(int i=0; i<nx; ++i)
        T[i*ny]/=2;
    for(int j=0; j<ny; ++j)
        T[j]/=2;*/
    //multiplication(TT, T, nx, ny, ny, Fy);
}

void horscalar(double *res, double *x, double *T, int nx, int ny){
    for(int i=0; i<ny; ++i){
        res[i]=0;
        for(int j=0; j<nx; ++j){
            res[i]+=(x[j]*T[j*ny+i]);
        }
    }
}

void verscalar(double *res, double *T, double *y, int nx, int ny){
    for(int i=0; i<nx; ++i){
        res[i]=0;
        for(int j=0; j<ny; ++j){
            res[i]+=(T[i*ny+j]*y[j]);
        }
    }

}

double scalar(double *x, double *y, int n){
    double res=0;
    for(int i=0; i<n; ++i){
        res+=(x[i]*y[i]);
    }
    return res;
}

double chebyshevvalue(double x, double y, double *T, int nx, int ny, double a, double b, double c, double d){
    double cx0=1;
    double cx1=chebyshevtrans(a, b, x);
    double tx=2*cx1;
    double cx2=tx*cx1-cx0;
    double cy0=1;
    double cy1=chebyshevtrans(c, d, y);
    double ty=2*cy1;
    double cy2=ty*cy1-cy0;
    double res=T[0];
    res+=(cx1*T[ny]);
    res+=(cy1*T[1]);
    res+=(cx1*cy1*T[ny+1]);
    res+=(cy2*T[2]);
    res+=(cx2*T[2*ny]);
    res+=(cx1*cy2*T[ny+2]);
    res+=(cx2*cy1*T[2*ny+1]);
    res+=(cx2*cy2*T[2*ny+2]);
    for(int i=3; i<nx; ++i){
        if(i%3==0){
            res+=(T[i*ny]*cx0);
            res+=(T[i*ny+1]*cx0*cy1);
            res+=(T[i*ny+2]*cx0*cy2);
        }else if(i%3==1){
            res+=(T[i*ny]*cx1);
            res+=(T[i*ny+1]*cx1*cy1);
            res+=(T[i*ny+2]*cx1*cy2);
        }else{
            res+=(T[i*ny]*cx2);
            res+=(T[i*ny+1]*cx2*cy1);
            res+=(T[i*ny+2]*cx2*cy2);
            cx0=cx2*tx-cx1;
            cx1=cx0*tx-cx2;
            cx2=cx1*tx-cx0;
        }
    }
    cx0=1;
    cx1=chebyshevtrans(a, b, x);
    cx2=tx*cx1-cx0;
    for(int j=3; j<ny; ++j){
        if(j%3==0){
            res+=(T[j]*cy0);
            res+=(T[ny+j]*cy0*cx1);
            res+=(T[2*ny+j]*cy0*cx2);
        }else if(j%3==1){
            res+=(T[j]*cy1);
            res+=(T[ny+j]*cy1*cx1);
            res+=(T[2*ny+j]*cy1*cx2);
        }else{
            res+=(T[j]*cy2);
            res+=(T[ny+j]*cy2*cx1);
            res+=(T[2*ny+j]*cy2*cx2);
            cy0=cy2*ty-cy1;
            cy1=cy0*ty-cy2;
            cy2=cy1*ty-cy0;
        }
    }
    cy0=1;
    cy1=chebyshevtrans(c, d, y);
    cy2=ty*cy1-cy0;
    for(int i=3; i<nx; ++i){
        double comp;
        if(i%3==0){
            comp=cx0;
        }else if(i%3==1){
            comp=cx1;
        }else{
            comp=cx2;
            cx0=cx2*tx-cx1;
            cx1=cx0*tx-cx2;
            cx2=cx1*tx-cx0;
        }
        for(int j=3; j<ny; ++j){
            if(j%3==0){
                comp*=cy0;
            }else if(i%3==1){
                comp*=cy1;
            }else{
                comp*=cy2;
                cy0=cx2*ty-cy1;
                cy1=cx0*ty-cy2;
                cy2=cy1*ty-cy0;
            }
            res+=(T[i*ny+j]*comp);
        }
    }
    return res;
}

void copy(double *a1, double *a2, int n){
    for(int i=0; i<n; ++i)
        a1[i]=a2[i];
}

void copy(int i, double *a, double *T, int n){
    for(int j=0; j<n; ++j)
        a[j]=T[j*n+i];
}

void copy(double *a, double *T, int n, int i){
    for(int j=0; j<n; ++j)
        a[j]=T[i*n+j];
}

double scalar(double *x, double *y, double *T, int nx, int ny){
    double res=0;
    for(int i=0; i<nx; ++i){
        for(int j=0; j<ny; ++j)
        res+=(x[i]*T[i*ny+j]*y[j]);
    }
    return res;
}
