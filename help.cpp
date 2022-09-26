#include "help.hpp"

double f0(double x, double y){
    return 1;
}

double f1(double x, double y){
    return x;
}

double f2(double x, double y){
    return y;
}

double f3(double x, double y){
    return x+y;
}

double f4(double x, double y){
    return sqrt(x*x+y*y);
}

double f5(double x, double y){
    return x*x+y*y;
}

double f6(double x, double y){
    return exp(x*x-y*y);
}

double f7(double x, double y){
    return 1/(25*(x*x+y*y)+1);
}

double max4(double a, double b, double c, double d){
    double res=a;
    if(res<b)
        res=b;
    if(res<c)
        res=c;
    if(res<d)
        res=d;
    return res;
}

double min4(double a, double b, double c, double d){
    double res=a;
    if(res>b)
        res=b;
    if(res>c)
        res=c;
    if(res>d)
        res=d;
    return res;
}

double max(double a, double b){
    if(a>b)
        return a;
    return b;
}

bool zeroinsquare(double a, double b, double c, double d){
    if(b>1e-6 && a<-(1e-6) && d>1e-6 && c<-(1e-6))
        return true;
    return false;
}

double** allocate_matrix(int n, int m){ //n rows, m columns
    double** a=new double*[n];
    for(int i=0; i<n; ++i)
        a[i]=new double[m];
    return a;
}

void delete_matrix(double ** a, int n){
    for(int i=0; i<n; ++i){
        delete[] a[i];
    }
    delete[] a;
}



