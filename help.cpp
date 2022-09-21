#include <stdio.h>
#include <iostream>
#include "help.hpp"

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



