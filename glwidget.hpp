#pragma once
#include <qgl.h>
#include <QKeyEvent>
#include <qnamespace.h>
#include <QMainWindow>
#include <QPainter>
#include <QDebug>
#include <QtOpenGL>
//#include <gl/GL.h>
//#include <gl/GLU.h>
#include <QtWidgets/QWidget>
#include "chebyshev.hpp"
#include "help.hpp"

class MainWindow :public QGLWidget{
    Q_OBJECT
private:
    int nx; // number of points by X
    int ny; //number of points by Y
    int k; //id of the approxiamted function
    double a; //left end by X
    double b; //right end by X
    double c; //left end by Y
    double d; //right end by Y
    double min_z; //function's minimum
    double max_z; //function's maximum
    double absmax; //function's absolute maximum
    const char *f_name;
    double(*f)(double, double); //function
    int view_id=0; //what to be viewed
    int s=0; //squeeze-stretch
    int p=0; //disturbance
    double *F; //matrix of values
    double *Fx; //matrix of Chebyshev values by X
    double *Fy; //matrix of Chebyshev values by Y
    double *T; //matrix of interpolation coefficients
    double *TT; //intermediate matrix
    double *cx; //array of points by X
    double *cy; //array of points by Y
    double *cx1; //vector of chebyshev values in x1
    double *cy1; //vector of chebyshev values in y1
    double *cx3; //vector of chebyshev values in x3
    double *cy3; //vector of chebyshev values in y3
    double extr[2]; //factic extrema
    double B[4];
public:
    //MainWindow(QGLWidget *parent=nullptr);
    int parse_command_line(int argc, char *argv[]);
    void allocate();
    void print_console();
    void change_func();
    void extrema_hunt();
    void func_graph();
    void appr_graph();
    void err_graph();
    void printwindow();
    ~MainWindow();
protected:
    virtual void paintGL();
    virtual void initializeGL();
    virtual void resizeGL(int nWidth, int nHeight);
    virtual void keyPressEvent(QKeyEvent* e);
    void setProjectionMatrix();
    void setDefaultCamera();
    float angle_h, angle_v;
    float camera_p;
    float aspect;
};
