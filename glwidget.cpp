#include "glwidget.hpp"
#include "ui_mainwindow.hpp"
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <QtWidgets/QMessageBox>

static double f0(double x, double y){
    return 1;
}

static double f1(double x, double y){
    return x;
}

static double f2(double x, double y){
    return y;
}

static double f3(double x, double y){
    return x+y;
}

static double f4(double x, double y){
    return sqrt(x*x+y*y);
}

static double f5(double x, double y){
    return x*x+y*y;
}

static double f6(double x, double y){
    return exp(x*x-y*y);
}

static double f7(double x, double y){
    return 1/(25*(x*x+y*y)+1);
}

int MainWindow::parse_command_line(int argc, char *argv[]){
    FILE *f;
    if(argc!=5)
        return -3;
    f=fopen(argv[1], "r");
    if(!f)
        return -1;
    fscanf(f, "%lf", &a);
    fscanf(f, "%lf", &c);
    fscanf(f, "%lf", &b);
    fscanf(f, "%lf", &d);
    if(fabs(a-b)<1e-6 || fabs(c-d)<1e-6)
        return -2;
    fclose(f);
    nx=atoi(argv[2]);
    if(nx<3)
        return -2;
    ny=atoi(argv[3]);
    if(ny<3)
        return -2;
    k=atoi(argv[4]);
    if(k<0 || k>7)
        return -2;
    return 0;
}

void MainWindow::allocate(){
    F=(double*)malloc(nx*ny*sizeof(double));
    Fx=(double*)malloc(nx*nx*sizeof(double));
    Fy=(double*)malloc(nx*ny*sizeof(double));
    T=(double*)malloc(nx*ny*sizeof(double));
    TT=(double*)malloc(nx*ny*sizeof(double));
    cx=(double*)malloc(nx*sizeof(double));
    cy=(double*)malloc(ny*sizeof(double));
    cx1=(double*)malloc(nx*sizeof(double));
    cy1=(double*)malloc(ny*sizeof(double));
    cx3=(double*)malloc(nx*sizeof(double));
    cy3=(double*)malloc(ny*sizeof(double));
}

void MainWindow::print_console(){
    printf("format: %d\n", view_id);
    printf("segment: [%lf;%lf]x[%lf;%lf]\n", a,b,c,d);
    printf("squeeze-stretch: %d\n", s);
    printf("points: %d %d\n", nx,ny);
    printf("disturbance: %d\n", p);
    printf("function absmax: %lf\n", absmax);
    printf("factic absmax: %lf\n", absmax);
    printf("extrema: %lf %lf\n\n", extr[0], extr[1]);
}

void MainWindow::change_func(){
    switch(k){
        case 0:
            f_name="k=0 f(x,y)=1";
            f=f0;
            min_z=1;
            max_z=1;
            absmax=1;
            break;
        case 1:
            f_name="k=1 f(x,y)=x";
            f=f1;
            min_z=a;
            max_z=b;
            absmax=max(fabs(min_z), fabs(max_z));
            break;
        case 2:
            f_name="k=2 f(x,y)=y";
            f=f2;
            min_z=c;
            max_z=d;
            absmax=max(fabs(min_z), fabs(max_z));
            break;
        case 3:
            f_name="k=3 f(x,y)=x+y";
            f=f3;
            min_z=a+c;
            max_z=b+d;
            absmax=max(fabs(min_z), fabs(max_z));
            break;
        case 4:
            f_name="k=4 f(x,y)=sqrt(x^2+y^2)";
            f=f4;
            min_z=min4(f4(a,c), f4(a,d), f4(b,c), f4(b,d));
            max_z=max4(f4(a,c), f4(a,d), f4(b,c), f4(b,d));
            if(zeroinsquare(a,b,c,d))
                min_z=0;
            absmax=max_z;
            break;
        case 5:
            f_name="k=5 f(x,y)=x^2+y^2";
            f=f5;
            min_z=min4(f5(a,c), f5(a,d), f5(b,c), f5(b,d));
            max_z=max4(f5(a,c), f5(a,d), f5(b,c), f5(b,d));
            if(zeroinsquare(a,b,c,d))
                min_z=0;
            absmax=max_z;
            break;
        case 6:
            f_name="k=6 f(x,y)=exp(x^2-y^2)";
            f=f6;
            min_z=min4(f6(a,c), f6(a,d), f6(b,c), f6(b,d));
            max_z=max4(f6(a,c), f6(a,d), f6(b,c), f6(b,d));
            absmax=max_z;
            break;
        case 7:
            f_name="k=7 f(x,y)=1/(25*(x^2+y^2)+1)";
            f=f7;
            min_z=min4(f7(a,c), f7(a,d), f7(b,c), f7(b,d));
            max_z=max4(f7(a,c), f7(a,d), f7(b,c), f7(b,d));
            if(zeroinsquare(a,b,c,d))
                max_z=1;
            absmax=max_z;
            break;
    }
}

void MainWindow::extrema_hunt(){
    double x1, x2, x3, y1, y2, y3, z11, z13, z31, z33;
    double delta_x=0.0001*(b-a)/nx;
    double delta_y=0.0001*(d-c)/ny;
    if(view_id==1){
        x1=a;
        fill_chebysheva(cx1, nx);
        for(int i=0; i<=nx; ++i){
            if(i<nx)
                x2=cx[i];
            else
                x2=b;
            for(int j=0; j<=ny; ++j){
                if(j<ny)
                    y2=cy[j];
                else
                    y2=d;
                for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                    if(j==0){
                        y1=c;
                        fill_chebysheva(cy1, ny);
                    }else{
                        y1=cy[j-1];
                        copy(cy1,T,ny,j-1);
                    }
                    z11=scalar(cx1, cy1, T, nx, ny);
                    if(z11>extr[1])
                        extr[1]=z11;
                    if(z11<extr[0])
                        extr[0]=z11;
                    fill_chebyshev(cx3, nx, a, b, x3);
                    z31=scalar(cx3, cy1, T, nx, ny);
                    if(z31>extr[1])
                        extr[1]=z31;
                    if(z31<extr[0])
                        extr[0]=z31;
                    for(y3=y1+delta_y; y3-y2<1e-6; y3+=delta_y){
                        fill_chebyshev(cy3, ny, c, d, y3);
                        z13=scalar(cx1, cy3, T, nx, ny);
                        z33=scalar(cx3, cy3, T, nx, ny);
                        if(z13>extr[1])
                            extr[1]=z13;
                        if(z13<extr[0])
                            extr[0]=z13;
                        if(z33>extr[1])
                            extr[1]=z13;
                        if(z33<extr[0])
                            extr[0]=z33;
                        y1=y3;
                        copy(cy1,cy3,ny);
                        z11=z13;
                        z31=z33;
                    }
                    x1=x3;
                    copy(cx1,cx3,nx);
                }
            }
        }
    }
    if(view_id==2){
        extr[0]=0;
        extr[1]=0;
        x1=a;
        fill_chebysheva(cx1, nx);
        for(int i=0; i<=nx; ++i){
            if(i<nx)
                x2=cx[i];
            else
                x2=b;
            for(int j=0; j<=ny; ++j){
                if(j<ny)
                    y2=cy[j];
                else
                    y2=d;
                for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                    if(j==0){
                        y1=c;
                        fill_chebysheva(cy1, ny);
                    }else{
                        y1=cy[j-1];
                        copy(cy1,T,ny,j-1);
                    }
                    z11=f(x1,y1)-scalar(cx1, cy1, T, nx, ny);
                    if(i==nx/2+1 && j==ny/2+1)
                        z11=F[ny*(nx/2)+(ny/2)]-scalar(cx1, cy1, T, nx, ny);
                    if(z11>extr[1])
                        extr[1]=z11;
                    if(z11<extr[0])
                        extr[0]=z11;
                    fill_chebyshev(cx3, nx, a, b, x3);
                    z31=f(x3,y1)-scalar(cx3, cy1, T, nx, ny);
                    if(z31>extr[1])
                        extr[1]=z31;
                    if(z31<extr[0])
                        extr[0]=z31;
                    for(y3=y1+delta_y; y3-y2<1e-6; y3+=delta_y){
                        fill_chebyshev(cy3, ny, c, d, y3);
                        z13=f(x1,y3)-scalar(cx1, cy3, T, nx, ny);
                        z33=f(x3,y3)-scalar(cx3, cy3, T, nx, ny);
                        if(z13>extr[1])
                            extr[1]=z13;
                        if(z13<extr[0])
                            extr[0]=z13;
                        if(z33>extr[1])
                            extr[1]=z13;
                        if(z33<extr[0])
                            extr[0]=z33;
                        y1=y3;
                        copy(cy1,cy3,ny);
                        z11=z13;
                        z31=z33;
                    }
                    x1=x3;
                    copy(cx1,cx3,nx);
                }
            }
        }
    }
}

void MainWindow::func_graph(){
    double x1, x2, x3, y1, y2, y3, z11, z13, z33, z31;
    double delta_x=0.0001*(b-a)/nx;
    double delta_y=0.0001*(d-c)/ny;
    x1=a;
    glBegin(GL_QUADS);
    glColor3d(1.0,0.0,0.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    printf("function started");
    for(int i=0; i<=nx; ++i){
        if(i<nx)
            x2=cx[i];
        else
            x2=b;
        for(int j=0; j<=ny; ++j){
            if(j<ny)
                y2=cy[j];
            else
                y2=d;
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                if(j==0)
                    y1=c;
                else
                    y1=cy[j-1];
                z11=f(x1, y1);
                if(i==nx/2+1 && j==ny/2+1)
                    z11=F[ny*(nx/2)+(ny/2)];
                z31=f(x3, y1);
                for(y3=y1+delta_y; y3-y2<1e-6; y3+=delta_y){
                    z13=f(x1, y3);
                    z33=f(x3, y3);
                    printf("(%lf,%lf,%lf)\n",x1,y1,z11);
                    printf("(%lf,%lf,%lf)\n",x1,y3,z13);
                    printf("(%lf,%lf,%lf)\n",x3,y1,z31);
                    printf("(%lf,%lf,%lf)\n\n",x3,y3,z33);
                    /*glVertex3d(x1, y1, z11);
                    glVertex3d(x1, y3, z13);
                    glVertex3d(x3, y1, z31);
                    glVertex3d(x3, y3, z33);*/
                    y1=y3;
                    z11=z13;
                    z31=z33;
                }
                x1=x3;
            }
        }
    }
    printf("function ended\n");
    glEnd();
}

void MainWindow::appr_graph(){
    double x1, x2, x3, y1, y2, y3, z11, z13, z33, z31;
    double delta_x=0.0001*(b-a)/nx;
    double delta_y=0.0001*(d-c)/ny;
    x1=a;
    fill_chebysheva(cx1, nx);
    glBegin(GL_QUADS);
    glColor3d(0.0,1.0,0.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    printf("appr started\n");
    for(int i=0; i<=nx; ++i){
        if(i<nx)
            x2=cx[i];
        else
            x2=b;
        for(int j=0; j<=ny; ++j){
            if(j<ny)
                y2=cy[j];
            else
                y2=d;
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                if(j==0){
                    y1=c;
                    fill_chebysheva(cy1, ny);
                }else{
                    y1=cy[j-1];
                    copy(cy1,T,ny,j-1);
                }
                z11=scalar(cx1, cy1, T, nx, ny);
                fill_chebyshev(cx3, nx, a, b, x3);
                z31=scalar(cx3, cy1, T, nx, ny);
                for(y3=y1+delta_y; y3-y2<1e-6; y3+=delta_y){
                    fill_chebyshev(cy3, ny, c, d, y3);
                    z13=scalar(cx1, cy3, T, nx, ny);
                    z33=scalar(cx3, cy3, T, nx, ny);
                    printf("(%lf,%lf,%lf)\n",x1,y1,z11);
                    printf("(%lf,%lf,%lf)\n",x1,y3,z13);
                    printf("(%lf,%lf,%lf)\n",x3,y1,z31);
                    printf("(%lf,%lf,%lf)\n\n",x3,y3,z33);
                    /*glVertex3d(x1, y1, z11);
                    glVertex3d(x1, y3, z13);
                    glVertex3d(x3, y1, z31);
                    glVertex3d(x3, y3, z33);*/
                    y1=y3;
                    copy(cy1,cy3,ny);
                    z11=z13;
                    z31=z33;
                }
                x1=x3;
                copy(cx1,cx3,nx);
            }
        }
    }
    printf("appr ended\n");
    glEnd();
}

void MainWindow::err_graph(){
    double x1, x2, x3, y1, y2, y3, z11, z13, z33, z31;
    double delta_x=0.0001*(b-a)/nx;
    double delta_y=0.0001*(d-c)/ny;
    x1=a;
    fill_chebysheva(cx1, nx);
    glBegin(GL_QUADS);
    glColor3d(0.0,0.0,1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    printf("err started\n");
    for(int i=0; i<=nx; ++i){
        if(i<nx)
            x2=cx[i];
        else
            x2=b;
        for(int j=0; j<=ny; ++j){
            if(j<ny)
                y2=cy[j];
            else
                y2=d;
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                if(j==0){
                    y1=c;
                    fill_chebysheva(cy1, ny);
                }else{
                    y1=cy[j-1];
                    copy(cy1,T,ny,j-1);
                }
                z11=f(x1,y1)-scalar(cx1, cy1, T, nx, ny);
                if(i==nx/2+1 && j==ny/2+1)
                    z11=F[ny*(nx/2)+(ny/2)]-scalar(cx1, cy1, T, nx, ny);
                fill_chebyshev(cx3, nx, a, b, x3);
                z31=f(x3,y1)-scalar(cx3, cy1, T, nx, ny);
                for(y3=y1+delta_y; y3-y2<1e-6; y3+=delta_y){
                    fill_chebyshev(cy3, ny, c, d, y3);
                    z13=f(x1,y3)-scalar(cx1, cy3, T, nx, ny);
                    z33=f(x3,y3)-scalar(cx3, cy3, T, nx, ny);
                    printf("(%lf,%lf,%lf)\n",x1,y1,z11);
                    printf("(%lf,%lf,%lf)\n",x1,y3,z13);
                    printf("(%lf,%lf,%lf)\n",x3,y1,z31);
                    printf("(%lf,%lf,%lf)\n\n",x3,y3,z33);
                    /*glVertex3d(x1, y1, z11);
                    glVertex3d(x1, y3, z13);
                    glVertex3d(x3, y1, z31);
                    glVertex3d(x3, y3, z33);*/
                    y1=y3;
                    copy(cy1,cy3,ny);
                    z11=z13;
                    z31=z33;
                }
                x1=x3;
                copy(cx1,cx3,nx);
            }
        }
    }
    printf("err ended\n");
    glEnd();
}

MainWindow::~MainWindow(){
    free(F);
    free(T);
    free(Fx);
    free(Fy);
    free(TT);
    free(cx);
    free(cy);
    free(cx1);
    free(cy1);
    free(cx3);
    free(cy3);
}

void MainWindow::paintGL(){
    setProjectionMatrix();
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    allocate();
    printf("segment: [%lf;%lf]x[%lf;%lf]\n", a,b,c,d);
    printf("points: %d %d\n", nx,ny);
    change_func();
    chebyshevpoints(cx, cy, nx, ny, F, a, b, c, d, f);
    printf("ok\n");
    extr[0]=min_z;
    extr[1]=max_z;
    if(p!=0){
        F[ny*(nx/2)+(ny/2)]+=(p*0.1*absmax);
        if(p>0 && F[ny*(nx/2)+(ny/2)]>max_z){
            max_z=F[ny*(nx/2)+(ny/2)];
            extr[1]=F[ny*(nx/2)+(ny/2)];
        }
        if(p<0 && F[ny*(nx/2)+(ny/2)]<min_z){
            min_z=F[ny*(nx/2)+(ny/2)];
            extr[0]=F[ny*(nx/2)+(ny/2)];
        }
        absmax=max(fabs(min_z), fabs(max_z));
    }
    fill_Fx(Fx, cx, nx, a, b);
    fill_Fy(Fy, cy, ny, c, d);
    interpolation_tensor(T, TT, Fx, F, Fy, nx, ny);
    extrema_hunt();
    printf("%lf;%lf",extr[0], extr[1]);
    if(view_id==0){
        printf("func\n");
        func_graph();
    }
    if(view_id==1){
        printf("appr\n");
        appr_graph();
    }
    if(view_id==2){
        printf("err\n");
        err_graph();
    }
    print_console();
    glDisable(GL_DEPTH_TEST);
}

void MainWindow::initializeGL(){
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    setDefaultCamera();
    change_func();
}

void MainWindow::resizeGL(int nWidth, int nHeight){
    glViewport(0, 0, nWidth, nHeight);
    aspect=1.0*nWidth/nHeight;
    update();
}

void MainWindow::keyPressEvent(QKeyEvent* e){
    switch (e->key()){
        case Qt::Key_0:
            k=(k+1)%8;
            break;
        case Qt::Key_1:
            view_id=(view_id+1)%3;
            break;
        case Qt::Key_2:
            a=a-(b-a)/2;
            b=b+(b-a)/2;
            c=c-(d-c)/2;
            d=d+(d-c)/2;
            ++s;
            break;
        case Qt::Key_3:
            a=a+(b-a)/4;
            b=b-(b-a)/4;
            c=c+(d-c)/4;
            d=d-(d-c)/4;
            --s;
            break;
        case Qt::Key_4:
            nx*=2;
            ny*=2;
            break;
        case Qt::Key_5:
            nx/=2;
            ny/=2;
            break;
        case Qt::Key_6:
            ++p;
            break;
        case Qt::Key_7:
            --p;
            break;
        case Qt::Key_C:
            setDefaultCamera();
            break;
        case Qt::Key_Up:
        angle_v += 5.0;
        if (angle_v == 360)
          angle_v = 0;
        break;
    case Qt::Key_Down:
        angle_v -= 5.0;
        if (angle_v == 360)
          angle_v = 0;
        break;
    case Qt::Key_Left:
        angle_h -= 5.0;
        break;
    case Qt::Key_Right:
        angle_h += 5.0;
        break;
    case Qt::Key_Plus:
        camera_p = max(camera_p - 0.1, 7);
        break;
    }
    update();
}

void MainWindow::setProjectionMatrix(){
    GLfloat view[16]={0}, projection[16]={0}, tmp[16]={0};
    GLfloat near=5, top=2, bottom, left, right;
    bottom=-top;
    right=top*aspect;
    left=-right;
    projection[0]=2*near/(right-left);
    projection[8]=(right+left)/(right-left);
    projection[5]=2*near/(top-bottom);
    projection[9]=(top+bottom)/(top-bottom);
    projection[10]=-1;
    projection[14]=-2*near;
    projection[11]=-1;
    GLfloat cam_x, cam_y, cam_z;
    cam_x=0;
    cam_y=0;
    cam_z=camera_p;
    view[0]=1;
    view[6]=-1;
    view[9]=1;
    view[15]=1;
    view[12]=-cam_x;
    view[13]=-cam_y;
    view[14]=-cam_z;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glRotated(angle_h, 0, 0, 1);
    glGetFloatv(GL_PROJECTION_MATRIX, tmp);
    glLoadTransposeMatrixf(projection);
    glMultMatrixf(view);
    glRotated(angle_h, 0, 0, 1);
    glRotated(angle_v, tmp[0], tmp[4], tmp[8]);
}

void MainWindow::setDefaultCamera(){
    camera_p=30;
    angle_h=45;
    angle_v=20;
    aspect=1.0*width()/height();
}
