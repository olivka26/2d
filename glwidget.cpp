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
    printf("%lf ", &a);
    printf("%lf ", &c);
    printf("%lf ", &b);
    printf("%lf \n", &d);
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
    k=(k-1)%8;
    change_func();
    return 0;
}

void MainWindow::allocate(){
    F=(double*)malloc((nx+2)*(ny+2)*sizeof(double));
    cx=(double*)malloc((nx+2)*sizeof(double));
    cy=(double*)malloc((ny+2)*sizeof(double));
    if(view_id){
        Fx=(double*)malloc(nx*nx*sizeof(double));
        Fy=(double*)malloc(nx*ny*sizeof(double));
        T=(double*)malloc((nx+2)*(ny+2)*sizeof(double));
        TT=(double*)malloc((nx+2)*(ny+2)*sizeof(double));
    }
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
    k=(k+1)%8;
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
    update();
}

void MainWindow::extrema_hunt(){
    if(view_id==1){
        extr[0]=TT[0];
        extr[1]=TT[0];
        for(int i=0; i<=nx+1; ++i){
            for(int j=0; j<=ny+1; ++j){
                if(TT[i*(ny+2)+j]>extr[1])
                    extr[1]=TT[i*(ny+2)+j];
                if(TT[i*(ny+2)+j]<extr[0])
                    extr[0]=TT[i*(ny+2)+j];
            }
        }
    }
    if(view_id==2){
        TT[0]=F[0]-T[0];
        TT[0]=F[0]-T[0];
        extr[0]=TT[0];
        extr[1]=TT[0];
        for(int i=0; i<=nx+1; ++i){
            for(int j=0; j<=ny+1; ++j){
                TT[i*(ny+2)+j]=F[i*(ny+2)+j]-T[i*(ny+2)+j];
                if(TT[i*(ny+2)+j]>extr[1])
                    extr[1]=TT[i*(ny+2)+j];
                if(TT[i*(ny+2)+j]<extr[0])
                    extr[0]=TT[i*(ny+2)+j];
            }
        }
    }
}

void MainWindow::func_graph(){
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0,0.0,0.0);
    printf("function started");
    for(int i=0; i<=nx+1; ++i){
        for(int j=0; j<=nx+1; ++j){
            float x= (float) cx[i];
            float y= (float) cy[j];
            float z= (float) F[i*(ny+2)*i+j];
            glVertex3f(x,y,z);
        }
    }
    for(int j=0; j<=ny+1; ++j){
        for(int i=0; i<=nx+1; ++i){
            float x= (float) cx[i];
            float y= (float) cy[j];
            float z= (float) F[i*(ny+2)*i+j];
            glVertex3f(x,y,z);
        }
    }
    printf("function ended\n");
    glEnd();
    glPopMatrix();
}

void MainWindow::appr_graph(){
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,1.0,0.0);
    printf("appr started");
    for(int i=0; i<=nx+1; ++i){
        for(int j=0; j<=nx+1; ++j){
            float x= (float) cx[i];
            float y= (float) cy[j];
            float z= (float) TT[i*(ny+2)*i+j];
            glVertex3f(x,y,z);
        }
    }
    for(int j=0; j<=ny+1; ++j){
        for(int i=0; i<=nx+1; ++i){
            float x= (float) cx[i];
            float y= (float) cy[j];
            float z= (float) TT[i*(ny+2)*i+j];
            glVertex3f(x,y,z);
        }
    }
    printf("appr ended\n");
    glEnd();
    glPopMatrix();
}

void MainWindow::err_graph(){
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,1.0,0.0);
    printf("err started");
    for(int i=0; i<=nx+1; ++i){
        for(int j=0; j<=nx+1; ++j){
            float x= (float) cx[i];
            float y= (float) cy[j];
            float z= (float) TT[i*(ny+2)*i+j];
            glVertex3f(x,y,z);
        }
    }
    for(int j=0; j<=ny+1; ++j){
        for(int i=0; i<=nx+1; ++i){
            float x= (float) cx[i];
            float y= (float) cy[j];
            float z= (float) TT[i*(ny+2)*i+j];
            glVertex3f(x,y,z);
        }
    }
    printf("err ended\n");
    glEnd();
    glPopMatrix();
}

void MainWindow::printwindow(){
    QPainter painter(this);
    painter.setPen("black");
    painter.drawText(0, 20, f_name);
    painter.drawText(10, 30, QString("format: %1").arg(view_id));
    painter.drawText(10, 45, QString("scale: %1 [%2;%3]x[%4;%5]").arg(s).arg(a).arg(b).arg(c).arg(d));
    painter.drawText(10, 60, QString("points: %1,%2").arg(nx).arg(ny));
    painter.drawText(10, 75, QString("p: %1").arg(p));
    if(k==0 && p==0 && view_id!=4)
        painter.drawText(10, 90, QString("absmax(fact): %1( %2)").arg(absmax).arg(absmax));
    else
        painter.drawText(10, 90, QString("absmax(fact): %1( %2)").arg(absmax).arg(max(fabs(extr[0]), fabs(extr[1]))));
}

MainWindow::~MainWindow(){
    free(F);
    free(cx);
    free(cy);
    if(view_id){
        free(T);
        free(Fx);
        free(Fy);
        free(TT);
    }
}

void MainWindow::paintGL(){
    setProjectionMatrix();
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    allocate();
    printf("segment: [%lf;%lf]x[%lf;%lf]\n", a,b,c,d);
    printf("points: %d %d\n", nx,ny);
    chebyshevpoints(cx, cy, nx, ny, F, a, b, c, d, f);
    printf("ok\n");
    if(p!=0){
        F[(ny+2)*(1+nx/2)+(1+ny/2)]+=(p*0.1*absmax);
    }
    if(view_id){
        fill_Fx(Fx, cx, nx, a, b);
        fill_Fy(Fy, cy, ny, c, d);
        interpolation_tensor(T, Fx, F, Fy, nx, ny);
        Fill_TT(Fx, nx, T, Fy, ny, TT);
        extrema_hunt();
    }
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
    printwindow();
    print_console();
    glDisable(GL_DEPTH_TEST);
}

void MainWindow::initializeGL(){
    glClearColor(1.0, 1.0, 1.0, 1.0);
   glEnable(GL_DEPTH_TEST);
   glShadeModel(GL_FLAT);
   glEnable(GL_CULL_FACE);
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
            if(angle_v == 360)
                angle_v = 0;
            break;
        case Qt::Key_Down:
            angle_v -= 5.0;
            if(angle_v == 360)
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
    //glLoadTransposeMatrixf(projection);
    glLoadMatrixf(projection);
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
