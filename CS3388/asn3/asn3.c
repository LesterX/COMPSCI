/* Computer Science 3388
   Assignment 3
   Yimin Xu
   250876566 */

#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

#define FALSE 0
#define TRUE  1
#define COLS  1024
#define ROWS  1024

#define POSX  0
#define POSY  0

#define Ex -20.0
#define Ey -20.0
#define Ez 20.0

#define Gx 0.0
#define Gy 0.0
#define Gz 0.0

#define UPx 0.0
#define UPy 0.0
#define UPz 1.0

#define NP 5.0
#define FP 50.0

#define THETA 90.0
#define ASPECT 1.0

#define W  1024
#define H  1024

#define LX 30.0
#define LY 30.0
#define LZ 30.0

dmatrix_t *build_camera_matrix(dmatrix_t *E, dmatrix_t *G) {
    
    dmatrix_t N ; /* Viewing axis */
    
    N = *dmat_normalize(dmat_sub(E,G)) ;
    N.l = 3 ;

    dmatrix_t UP ;
    dmat_alloc(&UP,4,1) ;
    UP.l = 3 ;
    
    UP.m[1][1] = UPx ;
    UP.m[2][1] = UPy ;
    UP.m[3][1] = UPz ;
    UP.m[4][1] = 1.0 ;
    
    dmatrix_t U ;
    
    U = *dmat_normalize(dcross_product(&UP,&N)) ;
    
    dmatrix_t V ;
    V = *dcross_product(&N,&U) ;
    
    dmatrix_t Mv ; /* Build matrix M_v */
    dmat_alloc(&Mv,4,4) ;
    
    Mv.m[1][1] = U.m[1][1] ; 
    Mv.m[1][2] = U.m[2][1] ; 
    Mv.m[1][3] = U.m[3][1] ; 
    Mv.m[1][4] = -1.0*((*E).m[1][1]*U.m[1][1] + (*E).m[2][1]*U.m[2][1] + (*E).m[3][1]*U.m[3][1]) ;
    
    Mv.m[2][1] = V.m[1][1] ; 
    Mv.m[2][2] = V.m[2][1] ; 
    Mv.m[2][3] = V.m[3][1] ; 
    Mv.m[2][4] = -1.0*((*E).m[1][1]*V.m[1][1] + (*E).m[2][1]*V.m[2][1] + (*E).m[3][1]*V.m[3][1]) ;
    
    Mv.m[3][1] = N.m[1][1] ; 
    Mv.m[3][2] = N.m[2][1] ; 
    Mv.m[3][3] = N.m[3][1] ; 
    Mv.m[3][4] = -1.0*((*E).m[1][1]*N.m[1][1] + (*E).m[2][1]*N.m[2][1] + (*E).m[3][1]*N.m[3][1]) ;
    
    Mv.m[4][1] = 0.0 ; 
    Mv.m[4][2] = 0.0 ; 
    Mv.m[4][3] = 0.0 ; 
    Mv.m[4][4] = 1.0 ;
    
    dmatrix_t Mp ; /* Build matrix Mp */
    dmat_alloc(&Mp,4,4) ;
    Mp = *dmat_identity(&Mp) ;
    
    float a = -1.0*(FP + NP)/(FP - NP) ;
    float b = -2.0*(FP*NP)/(FP - NP) ;
    
    Mp.m[1][1] = NP ;
    Mp.m[2][2] = NP ;
    Mp.m[3][3] = a ;
    Mp.m[3][4] = b ;
    Mp.m[4][3] = -1.0 ;
    Mp.m[4][4] = 0.0 ;
    
    /* Build matrices T_1 and S_1 */
    
    /* Work out coordinates of near plane corners */
    
    float top = NP*tan(M_PI/180.0*THETA/2.0) ;
    float right = ASPECT*top ;
    float bottom = -top ;
    float left = -right ;
   
    dmatrix_t T1 ;
    dmat_alloc(&T1,4,4) ;
    
    T1 = *dmat_identity(&T1) ;
    T1.m[1][4] = -(right + left)/2.0 ;
    T1.m[2][4] = -(top + bottom)/2.0 ;

    dmatrix_t S1 ;
    dmat_alloc(&S1,4,4) ;
    
    S1 = *dmat_identity(&S1) ;
    S1.m[1][1] = 2.0/(right - left) ;
    S1.m[2][2] = 2.0/(top - bottom) ;

    /* Build matrices T2, S2, and W2 */
    
    dmatrix_t T2 ;
    dmatrix_t S2 ;
    dmatrix_t W2 ;
    
    dmat_alloc(&T2,4,4) ;
    dmat_alloc(&S2,4,4) ;
    dmat_alloc(&W2,4,4) ;
    
    T2 = *dmat_identity(&T2) ;
    S2 = *dmat_identity(&S2) ;
    W2 = *dmat_identity(&W2) ;
    
    T2.m[1][4] = 1.0 ;
    T2.m[2][4] = 1.0 ;

    S2.m[1][1] = W/2.0 ;
    S2.m[2][2] = H/2.0 ;
    
    W2.m[2][2] = -1.0 ;
    W2.m[2][4] = (double)H ;
    
    return dmat_mult(&W2,dmat_mult(&S2,dmat_mult(&T2,dmat_mult(&S1,dmat_mult(&T1,dmat_mult(&Mp,&Mv)))))) ;
}

dmatrix_t *perspective_projection(dmatrix_t *P) {

    (*P).m[1][1] /= (*P).m[4][1] ;
    (*P).m[2][1] /= (*P).m[4][1] ;
    (*P).m[3][1] /= (*P).m[4][1] ;
    (*P).m[4][1] /= (*P).m[4][1] ;

    return P ;
}


typedef struct {
    int x, y ;
} intpoint ;

Display *InitX(Display *d, Window *w, int *s) {
    
    d = XOpenDisplay(NULL) ;
    if(d == NULL) {
        printf("Cannot open display\n") ;
        exit(1) ;
    }
    *s = DefaultScreen(d) ;
    *w = XCreateSimpleWindow(d, RootWindow(d, *s), POSX, POSY, COLS, ROWS, 1, BlackPixel(d, *s), WhitePixel(d, *s)) ;
    Atom delWindow = XInternAtom(d, "WM_DELETE_WINDOW", 0) ;
    XSetWMProtocols(d, *w, &delWindow, 1) ;
    XSelectInput(d, *w, ExposureMask | KeyPressMask) ;
    XMapWindow(d, *w) ;
    return(d) ;
}

void SetCurrentColorX(Display *d, GC *gc, unsigned int r, unsigned int g, unsigned int b) {
    
    XSetForeground(d, *gc, r << 16 | g << 8 | b) ;
}

void SetPixelX(Display *d, Window w, int s, int i, int j) {
    
    XDrawPoint(d, w, DefaultGC(d, s), i, j) ;
}

void QuitX(Display *d, Window w) {
    
    XDestroyWindow(d,w) ;
    XCloseDisplay(d) ;
}

void exchangeInt(int *a, int *b)

{ int t ;
    
    t = *a ;
    *a = *b ;
    *b = t ;
}

void set_xyz(dmatrix_t p, double x, double y, double z)
{
    p.m[1][1] = x;
    p.m[2][1] = y;
    p.m[3][1] = z;
    p.m[4][1] = 1;
}

dmatrix_t get_line(dmatrix_t p1, dmatrix_t p2)
{
    dmatrix_t p;
    dmat_alloc(&p,4,1);
    p.m[1][1] = p2.m[1][1] - p1.m[1][1];
    p.m[2][1] = p2.m[2][1] - p1.m[2][1];
    p.m[3][1] = p2.m[3][1] - p1.m[3][1];
    p.m[4][1] = 1;

    return p;
}

dmatrix_t find_normal(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
{
    double l1x = x2 - x1;
    double l1y = y2 - y1;
    double l1z = z2 - z1;

    double l2x = x3 - x2;
    double l2y = y3 - y2;
    double l2z = z3 - z2;

    dmatrix_t result;
    dmat_alloc(&result,4,1);
    
    result.m[1][1] = 10000 * (l1y * l2z - l1z * l2y);
    result.m[2][1] = 10000 * (l1z * l2x - l1x * l2z);
    result.m[3][1] = 10000 * (l1x * l2y - l1y * l2x);
    result.m[4][1] = 1;

    return result;
}

double find_angle(dmatrix_t l1, dmatrix_t l2)
{
    double x1 = l1.m[1][1];
    double y1 = l1.m[2][1];
    double z1 = l1.m[3][1];
    double x2 = l2.m[1][1];
    double y2 = l2.m[2][1];
    double z2 = l2.m[3][1];

    return (x1 * x2 + y1 * y2 + z1 * z2) / (sqrt(x1 * x1 + y1 * y1 + z1 * z1) * sqrt(x2 * x2 + y2 * y2 + z2 * z1));
}

int backside(dmatrix_t E, dmatrix_t normal)
{
    if (find_angle(E,normal) < 0)
        return 1;
    else
        return 0;
}

double min(double n1, double n2)
{
    if (n1 < n2) return n1;
    else return n2;
}

double max(double n1, double n2)
{
    if (n1 > n2) return n1;
    else return n2;
}

int intersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
    double slope1 = (y2 - y1) / (x2 - x1);
    double intercept1 = y1 - x1 * slope1;
    double slope2 = (y4 - y3) / (x4 - x3);
    double intercept2 = y3 - x3 * slope2;
    double intersect_x = (intercept2 - intercept1) / (slope1 - slope2);
    double intersect_y = ((intercept2 - intercept1) / (slope1 - slope2)) * slope1 + intercept1;

    if (intersect_x > min(x1,x2) && intersect_x > min(x3,x4) && intersect_x < max(x1,x2) && intersect_x < max(x3,x4) && intersect_y > min(y1,y2) && intersect_y > min (y3,y4) && intersect_y < max(y1,y2) && intersect_y < max(y3,y4))
        return 1;
    else
        return 0;
}

int inside_polygon(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3)
{
    double center_x = (x1 + x3) / 2;
    double center_y = (y1 + y3) / 2;

    if (x0 == center_x && y0 == center_y)
        return 1;

    double x4 = x1 + x3 - x2;
    double y4 = y1 + y3 - y2;
    int count = 0;
    double arr_x[4] = {x1,x2,x3,x4};
    double arr_y[4] = {y1,y2,y3,y4};

    for (int i = 0; i < 4; i ++)
    {
        int j;
        if (i == 3)
            j = 0;
        else
            j = i + 1;

        if (intersect(x0,y0,center_x,center_y,arr_x[i],arr_y[i],arr_x[j],arr_y[j]))
            count ++;
    }

    if (count % 2 == 1)
        return 1;
    else
        return 0;
}

double max_3(double n1, double n2, double n3)
{
    if (n1 > n2 && n1 > n3)
        return n1;
    else if (n2 > n3)
        return n2;
    else
        return n3;
}

double min_3(double n1, double n2, double n3)
{
    if (n1 < n2 && n1 < n3)
        return n1;
    else if (n2 < n3)
        return n2;
    else
        return n3;
}

void set_color(Display* d, Window w, int s, dmatrix_t C, unsigned int r, unsigned int g, unsigned int b, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double intensity)
{
    SetCurrentColorX(d, &(DefaultGC(d, s)), r * intensity, g * intensity, b * intensity);

    dmatrix_t p1;
    dmat_alloc(&p1,4,1);

    p1.m[1][1] = x1;
    p1.m[2][1] = y1;
    p1.m[3][1] = z1;
    p1.m[4][1] = 1.0;
    p1 = *perspective_projection(dmat_mult(&C,&p1));

    dmatrix_t p2;
    dmat_alloc(&p2,4,1);

    p2.m[1][1] = x2;
    p2.m[2][1] = y2;
    p2.m[3][1] = z2;
    p2.m[4][1] = 1.0;
    p2 = *perspective_projection(dmat_mult(&C,&p2));    

    dmatrix_t p3;
    dmat_alloc(&p3,4,1);

    p3.m[1][1] = x1;
    p3.m[2][1] = y1;
    p3.m[3][1] = z1;
    p3.m[4][1] = 1.0;
    p3 = *perspective_projection(dmat_mult(&C,&p3));

    double x_min = min_3(p1.m[1][1],p2.m[1][1],p3.m[1][1]);
    double x_max = max_3(p1.m[1][1],p2.m[1][1],p3.m[1][1]);
    double y_min = min_3(p1.m[2][1],p2.m[2][1],p3.m[2][1]);
    double y_max = max_3(p1.m[2][1],p2.m[2][1],p3.m[2][1]);

    //printf("x1 = %f, y1 = %f, x2 = %f, y2 = %f, x3 = %f, y3 = %f\n", p1.m[1][1],p1.m[2][1],p2.m[1][1],p2.m[2][1],p3.m[1][1],p3.m[2][1]);

    for (int x = (int) x_min; x < x_max; x ++)
    {
        for (int y = (int) y_min; y < y_max; y ++)
        {
            if (inside_polygon(x,y,p1.m[1][1],p1.m[2][1],p2.m[1][1],p2.m[2][1],p3.m[1][1],p3.m[2][1]))
            {
                //printf("Set: (%d,%d)\n",x,y);
                SetPixelX(d,w,s,x,y);
            }
        }
    }

    SetCurrentColorX(d, &(DefaultGC(d, s)), 0, 0, 0);
}


void Bresenham(Display *d, Window w, int s, int x1, int y1, int x2, int y2)
{ int Transx, Transy ;
    int Pi, Dx, Dy, Two_Dx, Two_Dy, i, Inc1stcoord, Inc2ndcoord, Exchange ;
    
    Exchange = FALSE ;
    Inc1stcoord = 1 ;
    Inc2ndcoord = 1 ;
    
    Transx = -x1 ;
    Transy = -y1 ;
    
    x1 = 0 ;
    y1 = 0 ;
    
    x2 += Transx ;
    y2 += Transy ;
    
    Dx = x2 ;
    Dy = y2 ;
    
    Two_Dx = 2*x2 ;
    Two_Dy = 2*y2 ;
    
    if (Dy < 0) {
        Inc2ndcoord = -1 ;
        Dy *= -1 ;
        Two_Dy *= -1 ;
    }
    
    if (Dx < 0) {
        Inc1stcoord = -1 ;
        Dx *= -1 ;
        Two_Dx *= -1 ;
    }
    
    
    if (Dy > Dx) {
        Exchange = TRUE ;
        exchangeInt(&Two_Dx,&Two_Dy) ;
        exchangeInt(&Dx,&Dy) ;
        exchangeInt(&Inc1stcoord,&Inc2ndcoord) ;
    }
    
    Pi = Two_Dy - Dx ;
    if (Exchange) {
        SetPixelX(d, w, s, y1 - Transx, x1 - Transy) ;
    }
    else {
        SetPixelX(d, w, s, x1 - Transx, y1 - Transy) ;
    }
    for (i = 0 ; i < Dx ; i++) {
        if (Pi < 0) {
            Pi += Two_Dy ;
        }
        else {
            Pi += Two_Dy - Two_Dx ;
            y1 += Inc2ndcoord ;
        }
        x1 += Inc1stcoord ;
        if (Exchange) {
            SetPixelX(d, w, s, y1 - Transx, x1 - Transy) ;
        }
        else {
            SetPixelX(d, w, s, x1 - Transx, y1 - Transy) ;
        }
    }
}

void draw_line_3d(Display* d, Window w, int s, dmatrix_t C, double x1, double y1, double z1, double x2, double y2, double z2)
{
    dmatrix_t p1;
    dmat_alloc(&p1,4,1);

    p1.m[1][1] = x1;
    p1.m[2][1] = y1;
    p1.m[3][1] = z1;
    p1.m[4][1] = 1.0;

    p1 = *perspective_projection(dmat_mult(&C,&p1));

    dmatrix_t p2;
    dmat_alloc(&p2,4,1);

    p2.m[1][1] = x2;
    p2.m[2][1] = y2;
    p2.m[3][1] = z2;
    p2.m[4][1] = 1.0;

    p2 = *perspective_projection(dmat_mult(&C,&p2));

    Bresenham(d,w,s,p1.m[1][1],p1.m[2][1], p2.m[1][1], p2.m[2][1]);
}

void draw_sphere(Display* d, Window w, int s,dmatrix_t C, dmatrix_t E, double radius, double dt, double dp)
{
    for (double t = 0; t < 2 * M_PI;)
        {
            for (double p = 0; p < M_PI;)
            {
                double x1 = radius * cos(t) * sin(p);
                double y1 = radius * sin(t) * sin(p);
                double z1= radius * cos(p);
                
                p += dp;

                double x2 = radius * cos(t) * sin(p);
                double y2 = radius * sin(t) * sin(p);
                double z2= radius * cos(p);

                double x3 = radius * cos(t + dt) * sin(p);
                double y3 = radius * sin(t + dt) * sin(p);
                double z3= radius * cos(p);

                dmatrix_t center;
                dmat_alloc(&center, 4, 1);
                set_xyz(center,0,0,0);

                dmatrix_t mid;
                dmat_alloc(&mid,4, 1);
                set_xyz(mid,(x1 + x3) / 2, (y1 + y3) / 2, (z1 + z3) / 2);

                dmatrix_t outward = get_line(center, mid);

                dmatrix_t normal = find_normal(x1,y1,z1,x2,y2,z2,x3,y3,z3);

                set_color(d,w,s,C,255,0,0,x1,y1,z1,x2,y2,z2,x3,y3,z3,1);

                if (find_angle(outward,normal) < 0)
                    set_xyz(normal, 0 - normal.m[1][1], 0 - normal.m[2][1], 0 - normal.m[3][1]);

                if (!backside(E,normal))
                {
                    draw_line_3d(d,w,s,C,x1,y1,z1,x2,y2,z2);
                    draw_line_3d(d,w,s,C,x2,y2,z2,x3,y3,z3);   
                }
            }

            t += dt;
        }
}

void draw_torus(Display* d, Window w, int s, dmatrix_t C, dmatrix_t E, double radius_a,double radius_c,double t_du,double t_dv)
{
    for (double t_u = 0; t_u < 2 * M_PI; t_u += t_du)
        {
            for (double t_v = 0; t_v < 2 * M_PI;)
            {
                double x1 = (radius_c + radius_a * cos(t_v)) * cos(t_u) + 30;
                double y1 = (radius_c + radius_a * cos(t_v)) * sin(t_u);
                double z1 = radius_a * sin(t_v);
                
                t_v += t_dv;

                double x2 = (radius_c + radius_a * cos(t_v)) * cos(t_u) + 30;
                double y2 = (radius_c + radius_a * cos(t_v)) * sin(t_u);
                double z2 = radius_a * sin(t_v);

                

                double x3 = (radius_c + radius_a * cos(t_v)) * cos(t_u + t_du) + 30;
                double y3 = (radius_c + radius_a * cos(t_v)) * sin(t_u + t_du);
                double z3 = radius_a * sin(t_v);

                /*
                dmatrix_t center;
                dmat_alloc(&center, 4, 1);
                set_xyz(center,30,0,0);

                dmatrix_t mid;
                dmat_alloc(&mid,4, 1);
                set_xyz(mid,(x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3, (z1 + z2 + z3) / 3);

                dmatrix_t outward = get_line(center, mid);


                dmatrix_t normal = find_normal(x1,y1,z1,x2,y2,z2,x3,y3,z3);
    
                if (find_angle(outward,normal) < 0)
                    set_xyz(normal, 0 - normal.m[1][1], 0 - normal.m[2][1], 0 - normal.m[3][1]);
                */

                dmatrix_t normal = find_normal(x3,y3,z3,x2,y2,z2,x1,y1,z1);

                set_color(d,w,s,C,0,255,0,x1,y1,z1,x2,y2,z2,x3,y3,z3,1);
                if (!backside(E,normal))
                {
                    draw_line_3d(d,w,s,C,x1,y1,z1,x2,y2,z2);
                    draw_line_3d(d,w,s,C,x2,y2,z2,x3,y3,z3);   
                }
            }
        }   
}

void draw_cone(Display* d, Window w, int s, dmatrix_t C, dmatrix_t E, double h, double r, double c_du, double c_dt)
{
    for (double u = 0; u < h; u += c_du)
    {
        for (double c_t = 0; c_t < M_PI * 2;)
        {
            double x1 = ((h - u) / h) * r * cos(c_t);
            double y1 = ((h - u) / h) * r * sin(c_t) + 30;
            double z1 = u;

            c_t += c_dt;

            double x2 = ((h - u) / h) * r * cos(c_t);
            double y2 = ((h - u) / h) * r * sin(c_t) + 30;
            double z2 = u;

            

            double x3 = ((h - u - c_du) / h) * r * cos(c_t);
            double y3 = ((h - u - c_du) / h) * r * sin(c_t) + 30;
            double z3 = u + c_du;
            

            dmatrix_t center;
            dmat_alloc(&center, 4, 1);
            set_xyz(center,0,30,0);

            dmatrix_t mid;
            dmat_alloc(&mid,4, 1);
            set_xyz(mid,(x1 + x3) / 2, (y1 + y3) / 2, (z1 + z3) / 2);

            dmatrix_t outward = get_line(center, mid);


            dmatrix_t normal = find_normal(x3,y3,z3,x2,y2,z2,x1,y1,z1);
                
                
            set_color(d,w,s,C,0,0,255,x1,y1,z1,x2,y2,z2,x3,y3,z3,1);

            if (find_angle(outward,normal) < 0)
                set_xyz(normal, 0 - normal.m[1][1], 0 - normal.m[2][1], 0 - normal.m[3][1]);

            //dmatrix_t normal = find_normal(x3,y3,z3,x2,y2,z2,x1,y1,z1);

            if (!backside(E,normal))
            {
                draw_line_3d(d,w,s,C,x1,y1,z1,x2,y2,z2);
                draw_line_3d(d,w,s,C,x2,y2,z2,x3,y3,z3);   
            }
        }
    }
}

int main() 
{
	dmatrix_t E ; /* The centre of projection for the camera */
    
    dmat_alloc(&E,4,1) ;
    
    E.m[1][1] = Ex ;
    E.m[2][1] = Ey ;
    E.m[3][1] = Ez ;
    E.m[4][1] = 1.0 ;
    
    dmatrix_t G ; /* Point gazed at by camera */
    
    dmat_alloc(&G,4,1) ;
    
    G.m[1][1] = Gx ;
    G.m[2][1] = Gy ;
    G.m[3][1] = Gz ;
    G.m[4][1] = 1.0 ;

    dmatrix_t C ; /* The camera matrix */

    dmat_alloc(&C,4,4) ;
    C = *build_camera_matrix(&E,&G) ;

    Display *d ;
    Window w ;
    XEvent e ;
    int s ;
    
    unsigned int r, g, b ;
    
    r = g = b = 0;

    d = InitX(d, &w, &s) ;
    
    SetCurrentColorX(d, &(DefaultGC(d, s)), r, g, b);

    printf("Inside?: %d\n",inside_polygon(0,0,1,0,0,1,-1,0));

    while (1) {
        XNextEvent(d, &e) ;

        draw_line_3d(d,w,s,C,0,0,0,50,0,0);
        draw_line_3d(d,w,s,C,0,0,0,0,50,0);

        double radius = 5.0; //Radius of the sphere
        double dt = 2.0 * M_PI / 40.0; //Change in Longitude
        double dp = M_PI / 20.0; //Change in Latitude

        draw_sphere(d,w,s,C,E,radius,dt,dp);

        double radius_a = 2.0; // Radius of the tube
        double radius_c = 8.0; // Radius from the center of the hole to the center of the tube
        double t_du = M_PI / 20; // Change around the center of the tube
        double t_dv = M_PI / 10; // Change around the center of the hole

        draw_torus(d,w,s,C,E,radius_a,radius_c,t_du,t_dv);

        double height = 12;
        double c_radious = 6;
        double c_dt = M_PI / 20;
        double c_du = 1;

        draw_cone(d,w,s,C,E,height,c_radious,c_du,c_dt);

        if(e.type == KeyPress)
            break ;
        if(e.type == ClientMessage)
            break ;
  }
  QuitX(d,w) ;
}
