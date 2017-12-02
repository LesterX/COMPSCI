/* Computer Science 3388
   Assignment 2
   Yimin Xu
   250876566 */

/* My work starts at line 317 ,
   before that is Brensenham 
   and camera.c from the website*/

#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

#define FALSE 0
#define TRUE  1
#define COLS  512
#define ROWS  512

#define POSX  0
#define POSY  0

#define Ex 15.0
#define Ey 15.0
#define Ez 15.0

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

#define W  512
#define H  512

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
    
    SetCurrentColorX(d, &(DefaultGC(d, s)), r, g, b) ;

    //Here's my work
	
	double radius = 10.0; //Radius of the sphere
    double dt = 2.0 * M_PI / 20.0; //Change in Longitude
    double dp = M_PI / 10.0; //Change in Latitude

    while (1) {
        XNextEvent(d, &e) ;
        
        for (double t = 0; t < 2 * M_PI;)
        {
        	for (double p = 0; p < M_PI;)
        	{
        		dmatrix_t p1;
        		dmat_alloc(&p1,4,1);

        		p1.m[1][1] = radius * cos(t) * sin(p);
        		p1.m[2][1] = radius * sin(t) * sin(p);
        		p1.m[3][1] = radius * cos(p);
        		p1.m[4][1] = 1.0;

        		p1 = *perspective_projection(dmat_mult(&C,&p1));
        		
        		p += dp;

        		dmatrix_t p2;
        		dmat_alloc(&p2,4,1);

        		p2.m[1][1] = radius * cos(t) * sin(p);
        		p2.m[2][1] = radius * sin(t) * sin(p);
        		p2.m[3][1] = radius * cos(p);
        		p2.m[4][1] = 1.0;

        		p2 = *perspective_projection(dmat_mult(&C,&p2));

        		Bresenham(d,w,s,p1.m[1][1],p1.m[2][1], p2.m[1][1], p2.m[2][1]);

        		dmatrix_t p3;
        		dmat_alloc(&p3,4,1);

        		p3.m[1][1] = radius * cos(t + dt) * sin(p);
        		p3.m[2][1] = radius * sin(t + dt) * sin(p);
        		p3.m[3][1] = radius * cos(p);
        		p3.m[4][1] = 1.0;
        		p3 = *perspective_projection(dmat_mult(&C,&p3));

        		Bresenham(d,w,s,p2.m[1][1],p2.m[2][1], p3.m[1][1], p3.m[2][1]);
        	}

        	t += dt;
        }

        double radius_a = 3.0; // Radius of the tube
        double radius_c = 10.0; // Radius from the center of the hole to the center of the tube
        double t_du = M_PI / 20; // Change around the center of the tube
        double t_dv = M_PI / 10; // Change around the center of the hole

        for (double t_u = 0; t_u < 2 * M_PI;)
        {
        	for (double t_v = 0; t_v < 2 * M_PI;)
        	{
        		dmatrix_t p4;
        		dmat_alloc(&p4,4,1);

        		p4.m[1][1] = (radius_c + radius_a * cos(t_v)) * cos(t_u);
        		p4.m[2][1] = (radius_c + radius_a * cos(t_v)) * sin(t_u);
        		p4.m[3][1] = radius_a * sin(t_v);
        		p4.m[4][1] = 1.0;

        		p4 = *perspective_projection(dmat_mult(&C,&p4));
        		
        		t_v += t_dv;

        		dmatrix_t p5;
        		dmat_alloc(&p5,4,1);

        		p5.m[1][1] = (radius_c + radius_a * cos(t_v)) * cos(t_u);
        		p5.m[2][1] = (radius_c + radius_a * cos(t_v)) * sin(t_u);
        		p5.m[3][1] = radius_a * sin(t_v);
        		p5.m[4][1] = 1.0;

        		p5 = *perspective_projection(dmat_mult(&C,&p5));

        		Bresenham(d,w,s,p4.m[1][1],p4.m[2][1], p5.m[1][1], p5.m[2][1]);

        		dmatrix_t p6;
        		dmat_alloc(&p6,4,1);

        		p6.m[1][1] = (radius_c + radius_a * cos(t_v)) * cos(t_u + t_du);
        		p6.m[2][1] = (radius_c + radius_a * cos(t_v)) * sin(t_u + t_du);
        		p6.m[3][1] = radius_a * sin(t_v);
        		p6.m[4][1] = 1.0;

        		p6 = *perspective_projection(dmat_mult(&C,&p6));

        		Bresenham(d,w,s,p5.m[1][1],p5.m[2][1], p6.m[1][1], p6.m[2][1]);
        	}

        	t_u += t_du;
        }

        if(e.type == KeyPress)
            break ;
        if(e.type == ClientMessage)
            break ;
  }
  QuitX(d,w) ;
}
