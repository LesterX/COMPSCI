/* Program which uses X11, implements the synthetic camera, and traces parametric functions.
   Compile with: cc -o "$1" "$1".c -I/opt/X11/include -L/opt/X11/lib -lX11 on Mac OS X or use XCode */

#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

#define FALSE           0
#define TRUE            1
 
#define X               1 
#define Y               2 
#define Z               3

#define Ex              100.0
#define Ey              100.0
#define Ez              100.0

#define Gx              0.0
#define Gy              0.0
#define Gz              30.0

#define Lx              0.0
#define Ly              0.0
#define Lz              1000.0

#define UPx             0.0
#define UPy             0.0
#define UPz             1.0

#define NP              10.0
#define FP              50.0

#define THETA           45.0

#define W               1400
#define H               840
#define ASPECT          (double)W/(double)H
 
#define POSX            0
#define POSY            0

#define SPHERE          0
#define CYLINDER        1
#define CONE            2
#define TORUS           3
#define PLANE           4
#define TOP_DISK        5
#define BOTTOM_DISK     6
#define SHELL           7

#define OBJECTS         64
#define SURFACES        1048576
#define SURFACE_POINTS  4
#define WIREMESH        1
#define SHADING         1

#define RED             1
#define GREEN           2
#define BLUE            3

#define MAX_INTENSITY   255.0 
#define LIGHT_INTENSITY 1.0
#define EXPONENT        10.0 

#define EPSILON         0.0001


typedef struct {
  int wiremesh, shading ;
} render_t ;

typedef struct {
  dmatrix_t P[SURFACE_POINTS], shading_color, wiremesh_color ;
  render_t render ;
  int backface ; 
} surface_t ;


typedef struct {
  int type ; /* Type of parametric opbject */
  double f ; /* specularity exponent */
  dmatrix_t T, p, u, v, c ; /* Placement matrix, surface constants, u, v parameters, color */
  render_t render ; /* Either shading or wiremesh, or both */
} parametric_object_t ;


typedef struct {
  double I ;          /* Light intensity */
  dmatrix_t P, r, c ; /* Position, coefficients, color */
} light_source_t ;


Display *XInit(Display *d, Window *w, int *s) {
    
  d = XOpenDisplay(NULL) ;
  if(d == NULL) {
    printf("Cannot open display\n") ;
    exit(1) ;
  }
  *s = DefaultScreen(d) ;
  *w = XCreateSimpleWindow(d,RootWindow(d,*s),POSX,POSY,W,H,1,BlackPixel(d,*s),WhitePixel(d,*s)) ;
  Atom delWindow = XInternAtom(d,"WM_DELETE_WINDOW",0) ;
  XSetWMProtocols(d,*w,&delWindow,1) ;
  XSelectInput(d,*w,ExposureMask | KeyPressMask) ;
  XMapWindow(d,*w) ;
  return(d) ;
}


void XSetCurrentColor(Display *d, GC *gc, dmatrix_t *C) {
    
  XSetForeground(d,*gc,(unsigned int)((*C).m[RED][1]) << 16 | (unsigned int)((*C).m[GREEN][1]) << 8 | (unsigned int)((*C).m[BLUE][1])) ;
}


void XSetPixel(Display *d, Window w, int s, int i, int j) {
    
  XDrawPoint(d,w,DefaultGC(d,s),i,j) ;
}


void XQuit(Display *d, Window w) {
    
  XDestroyWindow(d,w) ;
  XCloseDisplay(d) ;
}


void exchange_int(int *a, int *b) { 

  int t ;
    
  t = *a ;
  *a = *b ;
  *b = t ;
}


void XTraceLine(Display *d, Window w, int s, int x1, int y1, int x2, int y2)

{ int Transx, Transy, Pi, Dx, Dy, Two_Dx, Two_Dy, i, Inc1stcoord, Inc2ndcoord, Exchange ;
    
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
    exchange_int(&Two_Dx,&Two_Dy) ;
    exchange_int(&Dx,&Dy) ;
    exchange_int(&Inc1stcoord,&Inc2ndcoord) ;
  }
  Pi = Two_Dy - Dx ;
  if (Exchange) {
    XSetPixel(d,w,s,y1 - Transx,x1 - Transy) ;
  }
  else {
    XSetPixel(d,w,s,x1 - Transx,y1 - Transy) ;
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
      XSetPixel(d,w,s,y1 - Transx,x1 - Transy) ;
    }
    else {
      XSetPixel(d,w,s,x1 - Transx,y1 - Transy) ;
    }
  }
}


dmatrix_t *build_view_matrix(dmatrix_t *E, dmatrix_t *G) {

  dmatrix_t *Mv, UP, U, V, N ;

  Mv = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
  dmat_alloc(Mv,4,4) ;

  N = *dmat_normalize(dmat_sub(E,G)) ;

  dmat_alloc(&UP,4,1) ;
  UP.m[X][1] = UPx ;
  UP.m[Y][1] = UPy ;
  UP.m[Z][1] = UPz ;
  UP.m[4][1] = 0.0 ;

  U = *dmat_normalize(dcross_product(from_homogeneous(&UP),from_homogeneous(&N))) ;
  V = *dmat_normalize(dcross_product(from_homogeneous(&N),&U)) ;

  (*Mv).m[1][1] = U.m[X][1] ; 
  (*Mv).m[1][2] = U.m[Y][1] ; 
  (*Mv).m[1][3] = U.m[Z][1] ; 
  (*Mv).m[1][4] = -1.0*ddot_product(E,to_homogeneous(&U,0.0)) ;
    
  (*Mv).m[2][1] = V.m[X][1] ; 
  (*Mv).m[2][2] = V.m[Y][1] ; 
  (*Mv).m[2][3] = V.m[Z][1] ; 
  (*Mv).m[2][4] = -1.0*ddot_product(E,to_homogeneous(&V,0.0)) ;
    
  (*Mv).m[3][1] = N.m[X][1] ; 
  (*Mv).m[3][2] = N.m[Y][1] ; 
  (*Mv).m[3][3] = N.m[Z][1] ; 
  (*Mv).m[3][4] = -1.0*ddot_product(E,&N) ;
    
  (*Mv).m[4][1] = 0.0 ; 
  (*Mv).m[4][2] = 0.0 ; 
  (*Mv).m[4][3] = 0.0 ; 
  (*Mv).m[4][4] = 1.0 ;
  
  return Mv ;
}


dmatrix_t *build_camera_matrix(int h, int w, double np, double fp, double theta, double aspect) {
    
  dmatrix_t Mp, T1, S1, T2, S2, W2, V ; 
    
  dmat_alloc(&Mp,4,4) ;
  Mp = *dmat_identity(&Mp) ;
   
  double a = -1.0*(fp + np)/(fp - np) ;
  double b = -2.0*(fp*np)/(fp - np) ;
    
  Mp.m[1][1] = np ;
  Mp.m[2][2] = np ;
  Mp.m[3][3] = a ;
  Mp.m[3][4] = b ;
  Mp.m[4][3] = -1.0 ;
  Mp.m[4][4] = 0.0 ;
    
  /* Build matrices T_1 and S_1 */
    
  /* Work out coordinates of near plane corners */
    
  double top = np*tan(M_PI/180.0*theta/2.0) ;
  double right = aspect*top ;
  double bottom = -top ;
  double left = -right ;
   
  dmat_alloc(&T1,4,4) ;
  T1 = *dmat_identity(&T1) ;
  T1.m[1][4] = -(right + left)/2.0 ;
  T1.m[2][4] = -(top + bottom)/2.0 ;

  dmat_alloc(&S1,4,4) ;
  S1 = *dmat_identity(&S1) ;
  S1.m[1][1] = 2.0/(right - left) ;
  S1.m[2][2] = 2.0/(top - bottom) ;

  /* Build matrices T2, S2, W2, and V */
    
  dmat_alloc(&T2,4,4) ;
  dmat_alloc(&S2,4,4) ;
  dmat_alloc(&W2,4,4) ;
  dmat_alloc(&V,4,4) ;
    
  T2 = *dmat_identity(&T2) ;
  S2 = *dmat_identity(&S2) ;
  W2 = *dmat_identity(&W2) ;
  V = *dmat_identity(&V) ;
  
  T2.m[1][4] = 1.0 ;
  T2.m[2][4] = 1.0 ;

  S2.m[1][1] = (double)w/2.0 ;
  S2.m[2][2] = (double)h/2.0 ;
    
  W2.m[2][2] = -1.0 ;
  W2.m[2][4] = (double)h ;

  return dmat_mult(&W2,dmat_mult(&S2,dmat_mult(&T2,dmat_mult(&S1,dmat_mult(&T1,&Mp))))) ;
}


parametric_object_t *create_parametric_object(int type, dmatrix_t *p, dmatrix_t *u, dmatrix_t *v, dmatrix_t *c, dmatrix_t *T, dmatrix_t *L, render_t *render, double f) {

  int i ;
  parametric_object_t *object ;

  object = (parametric_object_t *)malloc(sizeof(parametric_object_t)) ;

  (*object).type = type ;
  (*object).T = *dmat_duplicate(T) ;
  (*object).p = *dmat_duplicate(p) ;
  (*object).u = *dmat_duplicate(u) ;
  (*object).v = *dmat_duplicate(v) ;
  (*object).c = *dmat_duplicate(c) ;
  for (i = 1 ; i < 4 ; i++) (*object).c.m[i][1] *= (*L).m[i][1] ; /* Compute color of reflected light, given incident light color and object color */
  (*object).render.wiremesh = (*render).wiremesh ;
  (*object).render.shading = (*render).shading ;
  (*object).f = f ;

  return object ;
}


void XTracePolygon(Display *d, Window w, int s, dmatrix_t P[], int n) {

  int i ;

  for (i = 0 ; i < n ; i++) {
    XTraceLine(d,w,s,(int)P[i].m[X][1],(int)P[i].m[Y][1],(int)P[(i+1)%n].m[X][1],(int)P[(i+1)%n].m[Y][1]) ;
  }
}


dmatrix_t *perspective_projection(dmatrix_t *P) {

  (*P).m[X][1] /= (*P).m[4][1] ;
  (*P).m[Y][1] /= (*P).m[4][1] ;
  (*P).m[Z][1] /= (*P).m[4][1] ;
  (*P).m[4][1] /= (*P).m[4][1] ;
  return P ;
}


dmatrix_t *parametric_point(parametric_object_t *object, double u, double v) {

  dmatrix_t *P ;

  P = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
  dmat_alloc(P,4,1) ;
  
  switch((*object).type) {
         case        SPHERE : { (*P).m[X][1] = (*object).p.m[1][1]*cos(u)*sin(v) ;
                                (*P).m[Y][1] = (*object).p.m[1][1]*sin(u)*sin(v) ;
                                (*P).m[Z][1] = (*object).p.m[1][1]*cos(v) ;
                                (*P).m[4][1] = 1.0 ;
                                break ;
                              }
         case      CYLINDER : { (*P).m[X][1] = (*object).p.m[1][1]*sin(u) ;
                                (*P).m[Y][1] = (*object).p.m[1][1]*cos(u) ;
                                (*P).m[Z][1] = v ;
                                (*P).m[4][1] = 1.0 ; 
                                break ; 
                              }
         case          CONE : { (*P).m[X][1] = (*object).p.m[1][1]*((*object).p.m[2][1] - v)/(*object).p.m[2][1]*sin(u) ;
                                (*P).m[Y][1] = (*object).p.m[1][1]*((*object).p.m[2][1] - v)/(*object).p.m[2][1]*cos(u) ;
                                (*P).m[Z][1] = v ;
                                (*P).m[4][1] = 1.0 ; 
                                break ; 
                              }
         case         TORUS : { (*P).m[X][1] = ((*object).p.m[2][1] + (*object).p.m[1][1]*cos(v))*sin(u) ;
                                (*P).m[Y][1] = ((*object).p.m[2][1] + (*object).p.m[1][1]*cos(v))*cos(u) ;
                                (*P).m[Z][1] = (*object).p.m[1][1]*sin(v) ;
                                (*P).m[4][1] = 1.0 ;
                                break ;
                              }
         case         PLANE : { (*P).m[X][1] = (*object).p.m[1][1] + v ;
                                (*P).m[Y][1] = (*object).p.m[2][1] + u ;
                                (*P).m[Z][1] = (*object).p.m[3][1] ;
                                (*P).m[4][1] = 1.0 ;
                                break ;
                              }
         case      TOP_DISK : { (*P).m[X][1] = u*sin(v) ;
                                (*P).m[Y][1] = u*cos(v) ;
                                (*P).m[Z][1] = (*object).p.m[1][1] ;
                                (*P).m[4][1] = 1.0 ;
                                break ;
                              }
         case   BOTTOM_DISK : { (*P).m[X][1] = u*cos(v) ;
                                (*P).m[Y][1] = u*sin(v) ;
                                (*P).m[Z][1] = (*object).p.m[1][1] ;
                                (*P).m[4][1] = 1.0 ;
                                break ;
                              }
                 case SHELL : { (*P).m[X][1] = ((*object).p.m[1][1]*u + (*object).p.m[3][1]*sin(v))*sin(u) ;
                                (*P).m[Y][1] = ((*object).p.m[1][1]*u + (*object).p.m[3][1]*sin(v))*cos(u) ;
                                (*P).m[Z][1] = (*object).p.m[2][1]*u + (*object).p.m[3][1]*(1.0 - cos(v)) ;
                                (*P).m[4][1] = 1.0 ;
                                break ;
                              }
  }
  return P ;
}


double minimum_coordinate(int coordinate, dmatrix_t P[], int n) {

  int i ;
  double min ;

  min = P[0].m[coordinate][1] ;

  for (i = 1 ; i < n ; i++) {
    if (P[i].m[coordinate][1] < min) {
      min = P[i].m[coordinate][1] ;
    }
  }
  return min ;
}


double maximum_coordinate(int coordinate, dmatrix_t P[], int n) {

  int i ;
  double max ;

  max = P[0].m[coordinate][1] ;

  for (i = 1 ; i < n ; i++) {
    if (P[i].m[coordinate][1] > max) {
      max = P[i].m[coordinate][1] ;
    }
  }
  return max ;
}


int surface_behind_camera(dmatrix_t P[], int n) {

  return minimum_coordinate(Z,P,n) > -NP ;
}


int maximum_intersection(int intersections[], int n) {

  int i, max ;
  
  max = intersections[0] ;
  for (i = 1 ; i < n ; i++) {
    if (intersections[i] > max) {
      max = intersections[i] ;
    }
  }
  return max ;
}


int minimum_intersection(int intersections[], int n) {

  int i, min ;
  
  min = intersections[0] ;
  for (i = 1 ; i < n ; i++) {
    if (intersections[i] < min) {
      min = intersections[i] ;
    }
  }
  return min ;
}


void XFillConvexPolygon(Display *d, Window w, int s, dmatrix_t P[], int n) {

  int i, j ;
  int y, y_min, y_max, min_int, max_int ;
  double m, b ;
  int *active, *horizontal, *intersections ;

  horizontal = (int *)malloc(n*sizeof(int)) ; /* Allocate horizontal segment table */
  active = (int *)malloc(n*sizeof(int)) ; /* Allocate active segment table */
  intersections = (int *)malloc(n*sizeof(int)) ; /* Allocate intersection table */

  y_min = (int)minimum_coordinate(Y,P,n) ; /* Determine number of scan lines */
  y_max = (int)maximum_coordinate(Y,P,n) ;

  for (i = 0 ; i < n ; i++) {
   horizontal[i] = (int)P[i].m[Y][1] == (int)P[(i+1)%n].m[Y][1] ; /* Find horizontal segments */
  }
 
  for (y = y_min ; y <= y_max ; y++) { /* For each scan line y */
    for (i = 0 ; i < n ; i++) {  /* Update segment table */
      if (!horizontal[i]) {
        active[i] = (y >= (int)P[i].m[Y][1] && y <= (int)P[(i+1)%n].m[Y][1]) || (y <= (int)P[i].m[Y][1] && y >= (int)P[(i+1)%n].m[Y][1]) ;
      }
    }
    j = 0 ; 
    for (i = 0 ; i < n ; i++) { /* find intersection x-value. The y-value is given by the scan line */
      if (active[i] && !horizontal[i]) {
        if ((int)P[i].m[X][1] == (int)P[(i+1)%n].m[X][1]) { /* Vertical segment */
          intersections[j++] = (int)P[i].m[X][1] ; 
        }
        else {
          m = (double)((int)P[(i+1)%n].m[Y][1] - (int)P[i].m[Y][1])/(double)((int)P[(i+1)%n].m[X][1] - (int)P[i].m[X][1]) ; /* Compute slope and intercept */
          b = (double)((int)P[i].m[Y][1]) - m*(double)((int)P[i].m[X][1]) ;
          intersections[j++] = (int)(((double)y - b)/m) ; /* Compute intersection */
        }  
      }
    }
    min_int = minimum_intersection(intersections,j) ;
    max_int = maximum_intersection(intersections,j) + 1 ;
    for ( i = min_int ; i < max_int ; i++) { /* Tracing from minimum to maximum intersection */
      XSetPixel(d,w,s,i,y) ;
    }
  }
  free(horizontal) ;
  free(active) ;
  free(intersections) ;
}


dmatrix_t *surface_centroid(dmatrix_t P[], int n) {

  int i ;
  double a ;
  dmatrix_t *A ;

  A = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
  dmat_alloc(A,4,1) ;

  A = dmat_init(A,0.0) ;
  for (i = 0 ; i < n ; i++) {
    A = dmat_add(A,&P[i]) ;
  }
  
  a = 1.0/(double)n ;
  A = dmat_scalar_mult(A,a) ;

  return A ;
}


dmatrix_t *normal_to_surface(dmatrix_t P[]) {

  dmatrix_t U, V ;

  U = *from_homogeneous(dmat_normalize(dmat_sub(&P[2],&P[0]))) ;
  V = *from_homogeneous(dmat_normalize(dmat_sub(&P[3],&P[1]))) ; 
  return to_homogeneous(dmat_normalize(dcross_product(&U,&V)),0.0) ; 
}


double color_index(light_source_t *light, dmatrix_t *N, dmatrix_t *S, dmatrix_t *R, dmatrix_t *V, double f) {

  double index ;
  double diffuse, specular ;

  if ((diffuse = ddot_product(N,S)) < 0.0) diffuse = 0.0 ; 
  if ((specular = ddot_product(R,V)) < 0.0) specular = 0.0 ;
  index = MAX_INTENSITY*(*light).I*((*light).r.m[1][1] + (*light).r.m[2][1]*diffuse + (*light).r.m[3][1]*pow(specular,f)) ; 

  return index ;
}


void quicksort(surface_t surface[], int index[], int length) {

  int i, j ;
  double pivot ;

  if (length < 2) return ;
  pivot = surface[index[length/2]].P[0].m[3][1] ;

  for (i = 0, j = length - 1 ; ; i++, j--) {
    while (surface[index[i]].P[0].m[3][1] > pivot) i++ ; /* Sorting on pseudo-depth values */
    while (surface[index[j]].P[0].m[3][1] < pivot) j-- ;
 
    if (i >= j) break ;
    exchange_int(&index[i],&index[j]) ; 
  }
  quicksort(surface,index,i) ;
  quicksort(surface,index+i,length-i) ; 
}


void object_to_surfaces(light_source_t *light, dmatrix_t *C, dmatrix_t *Mv, dmatrix_t *L, parametric_object_t *object, surface_t surface[], int *n_surfaces) {

  int i ;
  double u, v ;
  dmatrix_t N, S, R, V ;

  for (u = (*object).u.m[1][1] ; u + (*object).u.m[3][1] <= (*object).u.m[2][1] + EPSILON ; u += (*object).u.m[3][1]) {
    for (v = (*object).v.m[1][1] ; v + (*object).v.m[3][1] <= (*object).v.m[2][1] + EPSILON ; v += (*object).v.m[3][1]) {

      surface[*n_surfaces].P[0] = *dmat_mult(Mv,dmat_mult(&(*object).T,parametric_point(object,u,v))) ; /* Compute points in camera coordinates */
      surface[*n_surfaces].P[1] = *dmat_mult(Mv,dmat_mult(&(*object).T,parametric_point(object,u,v+(*object).v.m[3][1]))) ;
      surface[*n_surfaces].P[2] = *dmat_mult(Mv,dmat_mult(&(*object).T,parametric_point(object,u+(*object).u.m[3][1],v+(*object).v.m[3][1]))) ;
      surface[*n_surfaces].P[3] = *dmat_mult(Mv,dmat_mult(&(*object).T,parametric_point(object,u+(*object).u.m[3][1],v))) ;

      if (!surface_behind_camera(surface[*n_surfaces].P,SURFACE_POINTS)) { /* Do not render any surface behind the near plane */

        N = *normal_to_surface(surface[*n_surfaces].P) ; /* Compute normal to surface in camera coordinates */
        S = *to_homogeneous(dmat_normalize(dmat_sub(from_homogeneous(dmat_mult(Mv,L)),from_homogeneous(surface_centroid(surface[*n_surfaces].P,SURFACE_POINTS)))),0.0) ; 
        R = *dmat_add(dmat_scalar_mult(&S,-1.0),dmat_scalar_mult(&N,2.0*ddot_product(&S,&N))) ;
        V = *to_homogeneous(dmat_normalize(from_homogeneous(dmat_scalar_mult(surface_centroid(surface[*n_surfaces].P,SURFACE_POINTS),-1.0))),0.0) ;

        surface[*n_surfaces].backface = ddot_product(surface_centroid(surface[*n_surfaces].P,SURFACE_POINTS),&N) > 0.0 ;
        surface[*n_surfaces].shading_color = *dmat_scalar_mult(&(*object).c,color_index(light,&N,&S,&R,&V,(*object).f)) ; /* Apply color model to color channels */
        surface[*n_surfaces].wiremesh_color = *dmat_scalar_mult(&(*object).c,MAX_INTENSITY-color_index(light,&N,&S,&R,&V,(*object).f)) ; /* Complementary color to object for wiremesh */
        surface[*n_surfaces].render.shading = (*object).render.shading ;
        surface[*n_surfaces].render.wiremesh = (*object).render.wiremesh ;

        for (i = 0 ; i < SURFACE_POINTS ; i++) {
          surface[*n_surfaces].P[i] = *perspective_projection(dmat_mult(C,&surface[*n_surfaces].P[i])) ; /* Compute image coordinates of surface */
        }
        (*n_surfaces)++ ;
      } 
    }
  }
}


void trace_surfaces(Display *d, Window w, int s, surface_t surface[], int index[], int n_surfaces) {

  int i ;
 
  for (i = 0 ; i < n_surfaces ; i++) {
    if (surface[index[i]].render.shading) {
      XSetCurrentColor(d,&(DefaultGC(d,s)),&surface[index[i]].shading_color) ;
      XFillConvexPolygon(d,w,s,surface[index[i]].P,SURFACE_POINTS) ; /* Fill the polygon */
    }
    if (surface[index[i]].render.wiremesh) {
      XSetCurrentColor(d,&(DefaultGC(d,s)),&surface[index[i]].wiremesh_color) ;
      XTracePolygon(d,w,s,surface[index[i]].P,SURFACE_POINTS) ; /* Trace the contour of the polygon */
    } 
  }
}


light_source_t *build_light_source(dmatrix_t *P, dmatrix_t *r, dmatrix_t *c, double intensity) {

  light_source_t *light ;

  light = (light_source_t *)malloc(sizeof(light_source_t)) ;

  (*light).I = intensity ;
  (*light).P = *dmat_duplicate(P) ;
  (*light).r = *dmat_duplicate(r) ;
  (*light).c = *dmat_duplicate(c) ;

  return light ;
} 


int main() {

  Display *d ;
  Window w ;
  XEvent e ;

  int i, s, n_objects, n_surfaces ;
  double light_intensity, f ;
  dmatrix_t E, G, Mv, C, L, T, r, lc, p, u, v, c ; 
  render_t render ;
  light_source_t light ;
  parametric_object_t object[OBJECTS] ;

  static surface_t surface[SURFACES] ;
  static int index[SURFACES] ;

  dmat_alloc(&E,4,1) ; /* The centre of projection for the camera */
  E.m[X][1] = Ex ;
  E.m[Y][1] = Ey ;
  E.m[Z][1] = Ez ;
  E.m[4][1] = 1.0 ;

  dmat_alloc(&G,4,1) ; /* Point gazed at by camera */
  G.m[X][1] = Gx ;
  G.m[Y][1] = Gy ;
  G.m[Z][1] = Gz ;
  G.m[4][1] = 1.0 ;

  dmat_alloc(&L,4,1) ; /* Position of light source */
  L.m[X][1] = Lx ;
  L.m[Y][1] = Ly ;
  L.m[Z][1] = Lz ;
  L.m[4][1] = 0.0 ;

  dmat_alloc(&r,4,1) ; /* Ambient, diffuse, and specular coefficients */
  r.m[1][1] = 0.2 ;
  r.m[2][1] = 0.4 ;
  r.m[3][1] = 0.4 ;
  r.m[4][1] = 0.0 ;

  dmat_alloc(&lc,4,1) ; /* Color of incident light */
  lc.m[1][1] = 1.0 ;
  lc.m[2][1] = 1.0 ;
  lc.m[3][1] = 1.0 ;
  lc.m[4][1] = 0.0 ;

  light_intensity = LIGHT_INTENSITY ;
  f = EXPONENT ;

  Mv = *build_view_matrix(&E,&G) ; /* The camera matrix */
  C = *build_camera_matrix((double)H,(double)W,NP,FP,THETA,ASPECT) ; /* The remaining part of the matrix pipeline */
  light = *build_light_source(&L,&r,&lc,light_intensity) ; /* Scene light with position, coefficients, color, and intensity */

  dmat_alloc(&T,4,4) ;
  dmat_alloc(&p,6,1) ;
  dmat_alloc(&u,3,1) ;
  dmat_alloc(&v,3,1) ;
  dmat_alloc(&c,3,1) ;

  n_objects = 0 ;
  n_surfaces = 0 ;

  /* Create plane */

  T = *dmat_identity(&T) ;
  
  u.m[1][1] = -100.0 ;
  u.m[2][1] = 100.0 ;
  u.m[3][1] = u.m[2][1]/200.0 ;

  v.m[1][1] = -100.0 ;
  v.m[2][1] = 100.0 ;
  v.m[3][1] = v.m[2][1]/200.0 ;

  p.m[1][1] = 0.0 ; /* Point in the plane */
  p.m[2][1] = 0.0 ; 
  p.m[3][1] = 0.0 ; 

  c.m[RED][1] = 0.1 ;
  c.m[GREEN][1] = 0.5 ;
  c.m[BLUE][1] = 0.5 ;

  render.wiremesh = !TRUE ;
  render.shading = TRUE ;

  object[n_objects++] = *create_parametric_object(PLANE,&p,&u,&v,&c,&T,&light.c,&render,f) ;

  /* Create torus */

  T = *dmat_identity(&T) ;
  T.m[1][4] = 10.0 ;
  T.m[2][4] = 10.0 ;
  T.m[3][4] = 20.0 ;

  u.m[1][1] = 0.0 ;
  u.m[2][1] = 2.0*M_PI ;
  u.m[3][1] = u.m[2][1]/360.0 ;

  v.m[1][1] = 0.0 ;
  v.m[2][1] = 2.0*M_PI ;
  v.m[3][1] = v.m[2][1]/360.0 ;

  p.m[1][1] = 3.0 ; /* Torus inner radius */
  p.m[2][1] = 9.0 ; /* Torus outer radius */

  c.m[RED][1] = 1.0 ;
  c.m[GREEN][1] = 0.0 ;
  c.m[BLUE][1] = 0.0 ;

  render.wiremesh = !TRUE ;
  render.shading = TRUE ;

  object[n_objects++] = *create_parametric_object(TORUS,&p,&u,&v,&c,&T,&light.c,&render,f) ; 

  /* Create sphere */

  T = *dmat_identity(&T) ;
  T.m[1][4] = -20.0 ; 
  T.m[3][4] = 12.0 ; 

  u.m[1][1] = 0.0 ;
  u.m[2][1] = 2.0*M_PI ;
  u.m[3][1] = u.m[2][1]/360.0 ;

  v.m[1][1] = 0.0 ;
  v.m[2][1] = M_PI ;
  v.m[3][1] = v.m[2][1]/180.0 ;

  p.m[1][1] = 10.0 ; /* Sphere radius */

  c.m[RED][1] = 1.0 ;
  c.m[GREEN][1] = 0.6 ;
  c.m[BLUE][1] = 0.1 ;

  render.wiremesh = !TRUE ;
  render.shading = TRUE ;

  object[n_objects++] = *create_parametric_object(SPHERE,&p,&u,&v,&c,&T,&light.c,&render,f) ;

  /* Create sphere */

  T = *dmat_identity(&T) ;
  T.m[1][4] = 30.0 ; 
  T.m[2][4] = 30.0 ; 
  T.m[3][4] = 18.0 ; 

  u.m[1][1] = 0.0 ;
  u.m[2][1] = 2.0*M_PI ;
  u.m[3][1] = u.m[2][1]/360.0 ;

  v.m[1][1] = 0.0 ;
  v.m[2][1] = M_PI ;
  v.m[3][1] = v.m[2][1]/180.0 ;

  p.m[1][1] = 10.0 ; /* Sphere radius */

  c.m[RED][1] = 0.6 ;
  c.m[GREEN][1] = 0.6 ;
  c.m[BLUE][1] = 0.1 ;

  render.wiremesh = !TRUE ;
  render.shading = TRUE ;

  object[n_objects++] = *create_parametric_object(SPHERE,&p,&u,&v,&c,&T,&light.c,&render,f) ;

  /* Create cone */

  T = *dmat_identity(&T) ;
  T.m[2][4] = -20.0 ; 
  T.m[3][4] = 12.0 ; 
  T.m[1][1] = cos(1.3) ; T.m[1][3] = sin(1.3) ; /* Rotate cone around the y-axis */
  T.m[3][1] = -sin(1.3) ; T.m[3][3] = cos(1.3) ; 

  u.m[1][1] = 0.0 ;
  u.m[2][1] = 2.0*M_PI ; 
  u.m[3][1] = u.m[2][1]/360.0 ; 

  v.m[1][1] = 0.0 ;
  v.m[2][1] = 20.0 ;
  v.m[3][1] = v.m[2][1]/100.0 ;

  p.m[1][1] = 10.0 ; /* Cone radius */
  p.m[2][1] = 20.0 ; /* Cone height */

  c.m[RED][1] = 0.0 ;
  c.m[GREEN][1] = 0.0 ;
  c.m[BLUE][1] = 1.0 ; 

  render.wiremesh = !TRUE ;
  render.shading = TRUE ;

  object[n_objects++] = *create_parametric_object(CONE,&p,&u,&v,&c,&T,&light.c,&render,f) ;

  /* Cap the cone base with a circle */ 

  T = *dmat_identity(&T) ;
  T.m[2][4] = -20.0 ;
  T.m[3][4] = 12.0 ;
  T.m[1][1] = cos(1.3) ; T.m[1][3] = sin(1.3) ; /* Rotate circle around the y-axis  */
  T.m[3][1] = -sin(1.3) ; T.m[3][3] = cos(1.3) ;

  u.m[1][1] = 0.0 ;
  u.m[2][1] = 10.0 ; /* Circle radius */
  u.m[3][1] = u.m[2][1]/100.0 ; 

  v.m[1][1] = 0.0 ;
  v.m[2][1] = 2.0*M_PI ;
  v.m[3][1] = v.m[2][1]/360.0 ;

  p.m[1][1] = 0.0 ;

  c.m[RED][1] = 0.0 ;
  c.m[GREEN][1] = 0.0 ;
  c.m[BLUE][1] = 1.0 ;

  render.wiremesh = !TRUE ;
  render.shading = TRUE ;

  object[n_objects++] = *create_parametric_object(BOTTOM_DISK,&p,&u,&v,&c,&T,&light.c,&render,f) ;

  /* Create cone */

  T = *dmat_identity(&T) ;
  T.m[1][4] = -30.0 ; 
  T.m[2][4] = -30.0 ; 
  T.m[3][4] = 30.0 ; 
  T.m[1][1] = cos(-0.8) ; T.m[1][3] = sin(-0.8) ; /* Rotate cone around the y-axis */
  T.m[3][1] = -sin(-0.8) ; T.m[3][3] = cos(-0.8) ; 

  u.m[1][1] = 0.0 ;
  u.m[2][1] = 2.0*M_PI ; 
  u.m[3][1] = u.m[2][1]/360.0 ; 

  v.m[1][1] = 0.0 ;
  v.m[2][1] = 20.0 ;
  v.m[3][1] = v.m[2][1]/100.0 ;

  p.m[1][1] = 10.0 ; 
  p.m[2][1] = 20.0 ; 

  c.m[RED][1] = 0.0 ;
  c.m[GREEN][1] = 1.0 ;
  c.m[BLUE][1] = 0.0 ;

  render.wiremesh = !TRUE ;
  render.shading = TRUE ;

  object[n_objects++] = *create_parametric_object(CONE,&p,&u,&v,&c,&T,&light.c,&render,f) ;

  /* Create cylinder */

  T = *dmat_identity(&T) ;
  T.m[1][4] = -60.0 ; 
  T.m[2][4] = -60.0 ; 
  T.m[3][4] = 30.0 ;
  T.m[1][1] = cos(0.9) ; T.m[1][3] = sin(0.9) ; /* Rotate cylinder around the y-axis  */
  T.m[3][1] = -sin(0.9) ; T.m[3][3] = cos(0.9) ;

  u.m[1][1] = 0.0 ;
  u.m[2][1] = 2.0*M_PI ; 
  u.m[3][1] = u.m[2][1]/360.0 ; 

  v.m[1][1] = 0.0 ;
  v.m[2][1] = 20.0 ; /* Cylinder height */
  v.m[3][1] = v.m[2][1]/100.0 ;

  p.m[1][1] = 10.0 ; /* Cylinder radius */

  c.m[RED][1] = 1.0 ;
  c.m[GREEN][1] = 1.0 ;
  c.m[BLUE][1] = 0.0 ;

  render.wiremesh = !TRUE ;
  render.shading = TRUE ;

  object[n_objects++] = *create_parametric_object(CYLINDER,&p,&u,&v,&c,&T,&light.c,&render,f) ;

  /* Create circle */

  T = *dmat_identity(&T) ;
  T.m[1][4] = -60.0 ; 
  T.m[2][4] = -60.0 ; 
  T.m[3][4] = 30.0 ;
  T.m[1][1] = cos(0.9) ; T.m[1][3] = sin(0.9) ; /* Rotate circle around the y-axis  */
  T.m[3][1] = -sin(0.9) ; T.m[3][3] = cos(0.9) ;

  u.m[1][1] = 0.0 ;
  u.m[2][1] = 10.0 ; /* Circle radius */
  u.m[3][1] = u.m[2][1]/100.0 ; 

  v.m[1][1] = 0.0 ;
  v.m[2][1] = 2.0*M_PI ;
  v.m[3][1] = v.m[2][1]/360.0 ;

  p.m[1][1] = 20.0 ;

  c.m[RED][1] = 1.0 ;
  c.m[GREEN][1] = 1.0 ;
  c.m[BLUE][1] = 0.0 ;

  render.wiremesh = !TRUE ;
  render.shading = TRUE ;

  object[n_objects++] = *create_parametric_object(TOP_DISK,&p,&u,&v,&c,&T,&light.c,&render,f) ;

  /* Create circle */

  T = *dmat_identity(&T) ;
  T.m[1][4] = -30.0 ; 
  T.m[2][4] = -30.0 ; 
  T.m[3][4] = 30.0 ;
  T.m[1][1] = cos(-0.8) ; T.m[1][3] = sin(-0.8) ; /* Rotate circle around the y-axis */
  T.m[3][1] = -sin(-0.8) ; T.m[3][3] = cos(-0.8) ; 

  u.m[1][1] = 0.0 ;
  u.m[2][1] = 10.0 ; /* Circle radius */
  u.m[3][1] = u.m[2][1]/100.0 ;

  v.m[1][1] = 0.0 ;
  v.m[2][1] = 2.0*M_PI ;
  v.m[3][1] = v.m[2][1]/360.0 ;

  p.m[1][1] = 0.0 ;

  c.m[RED][1] = 0.0 ;
  c.m[GREEN][1] = 1.0 ;
  c.m[BLUE][1] = 0.0 ;

  render.wiremesh = !TRUE ;
  render.shading = TRUE ;

  object[n_objects++] = *create_parametric_object(BOTTOM_DISK,&p,&u,&v,&c,&T,&light.c,&render,f) ;

/* Create shell */

  T = *dmat_identity(&T) ;
  T.m[1][4] = -70 ;

  u.m[1][1] = 0.0 ;
  u.m[2][1] = 4.0*M_PI ;
  u.m[3][1] = u.m[2][1]/36.0 ;

  v.m[1][1] = 0.0 ;
  v.m[2][1] = 4.0*M_PI ;
  v.m[3][1] = v.m[2][1]/36.0 ;

  p.m[1][1] = 2.3 ;
  p.m[2][1] = 4.3 ;
  p.m[3][1] = 15.0 ;

  c.m[RED][1] = 1.0 ;
  c.m[GREEN][1] = 0.0 ;
  c.m[BLUE][1] = 0.5 ;

  render.wiremesh = TRUE ;
  render.shading = TRUE ;

  object[n_objects++] = *create_parametric_object(SHELL,&p,&u,&v,&c,&T,&light.c,&render,f) ;

  for (i = 0 ; i < n_objects ; i++) {
    object_to_surfaces(&light,&C,&Mv,&L,&object[i],surface,&n_surfaces) ;
  }

  for (i = 0 ; i < n_surfaces ; i++) { /* Initialize index array */
    index[i] = i ;         
  }
  quicksort(surface,index,n_surfaces) ;

  d = XInit(d,&w,&s) ;
  XSetWindowBackground(d,w,0) ;
  XNextEvent(d,&e) ;

  while (1) {
    XNextEvent(d,&e) ;
    if (e.type == Expose) { 
      trace_surfaces(d,w,s,surface,index,n_surfaces) ;
    }
    if(e.type == KeyPress)
      break ;
    if(e.type == ClientMessage)
      break ;
  }

  XQuit(d,w) ;
  free_dmatrix(E.m,1,E.l,1,E.c) ;
  free_dmatrix(G.m,1,G.l,1,G.c) ;
  free_dmatrix(C.m,1,C.l,1,C.c) ;
  free_dmatrix(L.m,1,L.l,1,L.c) ;
  free_dmatrix(T.m,1,T.l,1,T.c) ;
  free_dmatrix(p.m,1,p.l,1,p.c) ;
  free_dmatrix(u.m,1,u.l,1,u.c) ;
  free_dmatrix(v.m,1,v.l,1,v.c) ;
  free_dmatrix(c.m,1,c.l,1,c.c) ;
  free_dmatrix(r.m,1,r.l,1,r.c) ;
  free_dmatrix(lc.m,1,lc.l,1,lc.c) ;
  free_dmatrix(Mv.m,1,Mv.l,1,Mv.c) ;
}