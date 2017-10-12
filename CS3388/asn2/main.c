/* Program which uses X11 and implements Bresenham's integer line drawing algorithm.
   Compile with: cc -o "$1" "$1".c -I/opt/X11/include -L/opt/X11/lib -lX11 on Mac OS X */

#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FALSE 0
#define TRUE  1
#define COLS  512
#define ROWS  512

#define POSX  0
#define POSY  0

// #define NUM 1024

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
int main() {

    Display *d ;
    Window w ;
    XEvent e ;
    int s ;
    
    unsigned int r, g, b ;
    int x1, y1, x2, y2 ;
    double t, dt ;
    
    r = g = b = 0 ;
    dt = 2.0*M_PI/200.0 ;

    d = InitX(d, &w, &s) ;
    
    SetCurrentColorX(d, &(DefaultGC(d, s)), r, g, b) ;
    
    while (1) {
        XNextEvent(d, &e) ;
        if (e.type == Expose)
            for (t = 0.0 ; t < 2.0*M_PI;) {
                x1 = 256 + (int)100.0*(1.5*cos(t) - cos(13.0*t)) ;
                y1 = 256 + (int)100.0*(1.5*sin(t) - sin(13.0*t)) ;
                t += dt ;
                x2 = 256 + (int)100.0*(1.5*cos(t) - cos(13.0*t)) ;
                y2 = 256 + (int)100.0*(1.5*sin(t) - sin(13.0*t)) ;
                Bresenham(d,w,s,x1,y1,x2,y2) ;
            }
        if(e.type == KeyPress)
            break ;
        if(e.type == ClientMessage)
            break ;
  }
  QuitX(d,w) ;
}
