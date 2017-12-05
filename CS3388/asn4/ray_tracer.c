/*            PURPOSE : Simple framework for ray-tracing
 
        PREREQUISITES : matrix.h
 */

#include <X11/Xlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

#define INFINITE_PLANE   0
#define PLANE            1
#define SPHERE           2

#define EPSILON          0.00001
#define INFTY            999999.9
#define N_OBJECTS        1024
#define MAX_INTENSITY    255.0

#define Ex               0.0
#define Ey               -10.0
#define Ez               6.0

#define Gx               0.0
#define Gy               0.0
#define Gz               0.0

#define UPx              0.0
#define UPy              0.0
#define UPz              1.0

#define Lx               0.0
#define Ly               0.0
#define Lz               8.0

#define Near             1.0
#define Far              25.0

#define THETA            45.0
#define ASPECT           1.5

#define H                300

typedef struct {
    int width, height ;
} window_t ;

typedef struct {
    dmatrix_t UP ;
    dmatrix_t E ;
    dmatrix_t G ;
    dmatrix_t u, v, n ;
} camera_t ;

typedef struct {
    double r, g, b ;
} color_t ;

typedef struct {
    int type ;
    double (*intersection_function)(dmatrix_t *,dmatrix_t *) ;
    dmatrix_t M, Minv ;
    color_t specular_color, diffuse_color, ambient_color ;
    double reflectivity, specular_coeff, diffuse_coeff, ambient_coeff, f ;
} object_t ;

typedef struct {
    dmatrix_t position ;
    color_t color ;
    color_t intensity ;
} light_t ;

object_t object[N_OBJECTS] ;
int nobjects = 0 ;

Display *InitX(Display *d, Window *w, int *s, window_t *Window) {
    
    d = XOpenDisplay(NULL) ;
    if(d == NULL) {
        printf("Cannot open display\n") ;
        exit(1) ;
    }
    *s = DefaultScreen(d) ;
    *w = XCreateSimpleWindow(d,RootWindow(d,*s),0,0,Window->width,Window->height,1,BlackPixel(d,*s),WhitePixel(d, *s)) ;
    Atom delWindow = XInternAtom(d,"WM_DELETE_WINDOW",0) ;
    XSetWMProtocols(d,*w,&delWindow,1) ;
    XSelectInput(d,*w,ExposureMask | KeyPressMask) ;
    XMapWindow(d,*w) ;
    return(d) ;
}

void SetCurrentColorX(Display *d, GC *gc, unsigned int r, unsigned int g, unsigned int b) {
    
    XSetForeground(d,*gc,r << 16 | g << 8 | b) ;
}

void SetPixelX(Display *d, Window w, int s, int i, int j) {
    
    XDrawPoint(d,w,DefaultGC(d,s),i,j) ;
}

void QuitX(Display *d, Window w) {
    
    XDestroyWindow(d,w) ;
    XCloseDisplay(d) ;
}

light_t *build_light(light_t *light, dmatrix_t *position, color_t color, color_t intensity) {
    
    dmat_alloc(&light->position,4,1) ;
    
    light->position = *position ;
    light->color.r = color.r ;
    light->color.g = color.g ;
    light->color.b = color.b ;
    light->intensity.r = intensity.r ;
    light->intensity.g = intensity.g ;
    light->intensity.b = intensity.b ;
    return light ;
}

window_t *build_window(window_t *Window, int height, float aspect) {
    
    Window->height = height ;
    Window->width =  aspect*height ;
    
    return(Window) ;
}

camera_t *build_camera(camera_t *Camera, window_t *Window) {
    
    dmat_alloc(&Camera->E,4,1) ;

    Camera->E.m[1][1] = Ex ;
    Camera->E.m[2][1] = Ey ;
    Camera->E.m[3][1] = Ez ;
    Camera->E.m[4][1] = 1.0 ;
    
    dmat_alloc(&Camera->G,4,1) ;
    
    Camera->G.m[1][1] = Gx ;
    Camera->G.m[2][1] = Gy ;
    Camera->G.m[3][1] = Gz ;
    Camera->G.m[4][1] = 1.0 ;

    dmat_alloc(&Camera->n,4,1) ;
    Camera->n = *dmat_normalize(dmat_sub(&Camera->E,&Camera->G)) ;
    Camera->n.l = 3 ;
    
    dmat_alloc(&Camera->UP,4,1) ;
    
    Camera->UP.l = 3 ;
    
    Camera->UP.m[1][1] = UPx ;
    Camera->UP.m[2][1] = UPy ;
    Camera->UP.m[3][1] = UPz ;
    Camera->UP.m[4][1] = 1.0 ;
    
    dmat_alloc(&Camera->u,4,1) ;
    
    Camera->u = *dmat_normalize(dcross_product(&Camera->UP,&Camera->n)) ;
    Camera->v = *dmat_normalize(dcross_product(&Camera->n,&Camera->u)) ;
    
    return(Camera) ;
}

dmatrix_t *intersection_coordinates(dmatrix_t *e, dmatrix_t *direction, double t) {
    
    dmatrix_t *intersection ;
    
    intersection = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(intersection,4,1) ;
    
    intersection->m[1][1] = e->m[1][1] + direction->m[1][1]*t ;
    intersection->m[2][1] = e->m[2][1] + direction->m[2][1]*t ;
    intersection->m[3][1] = e->m[3][1] + direction->m[3][1]*t ;
    intersection->m[4][1] = 1.0 ;
    
    return intersection ;
}

double infinite_plane_intersection(dmatrix_t *e, dmatrix_t *d) {
    
    if (d->m[3][1] >= 0.0) return -1.0 ; else return -1.0*e->m[3][1]/d->m[3][1] ;
}

double plane_intersection(dmatrix_t *e, dmatrix_t *d) {
    
    double t ;
    dmatrix_t *intersection ;
    
    if (d->m[3][1] >= 0.0) {
        return -1.0 ;
    }
    else {
        t = -1.0*e->m[3][1]/d->m[3][1] ;
        intersection = intersection_coordinates(e,d,t) ;
        if ((fabs(intersection->m[1][1]) < 3.0) && (fabs(intersection->m[2][1]) < 3.0)) {
            delete_dmatrix(intersection) ;
            return t ;
        }
        else {
            delete_dmatrix(intersection) ;
            return -1.0 ;
        }
    }
}

double sphere_intersection(dmatrix_t *e, dmatrix_t *d) {
    
    double a = ddot_product(from_homogeneous(d),from_homogeneous(d)) ;
    double b = ddot_product(from_homogeneous(e),from_homogeneous(d)) ;
    double c = ddot_product(from_homogeneous(e),from_homogeneous(e)) - 1.0 ;
    
    double discriminant = b*b - a*c ;

    if (discriminant < 0.0) {
        return -1.0 ;
    }
    else {
        if (discriminant < EPSILON) {
            return -b/a ;
        }
        else {
            double t1 = -b/a - sqrtl(discriminant)/a ;
            double t2 = -b/a + sqrtl(discriminant)/a ;
            if (t1 < t2) {
                if (t1 > EPSILON) {
                    return t1 ;
                }
                else {
                    return -1.0 ;
                }
            }
            else {
                return t2 ;
            }
        }
    }
}


dmatrix_t *normal_to_surface(object_t *object, dmatrix_t *intersection) {
    
    dmatrix_t *normal ;
    
    normal = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(normal,4,1) ;
    
    switch ((*object).type) {
            
        case PLANE          :   normal->m[1][1] = 0.0 ;
                                normal->m[2][1] = 0.0 ;
                                normal->m[3][1] = 1.0 ;
                                break ;
            
        case INFINITE_PLANE :   normal->m[1][1] = 0.0 ;
                                normal->m[2][1] = 0.0 ;
                                normal->m[3][1] = 1.0 ;
                                normal->m[4][1] = 0.0 ;
                                break ;
            
        case SPHERE         :   normal->m[1][1] = intersection->m[1][1] ;
                                normal->m[2][1] = intersection->m[2][1] ;
                                normal->m[3][1] = intersection->m[3][1] ;
                                normal->m[4][1] = 0.0 ;
                                break ;
            
        default: printf("No such object type\n") ;
            
    }
    return normal ;
}

int find_min_hit_time(double t0[N_OBJECTS]) {
    
    double min_t = INFTY ;
    int position = -1 ;
    
    for(int i = 0 ; i < nobjects ; i++) {
        if (t0[i] != -1.0) {
            if (t0[i] < min_t) {
                min_t = t0[i] ;
                position = i ;
            }
        }
    }
    return position ;
}

dmatrix_t *ray_direction(camera_t *Camera, window_t *Window, double height, double width, double i, double j) {
    
    dmatrix_t *d ;
    
    d = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(d,4,1) ;

    d->m[1][1] = -1.0*Near*Camera->n.m[1][1] +
    width*(2.0*i/Window->width - 1.0)*Camera->u.m[1][1] +
    height*(2.0*j/Window->height - 1.0)*Camera->v.m[1][1] ;

    d->m[2][1] = -1.0*Near*Camera->n.m[2][1] +
    width*(2.0*i/Window->width - 1.0)*Camera->u.m[2][1] +
    height*(2.0*j/Window->height - 1.0)*Camera->v.m[2][1] ;

    d->m[3][1] = -1.0*Near*Camera->n.m[3][1] +
    width*(2.0*i/Window->width - 1.0)*Camera->u.m[3][1] +
    height*(2.0*j/Window->height - 1.0)*Camera->v.m[3][1] ;

    d->m[4][1] = 0.0 ;
    
    return(d) ;
}

dmatrix_t *vector_to_light_source(dmatrix_t *intersection, dmatrix_t *light_position) {
    
    dmatrix_t *s, *sn ;
    
    s = dmat_sub(light_position,intersection) ;
    sn = dmat_normalize(s) ;
    delete_dmatrix(s) ;

    return sn ;
}

dmatrix_t *vector_to_center_of_projection(dmatrix_t *intersection, dmatrix_t *e) {
    
    dmatrix_t *v, *vn ;

    v = dmat_sub(e,intersection) ;
    vn = dmat_normalize(v) ;
    delete_dmatrix(v) ;
    
    return vn ;
}

dmatrix_t *vector_to_specular_reflection(dmatrix_t *N, dmatrix_t *S) {
    
    dmatrix_t *r, *rn ;
    
    r = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(r,4,1) ;
    
    double sn = 2.0*ddot_product(N,S) ;
    
    r->m[1][1] = -1.0*S->m[1][1] + sn*N->m[1][1] ;
    r->m[2][1] = -1.0*S->m[2][1] + sn*N->m[2][1] ;
    r->m[3][1] = -1.0*S->m[3][1] + sn*N->m[3][1] ;
    r->m[4][1] = 0.0 ;
    
    rn = dmat_normalize(r) ;
    delete_dmatrix(r) ;
    
    return rn ;
}

color_t color_init(double r, double g, double b) {
    
    color_t s ;
    
    s.r = r ;
    s.g = g ;
    s.b = b ;
    
    return s ;
}

color_t color_mult(double a, color_t c) {
    
    color_t s ;
    
    s.r = a*c.r ;
    s.g = a*c.g ;
    s.b = a*c.b ;
    
    return s ;
}

color_t color_add(color_t c1, color_t c2) {
    
    color_t s ;
    
    s.r = c1.r + c2.r ;
    s.g = c1.g + c2.g ;
    s.b = c1.b + c2.b ;
    
    return s ;
}

color_t shade(light_t *light, object_t *object, dmatrix_t *e, dmatrix_t *d, color_t color, color_t background) {

    color.r = background.r ;
    color.g = background.g ;
    color.b = background.b ;

    return color ;
}

object_t *build_object(int object_type, dmatrix_t *M, color_t ambient_color, color_t diffuse_color, color_t specular_color, double ambient_coeff, double diffuse_coeff, double specular_coeff, double f, double reflectivity) {
    
    object_t *object ;
    
    object = malloc(sizeof(*object));
    
    object->type = object_type ;
    object->M = *M ;
    dmat_alloc(&object->Minv,4,4) ;
    object->Minv = *dmat_inverse(&object->M) ;

    object->reflectivity = reflectivity ;
    
    object->specular_color.r = specular_color.r ;
    object->specular_color.g = specular_color.g ;
    object->specular_color.b = specular_color.b ;
    object->specular_coeff = specular_coeff ;
    object->f = f ;
    
    object->diffuse_color.r = diffuse_color.r ;
    object->diffuse_color.g = diffuse_color.g ;
    object->diffuse_color.b = diffuse_color.b ;
    object->diffuse_coeff = diffuse_coeff ;
    
    object->ambient_color.r = ambient_color.r ;
    object->ambient_color.g = ambient_color.g ;
    object->ambient_color.b = ambient_color.b ;
    object->ambient_coeff = ambient_coeff ;
    
    switch (object_type) {
          
        case SPHERE         :   object->intersection_function = &sphere_intersection ;
                                break ;
            
        case PLANE          :   object->intersection_function = &plane_intersection ;
                                break ;
            
        case INFINITE_PLANE :   object->intersection_function = &infinite_plane_intersection ;
                                break ;
    }
    
    nobjects++ ;
    return(object) ;
}

dmatrix_t *init_matrix(double x, double y, double z)
{
    dmatrix_t *matrix;
    matrix = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(matrix,4,1);
    matrix->m[1][1] = x ;
    matrix->m[2][1] = y ;
    matrix->m[3][1] = z ;
    matrix->m[4][1] = 0.0 ;

    return matrix;
}

color_t add_intensity(object_t *shape, dmatrix_t *s, dmatrix_t *v,dmatrix_t *n, dmatrix_t *r, light_t light)
{
    color_t light_intensity = light.intensity;
    double intensity_r, intensity_g, intensity_b;
    double ai_r, ai_g, ai_b;
    double di_r, di_g, di_b;
    double si_r, si_g, si_b;

    ai_r = shape->ambient_color.r * shape->ambient_coeff * 255;
    ai_g = shape->ambient_color.g * shape->ambient_coeff * 255;
    ai_b = shape->ambient_color.b * shape->ambient_coeff * 255;

    double SxN = ddot_product(s,n);
    if (SxN < 0) SxN = 0;

    di_r = shape->diffuse_color.r * SxN * light_intensity.r * shape->diffuse_coeff * 255;
    di_g = shape->diffuse_color.g * SxN * light_intensity.g * shape->diffuse_coeff * 255;
    di_b = shape->diffuse_color.b * SxN * light_intensity.b * shape->diffuse_coeff * 255;


    double RxV = pow(ddot_product(r,v), shape->f);
    if (RxV < 0) RxV = 0;

    si_r = shape->specular_color.r * RxV * light_intensity.r * shape->specular_coeff * 255;
    si_g = shape->specular_color.g * RxV * light_intensity.g * shape->specular_coeff * 255;
    si_b = shape->specular_color.b * RxV * light_intensity.b * shape->specular_coeff * 255;

    intensity_r = ai_r + di_r + si_r;
    intensity_g = ai_g + di_g + si_g;
    intensity_b = ai_b + di_b + si_b;

    color_t color = color_init(intensity_r, intensity_g, intensity_b);

    //printf("R %f   G %f   B %f\n", color.r,color.g,color.b);

    return color;
}

int shadowed(dmatrix_t *s, dmatrix_t *intersect, object_t *object, int nobjects)
{
    for (int k = 0; k < nobjects; k ++)
    {
        object_t* shape = &object[k];
        dmatrix_t *M_inverse = dmat_inverse(&shape->M);

        if (shape->type == 0);
        else if (shape->type == 1)
        {
            double u = plane_intersection(dmat_mult(M_inverse,intersect), dmat_mult(M_inverse, s));
            if (u > 0 && u < 1)
            { 
                return 1;   
            }
        }else if (shape->type == 2)
        {
            double u = sphere_intersection(dmat_mult(M_inverse,intersect), dmat_mult(M_inverse, s));
            if (u > 0 && u < 1) 
            {
                return 1;
            }
        }
    }
    return 0;
}

color_t add_intensity_shadow(object_t *shape)
{
    double ai_r = shape->ambient_color.r * shape->ambient_coeff * 255;
    double ai_g = shape->ambient_color.g * shape->ambient_coeff * 255;
    double ai_b = shape->ambient_color.b * shape->ambient_coeff * 255;

    //if (shape->type == 1)
        //printf("AR: %f, AG: %f, AB: %f\n", ai_r, ai_g, ai_b);

    return color_init(ai_r,ai_g,ai_b);
}

color_t get_color(object_t *object, int nobjects, dmatrix_t *point, dmatrix_t *direction, light_t light)
{
    dmatrix_t light_position = light.position;
    color_t color = color_init(50.0,50.0,50.0);

    double t_min = INFTY;
    int smallest_index = 0, done = 0;

    for (int i = 0; i < nobjects; i ++)
    {
        object_t* shape = &object[i];
        dmatrix_t *M_inverse = dmat_inverse(&shape->M);
        if (shape->type == 0);
        else if (shape->type == 1)
        {
            //printf("Checking plane\n");
            double t = plane_intersection(dmat_mult(M_inverse,point),dmat_mult(M_inverse,direction));
            if (t < t_min && t > 0)
            {
                t_min = t;
                smallest_index = i;
            }

            if (done)
            {
                //dmatrix_t *intersect = dmat_mult(&shape->M, intersection_coordinates(dmat_mult(M_inverse, point),dmat_mult(M_inverse, direction),t_min));
                dmatrix_t *intersect = dmat_mult(M_inverse,intersection_coordinates(point,direction,t_min));
                dmatrix_t *source = dmat_mult(M_inverse,&light_position);

                dmatrix_t *s = vector_to_light_source(intersect,source);
                dmatrix_t *v = vector_to_center_of_projection(intersect,dmat_mult(M_inverse, point));
                dmatrix_t *n = normal_to_surface(shape,intersect);
                dmatrix_t *r = vector_to_specular_reflection(n,s);

                if (!shadowed(s,dmat_mult(&shape->M,intersect),object,nobjects))
                    color = add_intensity(shape,s,v,n,r,light);
                else
                    color = add_intensity_shadow(shape);

                break;
            }

            if (i == nobjects - 1 && t_min != INFTY)
            {
                i = smallest_index - 1;
                done = 1;
            }
        }else if (shape->type == 2)
        {
            double t = sphere_intersection(dmat_mult(M_inverse,point),dmat_mult(M_inverse,direction));
            if (t < t_min && t > 0)
            {
                t_min = t;
                smallest_index = i;
            }

            if (done)
            {
                //dmatrix_t *intersect = dmat_mult(&shape->M, intersection_coordinates(dmat_mult(M_inverse, point),dmat_mult(M_inverse, direction),t_min));
                dmatrix_t *intersect = dmat_mult(M_inverse, intersection_coordinates(point,direction,t_min));
                dmatrix_t *source = dmat_mult(M_inverse,&light_position);

                dmatrix_t *s = vector_to_light_source(intersect, source);
                dmatrix_t *v = vector_to_center_of_projection(intersect,dmat_mult(M_inverse, point));
                dmatrix_t *n = normal_to_surface(shape,intersect);
                dmatrix_t *r = vector_to_specular_reflection(n,s);

                if (!shadowed(s,dmat_mult(&shape->M,intersect),object,nobjects))
                    color = add_intensity(shape,s,v,n,r,light);
                else
                    color = add_intensity_shadow(shape);
                break;
            }

            if (i == nobjects - 1 && t_min != INFTY)
            {
                i = smallest_index - 1;
                done = 1;
            }
        }else
            return color;
    }

    if (color.r > 255) color.r = 255;
    if (color.g > 255) color.g = 255;
    if (color.b > 255) color.b = 255;

    return color;
}



int main() {
    
    Display *d ;
    Window w ;
    XEvent e ;
    
    int i, j, s ;
    
    camera_t Camera ;
    window_t Window ;
    light_t light ;
    dmatrix_t M, light_position ;
    color_t pixel, background, light_intensity, light_color, ambient_color, diffuse_color, specular_color ;
    double height, width, aspect, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity ;
    
    /* Set the background color */
    
    background.r = 0.0 ;
    background.g = 0.0 ;
    background.b = 0.0 ;
    
    /* Set up light position, intensity, and color */
    
    dmat_alloc(&light_position,4,1) ;
    
    light_position.m[1][1] = Lx ;
    light_position.m[2][1] = Ly ;
    light_position.m[3][1] = Lz ;
    light_position.m[4][1] = 1.0 ;
    
    light_intensity.r = 1.0 ;
    light_intensity.g = 1.0 ;
    light_intensity.b = 1.0 ;
    
    light_color.r = 1.0 ;
    light_color.g = 1.0 ;
    light_color.b = 1.0 ;
    
    light = *build_light(&light,&light_position,light_color,light_intensity) ;
    
    /* Build display window and synthetic camera */
    
    Window = *build_window(&Window,H,ASPECT) ;
    Camera = *build_camera(&Camera,&Window) ;
    
    /* Create a sphere */
    
    dmat_alloc(&M,4,4) ;
    M = *dmat_identity(&M) ;

    M.m[1][4] = 0.0 ;
    M.m[2][4] = 0.0 ;
    M.m[3][4] = -1.0;
    
    reflectivity = 0.2 ;
    
    specular_color.r = 1.0 ;
    specular_color.g = 1.0 ;
    specular_color.b = 1.0 ;
    specular_coeff = 0.4 ;
    f = 10.0 ;
    
    diffuse_color.r = 1.0 ;
    diffuse_color.g = 0.0 ;
    diffuse_color.b = 0.0 ;
    diffuse_coeff = 0.4 ;
    
    ambient_color.r = 1.0 ;
    ambient_color.g = 0.0 ;
    ambient_color.b = 0.0 ;
    ambient_coeff = 0.2 ;
    
    object[nobjects] = *build_object(SPHERE,&M,ambient_color,diffuse_color,specular_color,ambient_coeff,diffuse_coeff,specular_coeff,f,reflectivity) ;
    
    
    //Sphere 2
    
    dmat_alloc(&M,4,4) ;
    M = *dmat_identity(&M) ;
    
    M.m[1][1] = 2.0 ;
    M.m[2][2] = 2.0 ;
    M.m[1][4] = 4.0 ;
    M.m[2][4] = 2.0 ;
    M.m[3][3] = 2.0 ;
 
    
    reflectivity = 0.2 ;
    
    specular_color.r = 1.0 ;
    specular_color.g = 1.0 ;
    specular_color.b = 1.0 ;
    specular_coeff = 0.4 ;
    f = 10.0 ;
    
    diffuse_color.r = 0.4 ;
    diffuse_color.g = 1.0 ;
    diffuse_color.b = 1.0 ;
    diffuse_coeff = 0.4 ;
    
    ambient_color.r = 0.4 ;
    ambient_color.g = 1.0 ;
    ambient_color.b = 1.0 ;
    ambient_coeff = 0.2 ;
    
    object[nobjects] = *build_object(SPHERE,&M,ambient_color,diffuse_color,specular_color,ambient_coeff,diffuse_coeff,specular_coeff,f,reflectivity) ;

    //Sphere 3
    dmat_alloc(&M,4,4) ;
    M = *dmat_identity(&M) ;

    M.m[1][1] = 2.0 ;
    M.m[2][2] = 2.0 ;
    M.m[1][4] = -4.0 ;
    M.m[2][4] = 2.0  ;
    M.m[3][3] = 2.0  ;
    
    reflectivity = 0.2 ;
    
    specular_color.r = 1.0 ;
    specular_color.g = 1.0 ;
    specular_color.b = 1.0 ;
    specular_coeff = 0.4 ;
    f = 10.0 ;
    
    diffuse_color.r = 0.4 ;
    diffuse_color.g = 1.0 ;
    diffuse_color.b = 1.0 ;
    diffuse_coeff = 0.4 ;
    
    ambient_color.r = 0.4 ;
    ambient_color.g = 1.0 ;
    ambient_color.b = 1.0 ;
    ambient_coeff = 0.2 ;
    
    object[nobjects] = *build_object(SPHERE,&M,ambient_color,diffuse_color,specular_color,ambient_coeff,diffuse_coeff,specular_coeff,f,reflectivity) ;
    
    //Sphere 4
    dmat_alloc(&M,4,4) ;
    M = *dmat_identity(&M) ;

    M.m[1][4] = 0.0 ;
    M.m[1][1] = 4.0 ;
    M.m[2][2] = 0.5 ;
    M.m[2][4] = -3.0 ;
    M.m[3][4] = -1.0 ;
    
    reflectivity = 0.2 ;
    
    specular_color.r = 1.0 ;
    specular_color.g = 1.0 ;
    specular_color.b = 1.0 ;
    specular_coeff = 0.4 ;
    f = 10.0 ;
    
    diffuse_color.r = 0.8 ;
    diffuse_color.g = 0.4 ;
    diffuse_color.b = 0.0 ;
    diffuse_coeff = 0.4 ;
    
    ambient_color.r = 0.8 ;
    ambient_color.g = 0.4 ;
    ambient_color.b = 0.0 ;
    ambient_coeff = 0.2 ;
    
    object[nobjects] = *build_object(SPHERE,&M,ambient_color,diffuse_color,specular_color,ambient_coeff,diffuse_coeff,specular_coeff,f,reflectivity) ;

    
    //Create tiled floor with planes
    
    for (i = -4 ; i < 4 ; i++) {
        for (j = -4 ; j < 4 ; j++) {
     
            dmat_alloc(&M,4,4) ;
            M = *dmat_identity(&M) ;
    
            M.m[1][4] = (double)2*i ;
            M.m[2][4] = (double)2*j ;
            M.m[3][4] = -2.0 ;
    
            reflectivity = 0.5 ;
    
            specular_color.r = 1.0 ;
            specular_color.g = 1.0 ;
            specular_color.b = 1.0 ;
            specular_coeff = 0.2 ;
            f = 10.0 ;
    
            diffuse_color.r = 1.0 ;
            diffuse_color.g = 1.0 ;
            diffuse_color.b = 1.0 ;
            diffuse_coeff = 0.2 ;
    
            ambient_color.r = (double)(abs(i+j)%2) ;
            ambient_color.g = (double)(abs(i+j)%2) ;
            ambient_color.b = (double)(abs(i+j)%2) ;
            ambient_coeff = 0.6 ;
            
            object[nobjects] = *build_object(PLANE,&M,ambient_color,diffuse_color,specular_color,ambient_coeff,diffuse_coeff,specular_coeff,f,reflectivity) ;
        }
    }
    
    
    /*
    dmat_alloc(&M,4,4) ;
    M = *dmat_identity(&M) ;

    M.m[1][1] = 2.0;
    M.m[2][2] = 2.0; 
    M.m[3][4] = -2.0;

    reflectivity = 0.1 ;

    specular_color.r = 1.0 ;
    specular_color.g = 1.0 ;
    specular_color.b = 1.0 ;
    specular_coeff = 0.4 ;
    f = 10.0 ;

    diffuse_color.r = 1.0 ;
    diffuse_color.g = 1.0 ;
    diffuse_color.b = 1.0 ;
    diffuse_coeff = 0.4 ;

    ambient_color.r = 1.0 ;
    ambient_color.g = 1.0 ;
    ambient_color.b = 1.0 ;
    ambient_coeff = 0.2 ;
    
    object[nobjects] = *build_object(PLANE,&M,ambient_color,diffuse_color,specular_color,ambient_coeff,diffuse_coeff,specular_coeff,f,reflectivity) ;
    
    */
 
    /* Set near plane dimensions */
    
    aspect = ASPECT ;
    height = Near*tan(M_PI/180.0*THETA/2.0) ;
    width = height*aspect ;
    
    dmatrix_t *direction ;
    
    d = InitX(d,&w,&s,&Window) ;
    XNextEvent(d, &e) ;
    
    while (1) {
        XNextEvent(d,&e) ;
        if (e.type == Expose) {
            
            for (i = 0 ; i < Window.width ; i++) {
                for (j = 0  ; j < Window.height ; j++) {
                    direction = ray_direction(&Camera,&Window,height,width,(double)i,(double)j) ;
                    pixel = shade(&light,object,&Camera.E,direction,pixel,background) ;
                    double t = INFTY;

                    color_t color = get_color(object,nobjects,&Camera.E,direction,light);
                    SetCurrentColorX(d,&(DefaultGC(d,s)),color.r,color.g,color.b);
                    SetPixelX(d,w,s,i,Window.height - (j + 1)) ;

                    delete_dmatrix(direction) ;
                }

            }
        }
        if (e.type == KeyPress)
            break ;
        if (e.type == ClientMessage)
            break ;
    }
    QuitX(d,w) ;
}