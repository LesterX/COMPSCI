/*
Computer Science 3388
Assignment 1
Yimin Xu
250876566
*/

#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//Error handling
static void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}

//Type Q to close the window
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_Q && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}

//Set a pixel to black
void set_pixel(int x, int y)
{
    glBegin(GL_POINTS);
    glColor3i(0,0,0);
    glVertex2i(x, y);
    glEnd();
}

//Draw a line from (x1,y1) to (x2,y2), following Bresenham Algorithm
void Bresenham(int x1, int y1, int x2, int y2)
{
    //In case there is no slope
    if (x1 == x2)
    {
        int s;
        int e;
        if (y1 < y2)
        {
            s = y1;
            e = y2;
        }
        else
        {
            s = y2;
            e = y1;
        }
        for (s; s < e; s ++)
            set_pixel(x1,s);
    }

    //Make x1 < x2 and if the two points are swapped, change it back at the end
    bool swapped = false;
    int temp1, temp2;
    if (x1 > x2)
    {
        swapped = true;
        temp1 = x1;
        temp2 = y1;
        x1 = x2;
        y1 = y2;
        x2 = temp1;
        y2 = temp2;
    }

    //In case the slope is negative
    bool neg_slope = false;
    if (y1 > y2)
        neg_slope = true;

    int p;
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);

    //Abs of slope is smaller than 1
    if (abs(dy) <= abs(dx))
    {
        set_pixel(x1, y1);

        int i = x1;
        for (i ; i < x2; i ++)
        {
            if (i == x1)
                p = 2 * dy - dx;
            else
            {
                if (p < 0)
                {
                    p = p + 2 * dy;
                }
                else
                {
                    p = p + 2 * dy - 2 * dx;
                    //If the slope is positive, go up, otherwise, go down
                    if (!neg_slope)
                        y1 ++;
                    else
                        y1 --;
                }
                x1 ++;

                set_pixel(x1, y1);
            }
        }
    }else //Abs of slope greater than 1, swap x and y, basically the same
    {
        set_pixel(x1, y1);

        int j = y1;
        for (j;;)
        {

            if (j == y1)
                p = 2 * dx - dy;
            else
            {
                if (p < 0)
                    p = p + 2 * dx;
                else
                {
                    p = p + 2 * dx - 2 * dy;
                    x1 ++;
                }
                if (neg_slope)
                    y1 --;
                else
                    y1 ++;

                set_pixel(x1, y1);
            }
            if (neg_slope)
            {
                j --;
                if (j <= y2)
                    break;
            }
            else
            {
                j ++;
                if (j >= y2)
                    break;
            }
        }
    }

    if (swapped)
    {
        swapped = true;
        temp1 = x1;
        temp2 = y1;
        x1 = x2;
        y1 = y2;
        x2 = temp1;
        y2 = temp2;
    }
}

int main(void)
{
    //Create the window
    GLFWwindow* window;
    glfwSetErrorCallback(error_callback);
    if (!glfwInit())
        exit(EXIT_FAILURE);

    //Size: 512*512
    window = glfwCreateWindow(512, 512, "Assignment 1", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }
    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);
    while (!glfwWindowShouldClose(window))
    {
        glClear(GL_COLOR_BUFFER_BIT);
        //Set background color to white
        glClearColor(255,255,255,0);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        glOrtho(0, 512, 0, 512, -1, 1);


        //Start to draw
        double dt = 2.0 * M_PI / 200.0;
        for (double t = 0.0; t < 2.0 * M_PI;)
        {
            int x1 = 256 + (int) 100.0 * (1.5 * cos(t) - cos(13.0 * t));
            int y1 = 256 + (int) 100.0 * (1.5 * sin(t) - sin(13.0 * t));
            t += dt;
            int x2 = 256 + (int) 100.0 * (1.5 * cos(t) - cos(13.0 * t));
            int y2 = 256 + (int) 100.0 * (1.5 * sin(t) - sin(13.0 * t));
            Bresenham(x1,y1,x2,y2);
        }

        glFlush();

        glEnd();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}
