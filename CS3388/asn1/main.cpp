/*
Computer Science 3388
Assignment 1
Yimin Xu
250876566
*/

#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>

static void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}

//Function to set a pixel
void set_pixel(int x, int y)
{
    glBegin(GL_POINTS);
    glColor3i(0,0,0);
    glVertex2i(x, y);
    glEnd();
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
        float ratio;
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        ratio = width / (float) height;
        glViewport(0, 0, width, height);
        glClear(GL_COLOR_BUFFER_BIT);
        //Set background color to white
        glClearColor(255,255,255,0);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        int x0 = 0;
        int y0 = 0;
        while (x0 < 100 && y0 < 100)
        {
            set_pixel(x0, y0);
            x0 ++;
            y0 ++;
        }

        glFlush();


        glEnd();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}
