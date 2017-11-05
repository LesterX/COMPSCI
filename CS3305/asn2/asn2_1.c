/*
Computer Science 3305
Assignment 2 Part 1
Yimin Xu
250876566
*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <wait.h>
#include <pthread.h>

//Package the three parameters
typedef struct nums
{
    int x;
    int y;
    int z;
} num;

//Thread function
void *thread_func(void* ptr)
{
    num *xyz = (num*) ptr;
    (*xyz).z = (*xyz).x + (*xyz).y;
    printf("Thread completed\n");

    return (void *) 0;
}

int main()
{
    num xyz;

    xyz.x = 10;
    xyz.y = 20;
    xyz.z = 0;

    pid_t pid = fork();

    if (pid < 0)
    {
        perror("fork");
        printf("Error\n");
        exit(pid);
    }else if (pid > 0)
    {
        //Parent
        wait(NULL);
        //Print the result after child is completed
        printf("z = %d\n", xyz.z);

        pthread_t t;

        //Create thread
        if (pthread_create(&t,NULL,thread_func,&xyz))
        {
            printf("Error: Thread\n");
            return 1;
        }

        //Wait for thread to complete
        if (pthread_join(t,NULL))
        {
            printf("Error: Joining\n");
            return 1;
        }

        //Print the result after thread is completed
	    printf("z = %d\n", xyz.z);
    }else
    {
        //Child
        xyz.z = xyz.x + xyz.y;
	printf("Child completed\n");
    }
}
