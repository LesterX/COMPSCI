/*
Computer Science 3305
Assignment 1
Yimin Xu
250876566
*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <wait.h>

int main()
{
    pid_t pid = fork();

    if (pid < 0)
    {
        perror("fork");
        printf("Error\n");
        exit(pid);
    }else if (pid > 0)
    {
        //Parent
        printf("parent process (PID %d) created child_1 (PID %d)\n",getpid(),pid);
        printf("parent (PID %d) is waiting for child_1 (PID %d) to complete before creating child_2\n",getpid(),pid);
        wait(NULL);
        pid = fork();
    	if (pid < 0)
    	{
        	perror("fork");
        	printf("Error\n");
        	exit(pid);
    	}else if (pid == 0)
    	{
            //Child 2
        	printf("child_2 (PID %d) is calling an external program B.out and leaving child_2...\n",getpid());
        	printf("From the external program B: \n");
        	execl("./B.out","./B.out",NULL);
    	}else
    	{
            //Parent
            printf("parent (PID %d) created child_2 (PID %d)\n",getpid(),pid);

    	}
    }else
    {
        //Child 1
        pid = fork();
        if (pid < 0)
        {
            perror("fork");
            printf("Error\n");
            exit(pid);
        }
        else if (pid > 0)
        {
            printf("child_1 (PID %d) created child_1.1 (PID %d)\n",getpid(),pid);
            printf("child_1 (PID %d) is now completed\n",pid);
        }
    }
}
