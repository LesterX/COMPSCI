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
        wait(NULL);
	pid = fork();
    	if (pid < 0)
    	{
        	perror("fork");
        	printf("Error\n");
        	exit(pid);
    	}else if (pid == 0)
    	{
        	printf("parent (PID %d) created child_2 (PID %d)\n",getppid(),getpid());
        	printf("child_2 (PID %d) is calling an external program B.out and leaving child_2...\n",getpid());
        	printf("From the external program B: \n");
        	execl("./B.out","./B.out",NULL);
    	}
    }else
    {
        printf("parent process (PID %d) created child_1 (PID %d)\n",getppid(),getpid());
        printf("parent (PID %d) is waiting for child_1 (PID %d) to complete before creating child_2\n",getppid(),getpid());
        pid = fork();
        if (pid < 0)
        {
            perror("fork");
            printf("Error\n");
            exit(pid);
        }
        else if (pid == 0)
        {
            printf("child_1 (PID %d) created child_1.1 (PID %d)\n",getppid(),getpid());
            printf("child_1 (PID %d) is now completed\n",getpid());
        }
    }
}
