#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <wait.h>

int main()
{
	int fd[2];
	pid_t pid;
	int X,Y;
    char in_X[100], in_Y[100];

    if (pipe(fd) < 0)
        perror("pipe error");

    printf("A pipe is crafted for communication between parent (PID %d) and child\n",getpid());

	pid = fork();
	if (pid < 0)
		perror("fork error");
	else if (pid == 0)
	{
		//Child
        printf("parent (PID %d) created a child (PID %d)\n",getppid(),getpid());
        close(fd[0]);
		printf("child (PID %d) Reading Y = ",getpid());
        scanf(" %s", in_Y);
        printf(" from the user\n");
		write(fd[1],in_Y,strlen(in_Y));
        printf("child (PID %d) Writing Y into the pipe\n",getpid());
	}else
	{
		//Parent
        wait(NULL);
        close(fd[1]);
		printf("parent (PID %d) Reading X = ",getpid());
        scanf(" %s", in_X);
        printf(" from the user\n");
		X = atoi(in_X);
        read(fd[0],in_Y,100);
        printf("parent (PID %d) Reading Y from the pipe (Y = %s)\n",getpid(),in_Y);
        Y = atoi(in_Y);
		printf("parent (PID %d) adding X + Y = %d\n",getpid(),X + Y);
	}
}
