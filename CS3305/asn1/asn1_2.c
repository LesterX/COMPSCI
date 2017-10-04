#include <stdio.h>
#include <stdlib.h>
#include <sys/type.h>
#include <unistd.h>
#include <errno.h>
#include <wait.h>

int main()
{
	int fd[2];
	pid_t pid;
	int X,Y;

	pid = fork();
	if (pid < 0)
		perror("pipe error");
	else if (pid == 0)
	{
		//Child
		char in_Y[100];
		printf("Enter the value of Y: \n");
		scanf("%s", in_Y);
		write(fd[1],in_Y,100);
	}else if (pid > 0)
	{
		//Parent
		char* in_X[100];
		printf("Enter the value of X: \n");
		scanf("%s", in_X);
		X = atoi(in_X);
		printf("X = %d\n", X);
	}
}
