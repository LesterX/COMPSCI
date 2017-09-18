#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>

int main()
{
	pid_t pid1,pid2,pid3;
	int a = 0;
	pid1 = fork();
	pid2 = fork();
	pid3 = fork();
	a = getpid();
	printf("pid = %d\n", a);	
}
