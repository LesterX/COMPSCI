#include <stdio.h>
#include <unistd.h>
#include <errno.h>

int main()
{
    FILE *fd1, *fd2;
    char c;
    pid_t pid;

    fd1 = fopen("foobar.txt","r");

    pid = fork();
    if (pid > 0)
    {
        fscanf(fd1,"%c\n",&c);
        printf("parent : c = %c",c);
    }else if(pid == 0)
    {
        fscanf(fd1,"%c\n",&c);
        printf("child : c = %c",c);
    }
}
