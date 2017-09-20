#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>

int main()
{
    FILE  *fd1, *fd2;
    char c;
    pid_t pid;

    fd1 = fopen("foo.txt", "w");
    fd2 = fopen("foo.txt", "w"); 

    fprintf(fd1, "%s", "Anwar");
    fprintf(fd2, "%s", "Haque");
    
}

