#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>

/*  

   This program forks a process.  The child and parent processes print to terminal to identify themselves

*/

int main()
{ 
  pid_t pid;          	
  int i;
         
  pid = fork();
  
   if (pid < 0){
    perror("fork");
    printf("\n main function: Unable to fork...errno number is %d\n",errno);
    exit(pid);
  }
  
  if (pid > 0){
   printf("\n From Parent Process: I am parent\n"); 
   wait(NULL);
   
  }
  
  if (pid == 0) { 
   printf("\n From Child Process: I am child\n"); 
  }
  return 0;
}
