#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>

/*

   Parent process sums numbers from 1 to 9; 
   Child  process subtracts from 0 the sum from 1 to 9

*/

int main()
{ 
  pid_t pid;          
  int i;
  int sum = 0; 
         
  pid = fork(); 

  if (pid < 0){
    perror("fork");
    printf("main function: errno number is %d\n",errno);
    exit(pid);
  }

  if( pid > 0 ) 
  { 
    for( i=0; i < 10; i++) 
         sum = sum + i;
    printf("parent: sum is %d\n",sum);
    wait(NULL);
  }
  else
  {
    for( i=0; i < 10; i++ ) 
         sum = sum - i;
    printf("child: sum is %d\n",sum);
  }
  return 0;
}
