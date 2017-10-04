#include <stdio.h>
#include <stdlib.h>

int main()
{
	char x_in[100];
	int X;
	printf("Enter the value of X: \n");
	scanf(" %s",x_in);
	X = atoi(x_in);
	printf("X = %d\n", X);	
}
