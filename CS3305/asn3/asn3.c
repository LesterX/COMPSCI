#include<stdio.h>
#include<string.h>
#include<pthread.h>
#include<stdlib.h>
#include<unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

/*
#define thread_1_req 400
#define thread_2_req 300
#define thread_3_req 200
#define thread_4_req 150

void  *withdraw_money_thread_1();
void  *withdraw_money_thread_2();
void  *withdraw_money_thread_3();
void  *withdraw_money_thread_4();
*/

typedef struct account_request
{
	char request_type;
	int account;
	int target;
	int amount;
} request;

pthread_mutex_t lock;

int shared_balance = 0;

void *withdraw_money_thread(int amount)
{

	pthread_mutex_lock(&lock);

	printf("\n Processing withdraw request.. current balance %d", shared_balance);


	if (shared_balance >= amount)
		printf("\n Withdraw amount of %d is Authorized", amount);
	else
	{
		printf("\n Request denied");
		return 0;
	}

	shared_balance = shared_balance - amount;
	printf("\n Updated balance %d", shared_balance);

	pthread_mutex_unlock(&lock);

	return 0;
}

int main(void)
{
	int k = 2,n = 2,m = 10;

	request rqs[m + n];

	//Read file
    FILE *fp;
    char line[100];

    fp = fopen("./input.txt","r");
    if (fp == NULL)
    {
        printf("Input file cannot be found");
        return 0;
    }

    int count = 0;

    while(fgets(line,sizeof(line),fp))
    {
    	printf("Line: %s\n", line);
    	char *str = strtok(line," ");
    	str = strtok(NULL," ");
    	
    	int str_count = 0;
    	while (str != NULL && str[0] != ' ')
    	{	
    		printf("	Type: %s\n",str);
    		rqs[count].request_type = str[0];
    		str = strtok(NULL," ");
    		rqs[count].account = str[3];
    		printf("	Account: %s\n",str);

    		if (rqs[count].request_type == 't')
    		{
    			str = strtok(NULL," ");
    			printf("	Target: %s\n",str);
    			rqs[count].target = atoi(str);
    		}
    		
    		str = strtok(NULL," ");
    		printf("	Amount: %s\n\n",str);

    		if (str == NULL)
    			printf("Hi\n");

    		rqs[count].amount = atoi(str);

    		str = strtok(NULL," ");
    		count ++;

    	}
    }
	
	/*
	int i;
	int err_thread;

	pthread_t thread_1, thread_2, thread_3, thread_4;


	if (pthread_mutex_init(&lock, NULL) != 0)
	    {
	        printf("\n mutex init failed\n");
	        return 1;
	    }


	err_thread = pthread_create(&thread_1, NULL, &withdraw_money_thread_1, NULL); 

	if (err_thread != 0)
		printf("\n Error creating thread_1");



	err_thread = pthread_create(&thread_2, NULL, &withdraw_money_thread_2, NULL); 

	if (err_thread != 0)
		printf("\n Error creating thread_2");

	err_thread = pthread_create(&thread_3, NULL, &withdraw_money_thread_3, NULL); 

	if (err_thread != 0)
		printf("\n Error creating thread_3");

	err_thread = pthread_create(&thread_4, NULL, &withdraw_money_thread_4, NULL); 

	if (err_thread != 0)
		printf("\n Error creating thread_4");

	pthread_join(thread_1, NULL);
	pthread_join(thread_2, NULL);
	pthread_join(thread_3, NULL);
	pthread_join(thread_4, NULL);


	pthread_mutex_destroy(&lock);

	return 0;
	*/
}
