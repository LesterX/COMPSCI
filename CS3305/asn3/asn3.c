/*
Computer Science 3305
Assignment 3
Yimin Xu
250876566
*/

#include<stdio.h>
#include<string.h>
#include<pthread.h>
#include<stdlib.h>
#include<unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

typedef struct account_request
{
	char request_type;
	int account;
	int target;
	int amount;
} request;

pthread_mutex_t lock;

int balance_acc1 = 0;
int balance_acc2 = 0;
int balance_acc3 = 0;
int balance_acc4 = 0;
int balance_acc5 = 0;

int *find_account(int account)
{
	if (account == 1)
		return &balance_acc1;
	else if (account == 2)
		return &balance_acc2;
	else if (account == 3)
		return &balance_acc3;
	else if (account == 4)
		return &balance_acc4;
	else if (account == 5)
		return &balance_acc5;
}

//Thread function
void *bank_account_request(void* rqst)
{
	request *rq = (request*) rqst;

	pthread_mutex_lock(&lock);

	//printf("\n\nProcessing request... \n");
	//printf("Current balance of Account 1: %d\nCurrent balance of Account 2: %d\n", balance_acc1, balance_acc2);

	if (rq->request_type == 'd')
	{
		//printf("Deposit to %d: $%d\n", rq->account, rq->amount);
		*find_account(rq->account) += rq->amount;
	}else if (rq->request_type == 'w')
	{
		//printf("Withdraw from %d: $%d\n", rq->account, rq->amount);
		if (*find_account(rq->account) < rq->amount)
		{
			//printf("Request denied: balance is not enough\n");
			pthread_mutex_unlock(&lock);
			return 0;
		}
		else
		{
			*find_account(rq->account) -= rq->amount;
			//printf("Withdraw Authorized\n");
		}
	}else if (rq->request_type == 't')
	{
		//printf("Transfer from %d to %d: $%d\n", rq->account, rq->target, rq->amount);
		if (*find_account(rq->account) < rq->amount)
		{
			//printf("Request denied: balance is not enough\n");
			pthread_mutex_unlock(&lock);
			return 0;
		}
		else
		{
			*find_account(rq->account) -= rq->amount;
			*find_account(rq->target) += rq->amount;
			//printf("Transfer Authorized\n");
		}
	}

	//printf("Updated Balance:\n");
	//printf("Account 1: %d\nAccount 2: %d\n", balance_acc1, balance_acc2);

	pthread_mutex_unlock(&lock);

	return 0;
}

int main(void)
{
	int k = 2,n = 2,m = 10;

	request d_rq[n * 4];
	request c_rq[m * 4];

	//Read file
    FILE *fp;
    char line[100];

    fp = fopen("./input.txt","r");
    if (fp == NULL)
    {
        printf("Input file cannot be found");
        return 0;
    }

    int d_count = 0, c_count = 0;

    while(fgets(line,sizeof(line),fp))
    {
    	//Build a request, separate depositors and clients

    	char *str = strtok(line," ");
    	char user = str[0];

    	str = strtok(NULL," ");
    	
    	int str_count = 0;
    	while (str != NULL && str[0] != ' ')
    	{	
    		if (user == 'd')
    		{
    			d_rq[d_count].request_type = str[0];

	    		str = strtok(NULL," ");
	    		
	    		d_rq[d_count].account = str[3] - '0';

	    		if ((d_rq)[d_count].request_type == 't')
	    		{
	    			str = strtok(NULL," ");
	    			d_rq[d_count].target = str[3] - '0';
	    		}else
	    			d_rq[d_count].target = d_rq[d_count].account;
	    		
	    		str = strtok(NULL," ");

	    		d_rq[d_count].amount = atoi(str);

	    		str = strtok(NULL," ");
	    		d_count ++;
    		}else
    		{
    			c_rq[c_count].request_type = str[0];

	    		str = strtok(NULL," ");
	    		
	    		c_rq[c_count].account = str[3] - '0';
	    	
	    		if ((c_rq)[c_count].request_type == 't')
	    		{
	    			str = strtok(NULL," ");
	    			c_rq[c_count].target = str[3] - '0';
	    		}else
	    		{
					c_rq[c_count].target = c_rq[c_count].account;
					//printf("Client Count = %d!!!!!!!!\n\n\n",c_count);
			    }
	    		
	    		str = strtok(NULL," ");
	    		c_rq[c_count].amount = atoi(str);
	    		str = strtok(NULL," ");
	    		c_count ++;
    		}
    	}
    }

    fclose(fp);
    //Read completed

    int err_thread;
    
    pthread_t* ptr_d = malloc(sizeof(pthread_t) * d_count);
    pthread_t* ptr_c = malloc(sizeof(pthread_t) * c_count);


    //Set the lock
    if (pthread_mutex_init(&lock, NULL) != 0)
	{
	        printf("\n mutex init failed\n");
	        return 1;
	}

	//Create the threads, do the depositors first
    for (int i = 0; i < d_count; i ++)
    {
		err_thread = pthread_create(&ptr_d[i], NULL, bank_account_request, (void *) &d_rq[i]);
		if (err_thread != 0)
			printf("Error in creating thread\n");
	}

    for (int i = 0; i < d_count; i ++)
    	pthread_join(ptr_d[i], NULL);

    for (int i = 0; i < c_count; i ++)
    {
    	pthread_create(&ptr_c[i], NULL, bank_account_request, (void *) &c_rq[i]);
    	if (err_thread != 0)
			printf("Error in creating thread\n");
	}

    for (int i = 0; i < c_count; i ++)
    	pthread_join(ptr_c[i], NULL);

    pthread_mutex_destroy(&lock);

    //Print the output to the file 
    fp = fopen("Assignment_3_output_file.txt","w");
    if (fp == NULL)
    {
    	printf("Error opening file\n");
    	exit(1);
    }

    for (int i = 1; i < 5; i ++)
    {
    	if (*find_account(i) > 0)
    		fprintf(fp,"Account %d: $%d    ", i, *find_account(i));
    }
    
    fclose(fp);

	return 0;
}
