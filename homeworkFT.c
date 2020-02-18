#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <complex.h>


#define PI 3.14159265358979323846

int P;
int N;
int nrElemsPerThread;

FILE *inputFile;
FILE *outputFile;

double *Xi;
double complex *var;


void getArgs(int argc, char **argv)
{
	if(argc < 4) {
		printf("Not enough paramters: ./program input output P\n");
		exit(-1);
	}
	
	inputFile = fopen(argv[1], "r");
	if(inputFile == NULL)	exit(2);

	outputFile = fopen(argv[2], "w");
	if(outputFile == NULL)	exit(3);

	P = atoi(argv[3]);
}




void* threadFunction(void *args)
{
	int thread_id = *(int*)args;

	int k, n;

	//calculez limitele subvectorului de prelucrat 
	//pentru thread-ul curent
	int start = thread_id * nrElemsPerThread;
	int end = fmin(N, start + nrElemsPerThread);

	//ultimul thread calculeaza pana la finalul sirului 
	if(thread_id == (P - 1)){
		end = N;
	}


	for(k = start; k < end; k++){
		var[k] = 0;

		for(n = 0; n < N; n++){
			var[k] += Xi[n] * cexp(-2 * PI * I / N * k * n);
		}
	}
	
	return NULL;
}


int main(int argc, char * argv[]) {
	getArgs(argc, argv);
	
	pthread_t tid[P];
	int thread_id[P];

	int i;
	int error;

	error = fscanf(inputFile, "%d", &N);
	if(error == 0)	exit(-2);

	fprintf(outputFile ,"%d\n", N);

	Xi = malloc(N * sizeof(double));
	if(Xi == NULL)	exit(-3);

	var = malloc(N * sizeof(double complex));
	if(var == NULL)	exit(-4);

	for(i = 0; i < N; i++){
		error = fscanf(inputFile, "%lf", &Xi[i]);
		if(error == 0)	exit(-3);
	}

	
	for(i = 0;i < P; i++)
		thread_id[i] = i;


	//numarul de elemenete ce revin fiecarui thread de prelucrat
	nrElemsPerThread = fmin((double) (N - ceil(N/P)), (double) ceil(N/P));

	for(i = 0; i < P; i++) {
		pthread_create(&(tid[i]), NULL, threadFunction, &(thread_id[i]));
	}

	for(i = 0; i < P; i++) {
		pthread_join(tid[i], NULL);
	}


	for(i = 0; i < N; i++){
		fprintf(outputFile, "%lf %lf\n", creal(var[i]), cimag(var[i]));
	}


	fclose(inputFile);
	fclose(outputFile);

	free(Xi);
	free(var);

	return 0;
}
