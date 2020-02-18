#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <complex.h>

typedef double complex cplx;
cplx *Xi;

int P;
int N;

FILE *inputFile;
FILE *outputFile;

double PI;

pthread_t thread[10];


//structura utila pentru a transforma fft 
//intr-o functie cu un singur parametru
typedef struct {
	cplx *buf;
	cplx *out;
	int step;
} Parameter;


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

//functia de pe rosetta modificata astfel incat sa aiba un singur argument
void *fft(void* argc){

	Parameter* parameter = (Parameter*) argc;

	if ((parameter -> step) < N) {
		Parameter p1, p2;

		p1.buf = parameter -> out;
		p1.out = parameter -> buf;
		p1.step = (parameter -> step) * 2;

		p2.buf = (parameter -> out) + (parameter -> step);
		p2.out = (parameter -> buf) + (parameter -> step);
		p2.step = (parameter -> step) * 2;

		fft(&p1);
		fft(&p2);
 
		for (int i = 0; i < N; i += 2 * (parameter -> step)) {
			cplx t = cexp(-I * PI * i / N) * (parameter -> out)[i + (parameter -> step)];
			(parameter -> buf)[i / 2]     = (parameter -> out)[i] + t;
			(parameter -> buf)[(i + N)/2] = (parameter -> out)[i] - t;
		}
	}

	return NULL;
}


int main(int argc, char * argv[]) {
	getArgs(argc, argv);

	PI = atan2(1, 1) * 4;
	
 	int error;

	error = fscanf(inputFile, "%d", &N);
	if(error == 0)	exit(-2);

	error = fprintf(outputFile ,"%d\n", N);
	if(error == 0)	exit(-3);

	Xi = malloc(sizeof(cplx) * N);
	if(Xi == NULL) {
		exit(-4);
	}

	//citesc din input vectorul de prelucrat, Xi
	for(int i = 0; i < N; i++){
		double tmp;

		error = fscanf(inputFile, "%lf", &tmp);
		if(error == 0)	exit(-5);

		Xi[i] = tmp;
	}


	cplx *out = (cplx*) malloc(sizeof(cplx) * N);

	//structura initiala(parametrul initial al functiei)
	Parameter parameter;
	parameter.buf = Xi;
	parameter.out = out;
	parameter.step = 1;
 
	//daca am un singur core nu creez thread-uri
	if(P == 1){
		fft(&parameter);
	}


	//caz doua core-uri
	if(P == 2){
		//creez parametrii corespunzatori apelurilor functiilor de la nivel 1
		//numerotarea nivelelor incepand de la 0
		Parameter p1,p2;
		p1.step = parameter.step * 2;
		p1.buf = parameter.out;
		p1.out = parameter.buf;

		p2.out = parameter.buf + 1;
		p2.buf = parameter.out + 1;
		p2.step = parameter.step * 2;

		pthread_create(&thread[0], NULL, fft, &p1);
		pthread_create(&thread[1], NULL, fft, &p2);
		
		for(int i = 0; i< P; i++){
			pthread_join(thread[i], NULL);
		}

		//execut for-ul de la primul nivel
		//dupa terminarea executiilor de la nivelele urmatoare
		for (int i = 0; i < N; i += 2 * parameter.step) {
			cplx t = cexp(-I * PI * i / N) * parameter.out[i + parameter.step];
			parameter.buf[i / 2]     = parameter.out[i] + t;
			parameter.buf[(i + N)/2] = parameter.out[i] - t;
		}
	}


	//caz 4 core-uri
	if(P == 4){

		//creez parametrii corespunzatori apelurilor de functii de la nivel 2
		//numerotarea nivelurilor incepand de la 0
		Parameter p1, p2, p3, p4;
		p1.step = parameter.step * 4;
		p1.buf = parameter.buf;
		p1.out = parameter.out;

		p2.out = parameter.out + 1;
		p2.buf = parameter.buf + 1;
		p2.step = parameter.step * 4;

		p3.out = parameter.out + 2;
		p3.buf = parameter.buf + 2;
		p3.step = parameter.step * 4;

		p4.out = parameter.out + 3;
		p4.buf = parameter.buf + 3;
		p4.step = parameter.step * 4;

		pthread_create(&thread[0], NULL, fft, &p1);
		pthread_create(&thread[1], NULL, fft, &p2);
		pthread_create(&thread[2], NULL, fft, &p3);
		pthread_create(&thread[3], NULL, fft, &p4);
		

		for(int i = 0 ;i < P; i++){
			pthread_join(thread[i], NULL);
		}
		

		//execut cele doua for-uri corespunzatoare apelurilor de la nivel 1
		for (int i = 0; i < N; i += 2 * (parameter.step * 2)) {

			cplx t = cexp(-I * PI * i / N) * parameter.buf[i + (parameter.step * 2)];
			parameter.out[i / 2]     = parameter.buf[i] + t;
			parameter.out[(i + N)/2] = parameter.buf[i] - t;
		}

		for (int i = 0; i < N; i += 2 * (parameter.step * 2)) {

			cplx t = cexp(-I * PI * i / N) * 
						(parameter.buf + 1)[i + (parameter.step * 2)];
			(parameter.out + 1)[i / 2]     = (parameter.buf + 1)[i] + t;
			(parameter.out + 1)[(i + N)/2] = (parameter.buf + 1)[i] - t;
		}


		//execut for-ul corespunzator apelului de la nivel 0
		for (int i = 0; i < N; i += 2 * parameter.step) {
			cplx t = cexp(-I * PI * i / N) * parameter.out[i + parameter.step];
			parameter.buf[i / 2]     = parameter.out[i] + t;
			parameter.buf[(i + N)/2] = parameter.out[i] - t;
		}

	}



	for(int i = 0; i < N; i++){
		fprintf(outputFile, "%lf %lf\n", creal(Xi[i]), cimag(Xi[i]));
	}


	fclose(inputFile);
	fclose(outputFile);

	free(Xi);
	free(out);

	return 0;
}
