#include<stdio.h>
#include<math.h>
#include"calculation.h"

#define PI 3.1415926535897

int samp[] = {0, 0, 0, 0,     0,   0,   0, 0, 0, 0,
							0, 0, 0, 0,     0,   0,   0, 0, 0, 0,
							0, 0, 0, 0,     0,   0,   0, 0, 0, 0,
							0, 0, 0, 100, 100, 100, 100, 0, 0, 0,
							0, 0, 0, 100, 100, 100, 100, 0, 0, 0,
							0, 0, 0, 100, 100, 100, 100, 0, 0, 0,
							0, 0, 0, 100, 100, 100, 100, 0, 0, 0,
							0, 0, 0, 0,     0,   0,   0, 0, 0, 0,
							0, 0, 0, 0,     0,   0,   0, 0, 0, 0,
							0, 0, 0, 0,     0,   0,   0, 0, 0, 0};

int main()
{
	double re[100], im[100];
	int i, j;

	for(i=0; i<100; i++){
		re[i] = samp[i];
		im[i] = 0;
	}

	for(i=0; i<10; i++){
		for(j=0; j<10; j++){
			printf("%.3d ", (int)re[i*10 + j]);
		}
		printf("\n");
	}

	dft_swap2(re, im, 10, 10);

	printf("\n");

	dft_idft2(re, im, 10, 10, DFT);

	dft_swap2(re, im, 10, 10);

	for(i=0; i<10; i++){
		for(j=0; j<10; j++){
			printf("%.3d ", (int)re[i*10 + j]);
		}
		printf("\n");
	}

	printf("\n");

	dft_swap2(re, im, 10, 10);

	dft_idft2(re, im, 10, 10, IDFT);

	dft_swap2(re, im, 10, 10);

	for(i=0; i<10; i++){
		for(j=0; j<10; j++){
			printf("%.3d ", (int)re[i*10 + j]);
		}
		printf("\n");
	}
		printf("\n");

	return 0;
}

