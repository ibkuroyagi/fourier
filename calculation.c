#include<stdio.h>
#include<math.h>
#include"calculation.h"

int dft_idft(double *re, double *im, int num, int flag)
{
	int i, j;
	double *temp_re, *temp_im;
	double coefficient;

	if((temp_re = (double*)malloc(sizeof(double)*num)) == NULL){
		fprintf(stderr, "Allocationerror!\n");
		return 1;
	}

	if((temp_im = (double*)malloc(sizeof(double)*num)) == NULL){
		fprintf(stderr, "Allocationerror!\n");
		free(temp_re);
		return 1;
	}

	for(i=0; i<num; i++){
		temp_re[i] = temp_im[i] = 0.0;
	}

	if(flag==IDFT)coefficient=num;
	else coefficient=1;

	for(i=0; i<num; i++){
		for(j=0; j<num; j++){
			temp_re[i] += re[j]*cos(2*PI*i*j/num) + flag*im[j]*sin(2*PI*i*j/num);
			temp_im[i] += -flag*re[j]*sin(2*PI*i*j/num) + im[j]*cos(2*PI*i*j/num);
		}
		temp_re[i] /= coefficient;
		temp_im[i] /= coefficient;
	}

	for(i=0; i<num; i++){
		re[i] = temp_re[i];
		im[i] = temp_im[i];
	}

	free(temp_re);
	free(temp_im);

	return 0;
}

int dft_idft2(double *re, double *im, int width, int height, int flag)
{
	double *temp_re;
	double *temp_im;
	int i, j;

	if((temp_re = (double*)malloc(sizeof(double)*width)) == NULL){
		fprintf(stderr, "Allocationerror!\n");
		return 1;
	}

	if((temp_im = (double*)malloc(sizeof(double)*width)) == NULL){
		fprintf(stderr, "Allocationerror!\n");
		free(temp_re);
		return 1;
	}

	for(i=0; i<height; i++){
		for(j=0; j<width; j++){
			temp_re[j] = re[i*width + j];
			temp_im[j] = im[i*width + j];
		}

		dft_idft(temp_re, temp_im, width, flag);

		for(j=0; j<width; j++){
			re[i*width + j] = temp_re[j];
			im[i*width + j] = temp_im[j];
		}
	}

	free(temp_re);
	free(temp_im);

	if((temp_re = (double*)malloc(sizeof(double)*height)) == NULL){
		fprintf(stderr, "Allocationerror!\n");
		return 1;
	}

	if((temp_im = (double*)malloc(sizeof(double)*height)) == NULL){
		fprintf(stderr, "Allocationerror!\n");
		free(temp_re);
		return 1;
	}

	for(j=0; j<width; j++){
		for(i=0; i<height; i++){
			temp_re[i] = re[i*width + j];
			temp_im[i] = im[i*width + j];
		}

		dft_idft(temp_re, temp_im, height, flag);

		for(i=0; i<height; i++){
			re[i*width + j] = temp_re[i];
			im[i*width + j] = temp_im[i];
		}
	}

	free(temp_re);
	free(temp_im);

	return 0;
}

void dft_swap(double *re, double *im, int num)
{
	int i;

	for(i=0; i<num/2; i++){
		swap(&re[i], &re[num/2 + i]);
		swap(&im[i], &im[num/2 + i]);
	}
}

void dft_swap2(double *re, double *im, int width, int height)
{
	int i, j;

	for(i=0; i<height/2; i++){
		for(j=0; j<width/2; j++){
			swap(&re[width*i + j], &re[width*(height/2+i) + width/2 + j]);
			swap(&im[width*i + j], &im[width*(height/2+i) + width/2 + j]);
			swap(&re[width*i + width/2 + j], &re[width*(height/2+i) + j]);
			swap(&im[width*i + width/2 + j], &im[width*(height/2+i) + j]);
		}
	}
}

void swap(double *a, double *b)
{
	double temp;

	temp = *a;
	*a = *b;
	*b = temp;
}


