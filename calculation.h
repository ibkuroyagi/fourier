#ifndef __CALCULATION_H_INCLUDED__
#define __CALCULATION_H_INCLUDED__

#define PI 3.1415926535897

//$B%U!<%j%(JQ49$+5U%U!<%j%(JQ49$+$N;XDj(B
#define DFT 1
#define IDFT -1

int dft_idft(double *re, double *im, int num, int flag);
int dft_idft2(double *re, double *im, int width, int height, int flag);
//dft$B$r$9$kA08e$G$N8r49$r$9$k(B
void dft_swap(double *re, double *im, int num);
void dft_swap2(double *re, double *im, int width, int height);
void swap(double *a, double *b);

#endif /* __CALCULATION_H_INCLUDED__ */
