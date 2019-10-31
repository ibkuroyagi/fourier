#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#pragma warning(disable : 4996)
#include "ImageIO.h"
const double PI = 3.141592653589793;
#define SIZE 4096


//�O���[�X�P�[���摜�̓�l�����s��
//���̊֐��̈����́C
//�������FImage<GRAY>�^�̓��͉摜(&���t���Ă��Ȃ��Ƃ��́C�R�s�[���n����邽�߁C�㏑������Ȃ�)
//�������FImage<GRAY>�^�̏o�͉摜(&���t���Ă���Ƃ��́C�㏑�������)
//��O�����Fint �^�̕ϐ��@臒l
//GRA
void Binarization(Image<GRAY> src, Image<GRAY> &dst, int threshold)
{
        int i, j;
        for (i = 0; i < src.H; i++)
        {
                for (j = 0; j < src.W; j++)
                {
                        if (src.data[i][j] < threshold) //GRAY�^�̏ꍇ�� .r , .g , .b �͎g��Ȃ��D
                                dst.data[i][j] = 0;
                        else
                                dst.data[i][j] = 255;
                }
        }
}
//n��2�̗ݏ悩�ǂ������m�F����
//n��2�̗ݏ�łȂ���΂����2�̗ݏ�ɐ؂�グ�����l��Ԃ�
int isPowerOfTwo(int n) {
        int pow = 2;
        if ((n & (n - 1)) == 0){
                return 0;
        }else{
                while(n > pow){
                        pow *= 2;
                }
                return pow;
        }
}

//�ŏ���2�̗ݏ�ɂȂ�܂�padding����
int padding_arr(double *arr, int n1, int n2){
  for(int i=n1;i<n2;i++){
    arr[i] = 0.0;
  }
}

// �����t�[���G�ϊ�
// int n �F�f�[�^���i�Q�ׂ̂���j
// int flg �F���ϊ�:-1, �t�ϊ�:1
// double *ar �F�f�[�^�z��i�����G�v�f�� n)
// double *ai �F�f�[�^�z��i�����G�v�f�� n�j
int fft(int n, int flg, double *ar, double *ai)
{
        long m, mh, i, j, k, irev;
        double wr, wi, xr, xi;
        double theta;
        int Is2bit = isPowerOfTwo(n);

        if(Is2bit != 0){
                printf("If you want to use fft or ifft, resize the image.\n");
                padding_arr(ar,n,Is2bit);
                padding_arr(ai,n,Is2bit);
                n = Is2bit;
        }

        theta = flg * 2 * PI / n;

        i = 0;
        for (j = 1; j < n - 1; j++)
        {
                for (k = n >> 1; k > (i ^= k); k >>= 1)
                        ;
                if (j < i)
                {
                        xr = ar[j];
                        xi = ai[j];
                        ar[j] = ar[i];
                        ai[j] = ai[i];
                        ar[i] = xr;
                        ai[i] = xi;
                }
        }
        for (mh = 1; (m = mh << 1) <= n; mh = m)
        {
                irev = 0;
                for (i = 0; i < n; i += m)
                {
                        wr = cos(theta * irev);
                        wi = sin(theta * irev);
                        for (k = n >> 2; k > (irev ^= k); k >>= 1)
                                ;
                        for (j = i; j < mh + i; j++)
                        {
                                k = j + mh;
                                xr = ar[j] - ar[k];
                                xi = ai[j] - ai[k];
                                ar[j] += ar[k];
                                ai[j] += ai[k];
                                ar[k] = wr * xr - wi * xi;
                                ai[k] = wr * xi + wi * xr;
                        }
                }
        }

        if (flg == -1)
        {
                for (i = 0; i < n; i++)
                {
                        ar[i] /= n;
                        ai[i] /= n;
                }
        }
        return 0;
}

// �t�[���G�ϊ�
// int N �F�f�[�^��
// double *re �F���̓f�[�^�z��i�����G�v�f�� N)
// double *im �F���̓f�[�^�z��i�����G�v�f�� N�j
// double *Re �F�o�̓f�[�^�z��i�����G�v�f�� N)
// double *Im �F�o�̓f�[�^�z��i�����G�v�f�� N�j
double dft(double re[], double im[], double Re[], double Im[], int N)
{
        //���������Ƌ��������ɕ����ăt�[���G�ϊ�
        for (int k = 0; k < N; k++)
        {
                for (int n = 0; n < N; n++)
                {
                        Re[k] += (re[n] * cos(2 * PI * k * n / N) + im[n] * sin(2 * PI * k * n / N));
                        Im[k] += (-1 * re[n] * sin(2 * PI * k * n / N) + im[n] * cos(2 * PI * k * n / N));
                        // Im[k] = 0.0;
                }
        }
}

// �t�t�[���G�ϊ�
// int N �F�f�[�^��
// double *Re �F���̓f�[�^�z��i�����G�v�f�� N)
// double *Im �F���̓f�[�^�z��i�����G�v�f�� N�j
// double *re �F�o�̓f�[�^�z��i�����G�v�f�� N)
// double *im �F�o�̓f�[�^�z��i�����G�v�f�� N�j
double idft(double Re[], double Im[], double re[], double im[], int N)
{
        //���������Ƌ��������ɕ����ċt�t�[���G�ϊ�
        for (int k = 0; k < N; k++)
        {
                for (int n = 0; n < N; n++)
                {
                        re[k] += (Re[n] * cos(2 * PI * k * n / N) - Im[n] * sin(2 * PI * k * n / N));
                        im[k] += (Re[n] * sin(2 * PI * k * n / N) + Im[n] * cos(2 * PI * k * n / N));
                        // im[k] = 0.0;
                }
                re[k] /= N;
                im[k] /= N;
        }
}

// void two_dimension_fourier(Image<GRAY> src, Image<GRAY> &dst, char function[], double re[256][256], double im[256][256])
void two_dimension_fourier(Image<GRAY> src, Image<GRAY> &dst, char function[], Image<GRAY> &dst_im)
{
        // printf("src.data[0][255] = %d\n",src.data[0][255]);
        double re_arr[SIZE], im_arr[SIZE], Re_tmp[SIZE], Im_tmp[SIZE];
        double Re_arr[256][256], Im_arr[256][256]; //������΂��̂Œ���
        
        for (int i = 0; i < src.H; i++)
        {
                for (int j = 0; j < src.W; j++)
                {
                        re_arr[j] = (double)src.data[i][j];
                        im_arr[j] = 0.0;
                }
                if (strcmp(function, "dft") == 0)
                {
                        dft(re_arr, im_arr, Re_tmp, Im_tmp, src.W);
                        // printf("i = %d, Re_tmp[0] = %lf\n",i,Re_tmp[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr[i][j] = Re_tmp[j];
                                Im_arr[i][j] = Im_tmp[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.W, -1, re_arr, im_arr);
                        // printf("i = %d, re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr[i][j] = re_arr[j];
                                Im_arr[i][j] = im_arr[j];
                        }
                }
        }
        for (int i = 0; i < src.W; i++)
        {
                for (int j = 0; j < src.H; j++)
                {
                        re_arr[j] = Re_arr[j][i];
                        im_arr[j] = Im_arr[j][i];
                }
                if (strcmp(function, "dft") == 0)
                {
                        dft(re_arr, im_arr, Re_tmp, Im_tmp, src.H);
                        for (int j = 0; j < src.H; j++)
                        {
                                dst.data[j][i] = (int)Re_tmp[j];
                                dst_im.data[j][i] = (int)Im_tmp[j];
                                // re[j][i] = re_arr[j];
                                // im[j][i] = im_arr[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.H, -1, re_arr, im_arr);
                        // printf("i = %d, cross_re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                dst.data[j][i] = (int)re_arr[j];
                                dst_im.data[j][i] = (int)im_arr[j];
                                // re[j][i] = re_arr[j];
                                // im[j][i] = im_arr[j];
                        }
                }
        }
}

// void inverse_two_dimension_fourier(Image<GRAY> src, Image<GRAY> &dst, char function[], double re[][256], double im[][256])
void inverse_two_dimension_fourier(Image<GRAY> src, Image<GRAY> &dst, char function[], Image<GRAY> src_im)
{
        double re_arr[SIZE], im_arr[SIZE], Re_tmp[SIZE], Im_tmp[SIZE];
        double Re_arr[256][256], Im_arr[256][256];
        for (int i = 0; i < src.H; i++)
        {
                for (int j = 0; j < src.W; j++)
                {
                        re_arr[j] = (double)src.data[j][i];
                        im_arr[j] = (double)src_im.data[j][i];
                        // re_arr[j] = re[j][i];
                        // im_arr[j] = im[j][i];
                }
                if (strcmp(function, "idft") == 0)
                {
                        idft(re_arr, im_arr, Re_tmp, Im_tmp, src.W);
                        // printf("i = %d, Re_tmp[0] = %lf\n",i,Re_tmp[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr[j][i] = Re_tmp[j];
                                Im_arr[j][i] = Im_tmp[j];
                        }
                }
                else if (strcmp(function, "ifft") == 0)
                {
                        fft(src.W, 1, re_arr, im_arr);
                        // printf("i = %d, re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr[j][i] = re_arr[j];
                                Im_arr[j][i] = im_arr[j];
                        }
                }
        }
        for (int i = 0; i < src.W; i++)
        {
                for (int j = 0; j < src.H; j++)
                {
                        re_arr[j] = Re_arr[i][j];
                        im_arr[j] = Im_arr[i][j];
                }
                if (strcmp(function, "idft") == 0)
                {
                        idft(re_arr, im_arr, Re_tmp, Im_tmp, src.H);
                        for (int j = 0; j < src.H; j++)
                        {
                                dst.data[i][j] = (int)Re_tmp[j];
                                // re[i][j] = re_arr[j];
                                // im[i][j] = im_arr[j];
                        }
                }
                else if (strcmp(function, "ifft") == 0)
                {
                        fft(src.H, 1, re_arr, im_arr);
                        // printf("i = %d, re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                dst.data[i][j] = (int)re_arr[j];
                                // re[i][j] = re_arr[j];
                                // im[i][j] = im_arr[j];
                        }
                }
        }
}



void two_dimension_fourier_color(Image<COLOR> src, Image<COLOR> &dst, char function[])
{       
        // printf("hello\n");
        // printf("src.data[0][255].r = %d\n",src.data[0][255].r);
        // double re_arr_r[SIZE], im_arr_r[SIZE], Re_tmp_r[SIZE], Im_tmp_r[SIZE];
        // double re_arr_b[SIZE], im_arr_b[SIZE], Re_tmp_b[SIZE], Im_tmp_b[SIZE];
        // double re_arr_g[SIZE], im_arr_g[SIZE], Re_tmp_g[SIZE], Im_tmp_g[SIZE];
        double re_arr_r[256], im_arr_r[256], Re_tmp_r[256], Im_tmp_r[256];
        double re_arr_b[256], im_arr_b[256], Re_tmp_b[256], Im_tmp_b[256];
        double re_arr_g[256], im_arr_g[256], Re_tmp_g[256], Im_tmp_g[256];
        double Re_arr_r[256][256], Im_arr_r[256][256]; 
        double Re_arr_b[256][256], Im_arr_b[256][256];
        double Re_arr_g[256][256], Im_arr_g[256][256];
        // printf("src.data[0][255].r = %d\n",src.data[0][255].r);

        for (int i = 0; i < src.H; i++)
        {
                for (int j = 0; j < src.W; j++)
                {
                        re_arr_r[j] = (double)src.data[i][j].r;
                        im_arr_r[j] = 0.0;
                        re_arr_b[j] = (double)src.data[i][j].b;
                        im_arr_b[j] = 0.0;
                        re_arr_g[j] = (double)src.data[i][j].g;
                        im_arr_g[j] = 0.0;
                }
                if (strcmp(function, "dft") == 0)
                {
                        dft(re_arr_r, im_arr_r, Re_tmp_r, Im_tmp_r, src.W);
                        dft(re_arr_b, im_arr_b, Re_tmp_b, Im_tmp_b, src.W);
                        dft(re_arr_g, im_arr_g, Re_tmp_g, Im_tmp_g, src.W);
                        // printf("i = %d, Re_tmp_r[0] = %lf\n",i,Re_tmp_r[10]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr_r[i][j] = Re_tmp_r[j];
                                Im_arr_r[i][j] = Im_tmp_r[j];
                                Re_arr_b[i][j] = Re_tmp_b[j];
                                Im_arr_b[i][j] = Im_tmp_b[j];
                                Re_arr_g[i][j] = Re_tmp_g[j];
                                Im_arr_g[i][j] = Im_tmp_g[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.W, -1, re_arr_r, im_arr_r);
                        fft(src.W, -1, re_arr_b, im_arr_b);
                        fft(src.W, -1, re_arr_g, im_arr_g);
                        // printf("i = %d, re_arr_r[0] = %lf\n",i,re_arr_r[10]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr_r[i][j] = re_arr_r[j];
                                Im_arr_r[i][j] = im_arr_r[j];
                                Re_arr_b[i][j] = re_arr_b[j];
                                Im_arr_b[i][j] = im_arr_b[j];
                                Re_arr_g[i][j] = re_arr_g[j];
                                Im_arr_g[i][j] = im_arr_g[j];
                        }
                }
        }
        for (int i = 0; i < src.W; i++)
        {
                for (int j = 0; j < src.H; j++)
                {
                        re_arr_r[j] = Re_arr_r[j][i];
                        im_arr_r[j] = Im_arr_r[j][i];
                        re_arr_b[j] = Re_arr_b[j][i];
                        im_arr_b[j] = Im_arr_b[j][i];
                        re_arr_g[j] = Re_arr_g[j][i];
                        im_arr_g[j] = Im_arr_g[j][i];
                }
                if (strcmp(function, "dft") == 0)
                {
                        dft(re_arr_r, im_arr_r, Re_tmp_r, Im_tmp_r, src.H);
                        dft(re_arr_b, im_arr_b, Re_tmp_b, Im_tmp_b, src.H);
                        dft(re_arr_g, im_arr_g, Re_tmp_g, Im_tmp_g, src.H);
                        for (int j = 0; j < src.H; j++)
                        {
                                dst.data[j][i].r = (int)Re_tmp_r[j];
                                dst.data[j][i].b = (int)Re_tmp_b[j];
                                dst.data[j][i].g = (int)Re_tmp_g[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.H, -1, re_arr_r, im_arr_r);
                        fft(src.H, -1, re_arr_b, im_arr_b);
                        fft(src.H, -1, re_arr_g, im_arr_g);
                        // printf("i = %d, cross_re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                dst.data[j][i].r = (int)re_arr_r[j];
                                dst.data[j][i].b = (int)re_arr_b[j];
                                dst.data[j][i].g = (int)re_arr_g[j];
                        }
                }
        }
}

void inverse_two_dimension_fourier_color(Image<COLOR> src, Image<COLOR> &dst, char function[])
{
        double re_arr_r[SIZE], im_arr_r[SIZE], Re_tmp_r[SIZE], Im_tmp_r[SIZE];
        double re_arr_b[SIZE], im_arr_b[SIZE], Re_tmp_b[SIZE], Im_tmp_b[SIZE];
        double re_arr_g[SIZE], im_arr_g[SIZE], Re_tmp_g[SIZE], Im_tmp_g[SIZE];
        double Re_arr_r[256][256], Im_arr_r[256][256]; 
        double Re_arr_b[256][256], Im_arr_b[256][256];
        double Re_arr_g[256][256], Im_arr_g[256][256];

        for (int i = 0; i < src.H; i++)
        {
                for (int j = 0; j < src.W; j++)
                {
                        re_arr_r[j] = (double)src.data[j][i].r;
                        im_arr_r[j] = 0.0;
                        re_arr_b[j] = (double)src.data[j][i].b;
                        im_arr_b[j] = 0.0;
                        re_arr_g[j] = (double)src.data[j][i].g;
                        im_arr_g[j] = 0.0;
                }
                if (strcmp(function, "idft") == 0)
                {
                        idft(re_arr_r, im_arr_r, Re_tmp_r, Im_tmp_r, src.W);
                        idft(re_arr_b, im_arr_b, Re_tmp_b, Im_tmp_b, src.W);
                        idft(re_arr_g, im_arr_g, Re_tmp_g, Im_tmp_g, src.W);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr_r[j][i] = Re_tmp_r[j];
                                Im_arr_r[j][i] = Im_tmp_r[j];
                                Re_arr_b[j][i] = Re_tmp_b[j];
                                Im_arr_b[j][i] = Im_tmp_b[j];
                                Re_arr_g[j][i] = Re_tmp_g[j];
                                Im_arr_g[j][i] = Im_tmp_g[j];
                        }
                }
                else if (strcmp(function, "ifft") == 0)
                {
                        fft(src.W, 1, re_arr_r, im_arr_r);
                        fft(src.W, 1, re_arr_b, im_arr_b);
                        fft(src.W, 1, re_arr_g, im_arr_g);
                        // printf("i = %d, re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr_r[j][i] = re_arr_r[j];
                                Im_arr_r[j][i] = im_arr_r[j];
                                Re_arr_b[j][i] = re_arr_b[j];
                                Im_arr_b[j][i] = im_arr_b[j];
                                Re_arr_g[j][i] = re_arr_g[j];
                                Im_arr_g[j][i] = im_arr_g[j];
                        }
                }
        }
        for (int i = 0; i < src.W; i++)
        {
                for (int j = 0; j < src.H; j++)
                {
                        re_arr_r[j] = Re_arr_r[i][j];
                        im_arr_r[j] = Im_arr_r[i][j];
                        re_arr_b[j] = Re_arr_b[i][j];
                        im_arr_b[j] = Im_arr_b[i][j];
                        re_arr_g[j] = Re_arr_g[i][j];
                        im_arr_g[j] = Im_arr_g[i][j];
                }
                if (strcmp(function, "idft") == 0)
                {
                        idft(re_arr_r, im_arr_r, Re_tmp_r, Im_tmp_r, src.H);
                        idft(re_arr_b, im_arr_b, Re_tmp_b, Im_tmp_b, src.H);
                        idft(re_arr_g, im_arr_g, Re_tmp_g, Im_tmp_g, src.H);
                        for (int j = 0; j < src.H; j++)
                        {
                                dst.data[i][j].r = (int)Re_tmp_r[j];
                                dst.data[i][j].b = (int)Re_tmp_b[j];
                                dst.data[i][j].g = (int)Re_tmp_g[j];
                        }
                }
                else if (strcmp(function, "ifft") == 0)
                {
                        fft(src.H, 1, re_arr_r, im_arr_r);
                        fft(src.H, 1, re_arr_b, im_arr_b);
                        fft(src.H, 1, re_arr_g, im_arr_g);
                        // printf("i = %d, re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                dst.data[i][j].r = (int)re_arr_r[j];
                                dst.data[i][j].b = (int)re_arr_b[j];
                                dst.data[i][j].g = (int)re_arr_g[j];
                        }
                }
        }
}


int file_load_and_write(char filename1[], char filename2[], char function[])
{
        FILE *fp1, *fp2; // FILE�^�\����
        double f1, f2 = 0.0;
        double Re[SIZE], Im[SIZE]; //�ϊ��O�̔z��
        double re[SIZE], im[SIZE]; //�ϊ���̔z��
        int i, k, i1, N = 0;
        // printf("fname1 = %s, fname2 = %s, function = %s\n",filename1,filename2,function);
        fp1 = fopen(filename1, "r"); // �t�@�C�����J���B���s�����NULL��Ԃ��B
        if (fp1 == NULL)
        {
                printf("%s file not open!\n", filename1);
                exit(1);
        }

        while (fscanf(fp1, "%d,%lf,%lf", &i1, &f1, &f2) != EOF)
        {
                // printf("i1 = %d,f1 = %f, f2 = %f\n", i1, f1, f2);
                Re[N] = f1;
                Im[N] = f2;
                N++;
        }
        // printf("Re[%d] = %f, Im[%d] = %f\n", N - 1, Re[N - 1], N - 1, Im[N - 1]);
        fclose(fp1); // �t�@�C�������

        if ((fp2 = fopen(filename2, "w")) == NULL)
        {
                printf("FILE2 not open\n");
                return -1;
        }
        //���������Ƌ��������ɕ����ăt�[���G�ϊ��܂��͋t�t�[���G�ϊ�
        if (strcmp(function, "dft") == 0)
        {
                dft(Re, Im, re, im, N);
        }
        else if (strcmp(function, "idft") == 0)
        {
                idft(Re, Im, re, im, N);
        }
        else if (strcmp(function, "fft") == 0)
        {
                fft(N, -1, Re, Im);
                for (int i = 0; i < N; i++)
                {
                        re[i] = Re[i];
                        im[i] = Im[i];
                }
        }
        else if (strcmp(function, "ifft") == 0)
        {
                fft(N, 1, Re, Im);
                for (int i = 0; i < N; i++)
                {
                        re[i] = Re[i];
                        im[i] = Im[i];
                }
        }
        // printf("re[%d] = %f, im[%d] = %f\n", N - 1, re[N - 1], N - 1, im[N - 1]);
        for (k = 0; k < N; k++)
        {
                fprintf(fp2, "%d, %f, %f\n", k, re[k], im[k]);
        }
        fclose(fp2);
        return 0;
}

int fft_and_ifft(Image<GRAY> src, Image<GRAY> &dst, char function[]){
        // printf("src.data[0][255] = %d\n",src.data[0][255]);
        double re_arr[SIZE], im_arr[SIZE], Re_tmp[SIZE], Im_tmp[SIZE];
        double Re_arr[256][256], Im_arr[256][256]; //������΂��̂Œ���
        double Re_arr_tmp[256][256], Im_arr_tmp[256][256]; //������΂��̂Œ���
        
        for (int i = 0; i < src.H; i++)
        {
                for (int j = 0; j < src.W; j++)
                {
                        re_arr[j] = (double)src.data[i][j];
                        im_arr[j] = 0.0;
                }
                if (strcmp(function, "dft") == 0)
                {
                        dft(re_arr, im_arr, Re_tmp, Im_tmp, src.W);
                        // printf("i = %d, Re_tmp[0] = %lf\n",i,Re_tmp[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr[i][j] = Re_tmp[j];
                                Im_arr[i][j] = Im_tmp[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.W, -1, re_arr, im_arr);
                        // printf("i = %d, re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr[i][j] = re_arr[j];
                                Im_arr[i][j] = im_arr[j];
                        }
                }
        }
        for (int i = 0; i < src.W; i++)
        {
                for (int j = 0; j < src.H; j++)
                {
                        re_arr[j] = Re_arr[j][i];
                        im_arr[j] = Im_arr[j][i];
                }
                if (strcmp(function, "dft") == 0)
                {
                        dft(re_arr, im_arr, Re_tmp, Im_tmp, src.H);
                        for (int j = 0; j < src.H; j++)
                        {
                                Re_arr_tmp[j][i] = Re_tmp[j];
                                Im_arr_tmp[j][i] = Im_tmp[j];
                                // dst.data[j][i] = (int)Re_tmp[j];
                                // dst_im.data[j][i] = (int)Im_tmp[j];
                                // re[j][i] = re_arr[j];
                                // im[j][i] = im_arr[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.H, -1, re_arr, im_arr);
                        // printf("i = %d, cross_re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr_tmp[j][i] = re_arr[j];
                                Im_arr_tmp[j][i] = im_arr[j];
                                // dst.data[j][i] = (int)re_arr[j];
                                // dst_im.data[j][i] = (int)im_arr[j];
                                // re[j][i] = re_arr[j];
                                // im[j][i] = im_arr[j];
                        }
                }
        }


        // double re_arr[SIZE], im_arr[SIZE], Re_tmp[SIZE], Im_tmp[SIZE];
        // double Re_arr[256][256], Im_arr[256][256];
        for (int i = 0; i < src.H; i++)
        {
                for (int j = 0; j < src.W; j++)
                {
                        re_arr[j] = Re_arr_tmp[j][i];
                        im_arr[j] = Im_arr_tmp[j][i];
                        // re_arr[j] = (double)src.data[j][i];
                        // im_arr[j] = (double)src_im.data[j][i];
                        // re_arr[j] = re[j][i];
                        // im_arr[j] = im[j][i];
                }
                if (strcmp(function, "dft") == 0)
                {
                        idft(re_arr, im_arr, Re_tmp, Im_tmp, src.W);
                        // printf("i = %d, Re_tmp[0] = %lf\n",i,Re_tmp[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr[j][i] = Re_tmp[j];
                                Im_arr[j][i] = Im_tmp[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.W, 1, re_arr, im_arr);
                        // printf("i = %d, re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr[j][i] = re_arr[j];
                                Im_arr[j][i] = im_arr[j];
                        }
                }
        }
        for (int i = 0; i < src.W; i++)
        {
                for (int j = 0; j < src.H; j++)
                {
                        re_arr[j] = Re_arr[i][j];
                        im_arr[j] = Im_arr[i][j];
                }
                if (strcmp(function, "dft") == 0)
                {
                        idft(re_arr, im_arr, Re_tmp, Im_tmp, src.H);
                        for (int j = 0; j < src.H; j++)
                        {
                                dst.data[i][j] = (int)Re_tmp[j];
                                // re[i][j] = re_arr[j];
                                // im[i][j] = im_arr[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.H, 1, re_arr, im_arr);
                        // printf("i = %d, re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                dst.data[i][j] = (int)re_arr[j];
                                // re[i][j] = re_arr[j];
                                // im[i][j] = im_arr[j];
                        }
                }
        }
}

int fft_and_ifft_color(Image<COLOR> src, Image<COLOR> &dst, char function[]){
        double re_arr_r[256], im_arr_r[256], Re_tmp_r[256], Im_tmp_r[256];
        double re_arr_b[256], im_arr_b[256], Re_tmp_b[256], Im_tmp_b[256];
        double re_arr_g[256], im_arr_g[256], Re_tmp_g[256], Im_tmp_g[256];
        double Re_arr_r[256][256], Im_arr_r[256][256]; 
        double Re_arr_b[256][256], Im_arr_b[256][256];
        double Re_arr_g[256][256], Im_arr_g[256][256];
        double Re_arr_r_tmp[256][256], Im_arr_r_tmp[256][256]; 
        double Re_arr_b_tmp[256][256], Im_arr_b_tmp[256][256];
        double Re_arr_g_tmp[256][256], Im_arr_g_tmp[256][256];
        // printf("src.data[0][255].r = %d\n",src.data[0][255].r);

        for (int i = 0; i < src.H; i++)
        {
                for (int j = 0; j < src.W; j++)
                {
                        re_arr_r[j] = (double)src.data[i][j].r;
                        im_arr_r[j] = 0.0;
                        re_arr_b[j] = (double)src.data[i][j].b;
                        im_arr_b[j] = 0.0;
                        re_arr_g[j] = (double)src.data[i][j].g;
                        im_arr_g[j] = 0.0;
                }
                if (strcmp(function, "dft") == 0)
                {
                        dft(re_arr_r, im_arr_r, Re_tmp_r, Im_tmp_r, src.W);
                        dft(re_arr_b, im_arr_b, Re_tmp_b, Im_tmp_b, src.W);
                        dft(re_arr_g, im_arr_g, Re_tmp_g, Im_tmp_g, src.W);
                        // printf("i = %d, Re_tmp_r[0] = %lf\n",i,Re_tmp_r[10]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr_r[i][j] = Re_tmp_r[j];
                                Im_arr_r[i][j] = Im_tmp_r[j];
                                Re_arr_b[i][j] = Re_tmp_b[j];
                                Im_arr_b[i][j] = Im_tmp_b[j];
                                Re_arr_g[i][j] = Re_tmp_g[j];
                                Im_arr_g[i][j] = Im_tmp_g[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.W, -1, re_arr_r, im_arr_r);
                        fft(src.W, -1, re_arr_b, im_arr_b);
                        fft(src.W, -1, re_arr_g, im_arr_g);
                        // printf("i = %d, re_arr_r[0] = %lf\n",i,re_arr_r[10]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr_r[i][j] = re_arr_r[j];
                                Im_arr_r[i][j] = im_arr_r[j];
                                Re_arr_b[i][j] = re_arr_b[j];
                                Im_arr_b[i][j] = im_arr_b[j];
                                Re_arr_g[i][j] = re_arr_g[j];
                                Im_arr_g[i][j] = im_arr_g[j];
                        }
                }
        }
        for (int i = 0; i < src.W; i++)
        {
                for (int j = 0; j < src.H; j++)
                {
                        re_arr_r[j] = Re_arr_r[j][i];
                        im_arr_r[j] = Im_arr_r[j][i];
                        re_arr_b[j] = Re_arr_b[j][i];
                        im_arr_b[j] = Im_arr_b[j][i];
                        re_arr_g[j] = Re_arr_g[j][i];
                        im_arr_g[j] = Im_arr_g[j][i];
                }
                if (strcmp(function, "dft") == 0)
                {
                        dft(re_arr_r, im_arr_r, Re_tmp_r, Im_tmp_r, src.H);
                        dft(re_arr_b, im_arr_b, Re_tmp_b, Im_tmp_b, src.H);
                        dft(re_arr_g, im_arr_g, Re_tmp_g, Im_tmp_g, src.H);
                        for (int j = 0; j < src.H; j++)
                        {
                                // dst.data[j][i].r = (int)Re_tmp_r[j];
                                // dst.data[j][i].b = (int)Re_tmp_b[j];
                                // dst.data[j][i].g = (int)Re_tmp_g[j];
                                Re_arr_r_tmp[j][i] = Re_tmp_r[j];
                                Re_arr_b_tmp[j][i] = Re_tmp_b[j];
                                Re_arr_g_tmp[j][i] = Re_tmp_g[j];
                                Im_arr_r_tmp[j][i] = Im_tmp_r[j];
                                Im_arr_b_tmp[j][i] = Im_tmp_b[j];
                                Im_arr_g_tmp[j][i] = Im_tmp_g[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.H, -1, re_arr_r, im_arr_r);
                        fft(src.H, -1, re_arr_b, im_arr_b);
                        fft(src.H, -1, re_arr_g, im_arr_g);
                        // printf("i = %d, cross_re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                // dst.data[j][i].r = (int)re_arr_r[j];
                                // dst.data[j][i].b = (int)re_arr_b[j];
                                // dst.data[j][i].g = (int)re_arr_g[j];
                                Re_arr_r_tmp[j][i] = re_arr_r[j];
                                Re_arr_b_tmp[j][i] = re_arr_b[j];
                                Re_arr_g_tmp[j][i] = re_arr_g[j];
                                Im_arr_r_tmp[j][i] = im_arr_r[j];
                                Im_arr_b_tmp[j][i] = im_arr_b[j];
                                Im_arr_g_tmp[j][i] = im_arr_g[j];
                        }
                }
        }




        ///////////
        for (int i = 0; i < src.H; i++)
        {
                for (int j = 0; j < src.W; j++)
                {
                        re_arr_r[j] = Re_arr_r_tmp[j][i];
                        im_arr_r[j] = Im_arr_r_tmp[j][i];
                        re_arr_b[j] = Re_arr_b_tmp[j][i];
                        im_arr_b[j] = Im_arr_b_tmp[j][i];
                        re_arr_g[j] = Re_arr_g_tmp[j][i];
                        im_arr_g[j] = Im_arr_g_tmp[j][i];
                }
                if (strcmp(function, "dft") == 0)
                {
                        idft(re_arr_r, im_arr_r, Re_tmp_r, Im_tmp_r, src.W);
                        idft(re_arr_b, im_arr_b, Re_tmp_b, Im_tmp_b, src.W);
                        idft(re_arr_g, im_arr_g, Re_tmp_g, Im_tmp_g, src.W);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr_r[j][i] = Re_tmp_r[j];
                                Im_arr_r[j][i] = Im_tmp_r[j];
                                Re_arr_b[j][i] = Re_tmp_b[j];
                                Im_arr_b[j][i] = Im_tmp_b[j];
                                Re_arr_g[j][i] = Re_tmp_g[j];
                                Im_arr_g[j][i] = Im_tmp_g[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.W, 1, re_arr_r, im_arr_r);
                        fft(src.W, 1, re_arr_b, im_arr_b);
                        fft(src.W, 1, re_arr_g, im_arr_g);
                        // printf("i = %d, re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                Re_arr_r[j][i] = re_arr_r[j];
                                Im_arr_r[j][i] = im_arr_r[j];
                                Re_arr_b[j][i] = re_arr_b[j];
                                Im_arr_b[j][i] = im_arr_b[j];
                                Re_arr_g[j][i] = re_arr_g[j];
                                Im_arr_g[j][i] = im_arr_g[j];
                        }
                }
        }
        for (int i = 0; i < src.W; i++)
        {
                for (int j = 0; j < src.H; j++)
                {
                        re_arr_r[j] = Re_arr_r[i][j];
                        im_arr_r[j] = Im_arr_r[i][j];
                        re_arr_b[j] = Re_arr_b[i][j];
                        im_arr_b[j] = Im_arr_b[i][j];
                        re_arr_g[j] = Re_arr_g[i][j];
                        im_arr_g[j] = Im_arr_g[i][j];
                }
                if (strcmp(function, "dft") == 0)
                {
                        idft(re_arr_r, im_arr_r, Re_tmp_r, Im_tmp_r, src.H);
                        idft(re_arr_b, im_arr_b, Re_tmp_b, Im_tmp_b, src.H);
                        idft(re_arr_g, im_arr_g, Re_tmp_g, Im_tmp_g, src.H);
                        for (int j = 0; j < src.H; j++)
                        {
                                dst.data[i][j].r = (int)Re_tmp_r[j];
                                dst.data[i][j].b = (int)Re_tmp_b[j];
                                dst.data[i][j].g = (int)Re_tmp_g[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.H, 1, re_arr_r, im_arr_r);
                        fft(src.H, 1, re_arr_b, im_arr_b);
                        fft(src.H, 1, re_arr_g, im_arr_g);
                        // printf("i = %d, re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                dst.data[i][j].r = (int)re_arr_r[j];
                                dst.data[i][j].b = (int)re_arr_b[j];
                                dst.data[i][j].g = (int)re_arr_g[j];
                        }
                }
        }

}
int kernel_filter(Image<GRAY> src, Image<GRAY> &dst)
{

        double kernel[3][3] = {
            {1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0},
            {1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0},
            {1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0}}; //�ړ����σt�B���^

        int k_size = 1;
        // kernel��K���ł��Ȃ����͂��̂܂܂̒l��������
        for (int i = 0; i < src.H - 1; i++)
        { //��A�E�A���A��
                dst.data[0][i] = src.data[0][i];
                dst.data[i][src.W - 1] = src.data[i][src.W - 1];
                dst.data[i + 1][0] = src.data[i + 1][0];
                dst.data[src.H - 1][i + 1] = src.data[src.H - 1][i + 1];
        }
        //kernel�Ɖ摜�̘a��ۑ�����z��
        double sum_result;
        for (int i = 1; i < src.H - 1; i++)
        {
                for (int j = 1; j < src.W - 1; j++)
                {
                        // src.data[i][j]���w�肵���Ƃ��̏������ȉ��ɋL�q����
                        // src.data[i][j]�̎��͂��w�肵�Akernel�Ƃ̐ς�sum_result�ɕۑ�����
                        // sum_result��dst.data[i][j]�ɕۑ�����
                        for (int p = 0; p < k_size * 2 + 1; p++)
                        {
                                for (int q = 0; q < k_size * 2 + 1; q++)
                                {
                                        sum_result += kernel[p][q] * (double)src.data[i - 1 + p][j - 1 + q];
                                }
                        }
                        dst.data[i][j] = (int)sum_result;
                        // printf("dst.data[%d][%d] = %d\n",i,j,dst.data[i][j]);
                }
        }
        return 0;
}


int main(void)
{
        printf("ready...\n");
        // 1�������U�t�[���G�ϊ�
        // file_load_and_write("test.csv", "inv_fourier.csv","dft");
        // printf("1dim dft is completed...\n");

        // 1�����������U�t�[���G�ϊ�
        // =3*SIN(0.01*RC[-1])+2*SIN(0.1*RC[-1])+SIN(RC[-1])
        file_load_and_write("test_fft1.csv", "inv_fourier_fft.csv", "fft");
        printf("1dim fft is completed...\n");

        // 1�����t���U�t�[���G�ϊ�
        // file_load_and_write("inv_fourier.csv", "check.csv","idft");
        // printf("1dim idft is completed...\n");

        // 1�����t�������U�t�[���G�ϊ�
        file_load_and_write("inv_fourier_fft.csv", "check_fft.csv", "ifft");
        printf("1dim ifft is completed...\n");
        // excel�ȂǂŃO���t�������Ċm�F�����Ă݂Ă��������B


        // 2�������U�t�[���G�ϊ�
        //�O���[�摜
        //���z��̐錾
        Image<GRAY> gray, gout_dft, gout_idft,gout_dft_im,gout_mean;
        // �摜�̃p�X���e���̊��ɕύX�����Ă��������B
        // char path2[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\pictures\\lenna.pgm";
        // char path2[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\pictures\\mandrill.pgm";
        // char path2[] = "./pictures/lenna.pgm";
        char path2[] = "./pictures/mandrill.pgm";
        //���摜���[�h
        gray.load(path2); //lenna.pgm��gray�Ƀ��[�h����
        //�摜�̈ꕔ���o��
        // for(int i=0;i<16;i++){
        //         for(int j=0;j<16;j++){
        //         printf("%4d, ",gray.data[i][j]);
        //         }
        //         printf("\n");
        // }
        //���o�͔z��m��
        gout_dft.create(gray.W, gray.H); // gout�� img.W(�摜img�̉���) x img.H(�摜img�̏c��) �̑傫���̉�f�z���p�ӂ���
        gout_dft_im.create(gray.W, gray.H);
        gout_idft.create(gray.W, gray.H);
        gout_mean.create(gray.W,gray.H);
        //���摜����//////////////////////////////////////////////////////////////
        // Binarization(gray, gout, 128); //gray��臒l128�œ�l������gout�ɏo��

        // 2�������U�t�[���G�ϊ�
        // two_dimension_fourier(gray,gout_dft,"dft");
        // 2�����������U�t�[���G�ϊ�
        // two_dimension_fourier(gray, gout_dft, "fft", gout_dft_im);
        //�摜�̈ꕔ���o��
        // for(int i=0;i<16;i++){
        //         for(int j=0;j<16;j++){
        //         printf("%4d, ",gout_dft.data[i][j]);
        //         }
        //         printf("\n");
        // }
        // char gout_dft_path[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\results\\two_dft_lenna.pgm";
        // char gout_dft_path[] = "./results/two_dft_lenna.pgm";
        
        // gout_dft.save(gout_dft_path);
        // 2�����t���U�t�[���G�ϊ�
        // inverse_two_dimension_fourier(gout_dft,gout_idft,"idft");

        // 2�����t�������U�t�[���G�ϊ�
        // inverse_two_dimension_fourier(gout_dft, gout_idft, "ifft",gout_dft_im);
        //�摜�̈ꕔ���o��
        // for(int i=0;i<16;i++){
        //         for(int j=0;j<16;j++){
        //         // printf("%4d, ",gray.data[i][j] - gout_idft.data[i][j]);
        //         printf("%4d, ",gout_idft.data[i][j]);
        //         }
        //         printf("\n");
        // }

        // char gout_idft_path[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\results\\two_idft_mandrill.pgm";
        // char gout_idft_path[] = "./results/two_idft_lenna.pgm";
        char gout_idft_path[] = "./results/two_idft_mandrill.pgm";
        fft_and_ifft(gray,gout_idft,"fft");
        gout_idft.save(gout_idft_path);
        // char gout_mean_path[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\results\\mandrill_mean.pgm";
        char gout_mean_path[] = "./results/mandrill_mean.pgm";
        kernel_filter(gray,gout_mean);
        gout_mean.save(gout_mean_path);


        

        Image<COLOR> img,img_dft,img_idft;	//�J���[�摜
        // char path1[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\pictures\\lenna.ppm";
        char path1[] = "./pictures/lenna.ppm";
        img.load( path1 );	//lenna.ppm��img�Ƀ��[�h����
        img_dft.create( img.W, img.H ); 
        img_idft.create( img.W, img.H );

        // for(int i=0;i<16;i++){
        //         for(int j=0;j<16;j++){
        //         printf("%4d, ",img.data[i][j].r);
        //         }
        //         printf("\n");
        // }
        // //���摜����//////////////////////////////////////////////////////////////
        // // int i, j;
        // //lenna.ppm�̉摜�̐F�����[���ɂ���
        // // for( i = 0; i < img.H; i++ ){
        // // 	for( j =0; j < img.W; j++ ){
        // // 		dst.data[i][j].r = img.data[i][j].r;
        // // 		dst.data[i][j].g = img.data[i][j].g;
        // // 		dst.data[i][j].b = 0;
        // // 	}
        // // }
        // two_dimension_fourier_color(img,img_dft,"fft");
        fft_and_ifft_color(img,img_dft,"fft");
        // for(int i=0;i<16;i++){
        //         for(int j=0;j<16;j++){
        //         printf("%4d, ",img_dft.data[i][j].r);
        //         }
        //         printf("\n");
        // }
        // char color_dft_path[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\results\\two_dft_lenna.ppm";
        char color_dft_path[] = "./results/two_dft_lenna.ppm";
        img_dft.save( color_dft_path );

        // inverse_two_dimension_fourier_color(img_dft,img_idft,"ifft");
        // for(int i=0;i<16;i++){
        //         for(int j=0;j<16;j++){
        //         printf("%4d, ",img_idft.data[i][j].r);
        //         }
        //         printf("\n");
        // }
        // char color_idft_path[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\results\\two_idft_lenna.ppm";
        // img_idft.save( color_idft_path );
        return 0;
}
