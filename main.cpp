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
int isPowerOfTwo(int n) {
    return (n & (n - 1)) == 0;
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

        if(isPowerOfTwo(n) == 0){
                printf("If you want to use fft or ifft, resize the image.\n");
                printf("It must be (2^n, 2^n)\n");
                return 0;
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

void two_dimension_fourier(Image<GRAY> src, Image<GRAY> &dst, char function[])
{
        double re_arr[SIZE], im_arr[SIZE], Re_tmp[SIZE], Im_tmp[SIZE];
        double Re_arr[256][256], Im_arr[256][256]; //������΂�
        // printf("src.data[0][255] = %d\n",src.data[0][255]);

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
                                // dst.data[j][i] = Im_tmp[j];
                        }
                }
                else if (strcmp(function, "fft") == 0)
                {
                        fft(src.H, -1, re_arr, im_arr);
                        // printf("i = %d, cross_re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                dst.data[j][i] = (int)re_arr[j];
                        }
                }
        }
}

void inverse_two_dimension_fourier(Image<GRAY> src, Image<GRAY> &dst, char function[])
{
        double re_arr[SIZE], im_arr[SIZE], Re_tmp[SIZE], Im_tmp[SIZE];
        double Re_arr[256][256], Im_arr[256][256];
        for (int i = 0; i < src.H; i++)
        {
                for (int j = 0; j < src.W; j++)
                {
                        re_arr[j] = (double)src.data[j][i];
                        im_arr[j] = 0.0;
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
                                // dst.data[j][i] = Im_tmp[j];
                        }
                }
                else if (strcmp(function, "ifft") == 0)
                {
                        fft(src.H, 1, re_arr, im_arr);
                        // printf("i = %d, re_arr[0] = %lf\n",i,re_arr[0]);
                        for (int j = 0; j < src.W; j++)
                        {
                                dst.data[i][j] = (int)re_arr[j];
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

int main(void)
{
        printf("ready...\n");
        // 1�������U�t�[���G�ϊ�
        // file_load_and_write("test.csv", "inv_fourier.csv","dft");
        // printf("1dim dft is completed...\n");

        // 1�����������U�t�[���G�ϊ�
        file_load_and_write("test.csv", "inv_fourier.csv", "fft");
        printf("1dim fft is completed...\n");

        // 1�����t���U�t�[���G�ϊ�
        // file_load_and_write("inv_fourier.csv", "check.csv","idft");
        // printf("1dim idft is completed...\n");

        // 1�����t�������U�t�[���G�ϊ�
        file_load_and_write("inv_fourier.csv", "check.csv", "ifft");
        printf("1dim ifft is completed...\n");
        // excel�ȂǂŃO���t�������Ċm�F�����Ă݂Ă��������B

        // char path1[SIZE],blue_save_path[SIZE],gray_save_path[SIZE];

        // 2�������U�t�[���G�ϊ�
        //�O���[�摜
        //���z��̐錾
        Image<GRAY> gray, gout_dft, gout_idft;
        // �摜�̃p�X���e���̊��ɕύX�����Ă��������B
        char path2[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\pictures\\lenna.pgm";
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
        gout_idft.create(gray.W, gray.H);
        //���摜����//////////////////////////////////////////////////////////////
        // Binarization(gray, gout, 128); //gray��臒l128�œ�l������gout�ɏo��

        // 1�������U�t�[���G�ϊ�
        // two_dimension_fourier(gray,gout_dft,"dft");

        // 1�����������U�t�[���G�ϊ�
        two_dimension_fourier(gray, gout_dft, "fft");
        //�摜�̈ꕔ���o��
        // for(int i=0;i<16;i++){
        //         for(int j=0;j<16;j++){
        //         printf("%4d, ",gout_dft.data[i][j]);
        //         }
        //         printf("\n");
        // }
        char gout_dft_path[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\results\\two_dft_lenna.pgm";
        gout_dft.save(gout_dft_path);
        // 1�����t���U�t�[���G�ϊ�
        // inverse_two_dimension_fourier(gout_dft,gout_idft,"idft");

        // 1�����t�������U�t�[���G�ϊ�
        inverse_two_dimension_fourier(gout_dft, gout_idft, "ifft");
        //�摜�̈ꕔ���o��
        // for(int i=0;i<16;i++){
        //         for(int j=0;j<16;j++){
        //         // printf("%4d, ",gray.data[i][j] - gout_idft.data[i][j]);
        //         printf("%4d, ",gout_idft.data[i][j]);
        //         }
        //         printf("\n");
        // }

        char gout_idft_path[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\results\\two_idft_lenna.pgm";
        gout_idft.save(gout_idft_path);

        // Image<COLOR> img,dst,nega;	//�J���[�摜
        // char path1[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\pictures\\lenna.ppm";
        // img.load( path1 );	//lenna.ppm��img�Ƀ��[�h����

        // dst.create( img.W, img.H );  // dst�Ɂ@img.W(�摜img�̉���) x img.H(�摜img�̏c��) �̑傫���̉�f�z���p�ӂ���
        // nega.create( img.W, img.H ); // nega�� img.W(�摜img�̉���) x img.H(�摜img�̏c��) �̑傫���̉�f�z���p�ӂ���

        //���摜����//////////////////////////////////////////////////////////////

        // int i, j;
        //lenna.ppm�̉摜�̐F�����[���ɂ���

        // for( i = 0; i < img.H; i++ ){
        // 	for( j =0; j < img.W; j++ ){
        // 		dst.data[i][j].r = img.data[i][j].r;
        // 		dst.data[i][j].g = img.data[i][j].g;
        // 		dst.data[i][j].b = 0;
        // 	}
        // }

        // //���ۑ�
        // char blue_save_path[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\results\\blueZero.ppm";
        // char gray_save_path[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\results\\Binarization.pgm";
        // dst.save( blue_save_path ); // blueZero.ppm �Ƃ������O��dst��ۑ�
        return 0;
}
