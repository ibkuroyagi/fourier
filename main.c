#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include<iostream>
#include <math.h>
#pragma warning(disable : 4996)
#include"ImageIO.h"

#define pi 3.141592653

//�O���[�X�P�[���摜�̓�l�����s��
//���̊֐��̈����́C
//�������FImage<GRAY>�^�̓��͉摜(&���t���Ă��Ȃ��Ƃ��́C�R�s�[���n����邽�߁C�㏑������Ȃ�)
//�������FImage<GRAY>�^�̏o�͉摜(&���t���Ă���Ƃ��́C�㏑�������)
//��O�����Fint �^�̕ϐ��@臒l
//GRA
void Binarization( Image<GRAY> src, Image<GRAY> &dst ,int threshold)
{
	int i,j;
	for( i=0; i <src.H; i++ ){
		for( j = 0; j < src.W; j++ ){
			if( src.data[i][j] < threshold ) //GRAY�^�̏ꍇ�� .r , .g , .b �͎g��Ȃ��D
				dst.data[i][j] = 0;
			else
				dst.data[i][j] = 255;
		}
	}
}



void two_dimension_fourier(double **f_re, double **f_im, double **F_re, double **F_im, int M, int N)
{
        // write your code
}

void inverse_two_dimension_fourier(double **F_re, double **F_im, double **f_re, double **f_im, int M, int N)
{
        // write your code
}

// ���U�t�[���G�ϊ��i�ǂݍ��ރt�@�C����, �������ރt�@�C�����j
int one_dimension_fourier(char filename1[], char filename2[])
{
        int i = 0, k, n, N = 0,i1;
        int max = 100000; // �ǂݍ��ރf�[�^���̏��
        double Re[max + 1], Im[max + 1], re, im;
        FILE *fp1,*fp2; // FILE�^�\����
        float f1,f2 = 0.0;

        fp1 = fopen(filename1, "r"); // �t�@�C�����J���B���s�����NULL��Ԃ��B
        if (fp1 == NULL)
        {
                printf("%s file not open!\n", filename1);
                exit(1);
        }

        while (fscanf(fp1, "%d,%f,%f", &i1,&f1,&f2) != EOF)
        {
                printf("i1 = %d,f1 = %f, f2 = %f\n",i1,f1,f2);
                Re[N] = f1;
                Im[N] = f2;
                N++;
        }
        N = N - 1;
        printf("Re[%d] = %f, Im[%d] = %f\n", N, Re[N],N,Im[N]);
        fclose(fp1); // �t�@�C�������
        // �t�@�C���I�[�v��(�t�[���G�ϊ��������f�[�^�t�@�C��, �t�[���G�ϊ���̃f�[�^�ۑ��p�t�@�C��)
        if ((fp2 = fopen(filename2, "w")) == NULL)
        {
                printf("FILE2 not open\n");
                return -1;
        }
   
        //���������Ƌ��������ɕ����ăt�[���G�ϊ�
        for (k = 0; k <= N; k++)
        {
                double re = 0.0;
                double im = 0.0;
                for (n = 0; n < N; n++)
                {
                        re += (Re[n] * cos(2 * pi * k * n / N) + Im[n] * sin(2 * pi * k * n / N)) ;
                        im += (-1*Re[n] * sin(2 * pi * k * n / N) + Im[n] * cos(2 * pi * k * n / N)) ;
                }
                fprintf(fp2, "%d, %f, %f\n", k, re, im);
        }
        fclose(fp2);
        return 0;
}

int inverse_one_dimension_fourier(char filename1[], char filename2[])
{
        int i = 0, k, n, N = 0,i1;
        int max = 100000; // �ǂݍ��ރf�[�^���̏��
        double Re[max + 1], Im[max + 1], re, im;
        FILE *fp1,*fp2; // FILE�^�\����
        float f1,f2 = 0.0;

        fp1 = fopen(filename1, "r"); // �t�@�C�����J���B���s�����NULL��Ԃ��B
        if (fp1 == NULL)
        {
                printf("%s file not open!\n", filename1);
                exit(1);
        }

        while (fscanf(fp1, "%d,%f,%f", &i1,&f1,&f2) != EOF)
        {
                printf("i1 = %d,f1 = %f, f2 = %f\n",i1,f1,f2);
                Re[N] = f1;
                Im[N] = f2;
                N++;
        }
        N = N - 1;
        printf("Re[%d] = %f, Im[%d] = %f\n", N, Re[N],N,Im[N]);
        fclose(fp1); // �t�@�C�������
        // �t�@�C���I�[�v��(�t�[���G�ϊ��������f�[�^�t�@�C��, �t�[���G�ϊ���̃f�[�^�ۑ��p�t�@�C��)
        if ((fp2 = fopen(filename2, "w")) == NULL)
        {
                printf("FILE2 not open\n");
                return -1;
        }
   
        //���������Ƌ��������ɕ����ċt�t�[���G�ϊ�
        for (k = 0; k <= N; k++)
        {
                double re = 0.0;
                double im = 0.0;
                for (n = 0; n < N; n++)
                {
                        re += (Re[n] * cos(2 * pi * k * n / N) - Im[n] * sin(2 * pi * k * n / N)) ;
                        im += (Re[n] * sin(2 * pi * k * n / N) + Im[n] * cos(2 * pi * k * n / N)) ;
                }
                fprintf(fp2, "%d, %f, %f\n", k, re/N, im/N);
        }
        fclose(fp2);
        return 0;
}

int main()
{
        // 1�������U�t�[���G�ϊ�
        // test.csv�̒��g��
        // t = [1:1024]*0.01
        one_dimension_fourier("test.csv", "inv_fourier.csv");
        printf("1dim dft is completed...\n");
        inverse_one_dimension_fourier("inv_fourier.csv","check.csv");
        printf("1dim idft is completed...\n");
        // excel�ȂǂŃO���t�������Ċm�F�����Ă݂Ă��������B
        return 0;

        // Image<GRAY> gray, Amp,Phase;		//�O���[�摜
        // gray.load("lenna.pgm");	//lenna.pgm��gray�Ƀ��[�h����

        //	//���z��̐錾
        //	Image<COLOR> img,dst,nega;	//�J���[�摜
        //	Image<GRAY> gray, gout;		//�O���[�摜
        //
        //	//���摜���[�h
        //	img.load( "lenna.ppm" );	//lenna.ppm��img�Ƀ��[�h����
        //	gray.load( "lenna.pgm" );	//lenna.pgm��gray�Ƀ��[�h����
        //
        //	//���z��m��
        //	dst.create( img.W, img.H );  // dst�Ɂ@img.W(�摜img�̉���) x img.H(�摜img�̏c��) �̑傫���̉�f�z���p�ӂ���
        //	nega.create( img.W, img.H ); // nega�� img.W(�摜img�̉���) x img.H(�摜img�̏c��) �̑傫���̉�f�z���p�ӂ���
        //	gout.create( gray.W, gray.H );// gout�� img.W(�摜img�̉���) x img.H(�摜img�̏c��) �̑傫���̉�f�z���p�ӂ���
        //
        //
        //
        //
        //	//���摜����//////////////////////////////////////////////////////////////
        //
        //	int i,j;
        //
        //	//lenna.ppm�̉摜�̐F�����[���ɂ���
        //
        //	for( i = 0; i < img.H; i++ ){
        //		for( j =0; j < img.W; j++ ){
        //			dst.data[i][j].r = img.data[i][j].r;
        //			dst.data[i][j].g = img.data[i][j].g;
        //			dst.data[i][j].b = 0;
        //		}
        //	}
        //
        //	//���ۑ�
        //	dst.save( "blueZero.ppm" ); // blueZero.ppm �Ƃ������O��dst��ۑ�
        //
        //
        //	//���摜����2/////////////////////////////////////////////////////////////////////////
        //
        //	Binarization( gray, gout , 128 ); //gray��臒l128�œ�l������gout�ɏo��
        //
        //	//���ۑ�
        //	gout.save( "Binarization.pgm" ); // Binarization.pgm �Ƃ������O��gout��ۑ�
        //

        //	Image<GRAY> gray, Amp,Phase;		//�O���[�摜
        //	gray.load( "lenna.pgm" );	//lenna.pgm��gray�Ƀ��[�h����
}
