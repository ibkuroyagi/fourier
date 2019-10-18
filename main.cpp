#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#pragma warning (disable:4996)
// #include"ImageIO.h"

#define pi 3.141592653


//グレースケール画像の二値化を行う
//この関数の引数は，
//第一引数：Image<GRAY>型の入力画像(&が付いていないときは，コピーが渡されるため，上書きされない)
//第二引数：Image<GRAY>型の出力画像(&が付いているときは，上書きされる)
//第三引数：int 型の変数　閾値
//GRA
// void Binarization( Image<GRAY> src, Image<GRAY> &dst ,int threshold)
// {
// 	int i,j;
// 	for( i=0; i <src.H; i++ ){
// 		for( j = 0; j < src.W; j++ ){
// 			if( src.data[i][j] < threshold ) //GRAY型の場合は .r , .g , .b は使わない．
// 				dst.data[i][j] = 0;
// 			else
// 				dst.data[i][j] = 255;
// 		}
// 	}
// }

// void one_dimension_fourier(double *f_re, double *f_im, double *F_re,double *F_im, int N){
//   // write your code
// }

void inverse_one_dimension_fourier(double *F_re, double *F_im, double *f_re,double *f_im, int N){
  // write your code
}

void two_dimension_fourier(double **f_re, double **f_im, double **F_re,double **F_im, int M, int N){
  // write your code
}

void inverse_two_dimension_fourier(double **F_re, double **F_im, double **f_re,double **f_im, int M, int N){
  // write your code
}

 // 離散フーリエ変換（読み込むファイル名, 書き込むファイル名）
int one_dimension_fourier(char filename1[], char filename2[])
{  
        int i, k, n, N = 0;
        int max = 1024;  // 読み込むデータ数の上限
        double Re[max+1], Im[max+1], re, im;

        FILE *fp1, *fp2;
    // ファイルオープン(フーリエ変換したいデータファイル, フーリエ変換後のデータ保存用ファイル)
		printf("f1 open\n");
        if((fp1=fopen(filename1,"r"))==NULL){
                printf("FILE1 not open\n");
        return -1;
        }
		printf("f2 open\n");
        if((fp2=fopen(filename2,"w"))==NULL){
                printf("FILE2 not open\n");
        return -1;
        }
		printf("データを読み込む\n");
		int ret;
        //データの読み込み
				while( ( fscanf( fp1, "%d,%lf,%lf\n", &i, &Re[N], &Im[N])) == NULL ){
          printf( " Re = %lf, Im =  %lf \n", Re[N], Im[N] );
					N++;
        }

        // for(N=0;N<max;N++) {
        //     if(fscanf(fp1,"%d, %lf, %lf", &i, &Re[N], &Im[N]) == EOF) {
        //       //  N--;
				// 			printf("i = %d, Re = %lf\n",i,Re[N]);
        //       break;
        //     }
        // }
			printf("red data\n");			
    	printf("Re = %lf\n",Re[1000]);
			printf("im = %lf\n",Im[1000]);

        //実数部分と虚数部分に分けてフーリエ変換
        for(k=0; k<N; k++) {
                double re = 0.0;
                double im = 0.0;
                for(n=0; n<N; n++) {
                        re += (Re[n]*cos(2*pi*k*n/N)-Im[n]*sin(2*pi*k*n/N))/N;
                        im += (Re[n]*sin(2*pi*k*n/N)+Im[n]*cos(2*pi*k*n/N))/N; 
                }
                fprintf(fp2,"%d, %f, %f\n", k, re, im);
        }
        fclose(fp1);
        fclose(fp2);
        return 0;
}

int main()
{
  // 1次元離散フーリエ変換
	printf("prepair for 1dim dft\n");
  one_dimension_fourier("fourier.csv", "inv_fourier.csv");
	printf("1dim dft is completed...\n");
	return 0;

	// Image<GRAY> gray, Amp,Phase;		//グレー画像
	// gray.load("lenna.pgm");	//lenna.pgmをgrayにロードする
	

//	//■配列の宣言
//	Image<COLOR> img,dst,nega;	//カラー画像
//	Image<GRAY> gray, gout;		//グレー画像
//
//	//■画像ロード
//	img.load( "lenna.ppm" );	//lenna.ppmをimgにロードする
//	gray.load( "lenna.pgm" );	//lenna.pgmをgrayにロードする
//	
//	//■配列確保
//	dst.create( img.W, img.H );  // dstに　img.W(画像imgの横幅) x img.H(画像imgの縦幅) の大きさの画素配列を用意する
//	nega.create( img.W, img.H ); // negaに img.W(画像imgの横幅) x img.H(画像imgの縦幅) の大きさの画素配列を用意する
//	gout.create( gray.W, gray.H );// goutに img.W(画像imgの横幅) x img.H(画像imgの縦幅) の大きさの画素配列を用意する
//
//
//
//
//	//■画像処理//////////////////////////////////////////////////////////////
//
//	int i,j;
//
//	//lenna.ppmの画像の青色だけゼロにする
//
//	for( i = 0; i < img.H; i++ ){
//		for( j =0; j < img.W; j++ ){
//			dst.data[i][j].r = img.data[i][j].r;
//			dst.data[i][j].g = img.data[i][j].g;
//			dst.data[i][j].b = 0;
//		}
//	}
//	
//	//■保存
//	dst.save( "blueZero.ppm" ); // blueZero.ppm という名前でdstを保存
//
//	
//	//■画像処理2/////////////////////////////////////////////////////////////////////////
//	
//	Binarization( gray, gout , 128 ); //grayを閾値128で二値化してgoutに出力
//	
//	//■保存
//	gout.save( "Binarization.pgm" ); // Binarization.pgm という名前でgoutを保存
//

//	Image<GRAY> gray, Amp,Phase;		//グレー画像
//	gray.load( "lenna.pgm" );	//lenna.pgmをgrayにロードする
	
}

