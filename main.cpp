#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#pragma warning(disable : 4996)
#include "ImageIO.h"

#define pi 3.141592653
#define SIZE 4096

//グレースケール画像の二値化を行う
//この関数の引数は，
//第一引数：Image<GRAY>型の入力画像(&が付いていないときは，コピーが渡されるため，上書きされない)
//第二引数：Image<GRAY>型の出力画像(&が付いているときは，上書きされる)
//第三引数：int 型の変数　閾値
//GRA
void Binarization(Image<GRAY> src, Image<GRAY> &dst, int threshold)
{
        int i, j;
        for (i = 0; i < src.H; i++)
        {
                for (j = 0; j < src.W; j++)
                {
                        if (src.data[i][j] < threshold) //GRAY型の場合は .r , .g , .b は使わない．
                                dst.data[i][j] = 0;
                        else
                                dst.data[i][j] = 255;
                }
        }
}

// one_dimension_fourier
double dft(double re[], double im[], double Re[], double Im[], int N)
{
        //実数部分と虚数部分に分けてフーリエ変換
        for (int k = 0; k < N; k++)
        {
                for (int n = 0; n < N; n++)
                {
                        Re[k] += (re[n] * cos(2 * pi * k * n / N) + im[n] * sin(2 * pi * k * n / N));
                        Im[k] += (-1 * re[n] * sin(2 * pi * k * n / N) + im[n] * cos(2 * pi * k * n / N));
                }
        }
}

// invert_one_dimension_fourier
double idft(double Re[], double Im[], double re[], double im[], int N)
{
        //実数部分と虚数部分に分けて逆フーリエ変換
        for (int k = 0; k < N; k++)
        {
                for (int n = 0; n < N; n++)
                {
                        re[k] += (Re[n] * cos(2 * pi * k * n / N) - Im[n] * sin(2 * pi * k * n / N));
                        im[k] += (Re[n] * sin(2 * pi * k * n / N) + Im[n] * cos(2 * pi * k * n / N));
                }
                re[k] /= N;
                im[k] /= N;
        }
}

void two_dimension_fourier(Image<GRAY> src, Image<GRAY> &dst)
{
        double re_arr[SIZE], im_arr[SIZE], Re_tmp[SIZE], Im_tmp[SIZE];
        double Re_arr[SIZE][SIZE], Im_arr[SIZE][SIZE];//こいつやばい
        printf("src.data[0][255] = %d\n",src.data[0][255]);

        // for (int i = 0; i < src.H; i++)
        // {
        //         for (int j = 0; j < src.W; j++)
        //         {
        //                 re_arr[j] = (double)src.data[i][j];
        //                 im_arr[j] = 0.0;
        //         }
        //         dft(re_arr, im_arr, Re_tmp, Im_tmp, src.W);
        //         printf("i = %d, Re_tmp[0] = %d",i,Re_tmp[0]);
        //         for (int j = 0; j < src.W; j++)
        //         {
        //                 Re_arr[i][j] = Re_tmp[j];
        //                 Im_arr[i][j] = Im_tmp[j];
        //         }
        // }
        // for (int i = 0; i < src.W; i++)
        // {
        //         for (int j = 0; j < src.H; j++)
        //         {
        //                 re_arr[j] = Re_arr[j][i];
        //                 im_arr[j] = Im_arr[j][i];
        //         }
        //         dft(re_arr, im_arr, Re_tmp, Im_tmp, src.H);
        //         for (int j = 0; j < src.H; j++)
        //         {
        //                 dst.data[j][i] = (int)Re_tmp[j];
        //                 // dst.data[j][i] = Im_tmp[j];
        //         }
        // }
}

void inverse_two_dimension_fourier(Image<GRAY> src, Image<GRAY> &dst)
{
        double re_arr[SIZE], im_arr[SIZE], Re_arr[SIZE][SIZE], Im_arr[SIZE][SIZE], Re_tmp[SIZE], Im_tmp[SIZE];
        for (int i = 0; i < src.H; i++)
        {
                for (int j = 0; j < src.W; j++)
                {
                        re_arr[j] = (double)src.data[i][j];
                        im_arr[j] = 0.0;
                }
                idft(re_arr, im_arr, Re_tmp, Im_tmp, src.W);
                printf("i = %d, Re_tmp[0] = %d",i,Re_tmp[0]);
                for (int j = 0; j < src.W; j++)
                {
                        Re_arr[i][j] = Re_tmp[j];
                        Im_arr[i][j] = Im_tmp[j];
                }
        }
        for (int i = 0; i < src.W; i++)
        {
                for (int j = 0; j < src.H; j++)
                {
                        re_arr[j] = Re_arr[j][i];
                        im_arr[j] = Im_arr[j][i];
                }
                idft(re_arr, im_arr, Re_tmp, Im_tmp, src.H);
                for (int j = 0; j < src.H; j++)
                {
                        dst.data[j][i] = (int)Re_tmp[j];
                        // dst.data[j][i] = Im_tmp[j];
                }
        }
}

int file_load_and_write(char filename1[], char filename2[],char function[]){
        FILE *fp1, *fp2; // FILE型構造体
        double f1, f2 = 0.0;
        double Re[SIZE], Im[SIZE]; //変換前の配列
        double re[SIZE], im[SIZE]; //変換後の配列
        int i,k,i1,N = 0;
        printf("fname1 = %s, fname2 = %s, function = %s\n",filename1,filename2,function);
        fp1 = fopen(filename1, "r"); // ファイルを開く。失敗するとNULLを返す。
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
        printf("Re[%d] = %f, Im[%d] = %f\n", N - 1, Re[N - 1], N - 1, Im[N - 1]);
        fclose(fp1); // ファイルを閉じる

        if ((fp2 = fopen(filename2, "w")) == NULL)
        {
                printf("FILE2 not open\n");
                return -1;
        }
        //実数部分と虚数部分に分けてフーリエ変換または逆フーリエ変換
        if(strcmp(function,"dft")==0){
                dft(Re, Im, re, im, N);
        }else if(strcmp(function,"idft")==0){
                idft(Re, Im, re, im, N);
        }
        printf("re[%d] = %f, rm[%d] = %f\n", N - 1, re[N - 1], N - 1, im[N - 1]);
        for (k = 0; k < N; k++)
        {
                fprintf(fp2, "%d, %f, %f\n", k, re[k], im[k]);
        }
        fclose(fp2);
        return 0;
}

int main(void)
{
        // 1次元離散フーリエ変換
        printf("ready...\n");
        // file_load_and_write("test.csv", "inv_fourier.csv","dft");
        // printf("1dim dft is completed...\n");
        // file_load_and_write("inv_fourier.csv", "check.csv","idft");
        // printf("1dim idft is completed...\n");
        // excelなどでグラフを書いて確認をしてみてください。
        // char path1[SIZE],blue_save_path[SIZE],gray_save_path[SIZE];

        // 2次元離散フーリエ変換
        //グレー画像
        //■配列の宣言
        Image<GRAY> gray, gout;
        char path2[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\pictures\\lenna.pgm";
        //■画像ロード
        gray.load(path2); //lenna.pgmをgrayにロードする
        //■出力配列確保
        gout.create(gray.W, gray.H); // goutに img.W(画像imgの横幅) x img.H(画像imgの縦幅) の大きさの画素配列を用意する
        //■画像処理//////////////////////////////////////////////////////////////
        // Binarization(gray, gout, 128); //grayを閾値128で二値化してgoutに出力
        // for(int i=0;i<gray.H;i++){
        //         for(int j=0;j<gray.W;j++){
        //                 printf("")
        //         }
        // }
        for(int i=0;i<gray.H;i++){
                printf("gray.data[0][%d] = %d\n",i,gray.data[0][i]);
        }

        two_dimension_fourier(gray,gout);

        for(int i=0;i<gout.H;i++){
                printf("gout.data[0][%d] = %d\n",i,gout.data[0][i]);
        }
        char gout_path[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\results\\two_dft_lenna.pgm";
        //■保存
        gout.save(gout_path); // Binarization.pgm という名前でgoutを保存


        // Image<COLOR> img,dst,nega;	//カラー画像
        // char path1[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\pictures\\lenna.ppm";
        // img.load( path1 );	//lenna.ppmをimgにロードする

        // dst.create( img.W, img.H );  // dstに　img.W(画像imgの横幅) x img.H(画像imgの縦幅) の大きさの画素配列を用意する
        // nega.create( img.W, img.H ); // negaに img.W(画像imgの横幅) x img.H(画像imgの縦幅) の大きさの画素配列を用意する

        //■画像処理//////////////////////////////////////////////////////////////

        // int i, j;
        //lenna.ppmの画像の青色だけゼロにする

        // for( i = 0; i < img.H; i++ ){
        // 	for( j =0; j < img.W; j++ ){
        // 		dst.data[i][j].r = img.data[i][j].r;
        // 		dst.data[i][j].g = img.data[i][j].g;
        // 		dst.data[i][j].b = 0;
        // 	}
        // }

        // //■保存
        // char blue_save_path[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\results\\blueZero.ppm";
        // char gray_save_path[] = "C:\\Users\\ibuki\\program\\c\\ImageIO\\results\\Binarization.pgm";
        // dst.save( blue_save_path ); // blueZero.ppm という名前でdstを保存    
        return 0;
}
