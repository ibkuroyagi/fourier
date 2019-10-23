#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double arr_test(double arr1[],double arr2[],int N){
  for(int i=0;i<N;i++){
    arr2[i] = arr1[i]*2+0.01;
  }
}

double make_arr(double arr1[],double arr2[]){
  double a[16][16];
  for(int i=0;i<10;i++){
    for(int j=0;j<10;j++){
      a[i][j] = 20.0*(double)i + (double)j + 0.101;
    }
  }  
  for(int k=0;k<10;k++){
    arr2[k] = a[k][1];
  }
}


//nが2の累乗かどうかを確認する
//nが2の累乗でなければそれを2の累乗に切り上げた数値を返す
int isPowerOfTwo(int n) {
        int pow = 2;
        if ((n & (n - 1)) == 0){
                return 1;
        }else{
                while(n > pow){
                        pow *= 2;
                }
                return pow;
        }
}


int padding_arr(double *arr, int n1, int n2){
  for(int i=n1;i<n2;i++){
    arr[i] = 0.0;
  }
}

int fft(int n, int flg, double *ar, double *ai)
{ 
  int padding_n = isPowerOfTwo(n);
  padding_arr(ar,n, padding_n);
}

int main(void) {
  int n1 = 16;
  int n2 = 30;
  printf("isPowerOfTwo(%d) == %d\n",n1,isPowerOfTwo(n1));
  printf("isPowerOfTwo(%d) == %d\n",n2,isPowerOfTwo(n2));
  
  double arr1[64],arr2[64];
  // double a[16][16];

  for(int i=0;i<30;i++){
    arr1[i] = i;
  }
  fft(n2,1,arr1,arr2);
  // arr_test(arr1,arr2,16);
  for(int i=0;i<36;i++){
    printf("arr1[%d] = %lf, arr2[%d] = %lf\n",i,arr1[i],i,arr2[i]);
  }
  
//  for(int i=0;i<10;i++){
//     for(int j=0;j<10;j++){
//       a[i][j] = 10.0*(double)i + (double)j;
//       // a[i][j] = 1;
//     }    
//   }
//   for(int i=0;i<10;i++){
//     for(int j=0;j<10;j++){
//       printf("%lf, ",a[i][j]);
//     }
//     printf("\n");
//   }
//   make_arr(arr1,arr2);
//   printf("%lf\n",arr2[1]); 
  return 0;
}

